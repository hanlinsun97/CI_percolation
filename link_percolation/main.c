#include "graph.h"
#include "output.h"
#include "removal.h"
#include "rng.h"
#include "stats.h"
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MILESTONE 10
#define REPORT_TIME_INTERVAL 3600

typedef enum { GRAPH_INPUT, GRAPH_ER } GraphMode;

static void print_usage(const char *prog) {
  printf("Usage: %s -g <input|ER> [options]\n", prog);
  printf("Required:\n");
  printf("  -g TYPE       Graph type: input|ER\n");
  printf("If -g input:  -f FILE\n");
  printf("If -g ER:     -N NUM -c AVG_DEGREE\n");
  printf("Options:\n");
  printf("  -M_graphs NUM       Number of graphs (default 1)\n");
  printf("  -M_realizations NUM Realizations per graph (default 1)\n");
  printf("  -n_threads NUM      OpenMP threads (default auto)\n");
  printf("  -o PREFIX           Output prefix (default LR)\n");
  printf("  -r SEED             Random seed (default time)\n");
}

static const char *path_basename(const char *path) {
  const char *slash = strrchr(path, '/');
  return slash ? slash + 1 : path;
}

int main(int argc, char **argv) {
  GraphMode mode = GRAPH_INPUT;
  const char *input_file = NULL;
  int N = 0;
  double avg_degree = 0.0;
  int M_graphs = 1;
  int M_realizations = 1;
  char output_prefix[128] = "LR";
  unsigned long seed = (unsigned long)time(NULL);
  int num_threads = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      return 0;
    } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
      const char *type = argv[++i];
      if (strcmp(type, "input") == 0) {
        mode = GRAPH_INPUT;
      } else if (strcmp(type, "ER") == 0) {
        mode = GRAPH_ER;
      } else {
        fprintf(stderr, "Unknown graph type '%s'\n", type);
        return 1;
      }
    } else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
      input_file = argv[++i];
    } else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
      avg_degree = atof(argv[++i]);
    } else if (strcmp(argv[i], "-M_graphs") == 0 && i + 1 < argc) {
      M_graphs = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-M_realizations") == 0 && i + 1 < argc) {
      M_realizations = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      strncpy(output_prefix, argv[++i], sizeof(output_prefix) - 1);
      output_prefix[sizeof(output_prefix) - 1] = '\0';
    } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      seed = strtoul(argv[++i], NULL, 10);
    } else if (strcmp(argv[i], "-n_threads") == 0 && i + 1 < argc) {
      num_threads = atoi(argv[++i]);
    } else {
      fprintf(stderr, "Unknown/incomplete option '%s'\n", argv[i]);
      print_usage(argv[0]);
      return 1;
    }
  }

  if (mode == GRAPH_INPUT && !input_file) {
    fprintf(stderr, "Input graph required for -g input\n");
    return 1;
  }
  if (mode == GRAPH_ER && (N <= 0 || avg_degree <= 0.0)) {
    fprintf(stderr, "-N and -c must be >0 for ER\n");
    return 1;
  }
  if (M_graphs <= 0 || M_realizations <= 0) {
    fprintf(stderr, "M_graphs and M_realizations must be positive\n");
    return 1;
  }

#ifdef _OPENMP
  if (num_threads > 0) {
    omp_set_num_threads(num_threads);
  }
#else
  (void)num_threads;
#endif

  srand((unsigned int)seed);

  Graph *input_graph = NULL;
  int graph_nodes = 0;
  int graph_edges = 0;
  if (mode == GRAPH_INPUT) {
    input_graph = read_graph_from_file(input_file);
    if (!input_graph) {
      fprintf(stderr, "Unable to load input graph\n");
      return 1;
    }
    graph_nodes = input_graph->n;
    graph_edges = input_graph->m;
  } else {
    graph_nodes = N;
  }

  const int series_length = NR_SERIES_LENGTH;
  OnlineStats *gc_stats = create_online_stats(series_length);
  OnlineStats *small_stats = create_online_stats(series_length);
  OnlineStats *lambda_stats = create_online_stats(series_length);
  OnlineStats *comp_stats = create_online_stats(series_length);
  OnlineStats *cycles_stats = create_online_stats(series_length);
  OnlineStats *core_stats = create_online_stats(series_length);
  if (!gc_stats || !small_stats || !lambda_stats || !comp_stats ||
      !cycles_stats || !core_stats) {
    fprintf(stderr, "Unable to allocate statistics buffers\n");
    return 1;
  }

  int total_runs = M_graphs * M_realizations;
  char ebe_file[512];
  char averaged_file[512];
  const char *graph_label =
      (mode == GRAPH_ER) ? "ER" : path_basename(input_file);
  if (mode == GRAPH_ER) {
    snprintf(averaged_file, sizeof(averaged_file),
             "%s_LR_ER_N%d_c%.2f_L0_Mg%d_Mr%d.dat", output_prefix,
             graph_nodes, avg_degree, M_graphs, M_realizations);
    snprintf(ebe_file, sizeof(ebe_file),
             "%s_EBE_LR_ER_N%d_c%.2f_L0_Mg%d_Mr%d.dat", output_prefix,
             graph_nodes, avg_degree, M_graphs, M_realizations);
  } else {
    snprintf(averaged_file, sizeof(averaged_file),
             "%s_LR_INPUT_%s_L0_Mg%d_Mr%d.dat", output_prefix,
             graph_label, M_graphs, M_realizations);
    snprintf(ebe_file, sizeof(ebe_file),
             "%s_EBE_LR_INPUT_%s_L0_Mg%d_Mr%d.dat", output_prefix,
             graph_label, M_graphs, M_realizations);
  }

  EBEWriter *ebe_writer = ebe_writer_begin(ebe_file, graph_nodes, avg_degree,
                                           M_graphs, M_realizations,
                                           total_runs);
  if (!ebe_writer) {
    fprintf(stderr, "Unable to open EBE output\n");
    return 1;
  }

  volatile int error_flag = 0;
  int completed_runs = 0;

#ifdef _OPENMP
  double last_print_time = omp_get_wtime();
#else
  clock_t last_print_time = clock();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    double *local_gc = (double *)malloc((size_t)series_length * sizeof(double));
    double *local_small =
        (double *)malloc((size_t)series_length * sizeof(double));
    double *local_lambda =
        (double *)malloc((size_t)series_length * sizeof(double));
    double *local_comp =
        (double *)malloc((size_t)series_length * sizeof(double));
    double *local_cycles =
        (double *)malloc((size_t)series_length * sizeof(double));
    double *local_core =
        (double *)malloc((size_t)series_length * sizeof(double));
    LinkPercolationResult *thread_result = NULL;

    if (!local_gc || !local_small || !local_lambda || !local_comp ||
        !local_cycles || !local_core || !thread_result) {
#ifdef _OPENMP
#pragma omp critical
#endif
      { error_flag = 1; }
    } else {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (int graph_idx = 0; graph_idx < M_graphs; graph_idx++) {
        if (error_flag)
          continue;
        Graph *current_graph = NULL;
        Graph *owned_graph = NULL;
        if (mode == GRAPH_INPUT) {
          current_graph = input_graph;
        } else {
          rng_state_t graph_rng;
          uint64_t graph_seed =
              seed ^ (0xbf58476d1ce4e5b9ULL * (uint64_t)(graph_idx + 1));
          owned_graph = generate_erdos_renyi_rng(N, avg_degree, &graph_rng);
          current_graph = owned_graph;
          if (!current_graph) {
#ifdef _OPENMP
#pragma omp critical
#endif
            { error_flag = 1; }
            if (owned_graph)
              free_graph(owned_graph);
            continue;
          }
        }

        thread_result = create_link_percolation_result(current_graph->n,
                                                       current_graph->m);
        if (!thread_result) {
#ifdef _OPENMP
#pragma omp critical
#endif
          { error_flag = 1; }
        }

        for (int realization = 0;
             realization < M_realizations && !error_flag; realization++) {
          int run = graph_idx * M_realizations + realization;
          reset_link_percolation_result(thread_result);
          uint64_t run_seed =
              seed ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(run + 1));
          PseudoCriticalPoints gc_ecp;
          if (run_link_percolation(current_graph, run_seed, thread_result,
                                   &gc_ecp) != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
            { error_flag = 1; }
            continue;
          }

          int sample_length = thread_result->series_length;
          double norm_n = (double)current_graph->n;
          for (int i = 0; i < sample_length; i++) {
            local_gc[i] = thread_result->giant_size[i] / norm_n;
            local_small[i] = thread_result->avg_small_size[i];
            local_lambda[i] = thread_result->lambda1[i];
            local_comp[i] = (double)thread_result->num_components[i];
            local_cycles[i] = (double)thread_result->num_cycles[i];
            local_core[i] = thread_result->core_size[i] / norm_n;
          }

          EnhancedEBEData record = {0};
          record.p_plus = gc_ecp.p_plus;
          record.p_plus_jump = gc_ecp.p_plus_jump;
          record.q_plus = gc_ecp.q_plus;
          record.q_zero = gc_ecp.q_zero;
          record.q_minus = gc_ecp.q_minus;
          record.p1 = gc_ecp.p1;
          record.p2 = gc_ecp.p2;
          record.P_p_plus = gc_ecp.P_at_p_plus;
          record.s_p_plus = gc_ecp.s_at_p_plus;
          record.nc_p_plus = gc_ecp.nc_at_p_plus;
          record.nl_p_plus = gc_ecp.nl_at_p_plus;
          record.P_q_plus = gc_ecp.P_at_q_plus;
          record.s_q_plus = gc_ecp.s_at_q_plus;
          record.nc_q_plus = gc_ecp.nc_at_q_plus;
          record.nl_q_plus = gc_ecp.nl_at_q_plus;
          record.P_q_zero = gc_ecp.P_at_q_zero;
          record.s_q_zero = gc_ecp.s_at_q_zero;
          record.nc_q_zero = gc_ecp.nc_at_q_zero;
          record.nl_q_zero = gc_ecp.nl_at_q_zero;
          record.P_q_minus = gc_ecp.P_at_q_minus;
          record.s_q_minus = gc_ecp.s_at_q_minus;
          record.nc_q_minus = gc_ecp.nc_at_q_minus;
          record.nl_q_minus = gc_ecp.nl_at_q_minus;
          record.core_p_star = thread_result->core_p_star;
          record.core_drop = thread_result->core_drop;
          record.core_size_star = thread_result->core_size_at_star;
          record.core_zero_p = thread_result->core_zero_p;

#ifdef _OPENMP
#pragma omp critical
#endif
          {
            if (!error_flag) {
              online_stats_update(gc_stats, local_gc);
              online_stats_update(small_stats, local_small);
              online_stats_update(lambda_stats, local_lambda);
              online_stats_update(comp_stats, local_comp);
              online_stats_update(cycles_stats, local_cycles);
              online_stats_update(core_stats, local_core);
              ebe_writer_record(ebe_writer, run, &record);
              completed_runs++;

#ifdef _OPENMP
              double now = omp_get_wtime();
              double elapsed = now - last_print_time;
#else
              clock_t now = clock();
              double elapsed =
                  (double)(now - last_print_time) / CLOCKS_PER_SEC;
#endif
              int progress_interval =
                  total_runs >= MILESTONE ? total_runs / MILESTONE : 1;
              bool at_milestone = (completed_runs % progress_interval) == 0;
              bool is_final = (completed_runs == total_runs);
              if (at_milestone || is_final ||
                  elapsed >= REPORT_TIME_INTERVAL) {
                printf("  Progress: %d/%d runs (%.1f%%)\n", completed_runs,
                       total_runs, 100.0 * completed_runs / total_runs);
                fflush(stdout);
                last_print_time = now;
              }
            }
          }
        }

        if (thread_result) {
          free_link_percolation_result(thread_result);
          thread_result = NULL;
        }

        if (mode == GRAPH_ER && owned_graph) {
          free_graph(owned_graph);
        }
      }
    }

    free(local_gc);
    free(local_small);
    free(local_lambda);
    free(local_comp);
    free(local_cycles);
    free(local_core);
    free_link_percolation_result(thread_result);
  }

  if (error_flag) {
    if (input_graph)
      free_graph(input_graph);
    ebe_writer_abort(ebe_writer);
    free_online_stats(gc_stats);
    free_online_stats(small_stats);
    free_online_stats(lambda_stats);
    free_online_stats(comp_stats);
    free_online_stats(cycles_stats);
    free_online_stats(core_stats);
    return 1;
  }

  save_averaged_results(averaged_file, graph_nodes,
                        (mode == GRAPH_INPUT && input_graph) ? input_graph->m
                                                             : (int)(avg_degree *
                                                                     graph_nodes /
                                                                     2.0),
                        M_graphs, M_realizations, gc_stats, small_stats,
                        lambda_stats, comp_stats, cycles_stats, core_stats);
  ebe_writer_finish(ebe_writer);
  free_online_stats(gc_stats);
  free_online_stats(small_stats);
  free_online_stats(lambda_stats);
  free_online_stats(comp_stats);
  free_online_stats(cycles_stats);
  free_online_stats(core_stats);
  if (input_graph)
    free_graph(input_graph);
  return 0;
}
