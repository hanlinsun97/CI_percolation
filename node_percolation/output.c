#include "output.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct EBEWriter {
  FILE *file;
  double sum_p_plus;
  double sum_p_plus_jump;
  double sum_q_plus;
  double sum_q_zero;
  double sum_q_minus;
  int total_runs;
  int runs_written;
  char *path;
} EBEWriter;

void save_averaged_results(const char *filename, int n,
                           const char *strategy_name, int M_graphs,
                           int M_realizations, OnlineStats *gc_stats,
                           OnlineStats *small_stats, OnlineStats *lambda_stats,
                           OnlineStats *comp_stats, OnlineStats *cycles_stats) {
  FILE *f = fopen(filename, "w");
  if (!f) {
    fprintf(stderr, "Error: Cannot create output file %s\n", filename);
    return;
  }

  int series_length = NR_SERIES_LENGTH;
  double *gc_std = (double *)malloc(series_length * sizeof(double));
  double *small_std = (double *)malloc(series_length * sizeof(double));
  double *lambda_std = (double *)malloc(series_length * sizeof(double));
  double *comp_std = (double *)malloc(series_length * sizeof(double));
  double *cycles_std = (double *)malloc(series_length * sizeof(double));

  online_stats_get_std(gc_stats, gc_std);
  online_stats_get_std(small_stats, small_std);
  online_stats_get_std(lambda_stats, lambda_std);
  online_stats_get_std(comp_stats, comp_std);
  online_stats_get_std(cycles_stats, cycles_std);

  fprintf(f, "# Node Removal Percolation Results (Averaged)\n");
  fprintf(f, "# N=%d Strategy=%s\n", n, strategy_name);
  fprintf(f, "# M_graphs=%d M_realizations=%d Total_runs=%d\n", M_graphs,
          M_realizations, M_graphs * M_realizations);
  fprintf(f, "# Columns:\n");
  fprintf(f, "#(1) nodes_remaining (2) frac_removed\n");
  fprintf(f, "#(3) gc_mean (4) gc_std\n");
  fprintf(f, "#(5) small_mean (6) small_std\n");
  fprintf(f, "#(7) comp_mean (8) comp_std\n");
  fprintf(f, "#(9) cycles_mean (10) cycles_std\n");
  fprintf(f, "#(11) lambda1_mean (12) lambda1_std\n");

  const int samples = NR_OUTPUT_SAMPLES;
  for (int s = 0; s <= samples; s++) {
    double frac_removed = (double)s / samples;
    int removed = (int)lround(frac_removed * n);
    if (removed > n)
      removed = n;
    int nodes_remaining = n - removed;
    fprintf(f, "%d %.8f ", nodes_remaining, frac_removed);
    int sample_index = s;
    fprintf(f, "%.8f %.8f ", gc_stats->mean[sample_index],
            gc_std[sample_index]);
    fprintf(f, "%.8f %.8f ", small_stats->mean[sample_index],
            small_std[sample_index]);
    fprintf(f, "%.8f %.8f ", comp_stats->mean[sample_index],
            comp_std[sample_index]);
    fprintf(f, "%.8f %.8f ", cycles_stats->mean[sample_index],
            cycles_std[sample_index]);
    fprintf(f, "%.8f %.8f\n", lambda_stats->mean[sample_index],
            lambda_std[sample_index]);
  }

  free(gc_std);
  free(small_std);
  free(lambda_std);
  free(comp_std);
  free(cycles_std);
  fclose(f);
  printf("Averaged results saved to %s\n", filename);
}

static void ebe_writer_write_header(FILE *f, int n, const char *strategy_name,
                                    int M_graphs, int M_realizations,
                                    int total_runs) {
  fprintf(f, "# Event-Based Ensemble (EBE) Results\n");
  fprintf(f, "# N=%d Strategy=%s\n", n, strategy_name);
  fprintf(f, "# M_graphs=%d M_realizations=%d Total_runs=%d\n", M_graphs,
          M_realizations, total_runs);
  fprintf(f, "# Columns:\n");
  fprintf(f, "#(1) run\n");
  fprintf(f, "#(2) p+ (3) dP(p+) (4) P(p+) (5) s(p+) (6) nc(p+) (7) nl(p+)\n");
  fprintf(f, "#(8) q+ (9) P(q+) (10) s(q+) (11) nc(q+) (12) nl(q+)\n");
  fprintf(f, "#(13) q0 (14) P(q0) (15) s(q0) (16) nc(q0) (17) nl(q0)\n");
  fprintf(f, "#(18) q- (19) P(q-) (20) s(q-) (21) nc(q-) (22) nl(q-)\n");
  fprintf(f, "#(23) p2 (24) p1\n");
}

EBEWriter *ebe_writer_begin(const char *filename, int n,
                            const char *strategy_name, int M_graphs,
                            int M_realizations, int total_runs) {
  FILE *f = fopen(filename, "w");
  if (!f) {
    fprintf(stderr, "Error: Cannot create EBE output file %s\n", filename);
    return NULL;
  }
  ebe_writer_write_header(f, n, strategy_name, M_graphs, M_realizations,
                          total_runs);
  EBEWriter *writer = (EBEWriter *)calloc(1, sizeof(EBEWriter));
  if (!writer) {
    fprintf(stderr, "Error: Unable to allocate EBE writer\n");
    fclose(f);
    return NULL;
  }
  writer->file = f;
  writer->total_runs = total_runs;
  size_t len = strlen(filename) + 1;
  writer->path = (char *)malloc(len);
  if (writer->path) {
    memcpy(writer->path, filename, len);
  }
  return writer;
}

void ebe_writer_record(EBEWriter *writer, int run_index,
                       const EnhancedEBEData *d) {
  if (!writer || !writer->file || !d)
    return;
  fprintf(writer->file, "%d ", run_index + 1);
  fprintf(writer->file, "%.8f %.8f %d %.8f %d %d ", d->p_plus, d->p_plus_jump,
          d->P_p_plus, d->s_p_plus, d->nc_p_plus, d->nl_p_plus);
  fprintf(writer->file, "%.8f %d %.8f %d %d ", d->q_plus, d->P_q_plus,
          d->s_q_plus, d->nc_q_plus, d->nl_q_plus);
  fprintf(writer->file, "%.8f %d %.8f %d %d ", d->q_zero, d->P_q_zero,
          d->s_q_zero, d->nc_q_zero, d->nl_q_zero);
  fprintf(writer->file, "%.8f %d %.8f %d %d ", d->q_minus, d->P_q_minus,
          d->s_q_minus, d->nc_q_minus, d->nl_q_minus);
  fprintf(writer->file, "%.8f %.8f\n", d->p1, d->p2);

  writer->sum_p_plus += d->p_plus;
  writer->sum_p_plus_jump += d->p_plus_jump;
  writer->sum_q_plus += d->q_plus;
  writer->sum_q_zero += d->q_zero;
  writer->sum_q_minus += d->q_minus;
  writer->runs_written++;
}

static void ebe_writer_close(EBEWriter *writer) {
  if (!writer)
    return;
  if (writer->file) {
    fclose(writer->file);
    writer->file = NULL;
  }
  free(writer->path);
  free(writer);
}

void ebe_writer_finish(EBEWriter *writer) {
  if (!writer)
    return;
  if (writer->file) {
    fprintf(writer->file, "# STATISTICS:\n");
    if (writer->runs_written > 0) {
      double inv = 1.0 / (double)writer->runs_written;
      fprintf(writer->file, "# Mean p+ = %.8f\n",
              writer->sum_p_plus * inv);
      fprintf(writer->file, "# Mean dP(p+) = %.8f\n",
              writer->sum_p_plus_jump * inv);
      fprintf(writer->file, "# Mean q+ = %.8f\n",
              writer->sum_q_plus * inv);
      fprintf(writer->file, "# Mean q0 = %.8f\n",
              writer->sum_q_zero * inv);
      fprintf(writer->file, "# Mean q- = %.8f\n",
              writer->sum_q_minus * inv);
    }
  }
  ebe_writer_close(writer);
}

void ebe_writer_abort(EBEWriter *writer) {
  if (!writer)
    return;
  if (writer->file) {
    fclose(writer->file);
    if (writer->path) {
      remove(writer->path);
    }
  }
  ebe_writer_close(writer);
}

// Convenience wrapper if buffering in memory is used.
void save_enhanced_ebe_results(const char *filename, int n,
                               const char *strategy_name, int M_graphs,
                               int M_realizations,
                               const EnhancedEBEData *ebe_data,
                               int total_runs) {
  EBEWriter *writer =
      ebe_writer_begin(filename, n, strategy_name, M_graphs, M_realizations,
                       total_runs);
  if (!writer)
    return;
  for (int i = 0; i < total_runs; i++) {
    ebe_writer_record(writer, i, &ebe_data[i]);
  }
  ebe_writer_finish(writer);
  printf("EBE results saved to %s\n", filename);
}
