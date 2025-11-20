#include "removal.h"
#include "newman_ziff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

NodePercolationResult *create_node_percolation_result(int n) {
  NodePercolationResult *result =
      (NodePercolationResult *)malloc(sizeof(NodePercolationResult));
  if (!result)
    return NULL;
  int sample_length = NR_SERIES_LENGTH;
  result->n = n;
  result->series_length = sample_length;
  result->giant_size = (double *)calloc(sample_length, sizeof(double));
  result->avg_small_size = (double *)calloc(sample_length, sizeof(double));
  result->lambda1 = (double *)calloc(sample_length, sizeof(double));
  result->num_components = (int *)calloc(sample_length, sizeof(int));
  result->num_cycles = (int *)calloc(sample_length, sizeof(int));
  result->removal_order = (int *)malloc(n * sizeof(int));
  result->removal_length = 0;
  if (!result->giant_size || !result->avg_small_size || !result->lambda1 ||
      !result->num_components || !result->num_cycles ||
      !result->removal_order) {
    free(result->giant_size);
    free(result->avg_small_size);
    free(result->lambda1);
    free(result->num_components);
    free(result->num_cycles);
    free(result->removal_order);
    free(result);
    return NULL;
  }
  for (int i = 0; i < sample_length; i++) {
    result->lambda1[i] = -1.0;
  }
  return result;
}

void free_node_percolation_result(NodePercolationResult *result) {
  if (!result)
    return;
  free(result->giant_size);
  free(result->avg_small_size);
  free(result->lambda1);
  free(result->num_components);
  free(result->num_cycles);
  free(result->removal_order);
  free(result);
}

void reset_node_percolation_result(NodePercolationResult *result) {
  if (!result)
    return;
  size_t len = (size_t)result->series_length;
  memset(result->giant_size, 0, len * sizeof(double));
  memset(result->avg_small_size, 0, len * sizeof(double));
  for (int i = 0; i < result->series_length; i++) {
    result->lambda1[i] = -1.0;
  }
  memset(result->num_components, 0, len * sizeof(int));
  memset(result->num_cycles, 0, len * sizeof(int));
  result->removal_length = 0;
}

int run_node_percolation(Graph *g, const NodeRemovalConfig *config,
                         NodePercolationResult *result,
                         PseudoCriticalPoints *ecp,
                         uint64_t rng_seed) {
  int n = g->n;
  bool *removed = (bool *)calloc(n, sizeof(bool));
  int *active_degree = (int *)malloc(n * sizeof(int));
  if (!removed || !active_degree) {
    free(removed);
    free(active_degree);
    return -1;
  }
  for (int i = 0; i < n; i++) {
    active_degree[i] = graph_degree(g, i);
  }

  rng_state_t rng_local;
  rng_init(&rng_local, rng_seed);

  NodeRanker *ranker =
      node_ranker_create(g, &config->ranking, removed, active_degree,
                         &rng_local);
  if (!ranker) {
    free(removed);
    free(active_degree);
    return -1;
  }

  int removal_pos = 0;
  while (removal_pos < n) {
    int node = node_ranker_next(ranker);
    if (node < 0) {
      break;
    }
    if (removed[node]) {
      continue;
    }

    removed[node] = true;
    result->removal_order[removal_pos++] = node;
    const int *neighbors = graph_neighbors(g, node);
    int degree = graph_degree(g, node);
    for (int i = 0; i < degree; i++) {
      int neighbor = neighbors[i];
      if (!removed[neighbor] && active_degree[neighbor] > 0) {
        active_degree[neighbor]--;
      }
    }
    active_degree[node] = 0;

    node_ranker_handle_removal(ranker, node);
  }

  for (int node = 0; node < n && removal_pos < n; node++) {
    if (!removed[node]) {
      removed[node] = true;
      result->removal_order[removal_pos++] = node;
    }
  }

  result->removal_length = removal_pos;

  if (removal_pos != n) {
    node_ranker_free(ranker);
    free(removed);
    free(active_degree);
    return -1;
  }

  reconstruct_with_newman_ziff(g, result->removal_order, n, result, ecp);

  node_ranker_free(ranker);
  free(removed);
  free(active_degree);
  return 0;
}
