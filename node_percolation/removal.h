#ifndef NODE_PERCOLATION_REMOVAL_H
#define NODE_PERCOLATION_REMOVAL_H

#include <stdbool.h>
#include "graph.h"
#include "output.h"
#include "ranking.h"

typedef struct {
  RankingConfig ranking;
} NodeRemovalConfig;

typedef struct NodePercolationResult {
  int n;
  int series_length;
  double *giant_size;
  double *avg_small_size;
  double *lambda1;
  int *num_components;
  int *num_cycles;
  int *removal_order;
  int removal_length;
} NodePercolationResult;

NodePercolationResult *create_node_percolation_result(int n);
void free_node_percolation_result(NodePercolationResult *result);
void reset_node_percolation_result(NodePercolationResult *result);

int run_node_percolation(Graph *g, const NodeRemovalConfig *config,
                         NodePercolationResult *result,
                         PseudoCriticalPoints *ecp,
                         uint64_t rng_seed);

#endif // NODE_PERCOLATION_REMOVAL_H
