#ifndef LINK_PERCOLATION_REMOVAL_H
#define LINK_PERCOLATION_REMOVAL_H

#include "graph.h"
#include "output.h"
#include "rng.h"
#include <stdint.h>

typedef struct {
  int n;
  int m;
  int series_length;
  double *giant_size;
  double *avg_small_size;
  double *core_size;
  double *lambda1;
  int *num_components;
  int *num_cycles;
  int *removal_order;
  int removal_length;
  double core_p_star;
  int core_drop;
  int core_size_at_star;
  double core_zero_p;
} LinkPercolationResult;

LinkPercolationResult *create_link_percolation_result(int n, int m,
                                                      int series_length);
void free_link_percolation_result(LinkPercolationResult *result);
void reset_link_percolation_result(LinkPercolationResult *result);

int run_link_percolation(Graph *g, uint64_t rng_seed,
                         LinkPercolationResult *result,
                         PseudoCriticalPoints *gc_ecp,
                         int use_static_core, int use_random_all,
                         int use_biased_selection,
                         int use_random_biased_original,
                         int use_random_biased_adaptive);

#endif // LINK_PERCOLATION_REMOVAL_H
