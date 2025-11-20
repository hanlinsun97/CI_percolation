#ifndef NODE_PERCOLATION_RANKING_H
#define NODE_PERCOLATION_RANKING_H

#include <stdbool.h>
#include "graph.h"
#include "rng.h"

typedef struct {
  const char *strategy;
  int ci_radius;
} RankingConfig;

typedef struct NodeRanker NodeRanker;

NodeRanker *node_ranker_create(Graph *g, const RankingConfig *config,
                               bool *removed, int *active_degree,
                               rng_state_t *rng);
int node_ranker_next(NodeRanker *ranker);
void node_ranker_handle_removal(NodeRanker *ranker, int node);
void node_ranker_free(NodeRanker *ranker);

#endif // NODE_PERCOLATION_RANKING_H
