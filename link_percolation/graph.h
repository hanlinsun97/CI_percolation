#ifndef GRAPH_H
#define GRAPH_H

#include <stdbool.h>
#include "rng.h"

typedef struct {
  int u;
  int v;
} GraphEdge;

typedef struct {
  int n;
  int m;
  int *offsets;   // length n+1
  int *neighbors; // length 2m
} Graph;

// Graph creation and manipulation
Graph *create_graph_from_edges(int n, const GraphEdge *edges, int m);
void free_graph(Graph *g);

static inline int graph_degree(const Graph *g, int node) {
  return g->offsets[node + 1] - g->offsets[node];
}

static inline const int *graph_neighbors(const Graph *g, int node) {
  return &g->neighbors[g->offsets[node]];
}

// Graph I/O
Graph *read_graph_from_file(const char *filename);

// Graph generation with thread-safe RNG
Graph *generate_erdos_renyi_rng(int n, double c, rng_state_t *rng);

// Backward compatibility: uses global rand()
Graph *generate_erdos_renyi(int n, double c);

#endif // GRAPH_H
