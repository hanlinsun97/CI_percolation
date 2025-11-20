#include "giant_component.h"
#include "graph.h"
#include <stdlib.h>

typedef struct {
  int *queue;
  int *visited;
  int capacity;
} BFSState;

static BFSState *bfs_create(int n) {
  BFSState *b = (BFSState *)malloc(sizeof(BFSState));
  b->queue = (int *)malloc(n * sizeof(int));
  b->visited = (int *)calloc(n, sizeof(int));
  b->capacity = n;
  if (!b->queue || !b->visited) {
    free(b->queue);
    free(b->visited);
    free(b);
    return NULL;
  }
  return b;
}

static void bfs_free(BFSState *b) {
  if (!b)
    return;
  free(b->queue);
  free(b->visited);
  free(b);
}

int compute_giant_component(Graph *g) {
  int n = g->n;
  BFSState *b = bfs_create(n);
  if (!b)
    return -1;

  int gc = 0;
  int qhead = 0, qtail = 0;
  for (int start = 0; start < n; start++) {
    if (b->visited[start])
      continue;
    b->visited[start] = 1;
    b->queue[qtail++] = start;
    int size = 0;
    while (qhead < qtail) {
      int u = b->queue[qhead++];
      size++;
      const int *neighbors = graph_neighbors(g, u);
      int deg = graph_degree(g, u);
      for (int i = 0; i < deg; i++) {
        int v = neighbors[i];
        if (!b->visited[v]) {
          b->visited[v] = 1;
          b->queue[qtail++] = v;
        }
      }
    }
    if (size > gc)
      gc = size;
  }

  bfs_free(b);
  return gc;
}
