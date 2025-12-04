#include "graph.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Graph *allocate_graph(int n, int m) {
  Graph *g = (Graph *)calloc(1, sizeof(Graph));
  if (!g)
    return NULL;
  g->n = n;
  g->m = m;
  g->offsets = (int *)calloc((size_t)n + 1, sizeof(int));
  g->neighbors = (int *)malloc((size_t)2 * (size_t)m * sizeof(int));
  if (!g->offsets || (m > 0 && !g->neighbors)) {
    free(g->offsets);
    free(g->neighbors);
    free(g);
    return NULL;
  }
  return g;
}

Graph *create_graph_from_edges(int n, const GraphEdge *edges, int m) {
  if (n <= 0 || m < 0)
    return NULL;
  Graph *g = allocate_graph(n, m);
  if (!g)
    return NULL;

  for (int i = 0; i < m; i++) {
    int u = edges[i].u;
    int v = edges[i].v;
    if (u < 0 || u >= n || v < 0 || v >= n || u == v) {
      free_graph(g);
      return NULL;
    }
    g->offsets[u + 1]++;
    g->offsets[v + 1]++;
  }

  for (int i = 0; i < n; i++) {
    g->offsets[i + 1] += g->offsets[i];
  }

  int *cursor = (int *)malloc((size_t)n * sizeof(int));
  if (!cursor) {
    free_graph(g);
    return NULL;
  }
  memcpy(cursor, g->offsets, (size_t)n * sizeof(int));

  for (int i = 0; i < m; i++) {
    int u = edges[i].u;
    int v = edges[i].v;
    g->neighbors[cursor[u]++] = v;
    g->neighbors[cursor[v]++] = u;
  }

  free(cursor);
  return g;
}

void free_graph(Graph *g) {
  if (!g)
    return;
  free(g->offsets);
  free(g->neighbors);
  free(g);
}

static int compare_edgepair(const void *a, const void *b) {
  const GraphEdge *ea = (const GraphEdge *)a;
  const GraphEdge *eb = (const GraphEdge *)b;
  if (ea->u != eb->u)
    return (ea->u < eb->u) ? -1 : 1;
  if (ea->v != eb->v)
    return (ea->v < eb->v) ? -1 : 1;
  return 0;
}

Graph *read_graph_from_file(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (!f) {
    fprintf(stderr, "Error: Cannot open file %s\n", filename);
    return NULL;
  }

  int n, m;
  if (fscanf(f, "%d %d", &n, &m) != 2) {
    fprintf(stderr, "Error: Invalid file format\n");
    fclose(f);
    return NULL;
  }
  if (n <= 0 || m < 0) {
    fprintf(stderr, "Error: Invalid graph size in %s\n", filename);
    fclose(f);
    return NULL;
  }

  GraphEdge *edges = (GraphEdge *)malloc((size_t)m * sizeof(GraphEdge));
  if (!edges) {
    fprintf(stderr, "Error: Unable to allocate edge buffer\n");
    fclose(f);
    return NULL;
  }

  int stored = 0;
  for (int i = 0; i < m; i++) {
    int u, v;
    if (fscanf(f, "%d %d", &u, &v) != 2) {
      fprintf(stderr, "Error: Invalid edge format at line %d\n", i + 2);
      free(edges);
      fclose(f);
      return NULL;
    }
    if (u == v) {
      fprintf(stderr, "Warning: Ignoring self-loop (%d, %d)\n", u, v);
      continue;
    }
    if (v < u) {
      int tmp = u;
      u = v;
      v = tmp;
    }
    edges[stored].u = u;
    edges[stored].v = v;
    stored++;
  }

  // Sort + deduplicate to avoid double edges in input.
  qsort(edges, (size_t)stored, sizeof(GraphEdge), compare_edgepair);
  int unique = 0;
  for (int i = 0; i < stored; i++) {
    if (i == 0 ||
        edges[i].u != edges[i - 1].u || edges[i].v != edges[i - 1].v) {
      edges[unique++] = edges[i];
    }
  }

  Graph *g = create_graph_from_edges(n, edges, unique);
  free(edges);
  fclose(f);
  return g;
}

static inline size_t mix64(uint64_t z) {
  z = (z ^ (z >> 33)) * 0xff51afd7ed558ccdULL;
  z = (z ^ (z >> 33)) * 0xc4ceb9fe1a85ec53ULL;
  return (size_t)(z ^ (z >> 33));
}

Graph *generate_erdos_renyi_rng(int n, double c, rng_state_t *rng) {
  if (n <= 0)
    return NULL;
  long long target = (long long)(c * (double)n / 2.0);
  if (target <= 0)
    target = 0;
  GraphEdge *edges = (target > 0)
                        ? (GraphEdge *)malloc((size_t)target * sizeof(GraphEdge))
                        : NULL;
  if (target > 0 && !edges) {
    fprintf(stderr, "Error: Unable to allocate edge list for ER graph\n");
    return NULL;
  }

  size_t hash_cap = 1;
  while (hash_cap < (size_t)target * 4)
    hash_cap <<= 1;
  int64_t *edge_hash =
      hash_cap ? (int64_t *)malloc(hash_cap * sizeof(int64_t)) : NULL;
  if (hash_cap && !edge_hash) {
    free(edges);
    fprintf(stderr, "Error: Unable to allocate hash table for ER graph\n");
    return NULL;
  }
  if (edge_hash) {
    for (size_t i = 0; i < hash_cap; i++)
      edge_hash[i] = -1;
  }
  size_t hash_mask = hash_cap ? (hash_cap - 1) : 0;

  int edge_count = 0;
  while (edge_count < target) {
    int u = rng_range(rng, n);
    int v = rng_range(rng, n - 1);
    if (v >= u)
      v++;
    int a = (u < v) ? u : v;
    int b = (u < v) ? v : u;
    uint64_t key = ((uint64_t)a << 32) | (uint32_t)b;

    bool duplicate = false;
    if (edge_hash) {
      size_t idx = mix64(key) & hash_mask;
      while (edge_hash[idx] != -1) {
        if ((uint64_t)edge_hash[idx] == key) {
          duplicate = true;
          break;
        }
        idx = (idx + 1) & hash_mask;
      }
      if (duplicate)
        continue;
      edge_hash[idx] = (int64_t)key;
    }

    edges[edge_count].u = u;
    edges[edge_count].v = v;
    edge_count++;
  }

  free(edge_hash);
  Graph *g = create_graph_from_edges(n, edges, edge_count);
  free(edges);
  return g;
}

Graph *generate_erdos_renyi(int n, double c) {
  rng_state_t rng;
  rng_init(&rng, (uint64_t)rand());
  return generate_erdos_renyi_rng(n, c, &rng);
}
