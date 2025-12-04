#include "removal.h"
#include "newman_ziff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int u;
  int v;
  unsigned char active;
  unsigned char in_core;
} EdgeState;

static int *build_incident_offsets(const Graph *g) {
  int *offsets = (int *)calloc((size_t)g->n + 1, sizeof(int));
  if (!offsets)
    return NULL;
  for (int node = 0; node < g->n; node++) {
    offsets[node + 1] = graph_degree(g, node);
  }
  for (int i = 0; i < g->n; i++) {
    offsets[i + 1] += offsets[i];
  }
  return offsets;
}

static EdgeState *build_edge_list(const Graph *g, int *incident_edges,
                                  const int *incident_offsets) {
  EdgeState *edges = (EdgeState *)calloc((size_t)g->m, sizeof(EdgeState));
  if (!edges)
    return NULL;

  int *cursor = (int *)malloc((size_t)g->n * sizeof(int));
  if (!cursor) {
    free(edges);
    return NULL;
  }
  memcpy(cursor, incident_offsets, (size_t)g->n * sizeof(int));

  int edge_id = 0;
  for (int u = 0; u < g->n; u++) {
    const int *neighbors = graph_neighbors(g, u);
    int deg = graph_degree(g, u);
    for (int j = 0; j < deg; j++) {
      int v = neighbors[j];
      if (u < v) {
        if (edge_id >= g->m) {
          free(edges);
          free(cursor);
          return NULL;
        }
        edges[edge_id].u = u;
        edges[edge_id].v = v;
        edges[edge_id].active = 1;
        edges[edge_id].in_core = 0;
        incident_edges[cursor[u]++] = edge_id;
        incident_edges[cursor[v]++] = edge_id;
        edge_id++;
      }
    }
  }

  free(cursor);
  if (edge_id != g->m) {
    free(edges);
    return NULL;
  }
  return edges;
}

LinkPercolationResult *create_link_percolation_result(int n, int m) {
  LinkPercolationResult *res =
      (LinkPercolationResult *)calloc(1, sizeof(LinkPercolationResult));
  if (!res)
    return NULL;
  const int len = NR_SERIES_LENGTH;
  res->n = n;
  res->m = m;
  res->series_length = len;
  res->giant_size = (double *)calloc(len, sizeof(double));
  res->avg_small_size = (double *)calloc(len, sizeof(double));
  res->core_size = (double *)calloc(len, sizeof(double));
  res->lambda1 = (double *)calloc(len, sizeof(double));
  res->num_components = (int *)calloc(len, sizeof(int));
  res->num_cycles = (int *)calloc(len, sizeof(int));
  res->removal_order = (int *)malloc((size_t)m * sizeof(int));
  if (!res->giant_size || !res->avg_small_size || !res->core_size ||
      !res->lambda1 || !res->num_components || !res->num_cycles ||
      !res->removal_order) {
    free_link_percolation_result(res);
    return NULL;
  }
  for (int i = 0; i < len; i++) {
    res->lambda1[i] = -1.0;
  }
  res->core_zero_p = -1.0;
  return res;
}

void free_link_percolation_result(LinkPercolationResult *res) {
  if (!res)
    return;
  free(res->giant_size);
  free(res->avg_small_size);
  free(res->core_size);
  free(res->lambda1);
  free(res->num_components);
  free(res->num_cycles);
  free(res->removal_order);
  free(res);
}

void reset_link_percolation_result(LinkPercolationResult *res) {
  if (!res)
    return;
  const size_t len = (size_t)res->series_length;
  memset(res->giant_size, 0, len * sizeof(double));
  memset(res->avg_small_size, 0, len * sizeof(double));
  memset(res->core_size, 0, len * sizeof(double));
  for (int i = 0; i < res->series_length; i++) {
    res->lambda1[i] = -1.0;
  }
  memset(res->num_components, 0, len * sizeof(int));
  memset(res->num_cycles, 0, len * sizeof(int));
  res->removal_length = 0;
  res->core_p_star = 0.0;
  res->core_drop = 0;
  res->core_size_at_star = 0;
  res->core_zero_p = -1.0;
}

static inline int other_endpoint(const EdgeState *edges, int edge_id, int u) {
  return (edges[edge_id].u == u) ? edges[edge_id].v : edges[edge_id].u;
}

static void core_edge_remove(int edge_id, EdgeState *edges, int *core_edges,
                             int *core_pos, int *core_edge_count) {
  if (!edges[edge_id].in_core)
    return;
  edges[edge_id].in_core = 0;
  int last_idx = *core_edge_count - 1;
  int pos = core_pos[edge_id];
  if (pos < 0)
    return;
  if (pos != last_idx) {
    int swapped = core_edges[last_idx];
    core_edges[pos] = swapped;
    core_pos[swapped] = pos;
  }
  (*core_edge_count)--;
  core_pos[edge_id] = -1;
}

static void enqueue_if_needed(int node, int *queue, int *qtail,
                              unsigned char *in_queue, int *core_degree) {
  if (core_degree[node] < 2 && !in_queue[node]) {
    in_queue[node] = 1;
    queue[(*qtail)++] = node;
  }
}

static void remove_core_node(int node, EdgeState *edges, const int *inc_offsets,
                             const int *inc_edges, unsigned char *in_core,
                             int *core_degree, int *core_nodes,
                             int *core_edges_arr, int *core_pos,
                             int *core_edge_count, int *queue, int *qtail,
                             unsigned char *in_queue) {
  if (!in_core[node])
    return;
  in_core[node] = 0;
  (*core_nodes)--;
  for (int idx = inc_offsets[node]; idx < inc_offsets[node + 1]; idx++) {
    int e = inc_edges[idx];
    if (!edges[e].active)
      continue;
    if (edges[e].in_core) {
      core_edge_remove(e, edges, core_edges_arr, core_pos, core_edge_count);
    }
    int nbr = other_endpoint(edges, e, node);
    if (in_core[nbr]) {
      core_degree[nbr]--;
      enqueue_if_needed(nbr, queue, qtail, in_queue, core_degree);
    }
  }
  core_degree[node] = 0;
}

static void sample_core_series(LinkPercolationResult *res, int removed_edges,
                               int core_nodes) {
  const int samples = res->series_length - 1;
  if (samples <= 0 || res->m <= 0)
    return;
  int idx =
      (int)lround((double)removed_edges * (double)samples / (double)res->m);
  if (idx < 0)
    idx = 0;
  if (idx > samples)
    idx = samples;
  // Only fill if not already set (simple sparse sampling).
  if (res->core_size[idx] == 0.0 || removed_edges == 0) {
    res->core_size[idx] = (double)core_nodes;
  }
}

int run_link_percolation(Graph *g, uint64_t rng_seed,
                         LinkPercolationResult *result,
                         PseudoCriticalPoints *gc_ecp) {
  const int n = g->n;
  const int m = g->m;
  if (m <= 0)
    return -1;

  int *incident_offsets = build_incident_offsets(g);
  int *incident_edges = incident_offsets
                            ? (int *)malloc((size_t)incident_offsets[n] *
                                            sizeof(int))
                            : NULL;
  if (!incident_offsets || !incident_edges) {
    free(incident_offsets);
    free(incident_edges);
    return -1;
  }

  EdgeState *edges = build_edge_list(g, incident_edges, incident_offsets);
  if (!edges) {
    free(incident_offsets);
    free(incident_edges);
    return -1;
  }

  unsigned char *in_core =
      (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
  unsigned char *in_queue =
      (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
  int *core_degree = (int *)malloc((size_t)n * sizeof(int));
  int *queue = (int *)malloc((size_t)n * sizeof(int));
  int *core_edges_arr = (int *)malloc((size_t)m * sizeof(int));
  int *core_pos = (int *)malloc((size_t)m * sizeof(int));
  if (!in_core || !in_queue || !core_degree || !queue || !core_edges_arr ||
      !core_pos) {
    free(incident_offsets);
    free(incident_edges);
    free(edges);
    free(in_core);
    free(in_queue);
    free(core_degree);
    free(queue);
    free(core_edges_arr);
    free(core_pos);
    return -1;
  }

  for (int i = 0; i < n; i++) {
    in_core[i] = 1;
    core_degree[i] = graph_degree(g, i);
  }
  for (int i = 0; i < m; i++) {
    core_pos[i] = -1;
  }

  int qhead = 0, qtail = 0;
  for (int node = 0; node < n; node++) {
    enqueue_if_needed(node, queue, &qtail, in_queue, core_degree);
  }

  int core_nodes = n;
  int core_edge_count = 0;
  // Prune initial 2-core.
  while (qhead < qtail) {
    int node = queue[qhead++];
    in_queue[node] = 0;
    remove_core_node(node, edges, incident_offsets, incident_edges, in_core,
                     core_degree, &core_nodes, core_edges_arr, core_pos,
                     &core_edge_count, queue, &qtail, in_queue);
  }

  for (int e = 0; e < m; e++) {
    if (in_core[edges[e].u] && in_core[edges[e].v]) {
      edges[e].in_core = 1;
      core_pos[e] = core_edge_count;
      core_edges_arr[core_edge_count++] = e;
    }
  }

  rng_state_t rng_local;
  rng_init(&rng_local, rng_seed);

  reset_link_percolation_result(result);

  sample_core_series(result, 0, core_nodes);
  double prev_core_size = (double)core_nodes;
  int removed_edges = 0;

  while (core_edge_count > 0) {
    int idx = rng_range(&rng_local, (uint32_t)core_edge_count);
    int edge_id = core_edges_arr[idx];
    if (!edges[edge_id].active || !edges[edge_id].in_core) {
      core_edge_remove(edge_id, edges, core_edges_arr, core_pos,
                       &core_edge_count);
      continue;
    }

    edges[edge_id].active = 0;
    core_edge_remove(edge_id, edges, core_edges_arr, core_pos,
                     &core_edge_count);
    int u = edges[edge_id].u;
    int v = edges[edge_id].v;
    if (in_core[u]) {
      core_degree[u]--;
      enqueue_if_needed(u, queue, &qtail, in_queue, core_degree);
    }
    if (in_core[v]) {
      core_degree[v]--;
      enqueue_if_needed(v, queue, &qtail, in_queue, core_degree);
    }

    while (qhead < qtail) {
      int node = queue[qhead++];
      in_queue[node] = 0;
      remove_core_node(node, edges, incident_offsets, incident_edges, in_core,
                       core_degree, &core_nodes, core_edges_arr, core_pos,
                       &core_edge_count, queue, &qtail, in_queue);
    }

    result->removal_order[removed_edges++] = edge_id;
    sample_core_series(result, removed_edges, core_nodes);

    double core_size_now = (double)core_nodes;
    double drop = prev_core_size - core_size_now;
    if (drop > (double)result->core_drop) {
      result->core_drop = (int)drop;
      result->core_p_star = (double)removed_edges / (double)m;
      result->core_size_at_star = core_nodes;
    }
    prev_core_size = core_size_now;
    if (core_nodes == 0 && result->core_zero_p < 0.0) {
      result->core_zero_p = (double)removed_edges / (double)m;
    }
  }

  if (result->core_zero_p < 0.0) {
    result->core_zero_p = (double)removed_edges / (double)m;
  }

  // Remove remaining (tree) edges uniformly at random.
  int remaining = 0;
  for (int e = 0; e < m; e++) {
    if (edges[e].active) {
      core_edges_arr[remaining++] = e;
    }
  }
  for (int i = remaining - 1; i >= 0; i--) {
    int r = (int)rng_range(&rng_local, (uint32_t)(i + 1));
    int chosen = core_edges_arr[r];
    core_edges_arr[r] = core_edges_arr[i];
    core_edges_arr[i] = chosen;
    edges[chosen].active = 0;
    result->removal_order[removed_edges++] = chosen;
    sample_core_series(result, removed_edges, 0);
  }

  result->removal_length = removed_edges;
  if (removed_edges != m) {
    free(incident_offsets);
    free(incident_edges);
    free(edges);
    free(in_core);
    free(in_queue);
    free(core_degree);
    free(queue);
    free(core_edges_arr);
    free(core_pos);
    return -1;
  }

  // Forward-fill any unsampled core sizes to avoid artificial zeros.
  double last_core = result->core_size[0];
  for (int i = 0; i < result->series_length; i++) {
    if (result->core_size[i] <= 0.0 && i > 0) {
      result->core_size[i] = last_core;
    } else {
      last_core = result->core_size[i];
    }
  }

  reconstruct_with_newman_ziff(g, result->removal_order, m, result, gc_ecp);

  free(incident_offsets);
  free(incident_edges);
  free(edges);
  free(in_core);
  free(in_queue);
  free(core_degree);
  free(queue);
  free(core_edges_arr);
  free(core_pos);
  return 0;
}
