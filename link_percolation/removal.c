#include "removal.h"
#include "newman_ziff.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int n;
  int *parent;
  int *size;
  int num_components;
  int largest_component;
  int active_edges;
  long long sum_s_squared;
} UF;

typedef struct {
  double giant;
  double avg_small;
  int num_components;
  int num_cycles;
} SimpleSnap;

static UF *uf_create(int n) {
  UF *uf = (UF *)calloc(1, sizeof(UF));
  if (!uf)
    return NULL;
  uf->parent = (int *)malloc((size_t)n * sizeof(int));
  uf->size = (int *)malloc((size_t)n * sizeof(int));
  if (!uf->parent || !uf->size) {
    free(uf->parent);
    free(uf->size);
    free(uf);
    return NULL;
  }
  uf->n = n;
  uf->num_components = n;
  uf->largest_component = 1;
  uf->active_edges = 0;
  uf->sum_s_squared = (long long)n; // n * 1^2
  for (int i = 0; i < n; i++) {
    uf->parent[i] = i;
    uf->size[i] = 1;
  }
  return uf;
}

static void uf_free(UF *uf) {
  if (!uf)
    return;
  free(uf->parent);
  free(uf->size);
  free(uf);
}

static int uf_find(UF *uf, int x) {
  while (uf->parent[x] != x) {
    uf->parent[x] = uf->parent[uf->parent[x]];
    x = uf->parent[x];
  }
  return x;
}

static void uf_union(UF *uf, int a, int b) {
  int ra = uf_find(uf, a);
  int rb = uf_find(uf, b);
  if (ra == rb)
    return;
  if (uf->size[ra] < uf->size[rb]) {
    int tmp = ra;
    ra = rb;
    rb = tmp;
  }
  uf->parent[rb] = ra;
  long long sa = uf->size[ra];
  long long sb = uf->size[rb];
  uf->sum_s_squared -= sa * sa + sb * sb;
  uf->size[ra] += uf->size[rb];
  long long snew = uf->size[ra];
  uf->sum_s_squared += snew * snew;
  uf->num_components--;
  if (uf->size[ra] > uf->largest_component) {
    uf->largest_component = uf->size[ra];
  }
}

static double uf_avg_small(const UF *uf) {
  if (uf->n <= 0)
    return 0.0;
  int giant = uf->largest_component;
  int nodes_in_small = uf->n - giant;
  if (nodes_in_small <= 0)
    return 0.0;
  long long small_sq = uf->sum_s_squared - (long long)giant * (long long)giant;
  if (small_sq <= 0)
    return 0.0;
  return (double)small_sq / (double)nodes_in_small;
}

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

LinkPercolationResult *create_link_percolation_result(int n, int m,
                                                      int series_length) {
  if (n <= 0 || m < 0 || series_length < 2)
    return NULL;
  LinkPercolationResult *res =
      (LinkPercolationResult *)calloc(1, sizeof(LinkPercolationResult));
  if (!res)
    return NULL;
  const int len = series_length;
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
                             int *core_pos, int *core_edge_count,
                             int maintain_list) {
  if (!edges[edge_id].in_core)
    return;
  edges[edge_id].in_core = 0;
  if (!maintain_list)
    return;
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

static inline int heap_better(int score_a, int edge_a, int score_b, int edge_b) {
  if (score_a != score_b)
    return score_a > score_b;
  // Deterministic tie-break to avoid dependence on enumeration order.
  return edge_a < edge_b;
}

static void heap_swap(int *heap_edges, int *heap_scores, int i, int j) {
  int te = heap_edges[i];
  heap_edges[i] = heap_edges[j];
  heap_edges[j] = te;
  int ts = heap_scores[i];
  heap_scores[i] = heap_scores[j];
  heap_scores[j] = ts;
}

static void heap_sift_down(int *heap_edges, int *heap_scores, int size,
                           int idx) {
  while (1) {
    int left = 2 * idx + 1;
    int right = left + 1;
    int best = idx;
    if (left < size &&
        heap_better(heap_scores[left], heap_edges[left], heap_scores[best],
                    heap_edges[best])) {
      best = left;
    }
    if (right < size &&
        heap_better(heap_scores[right], heap_edges[right], heap_scores[best],
                    heap_edges[best])) {
      best = right;
    }
    if (best == idx)
      break;
    heap_swap(heap_edges, heap_scores, idx, best);
    idx = best;
  }
}

static void heap_sift_up(int *heap_edges, int *heap_scores, int idx) {
  while (idx > 0) {
    int parent = (idx - 1) / 2;
    if (heap_better(heap_scores[parent], heap_edges[parent], heap_scores[idx],
                    heap_edges[idx])) {
      break;
    }
    heap_swap(heap_edges, heap_scores, parent, idx);
    idx = parent;
  }
}

static void heap_build(int *heap_edges, int *heap_scores, int size) {
  for (int i = size / 2 - 1; i >= 0; i--) {
    heap_sift_down(heap_edges, heap_scores, size, i);
  }
}

static int heap_pop(int *heap_edges, int *heap_scores, int *size, int *edge_id,
                    int *score) {
  if (*size <= 0)
    return 0;
  *edge_id = heap_edges[0];
  *score = heap_scores[0];
  (*size)--;
  if (*size > 0) {
    heap_edges[0] = heap_edges[*size];
    heap_scores[0] = heap_scores[*size];
    heap_sift_down(heap_edges, heap_scores, *size, 0);
  }
  return 1;
}

static void heap_push(int *heap_edges, int *heap_scores, int *size, int edge_id,
                      int score) {
  int idx = (*size)++;
  heap_edges[idx] = edge_id;
  heap_scores[idx] = score;
  heap_sift_up(heap_edges, heap_scores, idx);
}

typedef struct {
  int edge_id;
  int score;
} EdgeScore;

static int compare_edge_score_desc(const void *a, const void *b) {
  const EdgeScore *ea = (const EdgeScore *)a;
  const EdgeScore *eb = (const EdgeScore *)b;
  if (ea->score != eb->score)
    return (ea->score < eb->score) ? 1 : -1; // descending
  if (ea->edge_id != eb->edge_id)
    return (ea->edge_id > eb->edge_id) ? 1 : -1; // ascending
  return 0;
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
                             int *core_edge_count, int maintain_core_list,
                             int *queue, int *qtail,
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
      core_edge_remove(e, edges, core_edges_arr, core_pos, core_edge_count,
                       maintain_core_list);
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
                         PseudoCriticalPoints *gc_ecp,
                         int use_static_core, int use_random_all,
                         int use_biased_selection,
                         int use_random_biased_original,
                         int use_random_biased_adaptive) {
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
  int maintain_core_list = 1;
  // Prune initial 2-core.
  long long safety_core_prune = 0;
  long long max_core_prune = 10LL * n + 100000;
  while (qhead < qtail) {
    if (++safety_core_prune > max_core_prune) {
      fprintf(stderr, "Core prune safety hit (N=%d, m=%d, q=%d)\n", n, m,
              qtail - qhead);
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
    int node = queue[qhead++];
    in_queue[node] = 0;
    remove_core_node(node, edges, incident_offsets, incident_edges, in_core,
                     core_degree, &core_nodes, core_edges_arr, core_pos,
                     &core_edge_count, 1, queue, &qtail, in_queue);
  }

  for (int e = 0; e < m; e++) {
    if (in_core[edges[e].u] && in_core[edges[e].v]) {
      edges[e].in_core = 1;
      core_pos[e] = core_edge_count;
      core_edges_arr[core_edge_count++] = e;
    }
  }

  // Selection pools:
  // - adaptive: use current core_edges_arr/core_edge_count
  // - static: frozen copy of initial core edges
  // - random_all: all edges at start
  int *static_core_edges = NULL;
  EdgeScore *static_scored = NULL;
  int static_remaining = 0;
  int static_total = 0;
  int static_cursor = 0;
  int static_mode = (use_static_core && core_edge_count > 0) ? 1 : 0;
  int biased_mode = use_biased_selection ? 1 : 0;
  int adaptive_biased_mode = (biased_mode && !static_mode) ? 1 : 0;
  int static_biased_mode = (biased_mode && static_mode) ? 1 : 0;
  if (static_mode) {
    static_remaining = core_edge_count;
    static_total = core_edge_count;
    if (static_biased_mode) {
      static_scored =
          (EdgeScore *)malloc((size_t)static_total * sizeof(EdgeScore));
      if (!static_scored) {
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
      for (int i = 0; i < static_total; i++) {
        int eid = core_edges_arr[i];
        int u = edges[eid].u;
        int v = edges[eid].v;
        static_scored[i].edge_id = eid;
        static_scored[i].score = core_degree[u] + core_degree[v];
      }
      qsort(static_scored, (size_t)static_total, sizeof(EdgeScore),
            compare_edge_score_desc);
    } else {
      static_core_edges = (int *)malloc((size_t)static_remaining * sizeof(int));
      if (!static_core_edges) {
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
      memcpy(static_core_edges, core_edges_arr,
             (size_t)static_remaining * sizeof(int));
    }
  }

  int *all_edges_pool = NULL;
  int all_remaining = 0;
  int random_mode = use_random_all ? 1 : 0;
  int random_biased_original_mode = use_random_biased_original ? 1 : 0;
  int random_biased_adaptive_mode = use_random_biased_adaptive ? 1 : 0;
  int random_biased_mode =
      (random_biased_original_mode || random_biased_adaptive_mode) ? 1 : 0;
  if (random_mode) {
    all_remaining = m;
    all_edges_pool = (int *)malloc((size_t)all_remaining * sizeof(int));
    if (!all_edges_pool) {
      free(incident_offsets);
      free(incident_edges);
      free(edges);
      free(in_core);
      free(in_queue);
      free(core_degree);
      free(queue);
      free(core_edges_arr);
      free(core_pos);
      free(static_core_edges);
      free(static_scored);
      return -1;
    }
    for (int i = 0; i < m; i++) {
      all_edges_pool[i] = i;
    }
  }

  int *active_degree = NULL;
  if (random_biased_adaptive_mode) {
    active_degree = (int *)malloc((size_t)n * sizeof(int));
    if (!active_degree) {
      free(incident_offsets);
      free(incident_edges);
      free(edges);
      free(in_core);
      free(in_queue);
      free(core_degree);
      free(queue);
      free(core_edges_arr);
      free(core_pos);
      free(static_core_edges);
      free(static_scored);
      free(all_edges_pool);
      return -1;
    }
    for (int i = 0; i < n; i++) {
      active_degree[i] = graph_degree(g, i);
    }
  }

  EdgeScore *all_scored = NULL;
  int all_scored_cursor = 0;
  if (random_biased_original_mode) {
    all_scored = (EdgeScore *)malloc((size_t)m * sizeof(EdgeScore));
    if (!all_scored) {
      free(incident_offsets);
      free(incident_edges);
      free(edges);
      free(in_core);
      free(in_queue);
      free(core_degree);
      free(queue);
      free(core_edges_arr);
      free(core_pos);
      free(static_core_edges);
      free(static_scored);
      free(all_edges_pool);
      free(active_degree);
      return -1;
    }
    for (int e = 0; e < m; e++) {
      int u = edges[e].u;
      int v = edges[e].v;
      all_scored[e].edge_id = e;
      all_scored[e].score = graph_degree(g, u) + graph_degree(g, v);
    }
    qsort(all_scored, (size_t)m, sizeof(EdgeScore), compare_edge_score_desc);
    all_scored_cursor = 0;
  }

  rng_state_t rng_local;
  rng_init(&rng_local, rng_seed);

  if (random_mode) {
    // Shuffle once to ensure random order independent of enumeration.
    for (int i = m - 1; i > 0; i--) {
      int r = (int)rng_range(&rng_local, (uint32_t)(i + 1));
      int tmp = all_edges_pool[i];
      all_edges_pool[i] = all_edges_pool[r];
      all_edges_pool[r] = tmp;
    }
  }

  reset_link_percolation_result(result);

  sample_core_series(result, 0, core_nodes);
  double prev_core_size = (double)core_nodes;
  int removed_edges = 0;
  int heap_size = 0;
  int *heap_edges = core_edges_arr;
  int *heap_scores = core_pos;
  if (adaptive_biased_mode && core_edge_count > 0) {
    maintain_core_list = 0;
    heap_size = core_edge_count;
    for (int i = 0; i < heap_size; i++) {
      int eid = core_edges_arr[i];
      int u = edges[eid].u;
      int v = edges[eid].v;
      heap_scores[i] = core_degree[u] + core_degree[v];
    }
    heap_build(heap_edges, heap_scores, heap_size);
  }
  int all_heap_size = 0;
  if (random_biased_adaptive_mode) {
    maintain_core_list = 0;
    all_heap_size = m;
    for (int i = 0; i < m; i++) {
      heap_edges[i] = i;
      int u = edges[i].u;
      int v = edges[i].v;
      heap_scores[i] = active_degree[u] + active_degree[v];
    }
    heap_build(heap_edges, heap_scores, all_heap_size);
  }

  long long safety_core_loop = 0;
  long long max_core_loop = 10LL * m + 100000;
  while (1) {
    int active_pool = 0;
    if (random_biased_original_mode) {
      active_pool = m - all_scored_cursor;
    } else if (random_biased_adaptive_mode) {
      active_pool = all_heap_size;
    } else if (random_mode) {
      active_pool = all_remaining;
    } else if (static_mode) {
      active_pool = static_biased_mode ? (static_total - static_cursor)
                                       : static_remaining;
    } else if (adaptive_biased_mode) {
      // Heap contains stale/inactive edges; core_nodes controls termination.
      active_pool = heap_size;
    } else {
      active_pool = core_edge_count;
    }
    if (!static_mode && !(random_mode || random_biased_mode) && core_nodes <= 0)
      break;
    if (active_pool <= 0)
      break;
    if (++safety_core_loop > max_core_loop) {
      fprintf(stderr, "Core loop safety hit (N=%d, m=%d, core_edges=%d)\n", n,
              m, core_edge_count);
      free(static_core_edges);
      free(static_scored);
      free(all_scored);
      free(all_edges_pool);
      free(active_degree);
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
    int edge_id = -1;
    if (random_biased_original_mode) {
      while (all_scored_cursor < m && edge_id < 0) {
        int candidate = all_scored[all_scored_cursor++].edge_id;
        if (!edges[candidate].active)
          continue;
        edge_id = candidate;
      }
      if (edge_id < 0)
        break;
    } else if (random_biased_adaptive_mode) {
      while (all_heap_size > 0 && edge_id < 0) {
        int cand = -1;
        int cand_score = 0;
        if (!heap_pop(heap_edges, heap_scores, &all_heap_size, &cand,
                      &cand_score)) {
          break;
        }
        if (!edges[cand].active) {
          continue;
        }
        int u = edges[cand].u;
        int v = edges[cand].v;
        int current = active_degree[u] + active_degree[v];
        if (current < cand_score) {
          heap_push(heap_edges, heap_scores, &all_heap_size, cand, current);
          continue;
        }
        edge_id = cand;
      }
      if (edge_id < 0)
        break;
    } else if (random_mode) {
      while (all_remaining > 0 && edge_id < 0) {
        int idx = (int)rng_range(&rng_local, (uint32_t)all_remaining);
        int candidate = all_edges_pool[idx];
        if (!edges[candidate].active) {
          all_edges_pool[idx] = all_edges_pool[all_remaining - 1];
          all_remaining--;
          continue;
        }
        edge_id = candidate;
        all_edges_pool[idx] = all_edges_pool[all_remaining - 1];
        all_remaining--;
      }
      if (edge_id < 0) {
        break;
      }
    } else if (static_mode) {
      // Pick from the original core edge list, skipping any already inactive.
      if (static_biased_mode) {
        while (static_cursor < static_total && edge_id < 0) {
          int candidate = static_scored[static_cursor++].edge_id;
          if (!edges[candidate].active)
            continue;
          edge_id = candidate;
        }
      } else {
        while (static_remaining > 0 && edge_id < 0) {
          int idx = (int)rng_range(&rng_local, (uint32_t)static_remaining);
          int candidate = static_core_edges[idx];
          if (!edges[candidate].active) {
            static_core_edges[idx] = static_core_edges[static_remaining - 1];
            static_remaining--;
            continue;
          }
          edge_id = candidate;
          static_core_edges[idx] = static_core_edges[static_remaining - 1];
          static_remaining--;
        }
      }
      if (edge_id < 0) {
        break;
      }
    } else if (adaptive_biased_mode) {
      int stored = 0;
      while (heap_size > 0 && edge_id < 0) {
        int cand = -1;
        int cand_score = 0;
        if (!heap_pop(heap_edges, heap_scores, &heap_size, &cand, &cand_score))
          break;
        if (!edges[cand].active || !edges[cand].in_core) {
          continue;
        }
        int u = edges[cand].u;
        int v = edges[cand].v;
        int current = core_degree[u] + core_degree[v];
        if (current < cand_score) {
          heap_push(heap_edges, heap_scores, &heap_size, cand, current);
          continue;
        }
        edge_id = cand;
        stored = current;
      }
      (void)stored;
      if (edge_id < 0)
        break;
    } else {
      int idx = rng_range(&rng_local, (uint32_t)core_edge_count);
      edge_id = core_edges_arr[idx];
    }
    if (!edges[edge_id].active) {
      core_edge_remove(edge_id, edges, core_edges_arr, core_pos,
                       &core_edge_count, maintain_core_list);
      continue;
    }

    // Physically remove only the chosen edge.
    int was_core_edge = edges[edge_id].in_core;
    edges[edge_id].active = 0;
    if (random_biased_adaptive_mode) {
      int u = edges[edge_id].u;
      int v = edges[edge_id].v;
      if (active_degree[u] > 0)
        active_degree[u]--;
      if (active_degree[v] > 0)
        active_degree[v]--;
    }
    core_edge_remove(edge_id, edges, core_edges_arr, core_pos,
                     &core_edge_count, maintain_core_list);
    int u = edges[edge_id].u;
    int v = edges[edge_id].v;
    if (was_core_edge && in_core[u]) {
      core_degree[u]--;
      enqueue_if_needed(u, queue, &qtail, in_queue, core_degree);
    }
    if (was_core_edge && in_core[v]) {
      core_degree[v]--;
      enqueue_if_needed(v, queue, &qtail, in_queue, core_degree);
    }

    while (qhead < qtail) {
      int node = queue[qhead++];
      in_queue[node] = 0;
      remove_core_node(node, edges, incident_offsets, incident_edges, in_core,
                       core_degree, &core_nodes, core_edges_arr, core_pos,
                       &core_edge_count, maintain_core_list, queue, &qtail,
                       in_queue);
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
  if (!(random_mode || random_biased_mode)) {
    long long safety_tree_loop = 0;
    int remaining = 0;
    for (int e = 0; e < m; e++) {
      if (edges[e].active) {
        core_edges_arr[remaining++] = e;
      }
    }
    long long max_tree_loop = 10LL * remaining + 100000;
    for (int i = remaining - 1; i >= 0; i--) {
      if (++safety_tree_loop > max_tree_loop) {
        fprintf(stderr,
                "Tree loop safety hit (N=%d, m=%d, remaining=%d, i=%d)\n", n,
                m, remaining, i);
        free(static_core_edges);
        free(static_scored);
        free(all_scored);
        free(all_edges_pool);
        free(active_degree);
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
      int r = (int)rng_range(&rng_local, (uint32_t)(i + 1));
      int chosen = core_edges_arr[r];
      core_edges_arr[r] = core_edges_arr[i];
      core_edges_arr[i] = chosen;
      edges[chosen].active = 0;
      result->removal_order[removed_edges++] = chosen;
      sample_core_series(result, removed_edges, 0);
    }
  }

  result->removal_length = removed_edges;
  if (removed_edges != m) {
    fprintf(stderr,
            "Error: removal order length mismatch (removed=%d, m=%d)\n",
            removed_edges, m);
    free(static_core_edges);
    free(static_scored);
    free(all_scored);
    free(all_edges_pool);
    free(active_degree);
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
  double zero_p = result->core_zero_p;
  for (int i = 0; i < result->series_length; i++) {
    double p_here = (result->series_length <= 1)
                        ? 1.0
                        : (double)i / (double)(result->series_length - 1);
    if (result->core_size[i] <= 0.0) {
      if (zero_p >= 0.0 && p_here >= zero_p) {
        result->core_size[i] = 0.0; // past core extinction, stay zero
        last_core = 0.0;
      } else if (i > 0) {
        result->core_size[i] = last_core; // fill early sparse samples
      }
    } else {
      last_core = result->core_size[i];
    }
  }
  // Ensure final sample at p=1 reflects empty 2-core.
  if (result->series_length > 0) {
    result->core_size[result->series_length - 1] = 0.0;
  }

  reconstruct_with_newman_ziff(g, result->removal_order, result->removal_length,
                               result, gc_ecp);

  free(static_core_edges);
  free(static_scored);
  free(all_scored);
  free(all_edges_pool);
  free(active_degree);
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
