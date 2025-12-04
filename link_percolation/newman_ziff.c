#include "newman_ziff.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int n;
  int *parent;
  int *size;
  long long sum_s_squared;
  int num_components;
  int largest_component;
  long long active_edges;
} NZUnionFind;

typedef struct {
  double giant;
  double avg_small;
  int num_components;
  int num_cycles;
} Snapshot;

static NZUnionFind *nzuf_create(int n) {
  NZUnionFind *uf = (NZUnionFind *)calloc(1, sizeof(NZUnionFind));
  if (!uf)
    return NULL;
  uf->n = n;
  uf->parent = (int *)malloc((size_t)n * sizeof(int));
  uf->size = (int *)malloc((size_t)n * sizeof(int));
  if (!uf->parent || !uf->size) {
    free(uf->parent);
    free(uf->size);
    free(uf);
    return NULL;
  }
  for (int i = 0; i < n; i++) {
    uf->parent[i] = i;
    uf->size[i] = 1;
  }
  uf->sum_s_squared = n;
  uf->num_components = n;
  uf->largest_component = n > 0 ? 1 : 0;
  uf->active_edges = 0;
  return uf;
}

static void nzuf_free(NZUnionFind *uf) {
  if (!uf)
    return;
  free(uf->parent);
  free(uf->size);
  free(uf);
}

static int nzuf_find(NZUnionFind *uf, int x) {
  int root = x;
  while (uf->parent[root] != root) {
    root = uf->parent[root];
  }
  while (x != root) {
    int next = uf->parent[x];
    uf->parent[x] = root;
    x = next;
  }
  return root;
}

static void nzuf_union(NZUnionFind *uf, int a, int b) {
  int root_a = nzuf_find(uf, a);
  int root_b = nzuf_find(uf, b);
  if (root_a == root_b)
    return;

  int size_a = uf->size[root_a];
  int size_b = uf->size[root_b];
  uf->sum_s_squared -= (long long)size_a * size_a;
  uf->sum_s_squared -= (long long)size_b * size_b;

  if (size_a < size_b) {
    int tmp = root_a;
    root_a = root_b;
    root_b = tmp;
    int tmp_size = size_a;
    size_a = size_b;
    size_b = tmp_size;
  }

  uf->parent[root_b] = root_a;
  uf->size[root_a] = size_a + size_b;
  long long new_size = uf->size[root_a];
  uf->sum_s_squared += new_size * new_size;
  if (uf->size[root_a] > uf->largest_component) {
    uf->largest_component = uf->size[root_a];
  }
  uf->num_components--;
}

static double nzuf_avg_small_size(const NZUnionFind *uf) {
  if (uf->largest_component == 0 || uf->num_components <= 0)
    return 0.0;
  int nodes_in_small = uf->n - uf->largest_component;
  if (nodes_in_small <= 0)
    return 0.0;
  long long small_sq =
      uf->sum_s_squared -
      (long long)uf->largest_component * uf->largest_component;
  if (small_sq <= 0)
    return 0.0;
  return (double)small_sq / (double)nodes_in_small;
}

static void fill_snapshot(Snapshot *snap, const NZUnionFind *uf,
                          double avg_small) {
  snap->giant = (double)uf->largest_component;
  snap->avg_small = avg_small;
  snap->num_components = uf->num_components;
  long long cycles = uf->active_edges - uf->n + uf->num_components;
  if (cycles < 0)
    cycles = 0;
  snap->num_cycles = (int)cycles;
}

static GraphEdge *enumerate_edges(const Graph *g) {
  GraphEdge *edges =
      (GraphEdge *)malloc((size_t)g->m * sizeof(GraphEdge));
  if (!edges)
    return NULL;
  int edge_id = 0;
  for (int u = 0; u < g->n; u++) {
    const int *neighbors = graph_neighbors(g, u);
    int deg = graph_degree(g, u);
    for (int i = 0; i < deg; i++) {
      int v = neighbors[i];
      if (u < v) {
        if (edge_id >= g->m) {
          free(edges);
          return NULL;
        }
        edges[edge_id].u = u;
        edges[edge_id].v = v;
        edge_id++;
      }
    }
  }
  if (edge_id != g->m) {
    free(edges);
    return NULL;
  }
  return edges;
}

void reconstruct_with_newman_ziff(Graph *g, const int *removal_order,
                                  int order_length,
                                  LinkPercolationResult *result,
                                  PseudoCriticalPoints *ecp) {
  int n = g->n;
  int m = order_length;
  NZUnionFind *uf = nzuf_create(n);
  GraphEdge *edges = enumerate_edges(g);
  if (!uf || !edges) {
    nzuf_free(uf);
    free(edges);
    return;
  }

  unsigned char sampled[NR_SERIES_LENGTH] = {0};

  Snapshot snap = {0};
  Snapshot p_plus_snap = {0}, q_plus_snap = {0}, q_zero_snap = {0},
           q_minus_snap = {0};
  int p_plus_idx = 0, q_plus_idx = 0, q_zero_idx = 0, q_minus_idx = 0;
  int p1_idx = -1, p2_idx = -1;
  const int threshold_p1 = n / 2;
  const double threshold_p2 = sqrt((double)n);

  double prev_gc = 0.0;
  double prev_small = 0.0;
  double max_gc_jump = 0.0;
  double max_positive_small = 0.0;
  double max_negative_small = 0.0;
  double max_small_size = 0.0;

  const int samples = result->series_length - 1;
  int sample_idx0 = samples;
  if (sample_idx0 >= 0 && sample_idx0 < result->series_length) {
    double avg_small = nzuf_avg_small_size(uf);
    result->giant_size[sample_idx0] = (double)uf->largest_component;
    result->avg_small_size[sample_idx0] = avg_small;
    result->num_components[sample_idx0] = uf->num_components;
    result->num_cycles[sample_idx0] =
        (int)(uf->active_edges - uf->n + uf->num_components);
    result->lambda1[sample_idx0] = -1.0;
    sampled[sample_idx0] = 1;
  }

  for (int idx = 1; idx <= m; idx++) {
    int edge_id = removal_order[m - idx];
    int u = edges[edge_id].u;
    int v = edges[edge_id].v;

    uf->active_edges++;
    nzuf_union(uf, u, v);

    double giant = (double)uf->largest_component;
    double avg_small = nzuf_avg_small_size(uf);
    fill_snapshot(&snap, uf, avg_small);

    double gc_jump = giant - prev_gc;
    if (gc_jump > max_gc_jump) {
      max_gc_jump = gc_jump;
      p_plus_idx = idx;
      p_plus_snap = snap;
    }

    double small_jump = avg_small - prev_small;
    if (small_jump > max_positive_small) {
      max_positive_small = small_jump;
      q_plus_idx = idx;
      q_plus_snap = snap;
    }
    if (small_jump < max_negative_small) {
      max_negative_small = small_jump;
      q_minus_idx = idx;
      q_minus_snap = snap;
    }
    if (avg_small > max_small_size) {
      max_small_size = avg_small;
      q_zero_idx = idx;
      q_zero_snap = snap;
    }

    if (giant < (double)threshold_p1) {
      p1_idx = idx;
    }
    if (giant < threshold_p2) {
      p2_idx = idx;
    }

    int sample_idx =
        (int)lround((double)(m - idx) * (double)samples / (double)m);
    if (sample_idx < 0)
      sample_idx = 0;
    if (sample_idx > samples)
      sample_idx = samples;
    if (!sampled[sample_idx]) {
      result->giant_size[sample_idx] = snap.giant;
      result->avg_small_size[sample_idx] = snap.avg_small;
      result->num_components[sample_idx] = snap.num_components;
      result->num_cycles[sample_idx] = snap.num_cycles;
      result->lambda1[sample_idx] = -1.0;
      sampled[sample_idx] = 1;
    }

    prev_gc = giant;
    prev_small = avg_small;
  }

  nzuf_free(uf);
  free(edges);

  if (!ecp)
    return;

  ecp->p_plus = 1.0 - (double)p_plus_idx / (double)m;
  ecp->p_plus_jump = max_gc_jump / (double)n;
  ecp->q_plus = 1.0 - (double)q_plus_idx / (double)m;
  ecp->q_zero = 1.0 - (double)q_zero_idx / (double)m;
  ecp->q_minus = 1.0 - (double)q_minus_idx / (double)m;
  ecp->p1 = (p1_idx >= 0) ? 1.0 - (double)(p1_idx + 1) / (double)m : 0.0;
  ecp->p2 = (p2_idx >= 0) ? 1.0 - (double)(p2_idx + 1) / (double)m : 0.0;

  ecp->P_at_p_plus = (int)p_plus_snap.giant;
  ecp->s_at_p_plus = p_plus_snap.avg_small;
  ecp->nc_at_p_plus = p_plus_snap.num_components;
  ecp->nl_at_p_plus = p_plus_snap.num_cycles;

  ecp->P_at_q_plus = (int)q_plus_snap.giant;
  ecp->s_at_q_plus = q_plus_snap.avg_small;
  ecp->nc_at_q_plus = q_plus_snap.num_components;
  ecp->nl_at_q_plus = q_plus_snap.num_cycles;

  ecp->P_at_q_zero = (int)q_zero_snap.giant;
  ecp->s_at_q_zero = q_zero_snap.avg_small;
  ecp->nc_at_q_zero = q_zero_snap.num_components;
  ecp->nl_at_q_zero = q_zero_snap.num_cycles;

  ecp->P_at_q_minus = (int)q_minus_snap.giant;
  ecp->s_at_q_minus = q_minus_snap.avg_small;
  ecp->nc_at_q_minus = q_minus_snap.num_components;
  ecp->nl_at_q_minus = q_minus_snap.num_cycles;
}
