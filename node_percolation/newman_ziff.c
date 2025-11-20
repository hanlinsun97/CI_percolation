#include "newman_ziff.h"
#include "removal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int n;
  int *parent;
  int *size;
  unsigned char *active;
  int active_nodes;
  long long active_edges;
  long long sum_s_squared;
  int num_components;
  int largest_component;
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
  uf->parent = (int *)malloc(n * sizeof(int));
  uf->size = (int *)malloc(n * sizeof(int));
  uf->active = (unsigned char *)calloc(n, sizeof(unsigned char));
  if (!uf->parent || !uf->size || !uf->active) {
    free(uf->parent);
    free(uf->size);
    free(uf->active);
    free(uf);
    return NULL;
  }
  return uf;
}

static void nzuf_free(NZUnionFind *uf) {
  if (!uf)
    return;
  free(uf->parent);
  free(uf->size);
  free(uf->active);
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

static void nzuf_activate(NZUnionFind *uf, int node) {
  uf->active[node] = 1;
  uf->parent[node] = node;
  uf->size[node] = 1;
  uf->active_nodes++;
  uf->num_components++;
  uf->sum_s_squared += 1;
  if (uf->largest_component < 1)
    uf->largest_component = 1;
}

static int nzuf_is_active(const NZUnionFind *uf, int node) {
  return uf->active[node] != 0;
}

static double nzuf_avg_small_size(const NZUnionFind *uf) {
  if (uf->largest_component == 0 || uf->num_components <= 0)
    return 0.0;
  int nodes_in_small = uf->active_nodes - uf->largest_component;
  if (nodes_in_small <= 0)
    return 0.0;
  long long small_sq =
      uf->sum_s_squared - (long long)uf->largest_component * uf->largest_component;
  if (small_sq <= 0)
    return 0.0;
  return (double)small_sq / (double)nodes_in_small;
}

static void fill_snapshot(Snapshot *snap, const NZUnionFind *uf,
                          double avg_small, long long cycles) {
  snap->giant = (double)uf->largest_component;
  snap->avg_small = avg_small;
  snap->num_components = uf->num_components;
  snap->num_cycles = (int)cycles;
}

void reconstruct_with_newman_ziff(Graph *g, const int *removal_order,
                                  int order_length,
                                  NodePercolationResult *result,
                                  PseudoCriticalPoints *ecp) {
  int n = order_length;
  NZUnionFind *uf = nzuf_create(n);
  if (!uf)
    return;

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

  int samples = result->series_length - 1;
  int sample_idx0 = samples; // removed fraction 1 corresponds to idx=0 active
  if (sample_idx0 >= 0 && sample_idx0 < result->series_length) {
    result->giant_size[sample_idx0] = 0.0;
    result->avg_small_size[sample_idx0] = 0.0;
    result->num_components[sample_idx0] = 0;
    result->num_cycles[sample_idx0] = 0;
    result->lambda1[sample_idx0] = -1.0;
    sampled[sample_idx0] = 1;
  }

  for (int idx = 1; idx <= n; idx++) {
    int node = removal_order[n - idx];
    nzuf_activate(uf, node);

    const int *neighbors = graph_neighbors(g, node);
    int degree = graph_degree(g, node);
    for (int j = 0; j < degree; j++) {
      int neighbor = neighbors[j];
      if (!nzuf_is_active(uf, neighbor))
        continue;
      uf->active_edges++;
      nzuf_union(uf, node, neighbor);
    }

    double giant = (double)uf->largest_component;
    double avg_small = nzuf_avg_small_size(uf);
    long long cycles =
        uf->active_edges - uf->active_nodes + uf->num_components;
    if (cycles < 0)
      cycles = 0;
    fill_snapshot(&snap, uf, avg_small, cycles);

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

    int sample_idx = (int)lround((double)(n - idx) * samples / (double)n);
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
  if (!ecp)
    return;

  ecp->p_plus = 1.0 - (double)p_plus_idx / (double)n;
  ecp->p_plus_jump = max_gc_jump / (double)n;
  ecp->q_plus = 1.0 - (double)q_plus_idx / (double)n;
  ecp->q_zero = 1.0 - (double)q_zero_idx / (double)n;
  ecp->q_minus = 1.0 - (double)q_minus_idx / (double)n;
  ecp->p1 = (p1_idx >= 0) ? 1.0 - (double)(p1_idx + 1) / (double)n : 0.0;
  ecp->p2 = (p2_idx >= 0) ? 1.0 - (double)(p2_idx + 1) / (double)n : 0.0;

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
