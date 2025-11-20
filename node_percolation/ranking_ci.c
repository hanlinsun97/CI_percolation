#include "ranking.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
  RANKER_CI,
  RANKER_RANDOM
} RankerType;

typedef struct NodeRanker {
  RankerType type;
  Graph *g;
  bool *removed;
  int *active_degree;
  int radius;
  int n;
  double *ci_values;
  int *heap;
  int *heap_pos;
  int heap_size;
  int *queue;
  int *level_counts;
  uint64_t *visited_bits;
  uint64_t *dirty_bits;
  int bitset_words;
} CIRanker;

typedef struct {
  RankerType type;
  Graph *g;
  bool *removed;
  int *active_nodes;
  int active_count;
  rng_state_t rng;
} RandomRanker;

static inline void bitset_set(uint64_t *bits, int idx) {
  bits[idx >> 6] |= (uint64_t)1 << (idx & 63);
}

static inline void bitset_clear(uint64_t *bits, int idx) {
  bits[idx >> 6] &= ~((uint64_t)1 << (idx & 63));
}

static inline int bitset_test(const uint64_t *bits, int idx) {
  return (int)((bits[idx >> 6] >> (idx & 63)) & 1U);
}

static inline void swap_heap(CIRanker *ranker, int i, int j) {
  int node_i = ranker->heap[i];
  int node_j = ranker->heap[j];
  ranker->heap[i] = node_j;
  ranker->heap[j] = node_i;
  ranker->heap_pos[node_i] = j;
  ranker->heap_pos[node_j] = i;
}

static void heapify_down(CIRanker *ranker, int idx) {
  while (1) {
    int left = 2 * idx + 1;
    int right = left + 1;
    int largest = idx;
    if (left < ranker->heap_size &&
        ranker->ci_values[ranker->heap[left]] >
            ranker->ci_values[ranker->heap[largest]]) {
      largest = left;
    }
    if (right < ranker->heap_size &&
        ranker->ci_values[ranker->heap[right]] >
            ranker->ci_values[ranker->heap[largest]]) {
      largest = right;
    }
    if (largest == idx)
      break;
    swap_heap(ranker, idx, largest);
    idx = largest;
  }
}

static void heap_build(CIRanker *ranker) {
  for (int i = ranker->heap_size / 2 - 1; i >= 0; i--) {
    heapify_down(ranker, i);
  }
}

static double compute_ci(CIRanker *ranker, int node) {
  if (ranker->removed[node])
    return 0.0;

  int deg = ranker->active_degree[node];
  if (ranker->radius <= 0) {
    return deg > 0 ? (double)deg : 0.0;
  }
  if (deg <= 1)
    return 0.0;

  int L = ranker->radius;
  int *queue = ranker->queue;
  int *levels = ranker->level_counts;

  int head = 0, tail = 0;
  queue[tail++] = node;
  bitset_set(ranker->visited_bits, node);
  memset(levels, 0, (L + 1) * sizeof(int));
  levels[0] = 1;
  int depth = 0;
  int processed_at_level = 0;
  int next_level_count = 0;

  while (head < tail && depth < L) {
    int current = queue[head++];
    processed_at_level++;
    const int *neighbors = graph_neighbors(ranker->g, current);
    int deg_adj = graph_degree(ranker->g, current);
    for (int i = 0; i < deg_adj; i++) {
      int neighbor = neighbors[i];
      if (ranker->removed[neighbor])
        continue;
      if (!bitset_test(ranker->visited_bits, neighbor)) {
        bitset_set(ranker->visited_bits, neighbor);
        queue[tail++] = neighbor;
        next_level_count++;
      }
    }
    if (processed_at_level == levels[depth]) {
      depth++;
      if (depth > L)
        break;
      levels[depth] = next_level_count;
      processed_at_level = 0;
      next_level_count = 0;
    }
  }

  int start_idx = 0;
  for (int d = 0; d < L; d++)
    start_idx += levels[d];

  double boundary_sum = 0.0;
  for (int i = 0; i < levels[L]; i++) {
    int idx = start_idx + i;
    int boundary_node = queue[idx];
    boundary_sum += (double)(ranker->active_degree[boundary_node] - 1);
  }

  double ci_value = (double)(ranker->active_degree[node] - 1) * boundary_sum;

  for (int i = 0; i < tail; i++) {
    bitset_clear(ranker->visited_bits, queue[i]);
  }
  return (ci_value > 0.0) ? ci_value : 0.0;
}

static int collect_nodes_to_update(CIRanker *ranker, int source) {
  int L = ranker->radius;
  int *queue = ranker->queue;
  int head = 0, tail = 0;
  queue[tail++] = source;
  bitset_set(ranker->visited_bits, source);
  int depth = 0;
  int nodes_remaining = 1;
  int next_level = 0;

  while (head < tail) {
    int current = queue[head++];
    nodes_remaining--;
    if (depth < L) {
      const int *neighbors = graph_neighbors(ranker->g, current);
      int deg_adj = graph_degree(ranker->g, current);
      for (int i = 0; i < deg_adj; i++) {
        int neighbor = neighbors[i];
        if (ranker->removed[neighbor] ||
            bitset_test(ranker->visited_bits, neighbor)) {
          continue;
        }
        bitset_set(ranker->visited_bits, neighbor);
        queue[tail++] = neighbor;
        next_level++;
      }
    }
    if (nodes_remaining == 0) {
      depth++;
      if (depth > L)
        break;
      nodes_remaining = next_level;
      next_level = 0;
    }
  }
  for (int i = 0; i < tail; i++) {
    bitset_clear(ranker->visited_bits, queue[i]);
  }
  return tail;
}

static int heap_pop(CIRanker *ranker) {
  if (ranker->heap_size == 0)
    return -1;
  int node = ranker->heap[0];
  ranker->heap[0] = ranker->heap[ranker->heap_size - 1];
  ranker->heap_pos[ranker->heap[0]] = 0;
  ranker->heap_size--;
  heapify_down(ranker, 0);
  ranker->heap_pos[node] = -1;
  return node;
}

static NodeRanker *ci_ranker_create(Graph *g, const RankingConfig *config,
                                    bool *removed, int *active_degree,
                                    rng_state_t *rng) {
  (void)rng;
  CIRanker *ranker = (CIRanker *)calloc(1, sizeof(CIRanker));
  if (!ranker)
    return NULL;
  ranker->g = g;
  ranker->removed = removed;
  ranker->active_degree = active_degree;
  ranker->type = RANKER_CI;
  ranker->radius = config->ci_radius >= 0 ? config->ci_radius : 0;
  ranker->n = g->n;
  ranker->heap_size = g->n;
  ranker->ci_values = (double *)calloc(g->n, sizeof(double));
  ranker->heap = (int *)malloc(g->n * sizeof(int));
  ranker->heap_pos = (int *)malloc(g->n * sizeof(int));
  ranker->queue = (int *)malloc(g->n * sizeof(int));
  ranker->level_counts = (int *)calloc(ranker->radius + 1, sizeof(int));
  ranker->bitset_words = (g->n + 63) / 64;
  if (ranker->bitset_words == 0)
    ranker->bitset_words = 1;
  ranker->visited_bits =
      (uint64_t *)calloc((size_t)ranker->bitset_words, sizeof(uint64_t));
  ranker->dirty_bits =
      (uint64_t *)calloc((size_t)ranker->bitset_words, sizeof(uint64_t));

  if (!ranker->ci_values || !ranker->heap || !ranker->heap_pos ||
      !ranker->queue || !ranker->level_counts || !ranker->visited_bits ||
      !ranker->dirty_bits) {
    fprintf(stderr, "Error: Unable to allocate CI ranker data.\n");
    free(ranker->ci_values);
    free(ranker->heap);
    free(ranker->heap_pos);
    free(ranker->queue);
    free(ranker->level_counts);
    free(ranker->visited_bits);
    free(ranker->dirty_bits);
    free(ranker);
    return NULL;
  }

  for (int i = 0; i < g->n; i++) {
    ranker->heap[i] = i;
    ranker->heap_pos[i] = i;
    ranker->ci_values[i] = compute_ci(ranker, i);
  }
  heap_build(ranker);
  return (NodeRanker *)ranker;
}

static NodeRanker *random_ranker_create(Graph *g, bool *removed,
                                        rng_state_t *rng) {
  RandomRanker *ranker = (RandomRanker *)calloc(1, sizeof(RandomRanker));
  if (!ranker)
    return NULL;
  ranker->type = RANKER_RANDOM;
  ranker->g = g;
  ranker->removed = removed;
  ranker->active_count = g->n;
  ranker->active_nodes = (int *)malloc(g->n * sizeof(int));
  if (!ranker->active_nodes) {
    free(ranker);
    return NULL;
  }
  for (int i = 0; i < g->n; i++) {
    ranker->active_nodes[i] = i;
  }
  ranker->rng = *rng;
  return (NodeRanker *)ranker;
}

NodeRanker *node_ranker_create(Graph *g, const RankingConfig *config,
                               bool *node_removed, int *active_degree,
                               rng_state_t *rng) {
  if (!config || !config->strategy) {
    fprintf(stderr, "Error: Ranking configuration is missing.\n");
    return NULL;
  }
  if (strcmp(config->strategy, "ci") == 0) {
    return ci_ranker_create(g, config, node_removed, active_degree, rng);
  } else if (strcmp(config->strategy, "random") == 0) {
    return random_ranker_create(g, node_removed, rng);
  }
  fprintf(stderr, "Error: Unknown ranking strategy '%s'.\n",
          config->strategy);
  return NULL;
}

int node_ranker_next(NodeRanker *base_ranker) {
  if (base_ranker->type == RANKER_CI) {
    CIRanker *ranker = (CIRanker *)base_ranker;
    while (ranker->heap_size > 0) {
      int node = ranker->heap[0];
      if (ranker->removed[node]) {
        heap_pop(ranker);
        continue;
      }
      if (bitset_test(ranker->dirty_bits, node)) {
        ranker->ci_values[node] = compute_ci(ranker, node);
        bitset_clear(ranker->dirty_bits, node);
        heapify_down(ranker, 0);
        continue;
      }
      heap_pop(ranker);
      return node;
    }
    return -1;
  } else {
    RandomRanker *ranker = (RandomRanker *)base_ranker;
    if (ranker->active_count == 0)
      return -1;
    int idx = rng_range(&ranker->rng, (uint32_t)ranker->active_count);
    int node = ranker->active_nodes[idx];
    ranker->active_nodes[idx] =
        ranker->active_nodes[ranker->active_count - 1];
    ranker->active_count--;
    return node;
  }
}

void node_ranker_handle_removal(NodeRanker *base_ranker, int removed_node) {
  if (base_ranker->type == RANKER_RANDOM)
    return;
  CIRanker *ranker = (CIRanker *)base_ranker;
  int count = collect_nodes_to_update(ranker, removed_node);
  int *queue = ranker->queue;
  for (int i = 0; i < count; i++) {
    int node = queue[i];
    if (ranker->removed[node]) {
      ranker->ci_values[node] = 0.0;
      bitset_clear(ranker->dirty_bits, node);
      continue;
    }
    bitset_set(ranker->dirty_bits, node);
  }
}

void node_ranker_free(NodeRanker *base_ranker) {
  if (!base_ranker)
    return;
  if (base_ranker->type == RANKER_CI) {
    CIRanker *ranker_ci = (CIRanker *)base_ranker;
    free(ranker_ci->ci_values);
    free(ranker_ci->heap);
    free(ranker_ci->heap_pos);
    free(ranker_ci->queue);
    free(ranker_ci->level_counts);
    free(ranker_ci->visited_bits);
    free(ranker_ci->dirty_bits);
    free(ranker_ci);
  } else {
    RandomRanker *ranker = (RandomRanker *)base_ranker;
    free(ranker->active_nodes);
    free(ranker);
  }
}
