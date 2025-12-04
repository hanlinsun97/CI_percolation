#include "stats.h"
#include <math.h>
#include <stdlib.h>

OnlineStats *create_online_stats(size_t length) {
  if (length == 0)
    return NULL;
  OnlineStats *stats = (OnlineStats *)calloc(1, sizeof(OnlineStats));
  if (!stats)
    return NULL;
  stats->mean = (double *)calloc(length, sizeof(double));
  stats->m2 = (double *)calloc(length, sizeof(double));
  if (!stats->mean || !stats->m2) {
    free(stats->mean);
    free(stats->m2);
    free(stats);
    return NULL;
  }
  stats->length = length;
  stats->count = 0;
  return stats;
}

void free_online_stats(OnlineStats *stats) {
  if (!stats)
    return;
  free(stats->mean);
  free(stats->m2);
  free(stats);
}

void online_stats_update(OnlineStats *stats, const double *values) {
  if (!stats || !values)
    return;
  stats->count++;
  double inv = 1.0 / (double)stats->count;
  for (size_t i = 0; i < stats->length; i++) {
    double delta = values[i] - stats->mean[i];
    stats->mean[i] += delta * inv;
    double delta2 = values[i] - stats->mean[i];
    stats->m2[i] += delta * delta2;
  }
}

void online_stats_get_std(const OnlineStats *stats, double *std_dev) {
  if (!stats || !std_dev)
    return;
  if (stats->count < 2) {
    for (size_t i = 0; i < stats->length; i++) {
      std_dev[i] = 0.0;
    }
    return;
  }
  double denom = (double)(stats->count - 1);
  for (size_t i = 0; i < stats->length; i++) {
    double variance = stats->m2[i] / denom;
    std_dev[i] = sqrt(variance > 0.0 ? variance : 0.0);
  }
}
