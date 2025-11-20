#ifndef NODE_PERCOLATION_STATS_H
#define NODE_PERCOLATION_STATS_H

#include <stddef.h>

typedef struct {
  double *mean;
  double *m2;
  size_t length;
  long long count;
} OnlineStats;

OnlineStats *create_online_stats(size_t length);
void free_online_stats(OnlineStats *stats);
void online_stats_update(OnlineStats *stats, const double *values);
void online_stats_get_std(const OnlineStats *stats, double *std_dev);

#endif // NODE_PERCOLATION_STATS_H
