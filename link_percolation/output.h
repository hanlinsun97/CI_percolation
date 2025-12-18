#ifndef LINK_PERCOLATION_OUTPUT_H
#define LINK_PERCOLATION_OUTPUT_H

#include "stats.h"

typedef struct {
  double p_plus;
  double p_plus_jump;
  double q_plus;
  double q_zero;
  double q_minus;
  double p1;
  double p2;

  int P_at_p_plus;
  double s_at_p_plus;
  int nc_at_p_plus;
  int nl_at_p_plus;

  int P_at_q_plus;
  double s_at_q_plus;
  int nc_at_q_plus;
  int nl_at_q_plus;

  int P_at_q_zero;
  double s_at_q_zero;
  int nc_at_q_zero;
  int nl_at_q_zero;

  int P_at_q_minus;
  double s_at_q_minus;
  int nc_at_q_minus;
  int nl_at_q_minus;
} PseudoCriticalPoints;

typedef struct {
  double p_plus, p_plus_jump, q_plus, q_zero, q_minus, p1, p2;
  int P_p_plus, nc_p_plus, nl_p_plus;
  double s_p_plus;
  int P_q_plus, nc_q_plus, nl_q_plus;
  double s_q_plus;
  int P_q_zero, nc_q_zero, nl_q_zero;
  double s_q_zero;
  int P_q_minus, nc_q_minus, nl_q_minus;
  double s_q_minus;

  double core_p_star;
  int core_drop;
  int core_size_star;
  double core_zero_p;
} EnhancedEBEData;

#define DEFAULT_OUTPUT_SAMPLES 100

typedef struct EBEWriter EBEWriter;

void save_averaged_results(const char *filename, int n, int m,
                           int M_graphs, int M_realizations,
                           OnlineStats *gc_stats, OnlineStats *small_stats,
                           OnlineStats *lambda_stats, OnlineStats *comp_stats,
                           OnlineStats *cycles_stats,
                           OnlineStats *core_stats);

EBEWriter *ebe_writer_begin(const char *filename, int n, double c,
                            int M_graphs, int M_realizations,
                            int total_runs);
void ebe_writer_record(EBEWriter *writer, int run_index,
                       const EnhancedEBEData *record);
void ebe_writer_finish(EBEWriter *writer);
void ebe_writer_abort(EBEWriter *writer);

void save_enhanced_ebe_results(const char *filename, int n, double c,
                               int M_graphs, int M_realizations,
                               const EnhancedEBEData *ebe_data,
                               int total_runs);

#endif // LINK_PERCOLATION_OUTPUT_H
