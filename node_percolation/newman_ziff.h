#ifndef NODE_PERCOLATION_NEWMAN_ZIFF_H
#define NODE_PERCOLATION_NEWMAN_ZIFF_H

#include "graph.h"
#include "output.h"

struct NodePercolationResult;

void reconstruct_with_newman_ziff(Graph *g, const int *removal_order,
                                  int order_length,
                                  struct NodePercolationResult *result,
                                  PseudoCriticalPoints *ecp);

#endif // NODE_PERCOLATION_NEWMAN_ZIFF_H
