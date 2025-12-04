#ifndef LINK_PERCOLATION_NEWMAN_ZIFF_H
#define LINK_PERCOLATION_NEWMAN_ZIFF_H

#include "graph.h"
#include "output.h"
#include "removal.h"

void reconstruct_with_newman_ziff(Graph *g, const int *removal_order,
                                  int order_length,
                                  LinkPercolationResult *result,
                                  PseudoCriticalPoints *ecp);

#endif // LINK_PERCOLATION_NEWMAN_ZIFF_H
