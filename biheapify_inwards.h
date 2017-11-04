/*
 * biheapify_inwards.h
 *
 *  Created on: Oct 14, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_INWARDS_H_
#define BIHEAPIFY_INWARDS_H_

#include <algorithm>

#include "biheapify.h"

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwards(RAI first, size_type total_num_nodes) {
  while (total_num_nodes > 10) {
    BiHeapify(first, total_num_nodes);
    size_type heap_size              = HeapSize(total_num_nodes);
    size_type first_extended_in_node = total_num_nodes - heap_size;//= FLIP(num_nodes_in_heap - 1);
    first                           += first_extended_in_node;
    total_num_nodes                  = heap_size - first_extended_in_node;
  }
  if (total_num_nodes > 1)
    std::sort(first, first + total_num_nodes);
  return ;
}

#endif /* BIHEAPIFY_INWARDS_H_ */
