/*
 * biheapify_inwards_symmetrically_alternating_towards_ends.h
 *
 * The function BiHeapifyInwardsSymmetricallyAlternatingTowardsEnds()
 *  appears to always emplace the medians in the middle of the list
 *  and arranges the elements so that all elements before (resp. after)
 *  the middle element are less than or equal to (resp. greater than
 *  or equal to) the median value(s).
 *
 *  Created on: Nov 22, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_INWARDS_SYMMETRICALLY_ALTERNATING_TOWARDS_ENDS_H_
#define BIHEAPIFY_INWARDS_SYMMETRICALLY_ALTERNATING_TOWARDS_ENDS_H_

#include <algorithm>

#include "biheapify.h"
#include "biheapify_lambda.h"
#include "biheapify_flip_ordered.h"

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyTowardsStartAndEndParents(RAI V, size_type N) {
  size_type heap_size     = HeapSize(N);
  size_type first_in_node = N - heap_size;

  BiHeapifyFlipsOrdered<RAI, size_type>(V, N);
  BiHeapifyFlipsOrdered<RAI, size_type>(V + first_in_node, heap_size - first_in_node);

  size_type start_hc = first_in_node;
  size_type one_past_end_hc = N / 2 + 1;
  while (one_past_end_hc > 2) {
    BiHeapifyFlipsOrdered<RAI, size_type>(V, one_past_end_hc);
    BiHeapifyFlipsOrdered<RAI, size_type>(V + (N - one_past_end_hc), one_past_end_hc);

    size_type one_past_mid_point = (one_past_end_hc + start_hc) / 2 + 1;
    BiHeapifyFlipsOrdered<RAI, size_type>(V, one_past_mid_point);
    BiHeapifyFlipsOrdered<RAI, size_type>(V + (N - one_past_mid_point), one_past_mid_point);

    BiHeapifyFlipsOrdered<RAI, size_type>(V, start_hc);
    BiHeapifyFlipsOrdered<RAI, size_type>(V + ((N - 1) - start_hc), start_hc);
    size_type parent = Parent<size_type>(one_past_end_hc);
    one_past_end_hc = start_hc;
    start_hc = parent;
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifySymmetricallyAlternatingTowardsEnds(RAI V, size_type N) {
  BiHeapifyFlipsOrdered<RAI, size_type>(V, N);
  if (N < 4)
    return ;
  BiHeapifyTowardsStartAndEndParents<RAI, size_type>(V, N);
  size_type heap_size     = HeapSize(N);
  size_type first_in_node = N - heap_size;
  size_type current_half_height = 0;
  auto lambda = [](size_type N, size_type pos_hc) { return (N - 1) - pos_hc; };

  size_type last_node_hc = (N / 2) - 1;//(N % 2 == 0);
  size_type first_node_hc = first_in_node + (N % 3 == 2);
  while (last_node_hc > first_node_hc) {
    bool is_last_node_right_child = last_node_hc % 2 == 0;
    size_type num_nodes_to_biheapify = (last_node_hc + 1) - first_node_hc;

    BiHeapifyFlipsOrdered<RAI, size_type>(V, first_node_hc + num_nodes_to_biheapify);
    BiHeapifyFlipsOrdered<RAI, size_type>(V + ((N - 1) - last_node_hc), first_node_hc + num_nodes_to_biheapify);
    if (current_half_height % 4 == 0 || current_half_height % 4 == 1) {
      BiHeapifyFlipsOrdered<RAI, size_type>(V + first_node_hc, num_nodes_to_biheapify);
      BiHeapifyFlipsOrdered<RAI, size_type>(V + ((N - 1) - last_node_hc), num_nodes_to_biheapify, lambda);
    } else {
      BiHeapifyFlipsOrdered<RAI, size_type>(V + first_node_hc, num_nodes_to_biheapify, lambda);
      BiHeapifyFlipsOrdered<RAI, size_type>(V + ((N - 1) - last_node_hc), num_nodes_to_biheapify);
    }
    size_type parent_of_last_node_hc = Parent<size_type>(last_node_hc);
    last_node_hc = first_node_hc - 1;
    first_node_hc = parent_of_last_node_hc + is_last_node_right_child;
    current_half_height++;
  }
  BiHeapifyFlipsOrdered<RAI, size_type>(V, N);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsSymmetricallyAlternatingTowardsEnds(RAI V, size_type N) {
  while (N > 30) {
    BiHeapifySymmetricallyAlternatingTowardsEnds(V, N);
    size_type heap_size     = HeapSize(N);
    size_type first_in_node = N - heap_size;
    V                      += first_in_node;
    N                       = heap_size - first_in_node;
  }
  if (N > 1)
    std::sort(V, V + N);
  return ;
}

#endif /* BIHEAPIFY_INWARDS_SYMMETRICALLY_ALTERNATING_TOWARDS_ENDS_H_ */
