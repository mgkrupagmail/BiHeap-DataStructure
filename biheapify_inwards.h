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
inline void BiHeapifyInwards(RAI first, size_type N) {
  while (N > 9) {
    BiHeapify(first, N);
    size_type heap_size              = HeapSize(N);
    size_type first_extended_in_node = N - heap_size;//= FLIP(heap_size - 1);
    first                           += first_extended_in_node;
    N                                = heap_size - first_extended_in_node;
  }
  if (N > 1)
    std::sort(first, first + N);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsNicerMath(RAI first, size_type N) {
  while (N > 9) {
    if (N % 3 == 1) {
      BiHeapify(first, N);
      first = first + 1;
      N = N - 2;
    }
    if (N % 3 == 2) {
      BiHeapify(first, N);
      first = first + 1;
      N = N - 2;
    }
    BiHeapify(first, N);
    size_type heap_size              = HeapSize(N);
    size_type first_extended_in_node = N - heap_size;
    first                           += first_extended_in_node;
    N                                = heap_size - first_extended_in_node;
  }
  if (N > 1)
    std::sort(first, first + N);
  return ;
}

#endif /* BIHEAPIFY_INWARDS_H_ */
