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

//Assumes that *first is the min value and *(first + (total_num_nodes - 1)) is the max value.
//Places ALL mins at the start and ALL maxs at the end with O(total_num_nodes) complexity.
//Returns a pair of integers (a, b) such that *(first + a) is the LAST element
// with minimal value and *(first + b) is the FIRST element of maximal value.
// Note that a < b if and only if the minimum and maximum values are distinct.
template<class RAI, typename size_type = std::size_t>
inline std::pair<size_type, size_type> EmplaceAllMinsAndMaxsAtEndsAssumingEndsAreMinAndMax(RAI first, size_type total_num_nodes) {
  size_type min_index = 0, max_index = total_num_nodes - 1;
  auto min_value = *(first + min_index);
  auto max_value = *(first + max_index);
  if (min_value == max_value)
    return std::pair<size_type, size_type>(max_index, 0);
  if (total_num_nodes < 3)
    return std::pair<size_type, size_type>(0, max_index);
  while (min_index + 1 < max_index && min_value == *(first + (min_index + 1)))
    min_index++;
  while (max_index > min_index + 1 && max_value == *(first + (max_index - 1)))
    max_index--;
  size_type i = min_index + 1;
  while (i < max_index) {
    auto current_value = *(first + i);
    if (*(first + i) == min_value) {
      min_index++;
      if (min_index < i) {
        std::iter_swap(first + min_index, first + i);
      } else {
        i++;
      }
    } else if (current_value == max_value) {
      max_index--;
      std::iter_swap(first + max_index, first + i);
    } else {
      i++;
    }
  }
  return std::pair<size_type, size_type>(min_index, max_index);
}


//Places ALL mins at the start and ALL maxs at the end with O(total_num_nodes) complexity.
//Returns a pair of integers (a, b) such that *(first + a) is the LAST element
// with minimal value and *(first + b) is the FIRST element of maximal value.
// Note that a < b if and only if the minimum and maximum values are distinct.
template<class RAI, typename size_type = std::size_t>
inline std::pair<size_type, size_type> EmplaceAllMinsAndMaxsAtEnds(RAI first, size_type N) {
  if (N < 2)
    return std::pair<size_type, size_type>(0, 0);
  BiHeapify(first, N);
  return EmplaceAllMinsAndMaxsAtEndsAssumingEndsAreMinAndMax(first, N);
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsNicerMath(RAI first, size_type N) {
  while (N > 9) {
    if (N < 30)
      std::sort(first, first + N);
    if (N < 10)
      return ;
    auto p = EmplaceAllMinsAndMaxsAtEnds(first, N);
    if (p.first >= p.second || p.first >= N / 2 || p.second <= (N - 1) / 2) //If we've emplaced the median.
      return ;
    if (N % 2 == 0 && (p.first == (N - 1) / 2 || p.second == N / 2)) {
      EmplaceAllMinsAndMaxsAtEnds(first + (p.first + 1), p.second - (p.first + 1));
      return ;
    }
    if (N % 3 == 2) {
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
