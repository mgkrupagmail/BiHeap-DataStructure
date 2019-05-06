/*
 * biheap_select.h
 *
 *  The main function in this file is BiHeapSelect(), which does the same thing as
 *   the Quickselect algorithm. However, unlike the Quickselect algorithm,
 *   BiHeapSelect() is O(N) (this is something that I'm still mathematically proving).
 *  Empirical testing by the functions in biheap_select_measure_complexity.h show that
 *   BiHeapSelect() emplaces the element in the desired position using at most 14 * N
 *   swap operations.
 *
 *  Created on: May 2, 2019
 *      Author: Matthew Gregory Krupa
 *      Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAP_SELECT_H_
#define BIHEAP_SELECT_H_

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <string>
#include <vector>

//#include "biheap_common.h"
#include "biheap_ostream.h"
#include "biheapify.h"


//Returns a pair (a, b) where a is the index that is 1 + (the index of last element < pivot_value)
// and b is the index of the first element that is > pivot_value.
//Note for any index i, *(first + i) == pivot value if and only if a <= i < b. Also, it's guaranteed that a <= b.
template<typename RAI, typename size_type, typename ValueType>
inline std::pair<size_type, size_type> DutchPartition(RAI first, size_type start, size_type one_past_last, ValueType pivot_value) {
  size_type e = start;  //= first_equal   = The first element in the string of elements equal to pivot_val
  size_type u = start;  //= first_unknown = The current unknown element
  size_type one_before_first_largest = one_past_last - 1; //= one_before_first_largest = The element directly before the string of elements > pivot_value

  while(u <= one_before_first_largest) {
    auto current_value = *(first + u);
    if(current_value < pivot_value) {
      std::iter_swap(first + u, first + e);
      u++;
      e++;
    } else if(current_value > pivot_value) {
      std::iter_swap(first + u, first + one_before_first_largest);
      one_before_first_largest--;
    } else { //if(arr[first_unknown] == pivot_val)
      u++;    //Move to the next unknown element
    }
  }
  return std::pair<size_type, size_type>(e, one_before_first_largest + 1);
}

//Returns a pair (a, b) where
// a is the index that is 1 + (the index of last element that is < pivot_value), and
// b is the index of the first element that is > pivot_value.
//Note for any index i, *(first + i) == pivot value if and only if a <= i < b. Also, it's guaranteed that a <= b.
template<typename RAI, typename size_type>
inline std::pair<size_type, size_type> DutchPartition(RAI first, size_type one_past_last, size_type pivot_index) {
  if (one_past_last <= 1)
    return std::pair<size_type, size_type>(0, 0);
  auto pivot_value = *(first + pivot_index);
  std::iter_swap(first, first + pivot_index);
  size_type e = 0;  //= first_equal   = The first element in the string of elements equal to pivot_val
  size_type u = 1;  //= first_unknown = The current unknown element
  size_type one_before_first_largest = one_past_last - 1; //= one_before_first_largest = The element directly before the string of elements > pivot_value

  while(u <= one_before_first_largest) {
    auto current_value = *(first + u);
    if(current_value < pivot_value) {
      std::iter_swap(first + u, first + e);
      u++;
      e++;
    } else if(current_value > pivot_value) {
      std::iter_swap(first + u, first + one_before_first_largest);
      one_before_first_largest--;
    } else {
      u++;
    }
  }
  return std::pair<size_type, size_type>(e, one_before_first_largest + 1);
}

//Places the minimum at the very left, and the maximum at the very right.
//You could also replace this entire function with the single call:
// BiHeapify(first, total_num_nodes);
template<class RAI, typename size_type = std::size_t>
inline void EmplaceMinAndMax(RAI first, size_type total_num_nodes) {
  if (total_num_nodes <= 2) {
    if (total_num_nodes == 2 && *first > *(first + 1))
      std::iter_swap(first, first + 1);
    return ;
  }
  size_type min_index = 0, max_index = 0;
  auto min_value = *first, max_value = *first;
  for (auto i = 0u; i < total_num_nodes; i++) {
    auto current_value = *(first + i);
    if (current_value < min_value) {
      min_value = current_value;
      min_index = i;
    } else if (current_value > max_value) {
      max_value = current_value;
      max_index = i;
    }
  }
  if (min_index == max_index) { //Then the sequence is constant.
assert(min_index == 0);
    return ;
  }
  if (min_index == 0) {
    std::iter_swap(first + max_index, first + (total_num_nodes - 1));
    return ;
  }
  if (max_index == total_num_nodes - 1) {
    std::iter_swap(first, first + min_index);
    return ;
  }
  if (max_index == 0 && min_index == total_num_nodes - 1) {
    std::iter_swap(first, first + min_index);
    return ;
  }
assert(min_value != max_value);
assert(min_index != max_index);
  if (min_index != total_num_nodes - 1) {
    std::iter_swap(first + max_index, first + (total_num_nodes - 1));
    std::iter_swap(first, first + min_index);
  } else {
    std::iter_swap(first, first + (total_num_nodes - 1));
    std::iter_swap(first + max_index, first + (total_num_nodes - 1));
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyWithAllMinsAndMaxsAtEnds(RAI first, size_type N) {
  BiHeapify<RAI, size_type>(first, N); //BiHeapify the whole thing
  EmplaceAllMinsAndMaxsAtEndsAssumingEndsAreMinAndMax(first, N);
  BiHeapify<RAI, size_type>(first, N);
  return ;
}




//Returns a pair (start_index, end_index) of the interval whose median should
// be found.
//The interval that is returned is the largest progenic interval such that
// its start index is <= desired_index.
// N              - The BiHeap size
template<typename size_type = std::size_t>
std::pair<size_type, size_type> BiHeapSelectGetRangeToFindMedianOf(size_type N, size_type desired_index) {
  bool is_desired_index_flipped = desired_index >= (N / 2);
  //The following transformation allows us to assume WLOG that desired_index < N / 2
  // so that in particular, it is in the min heap.
  if (is_desired_index_flipped) {
    desired_index = (N - 1) - desired_index;
  }
  size_type heap_size = HeapSize<size_type>(N);
  assert(desired_index < heap_size);

  size_type start_index = heap_size / 2;
  size_type end_index = heap_size - 1;

  //start_index = Parent<size_type>(start_index);
  //end_index = 2 * start_index;

  while (start_index > desired_index) {
    start_index = Parent<size_type>(start_index);
    end_index = 2 * start_index;
  }

  if (is_desired_index_flipped) {
    size_type new_end_index = (N - 1) - start_index;
    start_index = (N - 1) - end_index;
    end_index = new_end_index;
  }
  return std::pair<size_type, size_type>(start_index, end_index);
}

#define BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_THIS 12

//Emplaces the element V + desired index.
// That is, it rearranges V, V + 1, ..., V+ (N - 1) so that
// (1) *(V + i) <= *(V + desired_index) for all 0 <= i < desired_index, and
// (2) *(V + i) >= *(V + desired_index) for all desired_index < i < N.
template<class RAI, typename size_type = std::size_t>
void BiHeapSelect(RAI V, size_type N, size_type desired_index) {
  assert(desired_index < N || !(std::cout << "desired_index = " << desired_index << " \tN = " << N << std::endl));
  while(N > 0) {
    if (N <= BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_THIS) {
      std::sort<RAI>(V, V + N);
      break ;
    }
    assert(desired_index >= 0);
    assert(desired_index < N || !(std::cout << "desired_index = " << desired_index << " \tN = " << N << std::endl));

    BiHeapify<RAI, size_type>(V, N);

    std::pair<size_type, size_type> next_range_pair = BiHeapSelectGetRangeToFindMedianOf<size_type>(N, desired_index);
    size_type next_range_start = next_range_pair.first;
    size_type next_range_end = next_range_pair.second;
    assert(next_range_end >= next_range_start);
    size_type next_range_length = (next_range_end + 1) - next_range_start;
    size_type middle_of_next_range = next_range_start + (next_range_length / 2);
    //std::cout << "next_range_start = " << next_range_start << " \tnext_range_end = " << next_range_end << " \tnext_range_length = " << next_range_length << " \tmiddle_of_next_range = " << middle_of_next_range << std::endl;
    assert(next_range_start <= middle_of_next_range);
    assert(middle_of_next_range <= next_range_end);

    //Find the median of the range.
    BiHeapSelect<RAI, size_type>(V + next_range_start, next_range_length, middle_of_next_range - next_range_start);

    //Dutch flag partition all N elements about the median of the range that was found by the above call to BiHeapSelect().
    std::pair<size_type, size_type> pivot_range = DutchPartition<RAI, size_type>(V, N, middle_of_next_range);
    auto a = pivot_range.first;
    auto b = pivot_range.second;
    // a is the index that is 1 + (the index of last element that is < the pivot value), and
    // b is the index of the first element that is > the pivot value.
    if (desired_index < a) {
      N = a;
      assert(desired_index < N);
    } else if (desired_index >= b) {
      V = V + b;
      N -= b;
      desired_index -= b;
      assert(desired_index < N);
    } else { //a <= diesred_index < b so we've emplaced the desired index.
      break ;
    }
  }
  return ;
}











#endif /* BIHEAP_SELECT_H_ */
