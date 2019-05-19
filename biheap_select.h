/*
 * biheap_select.h
 *
 *  The main function in this file is BiHeapSelect(), which does the same thing as
 *   the Quickselect algorithm. However, unlike the Quickselect algorithm,
 *   BiHeapSelect() is O(N) (this is something that I'm still mathematically proving).
 *  Empirical testing by the functions in biheap_select_measure_complexity.h shows that
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
    } else if(pivot_value < current_value) {
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
inline std::pair<size_type, size_type> DutchPartition(RAI V, size_type N, size_type pivot_index) {
  if (N <= 1)
    return std::pair<size_type, size_type>(0, 0);
  auto pivot_value = *(V + pivot_index);
  std::iter_swap(V, V + pivot_index);
  size_type e = 0;  //= first_equal   = The first element in the string of elements equal to pivot_val
  size_type u = 1;  //= first_unknown = The current unknown element
  size_type one_before_first_largest = N - 1; //= one_before_first_largest = The element directly before the string of elements > pivot_value

  while(u <= one_before_first_largest) {
    auto current_value = *(V + u);
    if(current_value < pivot_value) {
      std::iter_swap(V + u, V + e);
      u++;
      e++;
    } else if(pivot_value < current_value) {
      std::iter_swap(V + u, V + one_before_first_largest);
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
inline void EmplaceMinAndMax(RAI V, size_type N) {
  if (N <= 2) {
    if (N == 2 && *(V + 1) < *V)
      std::iter_swap(V, V + 1);
    return ;
  }
  size_type min_index = 0, max_index = 0;
  auto min_value = *V, max_value = *V;
  for (auto i = 0u; i < N; i++) {
    auto current_value = *(V + i);
    if (current_value < min_value) {
      min_value = current_value;
      min_index = i;
    } else if (max_value < current_value) {
      max_value = current_value;
      max_index = i;
    }
  }
  if (min_index == max_index) { //Then the sequence is constant.
assert(min_index == 0);
    return ;
  }
  if (min_index == 0) {
    std::iter_swap(V + max_index, V + (N - 1));
    return ;
  }
  if (max_index == N - 1) {
    std::iter_swap(V, V + min_index);
    return ;
  }
  if (max_index == 0 && min_index == N - 1) {
    std::iter_swap(V, V + min_index);
    return ;
  }
assert(min_value != max_value);
assert(min_index != max_index);
  if (min_index != N - 1) {
    std::iter_swap(V + max_index, V + (N - 1));
    std::iter_swap(V, V + min_index);
  } else {
    std::iter_swap(V, V + (N - 1));
    std::iter_swap(V + max_index, V + (N - 1));
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyWithAllMinsAndMaxsAtEnds(RAI V, size_type N) {
  BiHeapify<RAI, size_type>(V, N); //BiHeapify the whole thing
  EmplaceAllMinsAndMaxsAtEndsAssumingEndsAreMinAndMax(V, N);
  BiHeapify<RAI, size_type>(V, N);
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
  //The following transformation in the if statement's body allows us to assume WLOG
  // that k < N / 2 so that in particular, node k will be in the min heap.
  if (is_desired_index_flipped) {
    desired_index = (N - 1) - desired_index;
  }
  size_type heap_size = HeapSize<size_type>(N);
  assert(desired_index < heap_size);

  size_type s = heap_size / 2;
  size_type e = heap_size - 1;

  //start_index = Parent<size_type>(start_index);
  //end_index = 2 * start_index;

  while (s > desired_index) {
    s = Parent<size_type>(s);
    e = 2 * s;
  }

  if (is_desired_index_flipped) {
    size_type new_e = (N - 1) - s;
    s = (N - 1) - e;
    e = new_e;
  }
  return std::pair<size_type, size_type>(s, e);
}

#define BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_OR_EQUAL_TO_THIS 12

//Emplaces the element V + desired index.
// That is, it rearranges V, V + 1, ..., V+ (N - 1) so that
// (1) *(V + i) <= *(V + desired_index) for all 0 <= i < desired_index, and
// (2) *(V + i) >= *(V + desired_index) for all desired_index < i < N.
template<class RAI, typename size_type = std::size_t>
void BiHeapSelect(RAI V, size_type N, size_type k) {
  assert(k < N || !(std::cout << "k = " << k << " \tN = " << N << std::endl));
  while(N > 0) {
    if (N <= BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_OR_EQUAL_TO_THIS) {
      std::sort<RAI>(V, V + N);
      break ;
    }
    assert(k >= 0);
    assert(k < N || !(std::cout << "k = " << k << " \tN = " << N << std::endl));

    BiHeapify<RAI, size_type>(V, N);

    //Handle some trivial cases.
    if (k == 0 || k == N - 1)
      break ; //Then we're done
    if (k == 1 || k == N - 2) {
      if (*(V + 2) < *(V + 1)) //Sort nodes V + 1 and V + 2
          std::iter_swap(V + 1, V + 2);
      if (*(V + (N - 2)) < *(V + (N - 3))) //Sort nodes V + (N - 2) and V + (N - 3)
          std::iter_swap(V + (N - 2), V + (N - 3));
      break ; //We're now done
    }

    std::pair<size_type, size_type> next_range_pair = BiHeapSelectGetRangeToFindMedianOf<size_type>(N, k);
    size_type s = next_range_pair.first;
    size_type e = next_range_pair.second;
    assert(e >= s);
    size_type length = (e + 1) - s;
    size_type middle = s + (length / 2);
    //std::cout << "next_range_start = " << next_range_start << " \tnext_range_end = " << next_range_end << " \tnext_range_length = " << next_range_length << " \tmiddle_of_next_range = " << middle_of_next_range << std::endl;
    assert(s <= middle);
    assert(middle <= e);

    //Find the median of the range.
    BiHeapSelect<RAI, size_type>(V + s, length, middle - s);

    //Dutch flag partition all N elements about the median of the range that was found by the above call to BiHeapSelect().
    std::pair<size_type, size_type> pivot_range = DutchPartition<RAI, size_type>(V, N, middle);
    auto a = pivot_range.first;
    auto b = pivot_range.second;
    // a is the index that is 1 + (the index of last element that is < the pivot value), and
    // b is the index of the first element that is > the pivot value.
    if (k < a) {
      N = a;
      assert(k < N);
    } else if (k >= b) {
      V = V + b;
      N -= b;
      k -= b;
      assert(k < N);
    } else { //a <= desired_index < b so we've emplaced the desired index.
      break ;
    }
  }
  return ;
}











#endif /* BIHEAP_SELECT_H_ */
