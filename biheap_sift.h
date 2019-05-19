/*
 * biheap_sift.h
 *
 *  Created on: Jul 25, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 *
 * This header file defines the BiHeapSift() function. Given a biheap on
 *  N nodes defined by the iterator V (so that *V,
 *  *(V + 1), ..., *(V + (N - 1) are the nodes' values),
 *  and given 0 <= pos_hc < N, if *(V + pos_hc) is
 *  changed then these nodes may no longer form a biheap.
 * Calling BiHeapSift(V, N, pos_hc) will make it into a
 *  biheap once again and it will do this in O(log(N)) time.
 */

#ifndef BIHEAP_SIFT_H_
#define BIHEAP_SIFT_H_

#include <algorithm>

#include "biheapify.h"

#define FLIP(a) (N_minus1 - (a))

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftUpMaxHeapUnboundedMC(RAI V, SizeType N, SizeType N_minus1,
                            SizeType pos_mc, LambdaType lambda) {
  SizeType parent;
  auto pos_it = V + lambda(N, N_minus1 - pos_mc);
  auto pos_value = *pos_it;
  while (pos_mc > 0) {
    parent = ParentNotRoot<SizeType>(pos_mc);
    auto parent_it = V + lambda(N, N_minus1 - parent);
    if (*parent_it < pos_value) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftUpMinHeapUnboundedHC(RAI V, SizeType N, SizeType N_minus1,
                            SizeType pos_hc, LambdaType lambda) {
  SizeType parent;
  auto pos_it = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_hc > 0) {
    parent = ParentNotRoot<SizeType>(pos_hc);
    auto parent_it = V + lambda(N, parent);
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftFromMinToMaxUnbounded(RAI V, SizeType N, SizeType N_minus1,
                             SizeType heap_size, SizeType first_in_node,
                             SizeType pos_hc, LambdaType lambda) {
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_child_it  = V + lambda(N, left_child_hc);
    auto right_child_it = V + lambda(N, right_child_hc);

    //assert((left_child  < heap_size) && (left_child  < N) && (right_child < N));
    bool is_right_child_valid = right_child_hc < heap_size;
    RAI smaller_it;
    decltype(pos_value) smaller_value = *left_child_it, right_value;
    if (is_right_child_valid && (right_value = *right_child_it) < smaller_value) {
      smaller_value = right_value;
      smaller_it    = right_child_it;
      pos_hc        = right_child_hc;
    } else {
      smaller_it    = left_child_it;
      pos_hc        = left_child_hc;
    }
    if (smaller_value < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  }
  SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, N_minus1 - pos_hc, lambda);
  return ;
}

/* Assumes that N is odd, that the node pos_mc belongs to the
 *  max heap, and that FLIP(pos_mc) >= smallest_node_in_biheap_hc.
 */
template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftFromMaxToMinUnbounded(RAI V, SizeType N, SizeType N_minus1,
                              SizeType heap_size, SizeType first_in_node,
                              SizeType pos_mc, LambdaType lambda) {
  auto pos_hc    = N_minus1 - pos_mc;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = N_minus1 - left_child_mc;//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1;  //= FLIP(right_child_mc)
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    //assert((left_child_mc < heap_size) && (left_child_hc >= 0) && (right_child_hc >= 0));
    bool is_right_child_valid = right_child_mc < heap_size;
    decltype(pos_value) larger_value = *left_it, right_value;
    RAI larger_it;
    if (is_right_child_valid && larger_value < (right_value = *right_it)) {
      larger_value = right_value;
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else {
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < larger_value){
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    }
    else
      return ;
  }
  SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void BiHeapSift(RAI V, SizeType N, SizeType pos_hc, LambdaType lambda) {
  SizeType N_minus1       = N - 1;
  SizeType heap_size      = HeapSize<SizeType>(N);
  SizeType first_in_node  = N - heap_size;
  SizeType pos_mc         = N_minus1 - pos_hc;
  auto pos_value          = *(V + lambda(N, pos_hc));
  bool is_node_in_min_heap = pos_hc < heap_size;
  bool is_node_in_max_heap = pos_mc < heap_size;
  if (is_node_in_min_heap && (pos_hc == 0 ||
      *(V + lambda(N, ParentNotRoot(pos_hc))) <= pos_value))
    SiftFromMinToMaxUnbounded<RAI, SizeType, LambdaType>(V, N, N_minus1,
                               heap_size, first_in_node, pos_hc, lambda);
  else if (is_node_in_max_heap && (pos_mc == 0 ||
          *(V + lambda(N, N_minus1 - ParentNotRoot(pos_mc))) >= pos_value))
    SiftFromMaxToMinUnbounded<RAI, SizeType, LambdaType>(V, N, N_minus1,
                               heap_size, first_in_node, pos_mc, lambda);
  //At this point pos_hc != 0,  pos_mc != 0, and at least one
  // of the following is true:
  // (1) is_node_in_min_heap && *(V + Parent(pos_hc))       > pos_value
  // (2) is_node_in_max_heap && *(V + FLIP(Parent(pos_mc))) < pos_value
  else if (is_node_in_min_heap && pos_value < *(V + lambda(N, Parent(pos_hc))))
    SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, lambda);
  else
    SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_mc, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void BiHeapSiftMC(RAI V, SizeType N, SizeType pos_mc, LambdaType lambda) {
  BiHeapSift<RAI, SizeType, LambdaType>(V, N, (N - 1) - (pos_mc), lambda);
}

/*
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapSift(RAI V, RAI one_past_last, RAI pos_hc) {
  BiHeapSift<RAI, SizeType>(V, std::distance(V, one_past_last),
                             std::distance(V, pos_hc));
}
*/

template<class RAI, typename SizeType = std::size_t>
inline void BiHeapSift(RAI V, SizeType N, SizeType pos_hc) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  BiHeapSift<RAI, SizeType, decltype(trivial_lambda)>(V, N, pos_hc,
                                                      trivial_lambda);
  return ;
}

//Assumes that pos_hc is an In node.
template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void BiHeapSiftInNode(RAI V, SizeType N, SizeType pos_hc, LambdaType lambda) {
  SizeType N_minus1              = N - 1;
  SizeType minh_parent_of_pos_hc = Parent<SizeType>(pos_hc);
  auto pos_it                    = V + lambda(N, pos_hc);
  auto minh_parent_of_pos_it     = V + lambda(N, minh_parent_of_pos_hc);
  auto pos_value                 = *pos_it;
  if (pos_value < *minh_parent_of_pos_it) {
    SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, lambda);
  } else {
    SizeType pos_mc                = N_minus1 - pos_hc;
    SizeType maxh_parent_of_pos_hc = N_minus1 - Parent<SizeType>(pos_mc);
    auto maxh_parent_of_pos_it     = V + lambda(N, maxh_parent_of_pos_hc);
    if (*maxh_parent_of_pos_it < pos_value)
      SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_mc, lambda);
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t>
inline void BiHeapSiftInNode(RAI V, SizeType N, SizeType pos_hc) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  BiHeapSiftInNode<RAI, SizeType, decltype(trivial_lambda)>(V, N, pos_hc,
                                                      trivial_lambda);
  return ;
}

#undef FLIP

#endif /* BIHEAP_SIFT_H_ */
