/*
 * biheapify_flip_ordered.h
 * 
 * This file implements BiHeapifyFlipOrdered() which forms a
 *  flip-ordered BiHeap, that is, a BiHeap with the additional
 *  property that for all i = 0, ..., CEIL(N/2), 
 *  *(V + i) <= *(V + Flip(i)).
 * Note that a call to BiHeapifyFlipOrdered() is slightly more 
 *  expensive a call to BiHeapify().
 * Note that RAI = Random Access Iterator
 * 
 *  Created on: Nov 21, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_FLIP_ORDERED_H_
#define BIHEAPIFY_FLIP_ORDERED_H_

#include <algorithm>

#include "biheapify.h"
#include "biheapify_lambda.h"

//#define FLIP(a) ((N - 1) - (a))

/*
 * ================== START: Definition of lambda version of BiHeapifyFlipOrdered ====================
 */

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrderedSiftFromMinToMax(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc,
                                      LambdaType lambda) {
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    auto pos_it         = V + lambda(N, pos_hc);

    {
      auto flip_pos_node = V + lambda(N, N_minus1 - pos_hc);
      if (*pos_it > *flip_pos_node)
        std::iter_swap(pos_it, flip_pos_node);
    }

    //assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    RAI smaller_it;
    if (is_right_child_valid && *right_it < *left_it) {
      smaller_it = right_it;
      pos_hc     = right_child_hc;
    } else {
      smaller_it = left_it;
      pos_hc     = left_child_hc;
    }
    if (*pos_it > *smaller_it)
      std::iter_swap(pos_it, smaller_it);
    else
      return ;
  }
  if (N % 3 != 2 || pos_hc != (N - 2) / 3) { //If the node is not the end of a double arrow.
    auto pos_it         = V + lambda(N, pos_hc);
    auto flip_of_pos_it = V + lambda(N, N_minus1 - pos_hc);
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMaxHeapMC<RAI, SizeType, LambdaType>(V, N, N_minus1, N_minus1 - pos_hc,
                                     N_minus1 - last_node_in_biheap_hc, lambda);
  return ;
}

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrderedSiftFromMaxToMin(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc = N_minus1 - pos_mc;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = N_minus1 - left_child_mc;//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = V + lambda(N, pos_hc);
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    {
      auto flip_pos_node = V + lambda(N, pos_mc);
      if (*pos_it < *flip_pos_node)
        std::iter_swap(pos_it, flip_pos_node);
    }

    //assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
    bool is_right_child_valid = right_child_mc < heap_size;
    RAI larger_it;
    if (is_right_child_valid && *right_it > *left_it) {
      larger_it = right_it;
      pos_hc    = right_child_hc;
      pos_mc    = right_child_mc;
    } else {
      larger_it = left_it;
      pos_hc    = left_child_hc;
      pos_mc    = left_child_mc;
    }
    if (*pos_it < *larger_it)
      std::iter_swap(pos_it, larger_it);
    else
      return ;
  }
  if (N % 3 != 2 || pos_hc != (N - 2) / 3) { //If the node is not the end of a double arrow.
    auto pos_it         = V + lambda(N, pos_hc);
    auto flip_of_pos_it = V + lambda(N, pos_mc);
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMinHeapHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, first_node_in_biheap_hc, lambda);
  return ;
}

/* This will BiHeapify all nodes in [0, N).
 * Assumes that N is odd.
 */
/*
 * Remark:
 *  This algorithm has complexity O(N). To see why, recall the
 *   argument showing that the heapify operation has O(n) complexity (e.g. as
 *   found on pp. 115 - 116 of "The Algorithm Design Manual" 2nd edition); this
 *   argument generalizes to prove that this algorithm also runs in O(n) times.
 *  The key observation is that the biheap is constructed by adding one node
 *   at a time with this node alternating between a node in a min heap and a
 *   node in the max heap. The complexity of this biheapification is easily
 *   seen to be twice the complexity of the above mentioned heapification
 *   operation plus a constant.
 */
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyFlipOrdered(RAI V, SizeType N, LambdaType lambda) {
  if(N < 2)
    return ;
  SizeType heap_size     = HeapSize(N);
  SizeType first_in_node = N - heap_size;
  SizeType N_minus1      = N - 1;
  {
    SizeType one_past_last_node_to_order = N / 2;
    for (SizeType i = first_in_node; i < one_past_last_node_to_order; i++) {
      auto i_node      = V + lambda(N, i);
      auto flip_i_node = V + lambda(N, (N - 1) - i);
      if (*i_node > *flip_i_node)
        std::iter_swap(i_node, flip_i_node);
    }
  }
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  SizeType last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                     - (N % 3 == 2);
  SizeType first_node_in_biheap_hc = N_minus1 - last_node_in_biheap_hc;
  while (first_node_in_biheap_hc > 0) {
    --first_node_in_biheap_hc;
    {
      auto lambda_first_hc = lambda(N, first_node_in_biheap_hc);
      auto lambda_last_hc  = lambda(N, last_node_in_biheap_hc + 1);
      if (*(V + lambda_first_hc) > *(V + lambda_last_hc))
        std::iter_swap(V + lambda_first_hc, V + lambda_last_hc);
    }
    BiHeapifyFlipsOrderedSiftFromMinToMax<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, first_node_in_biheap_hc,
        last_node_in_biheap_hc, lambda);
    BiHeapifyFlipsOrderedSiftFromMaxToMin<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, N_minus1 - (++last_node_in_biheap_hc),
        first_node_in_biheap_hc, lambda);
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyFlipOrdered(RAI V, RAI one_past_last, LambdaType lambda) {
  BiHeapifyFlipOrdered<RAI, SizeType, LambdaType>(V, std::distance(V, one_past_last), lambda);
}

/*
 * ================== END: Definition of lambda version of BiHeapifyFlipOrdered ====================
 */

/*
 * ================== START: Specializations of lambda version of BiHeapifyFlipOrdered ====================
 */


template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifyFlipOrdered(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  BiHeapifyFlipOrdered<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
  return ;
}

/*
 * ================== END: Specializations of lambda version of BiHeapifyFlipOrdered ====================
 */

/*
 * ================== START: Check if nodes form a Flip Ordered BiHeap ====================
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool ArePairwiseFlipOrdered(RAI V, SizeType N, LambdaType lambda) {
  if (N < 2)
    return true;
  SizeType num_nodes_to_check = N / 2;
  SizeType N_minus1           = N - 1;
  for (SizeType i = 0; i < num_nodes_to_check; i++) {
    if (*(V + lambda(N, i)) > *(V + lambda(N, (N_minus1 - i))))
      return false;
  }
  return true;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsBiHeapWithFlipsOrdered(RAI V, SizeType N, LambdaType lambda) {
  if (N <= 1)
    return true;
  if (!IsBiHeap<RAI, SizeType, LambdaType>(V, N, lambda))
    return false;
  return ArePairwiseFlipOrdered<RAI, SizeType, LambdaType>(V, N, lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool ArePairwiseFlipOrdered(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  return ArePairwiseFlipOrdered<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool IsBiHeapWithFlipsOrdered(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  return IsBiHeapWithFlipsOrdered<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
}

/*
 * ================== END: Check if nodes form a Flip Ordered BiHeap ====================
 */
//#undef FLIP

#endif /* BIHEAPIFY_FLIP_ORDERED_H_ */
