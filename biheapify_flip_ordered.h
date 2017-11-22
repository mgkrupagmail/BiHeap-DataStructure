/*
 * biheapify_flip_ordered.h
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

#define FLIP(a) (N - 1 - (a))

template<class RAI, typename size_type = std::size_t>
bool ArePairwiseFlipOrdered(RAI V, size_type N) {
  if (N < 2)
    return true;
  size_type num_nodes_to_check = N / 2;
  for (size_type i = 0; i < num_nodes_to_check; i++) {
    if (*(V + i) > *(V + ((N - 1) - i)))
      return false;
  }
  return true;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool ArePairwiseFlipOrdered(RAI V, size_type N, LambdaType lambda) {
  if (N < 2)
    return true;
  size_type num_nodes_to_check = N / 2;
  for (size_type i = 0; i < num_nodes_to_check; i++) {
    if (*(V + lambda(N, i)) > *(V + lambda(N, ((N - 1) - i))))
      return false;
  }
  return true;
}

template<class RAI, typename size_type = std::size_t>
bool IsBiHeapWithFlipsOrdered(RAI V, size_type N) {
  if (N <= 1)
    return true;
  if (!IsBiHeap<RAI, size_type>(V, N))
    return false;
  return ArePairwiseFlipOrdered<RAI, size_type>(V, N);
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool IsBiHeapWithFlipsOrdered(RAI V, size_type N, LambdaType lambda) {
  if (N <= 1)
    return true;
  if (!IsBiHeap<RAI, size_type, LambdaType>(V, N, lambda))
    return false;
  return ArePairwiseFlipOrdered<RAI, size_type, LambdaType>(V, N, lambda);
}

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyFlipsOrderedSiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_hc,
                                      size_type last_node_in_biheap_hc) {
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;
    auto pos_it         = V + pos_hc;

    if (*pos_it > *(V + FLIP(pos_hc)))
      std::iter_swap(pos_it, V + FLIP(pos_hc));

    assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
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
  if (N % 3 != 2 || pos_hc != (N - 2) / 3) {
    auto pos_it         = V + pos_hc;
    auto flip_of_pos_it = V + FLIP(pos_hc);
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMaxHeapMC<RAI, size_type>(V, N, FLIP(pos_hc),
                                  FLIP(last_node_in_biheap_hc));
  return ;
}

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyFlipsOrderedSiftFromMaxToMin(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_mc,
                                      size_type first_node_in_biheap_hc) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = V + pos_hc;
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;

    if (*pos_it < *(V + pos_mc))
      std::iter_swap(pos_it, V + pos_mc);

    assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
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
  if (N % 3 != 2 || pos_mc != (N - 2) / 3) {
    auto pos_it         = V + pos_hc;
    auto flip_of_pos_it = V + pos_mc;
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMinHeapHC<RAI, size_type>(V, pos_hc, first_node_in_biheap_hc);
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
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyFlipsOrdered(RAI V, size_type N) {
  if(N < 2)
    return ;
  size_type heap_size                  = HeapSize(N);
  size_type first_node_in_mirror_heap  = N - heap_size;
  {
    size_type last_node_to_order = N / 2;
    for (size_type i = first_node_in_mirror_heap; i < last_node_to_order; i++) {
      size_type flip_i = (N - 1) - i;
      if (*(V + i) > *(V + flip_i))
        std::iter_swap(V + i, V + flip_i);
    }
  }
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  size_type last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                         - (N % 3 == 2);
  size_type first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    --first_node_in_biheap_hc;
    assert(last_node_in_biheap_hc + 1 == ((N - 1) - first_node_in_biheap_hc));
    if (*(V + first_node_in_biheap_hc) > *(V + (last_node_in_biheap_hc + 1)))
      std::iter_swap(V + first_node_in_biheap_hc, V + (last_node_in_biheap_hc + 1));

    BiHeapifyFlipsOrderedSiftFromMinToMax<RAI, size_type>(V, N, heap_size,
        first_node_in_mirror_heap,
        first_node_in_biheap_hc, last_node_in_biheap_hc);
    BiHeapifyFlipsOrderedSiftFromMaxToMin<RAI, size_type>(V, N, heap_size,
        first_node_in_mirror_heap,
        FLIP(++last_node_in_biheap_hc), first_node_in_biheap_hc);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyFlipsOrdered(RAI V, RAI one_past_last) {
  BiHeapifyFlipsOrdered<RAI, size_type>(V, std::distance(V, one_past_last));
}


//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrderedSiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_hc,
                                      size_type last_node_in_biheap_hc,
                                      LambdaType lambda) {
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    auto pos_it         = V + lambda(N, pos_hc);

    if (*pos_it > *(V + lambda(N, FLIP(pos_hc))))
      std::iter_swap(pos_it, V + lambda(N, FLIP(pos_hc)));

    assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
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
  if (N % 3 != 2 || pos_hc != (N - 2) / 3) {
    auto pos_it         = V + lambda(N, pos_hc);
    auto flip_of_pos_it = V + lambda(N, FLIP(pos_hc));
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMaxHeapMC<RAI, size_type, LambdaType>(V, N, FLIP(pos_hc),
                                  FLIP(last_node_in_biheap_hc), lambda);
  return ;
}

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrderedSiftFromMaxToMin(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_mc,
                                      size_type first_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = V + lambda(N, pos_hc);
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);

    if (*pos_it < *(V + lambda(N, pos_mc)))
      std::iter_swap(pos_it, V + lambda(N, pos_mc));

    assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
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
  if (N % 3 != 2 || pos_hc != (N - 2) / 3) {
    auto pos_it         = V + lambda(N, pos_hc);
    auto flip_of_pos_it = V + lambda(N, pos_mc);
    if ((pos_hc <  N / 2 && *pos_it > *flip_of_pos_it) ||
        (pos_hc >= N / 2 && *pos_it < *flip_of_pos_it))
      std::iter_swap(pos_it, flip_of_pos_it);
  }
  SiftUpMinHeapHC<RAI, size_type, LambdaType>(V, N, pos_hc, first_node_in_biheap_hc, lambda);
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
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrdered(RAI V, size_type N, LambdaType lambda) {
  if(N < 2)
    return ;
  size_type heap_size                  = HeapSize(N);
  size_type first_node_in_mirror_heap  = N - heap_size;
  {
    size_type last_node_to_order = N / 2;
    for (size_type i = first_node_in_mirror_heap; i < last_node_to_order; i++) {
      size_type flip_i = (N - 1) - i;
      if (*(V + lambda(N, i)) > *(V + lambda(N, flip_i)))
        std::iter_swap(V + lambda(N, i), V + lambda(N, flip_i));
    }
  }
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  size_type last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                         - (N % 3 == 2);
  size_type first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    --first_node_in_biheap_hc;
    if (*(V + first_node_in_biheap_hc) > *(V + (last_node_in_biheap_hc + 1)))
      std::iter_swap(V + lambda(N, first_node_in_biheap_hc), V + lambda(N, last_node_in_biheap_hc + 1));
    BiHeapifyFlipsOrderedSiftFromMinToMax<RAI, size_type, LambdaType>(V, N, heap_size,
        first_node_in_mirror_heap,
        first_node_in_biheap_hc, last_node_in_biheap_hc, lambda);
    BiHeapifyFlipsOrderedSiftFromMaxToMin<RAI, size_type, LambdaType>(V, N, heap_size,
        first_node_in_mirror_heap,
        FLIP(++last_node_in_biheap_hc), first_node_in_biheap_hc, lambda);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyFlipsOrdered(RAI V, RAI one_past_last, LambdaType lambda) {
  BiHeapifyFlipsOrdered<RAI, size_type, LambdaType>(V, std::distance(V, one_past_last), lambda);
}

#undef FLIP

#endif /* BIHEAPIFY_FLIP_ORDERED_H_ */
