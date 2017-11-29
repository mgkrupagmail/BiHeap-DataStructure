/*
 * almost_biheapify.h
 *
 * The AlmostBiHeapify() function performs the same operation as
 *  BiHeapify() except that it "skips over" all In nodes (which
 *  are those nodes that belong to both the min heap and the max heap).
 * The basic idea is this:
 * We construct a new graph from the BiHeap graph on N nodes (see
 *  "BiHeaps and Pivot Selection.pdf" for details on BiHeap graphs).
 * Note that every In node in the BiHeap graph is incident to exactly
 *  two edges. For every In node, "fuse" these two edges
 *  into a single edge and remove the In node from the graph.
 * Having done this to every In node, we now have a new graph; call it
 *  the almost BiHeap graph.
 * The AlmostBiHeapify() operation rearranges the values in this
 *  almost BiHeap graph so that if two node, say with min heap coordinates
 *  i and j where i < j, are incident to the same edge then the value
 *  of node i is <= the value of node j.
 * Note that at the end of the AlmostBiHeapify() operation, the values
 *  of the BiHeap's In nodes are unchanged.
 * To make an almost BiHeap into a BiHeap, we can call the BiHeapify()
 *  operation, the effect of which is equivalent to performing
 *  NOTHING BUT a sequence of sift up operations (with no sift downs).
 * This allows one to better track what happens to the values of the
 *  In nodes, which may be useful for forming inequalities.
 *
 * This algorithm was first discovered on Oct 24, 2017.
 *
 *  Created on: Nov 28, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef ALMOST_BIHEAPIFY_H_
#define ALMOST_BIHEAPIFY_H_

#include "biheapify.h"

#include <algorithm>

#define FLIP(a) ((N - 1) - (a))

template<typename size_type = std::size_t>
inline size_type IsInNode(size_type pos_hc, size_type pos_mc, size_type heap_size) {
  return (pos_hc < heap_size) && (pos_mc < heap_size);
}

template<typename size_type = std::size_t>
inline size_type IsInNodeHC(size_type N, size_type pos_hc, size_type heap_size) {
  return IsInNode(pos_hc, FLIP(pos_hc), heap_size);
}

template<typename size_type = std::size_t>
inline size_type IsInNodeMC(size_type N, size_type pos_mc, size_type heap_size) {
  return IsInNode(FLIP(pos_mc), pos_mc, heap_size);
}

/*
 * ================== START: Definition of AlmostBiheapify ====================
 */

template<class RAI, typename size_type = std::size_t>
bool IsAlmostBiheapCheckAlmostTripleConditionAtInNode(RAI V, size_type N, size_type in_node_hc) {
  size_type min_heap_parent_hc = Parent<size_type>(in_node_hc);
  size_type max_heap_parent_hc = FLIP(Parent<size_type>(FLIP(in_node_hc)));
  return *(V + min_heap_parent_hc) <= *(V + max_heap_parent_hc);
}

template<class RAI, typename size_type = std::size_t>
bool IsAlmostBiheapCheckAlmostQuadrupleCondition(RAI V, size_type N) {
  //assert(N % 3 == 2);
  size_type pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  size_type parent_of_pure_min_heap_double_arrow_node_hc = Parent<size_type>(pure_min_heap_double_arrow_node_hc);
  return *(V + parent_of_pure_min_heap_double_arrow_node_hc) <= *(V + FLIP(parent_of_pure_min_heap_double_arrow_node_hc));
}

template<class RAI, typename size_type = std::size_t>
void AlmostBiheapifyEnsureAlmostQuadrupleCondition(RAI V, size_type N) {
  //assert(N % 3 == 2);
  size_type pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  size_type parent_of_pure_min_heap_double_arrow_node_hc = Parent<size_type>(pure_min_heap_double_arrow_node_hc);
  if (*(V + parent_of_pure_min_heap_double_arrow_node_hc) > *(V + FLIP(parent_of_pure_min_heap_double_arrow_node_hc)))
    std::iter_swap(V + parent_of_pure_min_heap_double_arrow_node_hc, V + FLIP(parent_of_pure_min_heap_double_arrow_node_hc));
  return ;
}

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define an almost BiHeap.
 */
template<class RAI, typename size_type = std::size_t>
bool IsAlmostBiHeap(RAI V, size_type N) {
  if (N <= 3) {
    if(N <= 2)
      return true;
    else if (N == 3)
      return *V <= *(V + 2);
  }
  bool is_N_mod_3_equal_to_2 = N % 3 == 2;
  if (is_N_mod_3_equal_to_2 && !IsAlmostBiheapCheckAlmostQuadrupleCondition(V, N)) {
    return false;
  }
  size_type heap_size = HeapSize<size_type>(N);
  size_type first_in_node_hc = N - heap_size;
  {
    size_type one_past_last_node = heap_size - is_N_mod_3_equal_to_2;
    for (size_type in_hc = first_in_node_hc + is_N_mod_3_equal_to_2; in_hc < one_past_last_node; in_hc++) {
      if (!IsAlmostBiheapCheckAlmostTripleConditionAtInNode(V, N, in_hc))
        return false;
    }
  }

  //Check the min heap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < first_in_node_hc; i++) {
      auto parent_value = *(V + i);
      //Check that the parent and left child satisfy the min heap condition.
      if (parent_value > *(V + (right_child - 1)))
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (parent_value > *(V + right_child)) {
        if (!(is_N_mod_3_equal_to_2 && right_child == (N - 2) / 3))
          return false;
      }
    }
    //If the min heap's last non-In element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < first_in_node_hc
          && *(V + i) > *(V + left_child))
        return false;
    }
  }
  //Check the max heap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < first_in_node_hc; i++) {
      auto parent_value = *(V + FLIP(i));
      size_type mirror_left_child_hc = FLIP(right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_value < *(V + mirror_left_child_hc))
        return false;

      //Check that the parent and right child satisfy the max heap condition.
      if (parent_value < *(V + (mirror_left_child_hc - 1))) {
        if (!(is_N_mod_3_equal_to_2 && right_child == (N - 2) / 3))
          return false;
      }
    }
    //If the max heap's last non-In element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < first_in_node_hc
          && *(V + FLIP(i)) < *(V + FLIP(left_child)))
        return false;
    }
  }
  return true;
}

template<class RAI, typename size_type = std::size_t>
inline void AlmostBiHeapifySiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_in_node,
                                      size_type pos_hc,
                                      size_type last_node_in_biheap_hc) {
  bool is_total_num_nodes_mod_3_equal_to_2 = N % 3 == 2;
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    if (IsInNodeHC(N, left_child_hc, heap_size)) {
      size_type left_child_mc = Parent<size_type>(FLIP(left_child_hc));
      if (is_total_num_nodes_mod_3_equal_to_2 && left_child_mc == (N - 2) / 3)
        left_child_mc = Parent<size_type>(left_child_mc);
      left_child_hc = FLIP(left_child_mc);
    }
    if (IsInNodeHC(N, right_child_hc, heap_size)) {
      size_type right_child_mc = Parent<size_type>(FLIP(right_child_hc));
      if (is_total_num_nodes_mod_3_equal_to_2 && right_child_mc == (N - 2) / 3)
        right_child_mc = Parent<size_type>(right_child_mc);
      right_child_hc = FLIP(right_child_mc);
    }
    auto left_it   = V + left_child_hc;
    auto right_it  = V + right_child_hc;
    auto pos_it    = V + pos_hc;
    RAI smaller_it;

    is_right_child_valid = is_right_child_valid && right_child_hc <= last_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc <= last_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    if (!is_left_child_valid || (is_right_child_valid && *right_it < *left_it)) {
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
  SiftUpMaxHeapMC<RAI, size_type>(V, N, FLIP(pos_hc),
                                  FLIP(last_node_in_biheap_hc));
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void AlmostBiHeapifySiftFromMaxToMin(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_in_node,
                                      size_type pos_mc,
                                      size_type first_node_in_biheap_hc) {
  auto pos_hc = FLIP(pos_mc);
  bool is_total_num_nodes_mod_3_equal_to_2 = N % 3 == 2;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = FLIP(left_child_mc);
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size;
    if (IsInNode(left_child_hc, left_child_mc, heap_size)) {
      left_child_hc = Parent<size_type>(FLIP(left_child_mc));
      if (is_total_num_nodes_mod_3_equal_to_2 && left_child_hc == (N - 2) / 3)
        left_child_hc = Parent<size_type>(left_child_hc);
      left_child_mc = FLIP(left_child_hc);
    }
    if (IsInNode(right_child_hc, right_child_mc, heap_size)) {
      right_child_hc = Parent<size_type>(FLIP(right_child_mc));
      if (is_total_num_nodes_mod_3_equal_to_2 && right_child_hc == (N - 2) / 3)
        right_child_hc = Parent<size_type>(right_child_hc);
      right_child_mc = FLIP(right_child_hc);
    }
    auto pos_it    = V + pos_hc;
    auto left_it   = V + left_child_hc;
    auto right_it  = V + right_child_hc;
    RAI larger_it;

    is_right_child_valid = is_right_child_valid && right_child_hc >= first_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= first_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    if (!is_left_child_valid || (is_right_child_valid && *right_it > *left_it)) {
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
  SiftUpMinHeapHC<RAI, size_type>(V, pos_hc, first_node_in_biheap_hc);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void AlmostBiHeapify(RAI V, size_type N) {
  if(N < 3)
    return ;
  size_type heap_size               = HeapSize<size_type>(N);
  size_type first_in_node           = N - heap_size;
  size_type last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize<size_type>(heap_size);
  size_type first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    AlmostBiHeapifySiftFromMinToMax<RAI, size_type>(V, N,
        heap_size, first_in_node,
        --first_node_in_biheap_hc, last_node_in_biheap_hc);
    AlmostBiHeapifySiftFromMaxToMin<RAI, size_type>(V, N,
        heap_size, first_in_node,
        FLIP(++last_node_in_biheap_hc), first_node_in_biheap_hc);
  }
  return ;
}

/*
 * ================== END: Definition of AlmostBiheapify ====================
 */

/*
 * ================== START: Definition of lambda version of AlmostBiheapify ====================
 */


template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostTripleConditionAtInNode(RAI V, size_type N, size_type in_node_hc, LambdaType lambda) {
  size_type min_heap_parent_hc = Parent<size_type>(in_node_hc);
  size_type max_heap_parent_hc = FLIP(Parent<size_type>(FLIP(in_node_hc)));
  return *(V + lambda(N, min_heap_parent_hc)) <= *(V + lambda(N, max_heap_parent_hc));
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostQuadrupleCondition(RAI V, size_type N, LambdaType lambda) {
  //assert(N % 3 == 2);
  size_type pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  size_type parent_of_pure_min_heap_double_arrow_node_hc = Parent<size_type>(pure_min_heap_double_arrow_node_hc);
  return *(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc)) <= *(V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc)));
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
void AlmostBiheapifyEnsureAlmostQuadrupleCondition(RAI V, size_type N, LambdaType lambda) {
  //assert(N % 3 == 2);
  size_type pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  size_type parent_of_pure_min_heap_double_arrow_node_hc = Parent<size_type>(pure_min_heap_double_arrow_node_hc);
  if (*(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc)) > *(V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc))))
    std::iter_swap(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc), V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc)));
  return ;
}

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define an almost BiHeap.
 */
template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool IsAlmostBiHeap(RAI V, size_type N, LambdaType lambda) {
  if (N <= 3) {
    if(N <= 2)
      return true;
    else if (N == 3)
      return *(V + lambda(N, 0)) <= *(V + lambda(N, 2));
  }
  bool is_N_mod_3_equal_to_2 = N % 3 == 2;
  if (is_N_mod_3_equal_to_2 && !IsAlmostBiheapCheckAlmostQuadrupleCondition(V, N, lambda))
    return false;
  size_type heap_size = HeapSize<size_type>(N);
  size_type first_in_node_hc = N - heap_size;
  {
    size_type one_past_last_node = heap_size - is_N_mod_3_equal_to_2;
    for (size_type in_hc = first_in_node_hc + is_N_mod_3_equal_to_2; in_hc < one_past_last_node; in_hc++) {
      if (!IsAlmostBiheapCheckAlmostTripleConditionAtInNode(V, N, in_hc, lambda))
        return false;
    }
  }

  //Check the min heap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < first_in_node_hc; i++) {
      auto parent_value = *(V + lambda(N, i));
      //Check that the parent and left child satisfy the min heap condition.
      if (parent_value > *(V + lambda(N, (right_child - 1))))
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (parent_value > *(V + lambda(N, right_child))) {
        if (!(is_N_mod_3_equal_to_2 && right_child == (N - 2) / 3))
          return false;
      }
    }
    //If the min heap's last non-In element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < first_in_node_hc
          && *(V + lambda(N, i)) > *(V + lambda(N, left_child)))
        return false;
    }
  }
  //Check the max heap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < first_in_node_hc; i++) {
      auto parent_value = *(V + lambda(N, FLIP(i)));
      size_type mirror_left_child_hc = FLIP(right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, mirror_left_child_hc)))
        return false;

      //Check that the parent and right child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, (mirror_left_child_hc - 1)))) {
        if (!(is_N_mod_3_equal_to_2 && right_child == (N - 2) / 3))
          return false;
      }
    }
    //If the max heap's last non-In element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < first_in_node_hc
          && *(V + lambda(N, FLIP(i))) < *(V + lambda(N, FLIP(left_child))))
        return false;
    }
  }
  return true;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_hc,
                                      size_type largest_node_in_biheap_hc,
                                      LambdaType lambda) {
  bool is_total_num_nodes_mod_3_equal_to_2 = N % 3 == 2;
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= largest_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    bool is_left_child_valid = left_child_hc <= largest_node_in_biheap_hc &&
                               left_child_hc < heap_size;
    if (IsInNodeHC(N, left_child_hc, heap_size)) {
      size_type left_child_mc = Parent<size_type>(FLIP(left_child_hc));
      if (is_total_num_nodes_mod_3_equal_to_2 && left_child_mc == (N - 2) / 3)
        left_child_mc = Parent<size_type>(left_child_mc);
      left_child_hc = FLIP(left_child_mc);
    }
    if (IsInNodeHC(N, right_child_hc, heap_size)) {
      size_type right_child_mc = Parent<size_type>(FLIP(right_child_hc));
      if (is_total_num_nodes_mod_3_equal_to_2 && right_child_mc == (N - 2) / 3)
        right_child_mc = Parent<size_type>(right_child_mc);
      right_child_hc = FLIP(right_child_mc);
    }

    auto left_it   = V + lambda(N, left_child_hc);
    auto right_it  = V + lambda(N, right_child_hc);
    auto pos_it    = V + lambda(N, pos_hc);
    RAI smaller_it;

    is_right_child_valid = is_right_child_valid && right_child_hc <= largest_node_in_biheap_hc;
    is_left_child_valid  = is_left_child_valid && left_child_hc <= largest_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    if (!is_left_child_valid || (is_right_child_valid && *right_it < *left_it)) {
      smaller_it = right_it;
      pos_hc     = right_child_hc;
    } else { //Here, the left child is valid.
      smaller_it = left_it;
      pos_hc     = left_child_hc;
    }
    if (*pos_it > *smaller_it)
      std::iter_swap(pos_it, smaller_it);
    else
      return ;
  }
  SiftUpMaxHeapMC<RAI, size_type, LambdaType>(V, N, FLIP(pos_hc),
                                  FLIP(largest_node_in_biheap_hc), lambda);
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMaxToMin(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_mc,
                                      size_type smallest_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  bool is_total_num_nodes_mod_3_equal_to_2 = N % 3 == 2;
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = FLIP(left_child_mc);
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size && right_child_hc >= smallest_node_in_biheap_hc;
    if (IsInNode(left_child_hc, left_child_mc, heap_size)) {
      left_child_hc = Parent<size_type>(FLIP(left_child_mc));
      if (is_total_num_nodes_mod_3_equal_to_2 && left_child_hc == (N - 2) / 3)
        left_child_hc = Parent<size_type>(left_child_hc);
      left_child_mc = FLIP(left_child_hc);
    }
    if (IsInNode(right_child_hc, right_child_mc, heap_size)) {
      right_child_hc = Parent<size_type>(FLIP(right_child_mc));
      if (is_total_num_nodes_mod_3_equal_to_2 && right_child_hc == (N - 2) / 3)
        right_child_hc = Parent<size_type>(right_child_hc);
      right_child_mc = FLIP(right_child_hc);
    }
    auto pos_it    = V + lambda(N, pos_hc);
    auto left_it   = V + lambda(N, left_child_hc);
    auto right_it  = V + lambda(N, right_child_hc);
    RAI larger_it;

    is_right_child_valid = is_right_child_valid && right_child_hc >= smallest_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= smallest_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    if (!is_left_child_valid || (is_right_child_valid && *right_it > *left_it)) {
      larger_it = right_it;
      pos_hc    = right_child_hc;
      pos_mc    = right_child_mc;
    } else { //Here, the left child is valid.
      larger_it = left_it;
      pos_hc    = left_child_hc;
      pos_mc    = left_child_mc;
    }
    if (*pos_it < *larger_it)
      std::iter_swap(pos_it, larger_it);
    else
      return ;
  }
  SiftUpMinHeapHC<RAI, size_type, LambdaType>(V, N, pos_hc, smallest_node_in_biheap_hc, lambda);
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void AlmostBiHeapify(RAI V, size_type N, LambdaType lambda) {
  if(N < 3)
    return ;
  size_type heap_size          = HeapSize<size_type>(N);
  size_type first_node_in_mirror_heap  = N - heap_size;
  size_type largest_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize<size_type>(heap_size);
  size_type smallest_node_in_biheap_hc = FLIP(largest_node_in_biheap_hc);
  while (smallest_node_in_biheap_hc > 0) {
    AlmostBiHeapifySiftFromMinToMax<RAI, size_type, LambdaType>(V, N,
        heap_size, first_node_in_mirror_heap,
        --smallest_node_in_biheap_hc, largest_node_in_biheap_hc, lambda);
    AlmostBiHeapifySiftFromMaxToMin<RAI, size_type, LambdaType>(V, N,
        heap_size, first_node_in_mirror_heap,
        FLIP(++largest_node_in_biheap_hc), smallest_node_in_biheap_hc, lambda);
  }
  return ;
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify ====================
 */

#undef FLIP

#endif /* ALMOST_BIHEAPIFY_H_ */
