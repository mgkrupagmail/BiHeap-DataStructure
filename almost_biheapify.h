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
#include "biheapify_lambda.h"

#include <algorithm>

#define FLIP(a) ((N - 1) - (a))

template<typename SizeType = std::size_t>
inline SizeType IsInNode(SizeType pos_hc, SizeType pos_mc, SizeType heap_size) {
  return (pos_hc < heap_size) && (pos_mc < heap_size);
}

template<typename SizeType = std::size_t>
inline SizeType IsInNodeHC(SizeType N, SizeType pos_hc, SizeType heap_size) {
  return IsInNode(pos_hc, FLIP(pos_hc), heap_size);
}

template<typename SizeType = std::size_t>
inline SizeType IsInNodeMC(SizeType N, SizeType pos_mc, SizeType heap_size) {
  return IsInNode(FLIP(pos_mc), pos_mc, heap_size);
}

/*
 * ================== START: Definition of lambda version of IsAlmostBiHeap ====================
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostTripleConditionAtInNode(RAI V, SizeType N, SizeType in_node_hc, LambdaType lambda) {
  SizeType min_heap_parent_hc = Parent<SizeType>(in_node_hc);
  SizeType max_heap_parent_hc = FLIP(Parent<SizeType>(FLIP(in_node_hc)));
  return *(V + lambda(N, min_heap_parent_hc)) <= *(V + lambda(N, max_heap_parent_hc));
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostQuadrupleCondition(RAI V, SizeType N, LambdaType lambda) {
  SizeType pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  SizeType parent_of_pure_min_heap_double_arrow_node_hc = Parent<SizeType>(pure_min_heap_double_arrow_node_hc);
  return *(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc)) <= *(V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc)));
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
void AlmostBiheapifyEnsureAlmostQuadrupleCondition(RAI V, SizeType N, LambdaType lambda) {
  SizeType pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  SizeType parent_of_pure_min_heap_double_arrow_node_hc = Parent<SizeType>(pure_min_heap_double_arrow_node_hc);
  if (*(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc)) > *(V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc))))
    std::iter_swap(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc), V + lambda(N, FLIP(parent_of_pure_min_heap_double_arrow_node_hc)));
  return ;
}

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define an almost BiHeap.
 */
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiHeap(RAI V, SizeType N, LambdaType lambda) {
  if (N <= 3) {
    if(N <= 2)
      return true;
    else if (N == 3)
      return *(V + lambda(N, 0)) <= *(V + lambda(N, 2));
  }
  bool is_N_mod_3_equal_to_2 = N % 3 == 2;
  if (is_N_mod_3_equal_to_2 && !IsAlmostBiheapCheckAlmostQuadrupleCondition(V, N, lambda))
    return false;
  SizeType heap_size = HeapSize<SizeType>(N);
  SizeType first_in_node_hc = N - heap_size;
  {
    SizeType one_past_last_node = heap_size - is_N_mod_3_equal_to_2;
    for (SizeType in_hc = first_in_node_hc + is_N_mod_3_equal_to_2; in_hc < one_past_last_node; in_hc++) {
      if (!IsAlmostBiheapCheckAlmostTripleConditionAtInNode(V, N, in_hc, lambda))
        return false;
    }
  }

  //Check the min heap condition.
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
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
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < first_in_node_hc
          && *(V + lambda(N, i)) > *(V + lambda(N, left_child)))
        return false;
    }
  }
  //Check the max heap condition.
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
                                                    < first_in_node_hc; i++) {
      auto parent_value = *(V + lambda(N, FLIP(i)));
      SizeType mirror_left_child_hc = FLIP(right_child - 1);
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
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < first_in_node_hc
          && *(V + lambda(N, FLIP(i))) < *(V + lambda(N, FLIP(left_child))))
        return false;
    }
  }
  return true;
}


template<class RAI, typename SizeType = std::size_t>
bool IsAlmostBiHeap(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType { return i; };
  return IsAlmostBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
}


/*
 * ================== END: Definition of lambda version of IsAlmostBiHeap ====================
 */

/*
 * ================== START: Definition of lambda version of IsAlmostBiHeap with some permitted In nodes ====================
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiHeap(RAI V, SizeType N, SizeType fuse_first_hc,
                    SizeType fuse_last_hc, LambdaType lambda) {
  if(N < 2)
    return true;
  if (fuse_first_hc > fuse_last_hc) //If all nodes are permitted then it should be a BiHeap.
    return IsBiHeap<RAI, SizeType, LambdaType>(V, N, lambda);
  if(N == 2 || (fuse_last_hc - fuse_first_hc <= 0))
    return true;
  if (!IsAlmostBiHeap<RAI, SizeType, LambdaType>(V, N, lambda))
    return false ;
  if (N % 3 == 2) { //Then make sure that the double arrow satisfies the necessary conditions.
    SizeType pmin_double_arrow_node_hc = (N - 2) / 3;
    SizeType pmax_double_arrow_node_hc = (2 * N - 1) / 3;
    SizeType minH_parent_of_pmin_double_arrow_node_hc = Parent<SizeType>(pmin_double_arrow_node_hc);
    SizeType maxH_parent_of_pmax_double_arrow_node_hc = FLIP(minH_parent_of_pmin_double_arrow_node_hc);
    RAI pmin_double_arrow_node_it = V + lambda(N, pmin_double_arrow_node_hc);
    RAI pmax_double_arrow_node_it = V + lambda(N, pmax_double_arrow_node_hc);
    RAI minH_parent_of_pmin_double_arrow_node_it = V + lambda(N, minH_parent_of_pmin_double_arrow_node_hc);
    RAI maxH_parent_of_pmax_double_arrow_node_it = V + lambda(N, maxH_parent_of_pmax_double_arrow_node_hc);
    bool is_pmin_double_arrow_node_forbidden =
        fuse_first_hc <= pmin_double_arrow_node_hc
        && pmin_double_arrow_node_hc <= fuse_last_hc;
    bool is_pmax_double_arrow_node_forbidden =
        fuse_first_hc <= pmax_double_arrow_node_hc
        && pmax_double_arrow_node_hc <= fuse_last_hc;
    if (!is_pmin_double_arrow_node_forbidden && *pmin_double_arrow_node_it < *minH_parent_of_pmin_double_arrow_node_it)
      return false;
    if (!is_pmax_double_arrow_node_forbidden && *pmax_double_arrow_node_it > *maxH_parent_of_pmax_double_arrow_node_it)
      return false;
    if (!is_pmin_double_arrow_node_forbidden && !is_pmax_double_arrow_node_forbidden
        && *pmin_double_arrow_node_it > *pmax_double_arrow_node_it)
      return false;
    else if (is_pmin_double_arrow_node_forbidden && !is_pmax_double_arrow_node_forbidden
        && *minH_parent_of_pmin_double_arrow_node_it > *pmax_double_arrow_node_it)
      return false;
    else if (!is_pmin_double_arrow_node_forbidden && is_pmax_double_arrow_node_forbidden
        && *pmin_double_arrow_node_it > *maxH_parent_of_pmax_double_arrow_node_it)
      return false;
  }
  SizeType heap_size              = HeapSize(N);
  SizeType first_in_node          = N - heap_size;
  SizeType first_in_node_to_check = first_in_node + (N % 3 == 2); //Avoid the double arrow if it exists.
  SizeType last_in_node_to_check  = FLIP(first_in_node_to_check);
  for (SizeType pos_hc = first_in_node_to_check; pos_hc < fuse_first_hc; pos_hc++) {
    SizeType minh_parent_of_pos_hc = Parent<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_hc = FLIP(Parent<SizeType>(FLIP(pos_hc)));
    auto pos_value = *(V + lambda(N, pos_hc));
    if (pos_value < *(V + lambda(N, minh_parent_of_pos_hc)) ||
        pos_value > *(V + lambda(N, maxh_parent_of_pos_hc))) {
      return false;
    }
  }
  for (SizeType pos_hc = fuse_last_hc + 1; pos_hc <= last_in_node_to_check; pos_hc++) {
    SizeType minh_parent_of_pos_hc = Parent<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_hc = FLIP(Parent<SizeType>(FLIP(pos_hc)));
    auto pos_value = *(V + lambda(N, pos_hc));
    if (pos_value < *(V + lambda(N, minh_parent_of_pos_hc)) ||
        pos_value > *(V + lambda(N, maxh_parent_of_pos_hc))) {
      return false;
    }
  }
  return true;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiHeap(RAI V, SizeType N, SizeType num_permitted_in_nodes, LambdaType lambda) {
  if(N < 2 || (N == 2 && num_permitted_in_nodes <= 1))
    return true;
  if (num_permitted_in_nodes <= 0)
    return IsAlmostBiHeap<RAI, SizeType, LambdaType>(V, N, lambda);
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  SizeType num_in_nodes  = heap_size - first_in_node;
  if (num_permitted_in_nodes >= num_in_nodes)
    return IsBiHeap<RAI, SizeType, LambdaType>(V, N, lambda);
  if (N == 2 && num_permitted_in_nodes <= 1)
    return true;//Then it's trivially an almost BiHeap.
  //If there is an num_permitted_in_nodes is odd, then
  // allow there to be one more permitted pure min heap
  // In node than there are permitted pure max heap In nodes.
  SizeType fuse_first_hc = first_in_node + (num_permitted_in_nodes + 1) / 2;
  SizeType fuse_last_hc = FLIP(first_in_node + (num_permitted_in_nodes / 2));
  return IsAlmostBiHeap<RAI, SizeType, LambdaType>(V, N, fuse_first_hc,
                                                   fuse_last_hc, lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool IsAlmostBiHeap(RAI V, SizeType N,  SizeType num_permitted_in_nodes) {
  auto trivial_lambda = [](SizeType N, SizeType i) -> SizeType { return i; };
  return IsAlmostBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                       num_permitted_in_nodes, trivial_lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool IsAlmostBiHeap(RAI V, SizeType N, SizeType fuse_first_hc, SizeType fuse_last_hc) {
  auto trivial_lambda = [](SizeType N, SizeType i) -> SizeType { return i; };
  return IsAlmostBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                  fuse_first_hc, fuse_last_hc, trivial_lambda);
}

/*
 * ================== END: Definition of lambda version of IsAlmostBiHeap with some permitted In nodes ====================
 */

/*
 * ================== START: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 != 2 case =============
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMinToMaxIgnoreDoubledHeadedArrow(RAI V, SizeType N,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc,
                                      SizeType fuse_first_hc,
                                      SizeType fuse_last_hc,
                                      LambdaType lambda) {
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    bool is_left_child_valid = left_child_hc <= last_node_in_biheap_hc &&
                               left_child_hc < heap_size;
    if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc)
      left_child_hc = FLIP(Parent<SizeType>(FLIP(left_child_hc)));
    if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc)
      right_child_hc = FLIP(Parent<SizeType>(FLIP(right_child_hc)));
    is_right_child_valid = is_right_child_valid && right_child_hc <= last_node_in_biheap_hc;
    is_left_child_valid  = is_left_child_valid && left_child_hc <= last_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    auto pos_it   = V + lambda(N, pos_hc);
    RAI smaller_it;
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
  SiftUpMaxHeapMC<RAI, SizeType, LambdaType>(V, N, FLIP(pos_hc),
                                  FLIP(last_node_in_biheap_hc), lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMaxToMinIgnoreDoubledHeadedArrow(RAI V, SizeType N,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      SizeType fuse_first_hc,
                                      SizeType fuse_last_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = FLIP(left_child_mc);
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size && right_child_hc >= first_node_in_biheap_hc;
    if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc) {
      left_child_hc = Parent<SizeType>(left_child_hc);
      left_child_mc = FLIP(left_child_hc);
    }
    if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc) {
      right_child_hc = Parent<SizeType>(right_child_hc);
      right_child_mc = FLIP(right_child_hc);
    }
    is_right_child_valid = is_right_child_valid && right_child_hc >= first_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= first_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto pos_it   = V + lambda(N, pos_hc);
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI larger_it;
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
  SiftUpMinHeapHC<RAI, SizeType, LambdaType>(V, N, pos_hc, first_node_in_biheap_hc, lambda);
  return ;
}

//Assumes that N % 3 != 2 or that N % 3 == 2 but neither
// endpoint of the double headed arrow is in the interval
// [fuse_first_hc, fuse_last_hc].
//If fuse_last_hc < fuse_first_hc then it calls BiHeapify.
//If fuse_first_hc is not an In node then it is set to
//  2 * HeapSize(N) - N, the first min heap In node.
//If fust_last_hc is not an In node then it is set to
//  HeapSize(N) - 1, the last In node.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifyIgnoreDoubledHeadedArrow(RAI V, SizeType N,
    SizeType fuse_first_hc, SizeType fuse_last_hc, LambdaType lambda) {
  if (N < 2)
    return ;
  if (fuse_last_hc < fuse_first_hc) { //Then all nodes are permitted.
    BiHeapify<RAI, SizeType, LambdaType>(V, N, lambda);
    return ;
  }
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  if (fuse_first_hc < first_in_node)
    fuse_first_hc = first_in_node;
  if (fuse_last_hc >= heap_size)
    fuse_last_hc = heap_size - 1;
  SizeType last_node_in_biheap_hc  = (heap_size - 1) - (N % 3 == 2);
  SizeType first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    AlmostBiHeapifySiftFromMinToMaxIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N,
        heap_size, first_in_node, --first_node_in_biheap_hc,
        last_node_in_biheap_hc, fuse_first_hc, fuse_last_hc, lambda);
    AlmostBiHeapifySiftFromMaxToMinIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N,
        heap_size, first_in_node, FLIP(++last_node_in_biheap_hc),
        first_node_in_biheap_hc, fuse_first_hc, fuse_last_hc, lambda);
  }
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 != 2 case ===============
 */

/*
 * ================== START: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 == 2 case =============
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMinToMaxWithDoubleHeadedArrow(RAI V,
                                          SizeType N,
                                          SizeType heap_size,
                                          SizeType first_in_node,
                                          SizeType pos_hc,
                                          SizeType last_node_in_biheap_hc,
                                          SizeType fuse_first_hc,
                                          SizeType fuse_last_hc,
                                          SizeType pmin_double_arrow_end_hc,
                                          SizeType pmax_double_arrow_end_hc,
                                          SizeType maxh_parent_of_pmax_double_arrow_end_hc,
                                          LambdaType lambda) {
  if (fuse_first_hc <= pos_hc && pos_hc <= fuse_last_hc) //This can happen with a double arrow.
    return ;
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    bool is_left_child_valid = left_child_hc <= last_node_in_biheap_hc &&
                               left_child_hc < heap_size;
    if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc) {
      left_child_hc = FLIP(Parent<SizeType>(FLIP(left_child_hc)));
      //Note: The following is satisfied if and only if the following holds:
      // FLIP(fuse_first_hc) >= left_child_mc && left_child_mc >= FLIP(fuse_last_hc)
      if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc)
        left_child_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }
    if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc) {
      right_child_hc = FLIP(Parent<SizeType>(FLIP(right_child_hc)));
      if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc)
        right_child_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }

    is_right_child_valid = is_right_child_valid && right_child_hc <= last_node_in_biheap_hc;
    is_left_child_valid  = is_left_child_valid  && left_child_hc  <= last_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    auto pos_it   = V + lambda(N, pos_hc);
    RAI smaller_it;
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
  //At this point, it's not possible to be be simultaneously
  //  forbidden, in the pure max heap, and incident to the double arrow.
  if (pos_hc == pmin_double_arrow_end_hc) {
    //If the other end of the double arrow is forbidden.
    if (fuse_first_hc <= pmax_double_arrow_end_hc && pmax_double_arrow_end_hc <= fuse_last_hc) {
      auto max_heap_parent_of_other_end_it = V + lambda(N, maxh_parent_of_pmax_double_arrow_end_hc);
      auto pos_it                          = V + lambda(N, pos_hc);
      //Perform one iteration of sifting up the min heap while skipping the
      // forbidden node.
      if (maxh_parent_of_pmax_double_arrow_end_hc <= last_node_in_biheap_hc &&
          *max_heap_parent_of_other_end_it < *pos_it)
        std::iter_swap(pos_it, max_heap_parent_of_other_end_it);
      else
        return ;
      pos_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }
  }
  SiftUpMaxHeapMC<RAI, SizeType, LambdaType>(V, N, FLIP(pos_hc),
                                      FLIP(last_node_in_biheap_hc), lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifySiftFromMaxToMinWithDoubleHeadedArrow(RAI V,
                                      SizeType N,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      SizeType fuse_first_hc,
                                      SizeType fuse_last_hc,
                                      SizeType pmin_double_arrow_end_hc,
                                      SizeType minh_parent_of_pmin_double_arrow_end_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  if (fuse_first_hc <= pos_hc && pos_hc <= fuse_last_hc)
    return ;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = FLIP(left_child_mc);
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size && right_child_hc >= first_node_in_biheap_hc;
    if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc) {
      left_child_hc = Parent<SizeType>(left_child_hc);
      if (fuse_first_hc <= left_child_hc && left_child_hc <= fuse_last_hc)
        left_child_hc = minh_parent_of_pmin_double_arrow_end_hc;
      left_child_mc = FLIP(left_child_hc);
    }
    if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc) {
      right_child_hc = Parent<SizeType>(right_child_hc);
      if (fuse_first_hc <= right_child_hc && right_child_hc <= fuse_last_hc)
        right_child_hc = minh_parent_of_pmin_double_arrow_end_hc;
      right_child_mc = FLIP(right_child_hc);
    }
    is_right_child_valid = is_right_child_valid && right_child_hc >= first_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= first_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto pos_it   = V + lambda(N, pos_hc);
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI larger_it;
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
  //At this point, it's not possible to be be simultaneously
  //  forbidden, in the pure min heap, and incident to the double arrow.
  if (pos_mc == pmin_double_arrow_end_hc) { //if and only if pos_hc = pmax_double_arrow_end_hc
    //If the other end of the double arrow is forbidden.
    if (fuse_first_hc <= pmin_double_arrow_end_hc && pmin_double_arrow_end_hc <= fuse_last_hc) {
      auto min_heap_parent_of_other_end_it = V + lambda(N, minh_parent_of_pmin_double_arrow_end_hc);
      auto pos_it                          = V + lambda(N, pos_hc);
      //Perform one iteration of sifting up the min heap while skipping the
      // forbidden node.
      if (minh_parent_of_pmin_double_arrow_end_hc >= first_node_in_biheap_hc &&
          *min_heap_parent_of_other_end_it > *pos_it)
        std::iter_swap(pos_it, min_heap_parent_of_other_end_it);
      else
        return ;
      pos_hc = minh_parent_of_pmin_double_arrow_end_hc;
    }
  }
  SiftUpMinHeapHC<RAI, SizeType, LambdaType>(V, N, pos_hc, first_node_in_biheap_hc, lambda);
  return ;
}

//Assumes that N % 3 == 2.
// If fuse_last_hc < fuse_first_hc then it calls BiHeapify.
// If fuse_first_hc is not an In node then it is set to
//   2 * HeapSize(N) - N, the first min heap In node.
// If fust_last_hc is not an In node then it is set to
//   HeapSize(N) - 1, the last In node.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapifyWithDoubleHeadedArrow(RAI V, SizeType N,
    SizeType fuse_first_hc, SizeType fuse_last_hc, LambdaType lambda) {
  if(N <= 2) {
    if (N == 2 && fuse_last_hc < fuse_first_hc) {
      RAI V_0 = V + lambda(N, 0);
      RAI V_1 = V + lambda(N, 1);
      if (*V_0 > *V_1)
        std::iter_swap(V_0, V_1);
    }
    return ;
  }
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  if (fuse_first_hc > first_in_node && fuse_last_hc < heap_size - 1) {
    //If we don't have to worry about skipping over either one of the end
    // nodes of the double headed arrow, then we may as well use the more
    // efficient AlmostBiHeapify algorithm.
    AlmostBiHeapifyIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N,
                                          fuse_first_hc, fuse_last_hc, lambda);
    return ;
  }
  if (fuse_last_hc < fuse_first_hc) { //Then all nodes are permitted.
    BiHeapify<RAI, SizeType, LambdaType>(V, N, lambda);
    return ;
  }
  if (fuse_first_hc < first_in_node)
    fuse_first_hc = first_in_node;
  if (fuse_last_hc >= heap_size)
    fuse_last_hc = heap_size - 1; //The last In node.
  SizeType last_node_in_biheap_hc  = heap_size - 2;
  SizeType first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  //To increase efficiency, precompute the following values and pass them to the two calls
  // in the while loop. Any half descent optimizer will avoid actually allocating
  // additional space on the stack and copying these values into local variables since
  // these functions are both inlined and templates.
  SizeType pmin_double_arrow_end_hc                = (N - 2) / 3;
  SizeType pmax_double_arrow_end_hc                = 2 * pmin_double_arrow_end_hc + 1;
  SizeType minh_parent_of_pmin_double_arrow_end_hc = Parent<SizeType>(pmin_double_arrow_end_hc);
  SizeType maxh_parent_of_pmax_double_arrow_end_hc = FLIP(minh_parent_of_pmin_double_arrow_end_hc);
  while (first_node_in_biheap_hc > 0) {
    AlmostBiHeapifySiftFromMinToMaxWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N,
        heap_size, first_in_node, --first_node_in_biheap_hc,
        last_node_in_biheap_hc, fuse_first_hc, fuse_last_hc,
        pmin_double_arrow_end_hc, pmax_double_arrow_end_hc,
        maxh_parent_of_pmax_double_arrow_end_hc, lambda);
    AlmostBiHeapifySiftFromMaxToMinWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N,
        heap_size, first_in_node, FLIP(++last_node_in_biheap_hc),
        first_node_in_biheap_hc, fuse_first_hc, fuse_last_hc,
        pmin_double_arrow_end_hc,
        minh_parent_of_pmin_double_arrow_end_hc, lambda);
  }
  return ;
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 == 2 case ===============
 */

/*
 * ================== START: Definition of lambda version of AlmostBiheapify with some permitted In nodes and specializations ============
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapify(RAI V, SizeType N,
    SizeType fuse_first_hc,
    SizeType fuse_last_hc, LambdaType lambda) {
  if (N % 3 != 2)
    AlmostBiHeapifyIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, fuse_first_hc, fuse_last_hc, lambda);
  else
    AlmostBiHeapifyWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, fuse_first_hc, fuse_last_hc, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapify(RAI V, SizeType N,
    SizeType num_permitted_in_nodes, LambdaType lambda) {
  if (N < 2)
    return ;
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  SizeType num_in_nodes  = heap_size - first_in_node;
  if (num_permitted_in_nodes > num_in_nodes)
    num_permitted_in_nodes = num_in_nodes;
  if (N == 2 && num_permitted_in_nodes <= 1)
    return ;//Then there's nothing to do.
  if (num_permitted_in_nodes == num_in_nodes) {
    BiHeapify<RAI, SizeType, LambdaType>(V, N, lambda);
    return ;
  }
  //If num_permitted_in_nodes is odd, then allow
  // for there to be one more permitted pure min heap
  // In node than there are permitted pure max heap In nodes.
  SizeType fuse_first_hc = first_in_node + (num_permitted_in_nodes + 1) / 2;
  SizeType fuse_last_hc  = FLIP(first_in_node + (num_permitted_in_nodes / 2));
  AlmostBiHeapify<RAI, SizeType, LambdaType>(V, N, fuse_first_hc,
                                             fuse_last_hc, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t>
inline void AlmostBiHeapify(RAI V, SizeType N,
    SizeType num_permitted_in_nodes) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType { return i; };
  AlmostBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                       num_permitted_in_nodes, trivial_lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t>
inline void AlmostBiHeapify(RAI V, SizeType N,
    SizeType fuse_first_hc,
    SizeType fuse_last_hc) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType { return i; };
  AlmostBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N, fuse_first_hc,
                                                     fuse_last_hc, trivial_lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void AlmostBiHeapify(RAI V, SizeType N, LambdaType lambda) {
  AlmostBiHeapify<RAI, SizeType, LambdaType>(V, N, 0, N, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t>
inline void AlmostBiHeapify(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType { return i; };
  AlmostBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
  return ;
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify with some permitted In nodes and specializations ==========
 */

#undef FLIP

#endif /* ALMOST_BIHEAPIFY_H_ */
