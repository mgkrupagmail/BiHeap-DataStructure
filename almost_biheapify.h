/*
 * almost_biheapify.h
 *
 * The FusedBiHeapify() function performs the same operation as
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
 * The FusedBiHeapify() operation rearranges the values in this
 *  almost BiHeap graph so that if two node, say with min heap coordinates
 *  i and j where i < j, are incident to the same edge then the value
 *  of node i is <= the value of node j.
 * Note that at the end of the FusedBiHeapify() operation, the values
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

//#define FLIP(a) ((N - 1) - (a))

template<typename SizeType = std::size_t>
inline SizeType IsInNode(SizeType pos_hc, SizeType pos_mc, SizeType heap_size) {
  return (pos_hc < heap_size) && (pos_mc < heap_size);
}

template<typename SizeType = std::size_t>
inline SizeType IsInNodeHC(SizeType N, SizeType pos_hc, SizeType heap_size) {
  return IsInNode(pos_hc, (N - 1) - pos_hc, heap_size);
}

template<typename SizeType = std::size_t>
inline SizeType IsInNodeMC(SizeType N, SizeType pos_mc, SizeType heap_size) {
  return IsInNode((N - 1) - pos_mc, pos_mc, heap_size);
}

/*
 * ================== START: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 != 2 case =============
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifySiftFromMinToMaxIgnoreDoubledHeadedArrow(RAI V,
                                      SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc,
                                      SizeType F_first_hc,
                                      SizeType F_last_hc,
                                      LambdaType lambda) {
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    bool is_left_child_valid = left_child_hc <= last_node_in_biheap_hc &&
                               left_child_hc < heap_size;
    if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc)
      left_child_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - left_child_hc);
    if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc)
      right_child_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - right_child_hc);
    is_right_child_valid = is_right_child_valid && right_child_hc <= last_node_in_biheap_hc;
    is_left_child_valid  = is_left_child_valid && left_child_hc <= last_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI smaller_it;
    if (!is_left_child_valid || (is_right_child_valid && *right_it < *left_it)) {
      smaller_it    = right_it;
      pos_hc        = right_child_hc;
    } else { //Here, the left child is valid.
      smaller_it    = left_it;
      pos_hc        = left_child_hc;
    }
    if (*smaller_it < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  }
  //Start sifting up the max heap.
  //At this point pos is an In node.
  auto pos_mc                 = N_minus1 - pos_hc;
  auto last_node_in_biheap_mc = N_minus1 - last_node_in_biheap_hc;
  SizeType parent_mc          = ParentNotRoot<SizeType>(pos_mc);
  if (pos_mc == 0 || parent_mc < last_node_in_biheap_mc)
    return ;
  do {
    auto parent_it = V + lambda(N, N_minus1 - parent_mc);
    if (pos_value > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent_mc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) >= last_node_in_biheap_mc);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifySiftFromMaxToMinIgnoreDoubledHeadedArrow(RAI V,
                                      SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      SizeType F_first_hc,
                                      SizeType F_last_hc,
                                      LambdaType lambda) {
  auto pos_hc    = N_minus1 - pos_mc;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = N_minus1 - left_child_mc;
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size && right_child_hc >= first_node_in_biheap_hc;
    if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc) {
      left_child_hc = ParentNotRoot<SizeType>(left_child_hc);
      left_child_mc = N_minus1 - left_child_hc;
    }
    if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc) {
      right_child_hc = ParentNotRoot<SizeType>(right_child_hc);
      right_child_mc = N_minus1 - right_child_hc;
    }
    is_right_child_valid = is_right_child_valid && right_child_hc >= first_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= first_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI larger_it;
    if (!is_left_child_valid || (is_right_child_valid && *left_it < *right_it)) {
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else { //Here, the left child is valid.
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < *larger_it) {
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    } else
      return ;
  }
  //Start sifting up the min heap.
  //At this point pos is an In node.
  SizeType parent_hc = ParentNotRoot<SizeType>(pos_hc);
  if (pos_hc == 0 || parent_hc < first_node_in_biheap_hc)
    return ;
  do {
    auto parent_it = V + lambda(N, parent_hc);
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent_hc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 &&
     (parent_hc = ParentNotRoot<SizeType>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

//Assumes that N % 3 != 2 or that N % 3 == 2 but neither
// endpoint of the double headed arrow is in the interval
// [F_first_hc, F_last_hc].
//If F_last_hc < F_first_hc then it calls BiHeapify.
//If F_first_hc is not an In node then it is set to
//  2 * HeapSize(N) - N, the first min heap In node.
//If fust_last_hc is not an In node then it is set to
//  HeapSize(N) - 1, the last In node.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifyIgnoreDoubledHeadedArrow(RAI V, SizeType N,
    SizeType F_first_hc, SizeType F_last_hc, LambdaType lambda) {
  if (N < 2)
    return ;
  if (F_last_hc < F_first_hc) { //Then all nodes are permitted.
    BiHeapify<RAI, SizeType, LambdaType>(V, N, lambda);
    return ;
  }
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  if (F_first_hc < first_in_node)
    F_first_hc = first_in_node;
  if (F_last_hc >= heap_size)
    F_last_hc = heap_size - 1;
  SizeType last_node_in_biheap_hc  = (heap_size - 1) - (N % 3 == 2);
  SizeType N_minus1                = N - 1;
  SizeType first_node_in_biheap_hc = N_minus1 - last_node_in_biheap_hc;
  while (first_node_in_biheap_hc > 0) {
    FusedBiHeapifySiftFromMinToMaxIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, --first_node_in_biheap_hc,
        last_node_in_biheap_hc, F_first_hc, F_last_hc, lambda);
    FusedBiHeapifySiftFromMaxToMinIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, N_minus1 - (++last_node_in_biheap_hc),
        first_node_in_biheap_hc, F_first_hc, F_last_hc, lambda);
  }
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 != 2 case ===============
 */

/*
 * ================== START: Definition of lambda version of AlmostBiheapify with some permitted In nodes, N mod 3 == 2 case =============
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifySiftFromMinToMaxWithDoubleHeadedArrow(RAI V,
                                          SizeType N, SizeType N_minus1,
                                          SizeType heap_size,
                                          SizeType first_in_node,
                                          SizeType pos_hc,
                                          SizeType last_node_in_biheap_hc,
                                          SizeType F_first_hc,
                                          SizeType F_last_hc,
                                          SizeType pmin_double_arrow_end_hc,
                                          SizeType pmax_double_arrow_end_hc,
                                          SizeType maxh_parent_of_pmax_double_arrow_end_hc,
                                          LambdaType lambda) {
  if (F_first_hc <= pos_hc && pos_hc <= F_last_hc) //This can happen with a double arrow.
    return ;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    bool is_left_child_valid = left_child_hc <= last_node_in_biheap_hc &&
                               left_child_hc < heap_size;
    if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc) {
      left_child_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - left_child_hc);
      //Note: The following is satisfied if and only if the following holds:
      // FLIP(F_first_hc) >= left_child_mc && left_child_mc >= FLIP(F_last_hc)
      if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc)
        left_child_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }
    if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc) {
      right_child_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - right_child_hc);
      if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc)
        right_child_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }

    is_right_child_valid = is_right_child_valid && right_child_hc <= last_node_in_biheap_hc;
    is_left_child_valid  = is_left_child_valid  && left_child_hc  <= last_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI smaller_it;
    if (!is_left_child_valid || (is_right_child_valid && *right_it < *left_it)) {
      smaller_it    = right_it;
      pos_hc        = right_child_hc;
    } else { //Here, the left child is valid.
      smaller_it    = left_it;
      pos_hc        = left_child_hc;
    }
    if (*smaller_it < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  }
  //At this point, it's not possible to be be simultaneously
  //  forbidden, in the pure max heap, and incident to the double arrow.
  if (pos_hc == pmin_double_arrow_end_hc) {
    //If the other end of the double arrow is forbidden.
    if (F_first_hc <= pmax_double_arrow_end_hc && pmax_double_arrow_end_hc <= F_last_hc) {
      auto max_heap_parent_of_other_end_it = V + lambda(N, maxh_parent_of_pmax_double_arrow_end_hc);
      //Perform one iteration of sifting up the min heap while skipping the
      // fused node.
      if (maxh_parent_of_pmax_double_arrow_end_hc <= last_node_in_biheap_hc &&
          *max_heap_parent_of_other_end_it < pos_value) {
        std::iter_swap(pos_it, max_heap_parent_of_other_end_it);
        pos_it = max_heap_parent_of_other_end_it;
      } else
        return ;
      pos_hc = maxh_parent_of_pmax_double_arrow_end_hc;
    }
  }
  //Start sifting up the max heap.
  //At this point pos is an In node.
  auto pos_mc                 = N_minus1 - pos_hc;
  auto last_node_in_biheap_mc = N_minus1 - last_node_in_biheap_hc;
  SizeType parent_mc          = ParentNotRoot<SizeType>(pos_mc);
  if (pos_mc == 0 || parent_mc < last_node_in_biheap_mc)
    return ;
  do {
    auto parent_it = V + lambda(N, N_minus1 - parent_mc);
    if (pos_value > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent_mc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) >= last_node_in_biheap_mc);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifySiftFromMaxToMinWithDoubleHeadedArrow(RAI V,
                                      SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      SizeType F_first_hc,
                                      SizeType F_last_hc,
                                      SizeType pmin_double_arrow_end_hc,
                                      SizeType minh_parent_of_pmin_double_arrow_end_hc,
                                      LambdaType lambda) {
  auto pos_hc = N_minus1 - pos_mc;
  if (F_first_hc <= pos_hc && pos_hc <= F_last_hc)
    return ;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1;
    auto left_child_hc  = N_minus1 - left_child_mc;
    auto right_child_hc = left_child_hc - 1;
    bool is_right_child_valid = right_child_mc < heap_size && right_child_hc >= first_node_in_biheap_hc;
    if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc) {
      left_child_hc = ParentNotRoot<SizeType>(left_child_hc);
      if (F_first_hc <= left_child_hc && left_child_hc <= F_last_hc)
        left_child_hc = minh_parent_of_pmin_double_arrow_end_hc;
      left_child_mc = N_minus1 - left_child_hc;
    }
    if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc) {
      right_child_hc = ParentNotRoot<SizeType>(right_child_hc);
      if (F_first_hc <= right_child_hc && right_child_hc <= F_last_hc)
        right_child_hc = minh_parent_of_pmin_double_arrow_end_hc;
      right_child_mc = N_minus1 - right_child_hc;
    }
    is_right_child_valid = is_right_child_valid && right_child_hc >= first_node_in_biheap_hc;
    bool is_left_child_valid = left_child_hc >= first_node_in_biheap_hc;
    if (!is_left_child_valid && !is_right_child_valid)
      return ;
    auto left_it  = V + lambda(N, left_child_hc);
    auto right_it = V + lambda(N, right_child_hc);
    RAI larger_it;
    if (!is_left_child_valid || (is_right_child_valid && *left_it < *right_it)) {
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else { //Here, the left child is valid.
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < *larger_it) {
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    } else
      return ;
  }
  //At this point, it's not possible to be be simultaneously
  //  forbidden, in the pure min heap, and incident to the double arrow.
  if (pos_mc == pmin_double_arrow_end_hc) { //if and only if pos_hc = pmax_double_arrow_end_hc
    //If the other end of the double arrow is forbidden.
    if (F_first_hc <= pmin_double_arrow_end_hc && pmin_double_arrow_end_hc <= F_last_hc) {
      auto min_heap_parent_of_other_end_it = V + lambda(N, minh_parent_of_pmin_double_arrow_end_hc);
      //Perform one iteration of sifting up the min heap while skipping the
      // forbidden node.
      if (minh_parent_of_pmin_double_arrow_end_hc >= first_node_in_biheap_hc &&
          pos_value < *min_heap_parent_of_other_end_it) {
        std::iter_swap(pos_it, min_heap_parent_of_other_end_it);
        pos_it = min_heap_parent_of_other_end_it;
      } else
        return ;
      pos_hc = minh_parent_of_pmin_double_arrow_end_hc;
    }
  }
  //Start sifting up the min heap.
  //At this point pos is an In node.
  SizeType parent_hc = ParentNotRoot<SizeType>(pos_hc);
  if (pos_hc == 0 || parent_hc < first_node_in_biheap_hc)
    return ;
  do {
    auto parent_it = V + lambda(N, parent_hc);
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent_hc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 &&
     (parent_hc = ParentNotRoot<SizeType>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

//Assumes that N % 3 == 2.
// If F_last_hc < F_first_hc then it calls BiHeapify.
// If F_first_hc is not an In node then it is set to
//   2 * HeapSize(N) - N, the first min heap In node.
// If fust_last_hc is not an In node then it is set to
//   HeapSize(N) - 1, the last In node.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifyWithDoubleHeadedArrow(RAI V, SizeType N,
    SizeType F_first_hc, SizeType F_last_hc, LambdaType lambda) {
  if(N <= 2) {
    if (N == 2 && F_last_hc < F_first_hc) {
      RAI V_0 = V + lambda(N, 0);
      RAI V_1 = V + lambda(N, 1);
      if (*V_0 > *V_1)
        std::iter_swap(V_0, V_1);
    }
    return ;
  }
  SizeType N_minus1      = N - 1;
  SizeType heap_size     = HeapSize<SizeType>(N);
  SizeType first_in_node = N - heap_size;
  if (F_first_hc > first_in_node && F_last_hc < heap_size - 1) {
    //If we don't have to worry about skipping over either one of the end
    // nodes of the double headed arrow, then we may as well use the more
    // efficient FusedBiHeapify algorithm.
    FusedBiHeapifyIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N,
                                          F_first_hc, F_last_hc, lambda);
    return ;
  }
  if (F_last_hc < F_first_hc) { //Then all nodes are permitted.
    BiHeapify<RAI, SizeType, LambdaType>(V, N, lambda);
    return ;
  }
  if (F_first_hc < first_in_node)
    F_first_hc = first_in_node;
  if (F_last_hc >= heap_size)
    F_last_hc = heap_size - 1; //The last In node.
  SizeType last_node_in_biheap_hc  = heap_size - 2;
  SizeType first_node_in_biheap_hc = N_minus1 - last_node_in_biheap_hc;
  //To increase efficiency, we compute the following values once and pass them to the
  // two calls in the while loop. Any half descent optimizer will avoid actually
  // allocating additional space on the stack and copying these values there since
  // these functions are both inlined and templates.
  SizeType pmin_double_arrow_end_hc                = (N - 2) / 3;
  SizeType pmax_double_arrow_end_hc                = 2 * pmin_double_arrow_end_hc + 1;
  SizeType minh_parent_of_pmin_double_arrow_end_hc = ParentNotRoot<SizeType>(pmin_double_arrow_end_hc);
  SizeType maxh_parent_of_pmax_double_arrow_end_hc = N_minus1 - minh_parent_of_pmin_double_arrow_end_hc;
  while (first_node_in_biheap_hc > 0) {
    FusedBiHeapifySiftFromMinToMaxWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, --first_node_in_biheap_hc,
        last_node_in_biheap_hc, F_first_hc, F_last_hc,
        pmin_double_arrow_end_hc, pmax_double_arrow_end_hc,
        maxh_parent_of_pmax_double_arrow_end_hc, lambda);
    FusedBiHeapifySiftFromMaxToMinWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, N_minus1 - (++last_node_in_biheap_hc),
        first_node_in_biheap_hc, F_first_hc, F_last_hc,
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
inline void FusedBiHeapify(RAI V, SizeType N, SizeType F_first_hc,
                            SizeType F_last_hc, LambdaType lambda) {
  if (N % 3 != 2)
    FusedBiHeapifyIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, F_first_hc, F_last_hc, lambda);
  else
    FusedBiHeapifyWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, F_first_hc, F_last_hc, lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapify(RAI V, SizeType N,
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
  SizeType F_first_hc = first_in_node + (num_permitted_in_nodes + 1) / 2;
  SizeType F_last_hc  = (N - 1) - (first_in_node + (num_permitted_in_nodes / 2));
  FusedBiHeapify<RAI, SizeType, LambdaType>(V, N, F_first_hc,
                                             F_last_hc, lambda);
  return ;
}

//Specialize to the non-lambda case.
template<class RAI, typename SizeType = std::size_t>
inline void FusedBiHeapify(RAI V, SizeType N, SizeType num_permitted_in_nodes) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  FusedBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                       num_permitted_in_nodes, trivial_lambda);
  return ;
}

//Specialize to the non-lambda case.
template<class RAI, typename SizeType = std::size_t>
inline void FusedBiHeapify(RAI V, SizeType N,
                            SizeType F_first_hc, SizeType F_last_hc) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  FusedBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N, F_first_hc,
                                                 F_last_hc, trivial_lambda);
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapify(RAI V, SizeType N, LambdaType lambda) {
  FusedBiHeapify<RAI, SizeType, LambdaType>(V, N, 0, N, lambda);
  return ;
}

//Specialize to the non-lambda case.
template<class RAI, typename SizeType = std::size_t>
inline void FusedBiHeapify(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  FusedBiHeapify<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
  return ;
}

/*
 * ================== END: Definition of lambda version of AlmostBiheapify with some permitted In nodes and specializations ==========
 */


/*
 * ================== START: Definition of lambda version of IsFusedBiHeap ====================
 */

//Assumes that N > 2.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostTripleConditionAtInNode(RAI V, SizeType N, SizeType N_minus1, SizeType in_node_hc, LambdaType lambda) {
  SizeType min_heap_parent_hc = ParentNotRoot<SizeType>(in_node_hc);
  SizeType max_heap_parent_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - in_node_hc);
  return *(V + lambda(N, min_heap_parent_hc)) <= *(V + lambda(N, max_heap_parent_hc));
}

//Assumes that N > 3 and N mod 3 == 2.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsAlmostBiheapCheckAlmostQuadrupleCondition(RAI V, SizeType N, LambdaType lambda) {
  SizeType pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  SizeType parent_of_pure_min_heap_double_arrow_node_hc = ParentNotRoot<SizeType>(pure_min_heap_double_arrow_node_hc);
  return *(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc))
      <= *(V + lambda(N, (N - 1) - (parent_of_pure_min_heap_double_arrow_node_hc)));
}

//Assumes that N > 3 and N mod 3 = 2
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
void AlmostBiheapifyEnsureAlmostQuadrupleCondition(RAI V, SizeType N, LambdaType lambda) {
  SizeType pure_min_heap_double_arrow_node_hc = (N - 2) / 3;
  SizeType parent_of_pure_min_heap_double_arrow_node_hc = ParentNotRoot<SizeType>(pure_min_heap_double_arrow_node_hc);
  if (*(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc))
      > *(V + lambda(N, (N - 1) - parent_of_pure_min_heap_double_arrow_node_hc)))
    std::iter_swap(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc),
        V + lambda(N, (N - 1) - parent_of_pure_min_heap_double_arrow_node_hc));
  return ;
}

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define an almost BiHeap.
 */
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsFusedBiHeap(RAI V, SizeType N, LambdaType lambda) {
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
  SizeType N_minus1 = N - 1;
  {
    SizeType one_past_last_node = heap_size - is_N_mod_3_equal_to_2;
    for (SizeType in_hc = first_in_node_hc + is_N_mod_3_equal_to_2; in_hc < one_past_last_node; in_hc++) {
      if (!IsAlmostBiheapCheckAlmostTripleConditionAtInNode(V, N, N_minus1, in_hc, lambda))
        return false;
    }
  }
  //Check the min heap condition.
  SizeType pmin_node_of_double_arrow = (N - 2) / 3;
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
        if (!(is_N_mod_3_equal_to_2 && right_child == pmin_node_of_double_arrow))
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
      auto parent_value = *(V + lambda(N, N_minus1 - i));
      SizeType mirror_left_child_hc = N_minus1 - (right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, mirror_left_child_hc)))
        return false;
      //Check that the parent and right child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, (mirror_left_child_hc - 1)))) {
        if (!(is_N_mod_3_equal_to_2 && right_child == pmin_node_of_double_arrow))
          return false;
      }
    }
    //If the max heap's last non-In element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < first_in_node_hc
          && *(V + lambda(N, N_minus1 - i)) < *(V + lambda(N, N_minus1 - left_child)))
        return false;
    }
  }
  return true;
}

//Specialize to the non-lambda case.
template<class RAI, typename SizeType = std::size_t>
bool IsFusedBiHeap(RAI V, SizeType N) {
  auto trivial_lambda = [](SizeType local_N, SizeType i) -> SizeType {
    return i;
  };
  return IsFusedBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N, trivial_lambda);
}

/*
 * ================== END: Definition of lambda version of IsFusedBiHeap ====================
 */

/*
 * ================== START: Definition of lambda version of IsFusedBiHeap with some permitted In nodes ====================
 */

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsFusedBiHeap(RAI V, SizeType N, SizeType F_first_hc,
                    SizeType F_last_hc, LambdaType lambda) {
  if(N < 2)
    return true;
  if (F_first_hc > F_last_hc) //If all nodes are permitted then it should be a BiHeap.
    return IsBiHeap<RAI, SizeType, LambdaType>(V, N, lambda);
  if(N == 2 || (F_last_hc - F_first_hc <= 0))
    return true;
  if (!IsFusedBiHeap<RAI, SizeType, LambdaType>(V, N, lambda))
    return false ;
  SizeType N_minus1 = N - 1;
  if (N % 3 == 2) { //Then make sure that the double arrow satisfies the necessary conditions.
    SizeType pmin_double_arrow_node_hc = (N - 2) / 3;
    SizeType pmax_double_arrow_node_hc = (2 * N - 1) / 3;
    SizeType minH_parent_of_pmin_double_arrow_node_hc = ParentNotRoot<SizeType>(pmin_double_arrow_node_hc);
    SizeType maxH_parent_of_pmax_double_arrow_node_hc = N_minus1 - minH_parent_of_pmin_double_arrow_node_hc;
    RAI pmin_double_arrow_node_it = V + lambda(N, pmin_double_arrow_node_hc);
    RAI pmax_double_arrow_node_it = V + lambda(N, pmax_double_arrow_node_hc);
    RAI minH_parent_of_pmin_double_arrow_node_it = V + lambda(N, minH_parent_of_pmin_double_arrow_node_hc);
    RAI maxH_parent_of_pmax_double_arrow_node_it = V + lambda(N, maxH_parent_of_pmax_double_arrow_node_hc);
    bool is_pmin_double_arrow_node_forbidden =
        F_first_hc <= pmin_double_arrow_node_hc
        && pmin_double_arrow_node_hc <= F_last_hc;
    bool is_pmax_double_arrow_node_forbidden =
        F_first_hc <= pmax_double_arrow_node_hc
        && pmax_double_arrow_node_hc <= F_last_hc;
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
  SizeType last_in_node_to_check  = N_minus1 - first_in_node_to_check;
  for (SizeType pos_hc = first_in_node_to_check; pos_hc < F_first_hc; pos_hc++) {
    SizeType minh_parent_of_pos_hc = ParentNotRoot<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - pos_hc);
    auto pos_value = *(V + lambda(N, pos_hc));
    if (pos_value < *(V + lambda(N, minh_parent_of_pos_hc)) ||
        pos_value > *(V + lambda(N, maxh_parent_of_pos_hc))) {
      return false;
    }
  }
  for (SizeType pos_hc = F_last_hc + 1; pos_hc <= last_in_node_to_check; pos_hc++) {
    SizeType minh_parent_of_pos_hc = ParentNotRoot<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_hc = N_minus1 - ParentNotRoot<SizeType>(N_minus1 - pos_hc);
    auto pos_value = *(V + lambda(N, pos_hc));
    if (pos_value < *(V + lambda(N, minh_parent_of_pos_hc)) ||
        pos_value > *(V + lambda(N, maxh_parent_of_pos_hc))) {
      return false;
    }
  }
  return true;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsFusedBiHeap(RAI V, SizeType N, SizeType num_permitted_in_nodes, LambdaType lambda) {
  if(N < 2 || (N == 2 && num_permitted_in_nodes <= 1))
    return true;
  if (num_permitted_in_nodes <= 0)
    return IsFusedBiHeap<RAI, SizeType, LambdaType>(V, N, lambda);
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
  SizeType F_first_hc = first_in_node + (num_permitted_in_nodes + 1) / 2;
  SizeType F_last_hc = (N - 1) - (first_in_node + (num_permitted_in_nodes / 2));
  return IsFusedBiHeap<RAI, SizeType, LambdaType>(V, N, F_first_hc,
                                                   F_last_hc, lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool IsFusedBiHeap(RAI V, SizeType N,  SizeType num_permitted_in_nodes) {
  auto trivial_lambda = [](SizeType N, SizeType i) -> SizeType {
    return i;
  };
  return IsFusedBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                       num_permitted_in_nodes, trivial_lambda);
}

template<class RAI, typename SizeType = std::size_t>
bool IsFusedBiHeap(RAI V, SizeType N, SizeType F_first_hc, SizeType F_last_hc) {
  auto trivial_lambda = [](SizeType N, SizeType i) -> SizeType {
    return i;
  };
  return IsFusedBiHeap<RAI, SizeType, decltype(trivial_lambda)>(V, N,
                                                F_first_hc, F_last_hc, trivial_lambda);
}

/*
 * ================== END: Definition of lambda version of IsFusedBiHeap with some permitted In nodes ====================
 */


/*
 * ================== START: Definition of FusedBiHeapifyJumpMiddle ====================
 */


/* In short, given a list of distance nodes: V, ..., V + (distance - 1)
 *  this takes the first h := num_nodes_to_biheapify / 2 nodes
 *  and the last h nodes and makes them into a BiHeap.
 * Specifically, if num_nodes_to_biheapify is even then it makes the nodes
 *  V, ..., V + (h-1), V + Flip(h-1), ..., V + (distance - 1)
 *  into a BiHeap while if num_nodes_to_biheapify is odd then it makes
 *  V, ..., V + (h-1), V + (distance / 2), V + Flip(h-1), ..., V + (distance - 1)
 *  into a BiHeap (with minimum at V and maximum at V + (distance - 1)).
 * If num_nodes_to_biheapify is odd then it is assumed that distance is odd
 *  (since otherwise there is no natural way to distinguish a middle node)
 *  and the middle node V + (distance / 2) will become the
 *  middle node this BiHeap of size num_nodes_to_biheapify.
 * Assumes that num_nodes_to_biheapify <= distance.
 */
template<class RAI, typename SizeType = std::size_t>
inline void FusedBiHeapifyJumpMiddle(RAI V, SizeType distance, SizeType num_nodes_to_biheapify) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (distance == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    FusedBiHeapify<RAI, SizeType>(V, distance);
    return ;
  }
  SizeType half_of_num_nodes = num_nodes_to_biheapify / 2;
  SizeType distance_minus_num_nodes = distance - num_nodes_to_biheapify;
  if (num_nodes_to_biheapify % 2 == 1 && distance % 2 == 1) {
    SizeType half_of_distance = distance / 2;
    //Note that local_N will equal num_nodes_to_biheapify in the subsequent call to BiHeapify.
    auto index_lambda = [half_of_distance, half_of_num_nodes, distance_minus_num_nodes]
         (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return i;
      else if (i != half_of_num_nodes) //So i > half_of_num_nodes
        return distance_minus_num_nodes + i; // = (distance - 1) - [(local_N - 1) - i]
      else //i == half_of_num_nodes
        return half_of_distance;
    };
    FusedBiHeapify<RAI, SizeType>(V, num_nodes_to_biheapify, index_lambda);
  } else {
    auto index_lambda = [half_of_num_nodes, distance_minus_num_nodes]
                   (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return i;
      else
        return distance_minus_num_nodes + i; // = (distance - 1) - [(local_N - 1) - i]
    };
    FusedBiHeapify<RAI, SizeType>(V, num_nodes_to_biheapify, index_lambda);
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void FusedBiHeapifyJumpMiddle(RAI V, SizeType distance, SizeType num_nodes_to_biheapify, LambdaType lambda) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (distance == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    FusedBiHeapify<RAI, SizeType, LambdaType>(V, distance, lambda);
    return ;
  }
  SizeType half_of_num_nodes = num_nodes_to_biheapify / 2;
  SizeType distance_minus_num_nodes = distance - num_nodes_to_biheapify;
  if (num_nodes_to_biheapify % 2 == 1 && distance % 2 == 1) {
    SizeType lambda_of_half_of_distance = lambda(num_nodes_to_biheapify, distance / 2);
    //Note that local_N will equal num_nodes_to_biheapify in the subsequent call to BiHeapify.
    auto index_lambda =
        [half_of_num_nodes, distance_minus_num_nodes, num_nodes_to_biheapify, lambda_of_half_of_distance, lambda]
         (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return lambda(num_nodes_to_biheapify, i);
      else if (i != half_of_num_nodes) //So i > half_of_num_nodes
        return lambda(num_nodes_to_biheapify, distance_minus_num_nodes + i); // = (distance - 1) - [(local_N - 1) - i]
      else //i == half_of_num_nodes
        return lambda_of_half_of_distance;
    };
    FusedBiHeapify<RAI, SizeType, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  } else {
    auto index_lambda =
        [half_of_num_nodes, distance_minus_num_nodes, num_nodes_to_biheapify, lambda]
         (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return lambda(num_nodes_to_biheapify, i);
      else
        return lambda(num_nodes_to_biheapify, distance_minus_num_nodes + i); // = (distance - 1) - [(local_N - 1) - i]
    };
    FusedBiHeapify<RAI, SizeType, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  }
  return ;
}

/*
 * ================== END: Definition of FusedBiHeapifyJumpMiddle ====================
 */

/*
 * ================== START: Definition of FusedBiHeapifySift ====================
 */

#include "biheap_sift.h"

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void FusedBiHeapifySiftWithDoubleHeadedArrow(RAI V, SizeType N, SizeType pos_hc,
                          SizeType F_first_hc, SizeType F_last_hc, LambdaType lambda) {
  SizeType N_minus1       = N - 1;
  SizeType heap_size      = HeapSize<SizeType>(N);
  SizeType first_in_node  = N - heap_size;
  SizeType pmin_double_arrow_end_hc                = (N - 2) / 3;
  SizeType pmax_double_arrow_end_hc                = 2 * pmin_double_arrow_end_hc + 1;
  SizeType minh_parent_of_pmin_double_arrow_end_hc = ParentNotRoot<SizeType>(pmin_double_arrow_end_hc);
  SizeType maxh_parent_of_pmax_double_arrow_end_hc = N_minus1 - minh_parent_of_pmin_double_arrow_end_hc;
  SizeType pos_mc         = N_minus1 - pos_hc;
  auto pos_it             = V + lambda(N, pos_hc);
  auto pos_value          = *pos_it;
  bool is_node_in_min_heap = pos_hc < heap_size;
  bool is_node_in_max_heap = pos_mc < heap_size;
  bool is_pmin_double_arrow_end_in_F = (F_first_hc <= pmin_double_arrow_end_hc) && (pmin_double_arrow_end_hc <= F_last_hc);
  bool is_pmax_double_arrow_end_in_F = (F_first_hc <= pmax_double_arrow_end_hc) && (pmax_double_arrow_end_hc <= F_last_hc);
  if ((N % 3 == 2) && (pos_hc == pmin_double_arrow_end_hc || pos_hc == pmax_double_arrow_end_hc) //if pos_hc is the end of a double arrow.
      && (is_pmin_double_arrow_end_in_F || is_pmax_double_arrow_end_in_F)) { //and at least one end of a double arrow is in F
    auto minh_parent_of_pmax_double_arrow_end_it = V + lambda(N, minh_parent_of_pmin_double_arrow_end_hc);
    auto maxh_parent_of_pmax_double_arrow_end_it = V + lambda(N, maxh_parent_of_pmax_double_arrow_end_hc);
    if (pos_hc == pmin_double_arrow_end_hc && is_pmax_double_arrow_end_in_F) { //If !is_pmax_double_arrow_end_in_F then you can just sift up.
      if (pos_value > *maxh_parent_of_pmax_double_arrow_end_it) {
        std::iter_swap(pos_it, maxh_parent_of_pmax_double_arrow_end_it);
        SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      } else if (pos_value < *minh_parent_of_pmax_double_arrow_end_it) {
        std::iter_swap(pos_it, minh_parent_of_pmax_double_arrow_end_it);
        SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      }
      return ; //Return should be inside this if statement.
    } else if (pos_hc == pmax_double_arrow_end_hc && is_pmin_double_arrow_end_in_F) { //If !is_pmin_double_arrow_end_in_F then you can just sift up.
      if (pos_value > *maxh_parent_of_pmax_double_arrow_end_it) {
        std::iter_swap(pos_it, maxh_parent_of_pmax_double_arrow_end_it);
        SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      } else if (pos_value < *minh_parent_of_pmax_double_arrow_end_it) {
        std::iter_swap(pos_it, minh_parent_of_pmax_double_arrow_end_it);
        SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      }
      return ; //Return should be inside this if statement.
    }
  }
  if (is_node_in_min_heap && is_node_in_max_heap) {
    SizeType minh_parent_of_pos_hc = Parent<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_mc = Parent<SizeType>(pos_mc);
    SizeType maxh_parent_of_pos_hc = N_minus1 - maxh_parent_of_pos_mc;
    auto minh_parent_of_pos_it     = V + lambda(N, minh_parent_of_pos_hc);
    auto maxh_parent_of_pos_it     = V + lambda(N, maxh_parent_of_pos_hc);
    if (pos_value > *maxh_parent_of_pos_it) {
      std::iter_swap(pos_it, maxh_parent_of_pos_it);
      SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, N_minus1 - maxh_parent_of_pos_hc, lambda);
    } else if (pos_value < *minh_parent_of_pos_it) {
      std::iter_swap(pos_it, minh_parent_of_pos_it);
      SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pos_hc, lambda);
    }
    return ;
  }
  //At this point, pos is not an In node.
  if (!is_node_in_max_heap) { //Then it's in the pure min heap and not in the max heap.
    if (pos_hc == 0 ||
        !(pos_value < *(V + lambda(N, ParentNotRoot<SizeType>(pos_hc))))) {
      FusedBiHeapifySiftFromMinToMaxWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
                                    heap_size, first_in_node, pos_hc, N_minus1, F_first_hc, F_last_hc,
                                    pmin_double_arrow_end_hc, pmax_double_arrow_end_hc,
                                    maxh_parent_of_pmax_double_arrow_end_hc, lambda);
    } else {
      SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, lambda);
    }
  } else { //Then it's in the pure max heap and not in the min heap.
    if (pos_mc == 0 ||
        !(*(V + lambda(N, N_minus1 - ParentNotRoot<SizeType>(pos_mc))) < pos_value)) {
      FusedBiHeapifySiftFromMaxToMinWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
                                    heap_size, first_in_node, pos_mc, 0, F_first_hc, F_last_hc,
                                    pmin_double_arrow_end_hc, minh_parent_of_pmin_double_arrow_end_hc,
                                    lambda);
    } else {
      SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_mc, lambda);
    }
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void FusedBiHeapifySiftIgnoreDoubledHeadedArrow(RAI V, SizeType N, SizeType pos_hc,
                          SizeType F_first_hc, SizeType F_last_hc, LambdaType lambda) {
  SizeType N_minus1       = N - 1;
  SizeType heap_size      = HeapSize<SizeType>(N);
  SizeType first_in_node  = N - heap_size;
  SizeType pos_mc         = N_minus1 - pos_hc;
  auto pos_it                    = V + lambda(N, pos_hc);
  auto pos_value                 = *pos_it;
  bool is_node_in_min_heap = pos_hc < heap_size;
  bool is_node_in_max_heap = pos_mc < heap_size;
  if (is_node_in_min_heap && is_node_in_max_heap) {
    SizeType minh_parent_of_pos_hc = Parent<SizeType>(pos_hc);
    SizeType maxh_parent_of_pos_mc = Parent<SizeType>(pos_mc);
    SizeType maxh_parent_of_pos_hc = N_minus1 - maxh_parent_of_pos_mc;
    auto minh_parent_of_pos_it     = V + lambda(N, minh_parent_of_pos_hc);
    auto maxh_parent_of_pos_it     = V + lambda(N, maxh_parent_of_pos_hc);
    if (pos_value > *maxh_parent_of_pos_it) {
      std::iter_swap(pos_it, maxh_parent_of_pos_it);
      SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, N_minus1 - maxh_parent_of_pos_hc, lambda);
    } else if (pos_value < *minh_parent_of_pos_it) {
      std::iter_swap(pos_it, minh_parent_of_pos_it);
      SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pos_hc, lambda);
    }
    return ;
  }
  //At this point, pos is not an In node.
  if (!is_node_in_max_heap) { //Then it's in the pure min heap and not in the max heap.
    if (pos_hc == 0 ||
        !(pos_value < *(V + lambda(N, ParentNotRoot<SizeType>(pos_hc))))) {
      FusedBiHeapifySiftFromMinToMaxIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
                                    heap_size, first_in_node, pos_hc, N_minus1, F_first_hc, F_last_hc,
                                    lambda);
    } else {
      SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_hc, lambda);
    }
  } else { //Then it's in the pure max heap and not in the min heap.
    if (pos_mc == 0 ||
        !(*(V + lambda(N, N_minus1 - ParentNotRoot<SizeType>(pos_mc))) < pos_value)) {
      FusedBiHeapifySiftFromMaxToMinIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, N_minus1,
                                    heap_size, first_in_node, pos_mc, 0, F_first_hc, F_last_hc,
                                    lambda);
    } else {
      SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, pos_mc, lambda);
    }
  }
  return ;
}


template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void FusedBiHeapifySift(RAI V, SizeType N, SizeType pos_hc,
                          SizeType F_first_hc, SizeType F_last_hc, LambdaType lambda) {
  if (N % 3 != 2) {
    FusedBiHeapifySiftIgnoreDoubledHeadedArrow<RAI, SizeType, LambdaType>(V, N, pos_hc, F_first_hc, F_last_hc, lambda);
  } else {
    FusedBiHeapifySiftWithDoubleHeadedArrow<RAI, SizeType, LambdaType>(V, N, pos_hc, F_first_hc, F_last_hc, lambda);
  }
}

template<class RAI, typename SizeType = std::size_t>
inline void FusedBiHeapifySift(RAI V, SizeType N, SizeType pos_hc,
                                SizeType F_first_hc, SizeType F_last_hc) {
  auto trivial_lambda = [](SizeType N, SizeType i) -> SizeType {
    return i;
  };
  FusedBiHeapifySift<RAI, SizeType, decltype(trivial_lambda)>(V, N, pos_hc,
                                              F_first_hc, F_last_hc, trivial_lambda);
}
/*
 * ================== END: Definition of FusedBiHeapifySift ====================
 */

//#undef FLIP

#endif /* ALMOST_BIHEAPIFY_H_ */
