/*
 * biheap_sift.h
 *
 *  Created on: Jul 25, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Copyright Matthew Gregory Krupa
 *
 * This header file defines the BiHeapSift() function. Given a biheap on
 *  total_num_nodes nodes defined by the iterator first (so that *first,
 *  *(first + 1), ..., *(first + (total_num_nodes - 1) are the nodes' values),
 *  and given 0 <= pos_hc < total_num_nodes, if *(first + pos_hc) is
 *  changed then these nodes may no longer form a biheap.
 * Calling BiHeapSift(first, total_num_nodes, pos_hc) will make it into a
 *  biheap once again and it will do this in O(log(total_num_nodes)) time.
 */

#ifndef BIHEAP_SIFT_H_
#define BIHEAP_SIFT_H_

#include <algorithm>

#include "biheapify.h"

#define FLIP(a) (total_num_nodes - 1 - (a))

template<class RAI, typename size_type = std::size_t>
inline void SiftUpMaxHeapMC(RAI first, size_type total_num_nodes,
                            size_type pos_mc) {
  size_type parent;
  auto pos_it = first + FLIP(pos_mc);
  while (pos_mc > 0) {
    parent = Parent<size_type>(pos_mc);
    auto parent_it = first + FLIP(parent);
    if (*pos_it > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void SiftUpMaxHeapHC(RAI first, size_type total_num_nodes,
                            size_type pos_hc) {
  SiftUpMaxHeapMC<RAI, size_type>(first, total_num_nodes, FLIP(pos_hc));
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename size_type = std::size_t>
inline void SiftUpMinHeapHC(RAI first, size_type pos_hc) {
  size_type parent;
  auto pos_it = first + pos_hc;
  while (pos_hc > 0) {
    parent = Parent<size_type>(pos_hc);
    auto parent_it = first + parent;
    if (*pos_it < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void SiftFromMinToMax(RAI first, size_type total_num_nodes,
                                size_type num_nodes_in_heap,
                                size_type first_node_in_mirror_heap,
                                size_type pos_hc) {
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child  = LeftChild<size_type>(pos_hc);
    auto right_child = left_child + 1;
    auto left_it  = first + left_child;
    auto right_it = first + right_child;
    auto pos_it   = first + pos_hc;

    //assert((left_child  < num_nodes_in_heap) && (left_child  < total_num_nodes) && (right_child < total_num_nodes));
    bool is_right_child_valid = right_child < num_nodes_in_heap;
    RAI smaller_it;
    if (is_right_child_valid && *right_it < *left_it) {
      smaller_it = right_it;
      pos_hc     = right_child;
    } else {
      smaller_it = left_it;
      pos_hc     = left_child;
    }
    if (*pos_it > *smaller_it)
      std::iter_swap(pos_it, smaller_it);
    else
      return ;
  }
  SiftUpMaxHeapHC<RAI, size_type>(first, total_num_nodes, pos_hc);
  return ;
}

/* Assumes that total_num_nodes is odd, that the node pos_mc belongs to the
 *  max heap, and that FLIP(pos_mc) >= smallest_node_in_biheap_hc.
 */
template<class RAI, typename size_type = std::size_t>
inline void SiftFromMaxToMin(RAI first, size_type total_num_nodes,
                                size_type num_nodes_in_heap,
                                size_type first_node_in_mirror_heap,
                                size_type pos_mc) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = first + pos_hc;
    auto left_it        = first + left_child_hc;
    auto right_it       = first + right_child_hc;

    //assert((left_child_mc < num_nodes_in_heap) && (left_child_hc >= 0) && (right_child_hc >= 0));
    bool is_right_child_valid = right_child_mc < num_nodes_in_heap;
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
  SiftUpMinHeapHC<RAI, size_type>(first, pos_hc);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapSift(RAI first, size_type total_num_nodes, size_type pos_hc) {
  auto num_nodes_in_heap = HeapSize(total_num_nodes);
  auto first_node_in_mirror_heap  = total_num_nodes - num_nodes_in_heap;

  if (pos_hc < num_nodes_in_heap) { //If the node is in the min heap.
    auto parent_hc = Parent<size_type>(pos_hc);
    if (pos_hc != 0 && *(first + parent_hc) > *(first + pos_hc)) {
      SiftUpMinHeapHC<RAI, size_type>(first, pos_hc);
    } else {
      SiftFromMinToMax<RAI, size_type>(first, total_num_nodes, num_nodes_in_heap,
                                       first_node_in_mirror_heap, pos_hc);
    }
    return ;
  }
  //Else the node is in the max heap.
  auto pos_mc = FLIP(pos_hc);
  auto parent_mc = Parent<size_type>(pos_mc);
  auto parent_hc = Parent<size_type>(parent_mc);
  if (pos_mc != 0 && *(first + parent_hc) < *(first + pos_hc))
    SiftUpMaxHeapMC<RAI, size_type>(first, total_num_nodes, pos_mc);
  SiftFromMaxToMin<RAI, size_type>(first, total_num_nodes, num_nodes_in_heap,
                                   first_node_in_mirror_heap, pos_mc);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapSiftMC(RAI first, size_type total_num_nodes, size_type pos_mc) {
  BiHeapSift<RAI, size_type>(first, total_num_nodes, FLIP(pos_mc));
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapSift(RAI first, RAI one_past_last, RAI pos_hc) {
  BiHeapSift<RAI, size_type>(first, std::distance(first, one_past_last),
                             std::distance(first, pos_hc));
}

#undef FLIP

#endif /* BIHEAP_SIFT_H_ */
