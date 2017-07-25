/*
 * biheap_sift.h
 *
 *  Created on: Jul 25, 2017
 *      Author: Matthew Gregory Krupa
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

#include "biheap_common.h"

/* Note that FlipCo(coord1) >= coord2 iff coord1 <= FlipCo(coord2)
 * (ditto for >, <=, and <).
 */
#define FLIP_COORDINATE(a) (total_num_nodes - 1 - (a))

namespace {

template<class RAI>
inline void SiftUpMaxHeapMC(RAI first, size_type total_num_nodes,
                            size_type pos_mc) {
  size_type parent;
  auto pos_it = first + FLIP_COORDINATE(pos_mc);
  while (pos_mc > 0) {
    parent = GetParentInBiheapNotRoot(pos_mc);
    auto parent_it = first + FLIP_COORDINATE(parent);
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

template<class RAI>
inline void SiftUpMaxHeapHC(RAI first, size_type total_num_nodes,
                            size_type pos_hc) {
  SiftUpMaxHeapMC(first, total_num_nodes, FLIP_COORDINATE(pos_hc));
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapHC(RAI first, size_type pos_hc) {
  size_type parent;
  auto pos_it = first + pos_hc;
  while (pos_hc > 0) {
    parent = GetParentInBiheapNotRoot(pos_hc);
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

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapMC(RAI first, size_type total_num_nodes,
                            size_type pos_mc) {
  SiftUpMinHeapHC(first, FLIP_COORDINATE(pos_mc));
}


/* Assumes that total_num_nodes is even and that the node pos_mc belongs to the
 *  min heap.
 */
template<class RAI>
inline void SiftFromMinToMaxEven(RAI first, size_type total_num_nodes,
                          size_type num_nodes_in_heap,
                          size_type first_node_in_mirror_heap,
                          size_type pos_hc) {
  while (pos_hc < first_node_in_mirror_heap) {
    if (pos_hc >= total_num_nodes) //If node pos_hc is not in the biheap.
      //Note that the inequality pos_hc >= smallest_node_in_biheap_hc necessarily holds.
      // which is why it suffices to check the above inequality.
      return ;

    auto left_child     = GetLeftChildInBiheap(pos_hc);
    auto right_child    = left_child + 1;
    bool is_right_child_valid = (right_child < total_num_nodes) && //Is the node in the biheap?
                                (right_child < num_nodes_in_heap);            //Is the node in the min heap?
    //Note that the equivalent of is_right_child_valid for the left child, which is
    // (left_child <= largest_node_in_biheap_hc) && left_child  < num_nodes_in_heap),
    // is always true when total_num_nodes is even.
    auto left_it  = first + left_child;
    auto right_it = first + right_child;
    auto pos_it   = first + pos_hc;

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
  SiftUpMaxHeapHC(first, total_num_nodes, pos_hc);
  return ;
}

/* Assumes that total_num_nodes is even and that the node pos_mc belongs to the
 *  max heap.
 */
template<class RAI>
inline void SiftFromMaxToMinEven(RAI first, size_type total_num_nodes,
                          size_type num_nodes_in_heap,
                          size_type first_node_in_mirror_heap,
                          size_type pos_mc) {
  auto pos_hc = FLIP_COORDINATE(pos_mc);
  //While node pos_mc is NOT in the min heap.
  while (pos_mc < first_node_in_mirror_heap) {
    if (pos_hc < 0) //If the node is not in the biheap.
      return ;

    auto left_child_mc  = GetLeftChildInBiheap(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP_COORDINATE(left_child_mc);
    auto right_child_hc = left_child_hc - 1; //= FLIP_COORDINATE(right_child_mc)
    auto pos_it   = first + pos_hc;
    auto left_it  = first + left_child_hc;
    auto right_it = first + right_child_hc;

    //Is the node in the max heap?
    bool is_right_child_valid = right_child_mc < num_nodes_in_heap;
    //Note that the equivalent of is_right_child_valid for the left child, which is
    // (left_child <= largest_node_in_biheap_hc) && left_child  < num_nodes_in_heap)
    // is always true when total_num_nodes is even.

    RAI larger_it;
    if (is_right_child_valid && *right_it > *left_it) {
      larger_it = right_it;
      pos_hc  = right_child_hc;
      pos_mc  = right_child_mc;
    } else {
      larger_it = left_it;
      pos_hc  = left_child_hc;
      pos_mc  = left_child_mc;
    }

    if (*pos_it < *larger_it)
      std::iter_swap(pos_it, larger_it);
    else
      return ;
  }
  SiftUpMinHeapHC(first, pos_hc);
  return ;
}


/* Assumes that total_num_nodes is odd, that the node pos_mc belongs to the
 *  min heap, and that pos_hc >= smallest_node_in_biheap_hc.
 */
template<class RAI>
inline void SiftFromMinToMaxOdd(RAI first, size_type total_num_nodes,
                         size_type num_nodes_in_heap,
                         size_type pos_hc) {
  while (pos_hc <= total_num_nodes / 2) {
    if (pos_hc >= total_num_nodes) //If the node is not in the biheap.
      return ;

    auto left_child     = GetLeftChildInBiheap(pos_hc);
    auto right_child    = left_child + 1;
    if (left_child >= num_nodes_in_heap) //If the node is not in the min heap.
      break ;
    bool is_right_child_valid = right_child < num_nodes_in_heap;

    auto left_it  = first + (left_child);
    auto right_it = first + (right_child);
    auto pos_it   = first + (pos_hc);

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
  SiftUpMaxHeapHC(first, total_num_nodes, pos_hc);
  return ;
}

/* Assumes that total_num_nodes is odd and that the node pos_mc belongs to the
 *  max heap.
 */
template<class RAI>
inline void SiftFromMaxToMinOdd(RAI first, size_type total_num_nodes,
                         size_type num_nodes_in_heap,
                         size_type pos_mc) {
  auto pos_hc = FLIP_COORDINATE(pos_mc);
  while (pos_mc <= total_num_nodes / 2) {
    if (pos_hc < 0) //If the node is not in the biheap.
      return ;

    auto left_child_mc  = GetLeftChildInBiheap(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP_COORDINATE(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP_COORDINATE(right_child_mc)
    auto pos_it   = first + pos_hc;
    auto left_it  = first + left_child_hc;
    auto right_it = first + right_child_hc;

    if (left_child_mc >= num_nodes_in_heap)//If the node is not in the max heap.
      break ;
    bool is_right_child_valid = right_child_mc < num_nodes_in_heap;

    RAI larger_it;
    if (is_right_child_valid && *right_it > *left_it) {
      larger_it = right_it;
      pos_hc  = right_child_hc;
      pos_mc  = right_child_mc;
    } else {
      larger_it = left_it;
      pos_hc  = left_child_hc;
      pos_mc  = left_child_mc;
    }

    if (*pos_it < *larger_it)
      std::iter_swap(pos_it, larger_it);
    else
      return ;
  }
  SiftUpMinHeapHC(first, pos_hc);
  return ;
}

} //End anonymous namespace

template<class RAI>
void BiHeapSiftEven(RAI first, size_type total_num_nodes, size_type pos_hc) {
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  auto first_node_in_mirror_heap  = total_num_nodes - num_nodes_in_heap;

  if (pos_hc < num_nodes_in_heap) { //If the node is in the min heap.
    auto parent_hc = GetParentInBiheapZero(pos_hc);
    if (pos_hc != 0 && *(first + parent_hc) > *(first + pos_hc)) {
      SiftUpMinHeapHC(first, pos_hc);
    } else {
      SiftFromMinToMaxEven<RAI>(first, total_num_nodes, num_nodes_in_heap,
                                first_node_in_mirror_heap, pos_hc);
    }
    return ;
  }

  auto pos_mc = FLIP_COORDINATE(pos_hc);
  auto parent_mc = GetParentInBiheapZero(pos_mc);
  if (pos_mc < num_nodes_in_heap) { //If the node is in the max heap.
    auto parent_hc = GetParentInBiheapZero(parent_mc);
    if (pos_hc != 0 && *(first + parent_hc) < *(first + pos_hc)) {
      SiftUpMaxHeapMC(first, total_num_nodes, pos_mc);
    } else {
      SiftFromMaxToMinEven<RAI>(first, total_num_nodes, num_nodes_in_heap,
                                first_node_in_mirror_heap, pos_mc);
    }
  }
  return ;
}

template<class RAI>
void BiHeapSiftOdd(RAI first, size_type total_num_nodes, size_type pos_hc) {
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);

  if (pos_hc < num_nodes_in_heap) { //If the node is in the min heap.
    auto parent_hc = GetParentInBiheapZero(pos_hc);
    if (pos_hc != 0 && *(first + parent_hc) > *(first + pos_hc)) {
      SiftUpMinHeapHC(first, pos_hc);
    } else {
      SiftFromMinToMaxOdd<RAI>(first, total_num_nodes, num_nodes_in_heap,
                               pos_hc);
    }
    return ;
  }

  auto pos_mc = FLIP_COORDINATE(pos_hc);
  auto parent_mc = GetParentInBiheapZero(pos_mc);
  if (pos_mc < num_nodes_in_heap) { //If the node is in the max heap.
    auto parent_hc = GetParentInBiheapZero(parent_mc);
    if (pos_hc != 0 && *(first + parent_hc) < *(first + pos_hc)) {
      SiftUpMaxHeapMC(first, total_num_nodes, pos_mc);
    } else {
      SiftFromMaxToMinOdd<RAI>(first, total_num_nodes, num_nodes_in_heap,
                               pos_mc);
    }
  }
  return ;
}

template<class RAI>
void BiHeapSift(RAI first, size_type total_num_nodes, size_type pos_hc) {
  if (total_num_nodes % 2 == 0) {
    BiHeapSiftEven(first, total_num_nodes, pos_hc);
  } else {
    BiHeapSiftOdd(first, total_num_nodes, pos_hc);
  }
  return ;
}

#undef FLIP_COORDINATE

#endif /* BIHEAP_SIFT_H_ */
