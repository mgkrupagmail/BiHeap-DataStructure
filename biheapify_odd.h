/*
 * biheapify_odd.h
 *
 *  Created on: Jun 22, 2017
 *      Author: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 */

#ifndef BIHEAPIFY_ODD_H_
#define BIHEAPIFY_ODD_H_

#include <algorithm>

#include "biheap_common.h"

/* Note that FlipCo(coord1) >= coord2 iff coord1 <= FlipCo(coord2)
 * (ditto for >, <=, and <).
 */
#define FLIP(a) (total_num_nodes - 1 - (a))

/* Assumes that total_num_nodes is odd, that the node pos_mc belongs to the
 *  min heap, and that pos_hc <= largest_node_in_biheap_hc.
 */
template<class RAI, typename size_type = std::size_t>
inline void SiftFromMinToMaxOdd(RAI first, size_type total_num_nodes,
                         size_type num_nodes_in_heap,
                         size_type pos_hc,
                         size_type largest_node_in_biheap_hc) {
  while (pos_hc <= total_num_nodes / 2) {
    auto left_child     = LeftChild<size_type>(pos_hc);
    auto right_child    = left_child + 1;

    //Note that removing the conditions:
    // left_child  <= largest_node_in_biheap_hc  and
    // right_child <= largest_node_in_biheap_hc
    // from the below definitions does not affect the correctness of the
    // biheapify algorithm.

    //Is the node in the biheap? && Is the node in the min heap?
    bool is_left_child_valid  = //(left_child <= largest_node_in_biheap_hc) &&
                                left_child  < num_nodes_in_heap;
    if (!is_left_child_valid)
      break ;
    //At this point left_child < num_nodes_in_heap.
    //If left_child > largest_node_in_biheap_hc (i.e. the left_child is not
    // in the biheap) then in the next iteration of the loop, pos_hc will be
    // > total_num_nodes / 2 and execution will leave the loop. Thus in this
    // case there is at most one additional iteration so that the omission of
    // checking whether or not left_child <= largest_node_in_biheap_hc) does
    // NOT affect this algorithm's O(n) complexity. In addition, this does
    // not affect this algorithm's correctness. The same argument applies to
    // right_child.
    
    bool is_right_child_valid = //(right_child <= largest_node_in_biheap_hc) &&
                                right_child < num_nodes_in_heap;
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
  SiftUpMaxHeapHC<RAI, size_type>(first, total_num_nodes, pos_hc,
                       FLIP(largest_node_in_biheap_hc));
  return ;
}

/* Assumes that total_num_nodes is odd, that the node pos_mc belongs to the
 *  max heap, and that FLIP(pos_mc) >= smallest_node_in_biheap_hc.
 */
template<class RAI, typename size_type = std::size_t>
inline void SiftFromMaxToMinOdd(RAI first, size_type total_num_nodes,
                         size_type num_nodes_in_heap,
                         size_type pos_mc,
                         size_type smallest_node_in_biheap_hc) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc <= total_num_nodes / 2) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it   = first + pos_hc;
    auto left_it  = first + left_child_hc;
    auto right_it = first + right_child_hc;

    //Note that removing the condition:
    // right_child <= largest_node_in_biheap_hc
    // from the definition of is_left_child_valid does not affect the
    // correctness of the biheapify algorithm, but the same is not true of
    // removing the condition: left_child_hc >= smallest_node_in_biheap_hc.

    //Is the node in the biheap? && Is the node in the min heap?
    //Note that unlike in SiftFromMinToMaxOdd() for left_child we must must
    // check that left_child_hc >= smallest_node_in_biheap_hc.
    //This is likely due to SiftFromMaxToMinOdd() being called after
    // SiftFromMinToMaxOdd() in BiHeapifyOdd().
    bool is_left_child_valid  = (left_child_hc >= smallest_node_in_biheap_hc) &&
                                (left_child_mc < num_nodes_in_heap);
    if (!is_left_child_valid)
      break ;
    bool is_right_child_valid =//right_child_hc >= smallest_node_in_biheap_hc &&
                                right_child_mc < num_nodes_in_heap;

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
  SiftUpMinHeapHC<RAI, size_type>(first, pos_hc, smallest_node_in_biheap_hc);
  return ;
}

/* This will BiHeapify all nodes in [0, total_num_nodes).
 * Assumes that total_num_nodes is odd.
 */
/*
 * Remark:
 *  (1) This algorithm has complexity O(total_num_nodes). To see why, recall the
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
inline void BiHeapifyOdd(RAI first, size_type total_num_nodes) {
  if(total_num_nodes < 2)
    return ;
  auto num_nodes_in_heap = HeapSize(total_num_nodes);

  size_type smallest_node_in_biheap_hc = (total_num_nodes / 2) + 1;
  size_type largest_node_in_biheap_hc  = smallest_node_in_biheap_hc - 1;

  largest_node_in_biheap_hc++;

  while (smallest_node_in_biheap_hc > 0) {
    --smallest_node_in_biheap_hc;
    SiftFromMinToMaxOdd<RAI, size_type>(first, total_num_nodes, num_nodes_in_heap,
                          smallest_node_in_biheap_hc,
                          largest_node_in_biheap_hc);
    if (largest_node_in_biheap_hc < total_num_nodes - 1) {
      ++largest_node_in_biheap_hc;
      SiftFromMaxToMinOdd<RAI, size_type>(first, total_num_nodes, num_nodes_in_heap,
          FLIP(largest_node_in_biheap_hc),
          smallest_node_in_biheap_hc);
    }
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyOdd(RAI first, RAI one_past_last) {
  BiHeapifyOdd<RAI, size_type>(first, std::distance(first, one_past_last));
}

#undef FLIP

#endif /* BIHEAPIFY_ODD_H_ */
