/*
 * biheapify_simple.h
 *
 *  Created on: Jun 23, 2017
 *      Author: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 * The BiHeapifySimple() function is an implementation of the biheapify
 *  combines the algorithms that constitute BiHeapifyEven() and
 *  BiHeapifyOdd(). The implementation of BiHeapifySimple() is half the size of
 *  BiHeapify() and largely consists of SiftFromMinToMaxSimple(), which consists
 *  of an implementation of sifting down a min heap and sifting up a max heap,
 *  and of SiftFromMaxToMinSimple(), which consists of an implementation of
 *  sifting down a max heap and sifting up a min heap.
 * In short, this means that the amount of work it takes to implement
 *  SiftFromMinToMaxSimple() is only marginally more than the work it takes
 *  to implement the basic sifting up and down operations for both a min heap
 *  and a max heap.
 * However, when total_num_nodes is even then BiHeapifySimple() does do slightly
 *  more computations than is necessary, which is a (small) part the reason why
 *  BiHeapify() is separated into the BiHeapifyEven() and BiHeapifyOdd() functions.
 * The most essential and non-trivial part of the implementation
 *  of this algorithm is knowing how to compute GetNumNodesInHeapContainedInBiheap(),
 *  which is defined in the biheap_common.h header.
 */

#ifndef BIHEAPIFY_SIMPLE_H_
#define BIHEAPIFY_SIMPLE_H_

#include <cassert>

#include <algorithm>

#include "biheap_common.h"

//Note that FlipCo(coord1) >= coord2 iff coord1 <= FlipCo(coord2) (ditto for >, <=, and <).
#define FLIP_COORDINATE(a) (total_num_nodes - 1 - (a))

typedef std::size_t size_type;

template<class RAI>
void SiftFromMinToMaxSimple(RAI first, size_type total_num_nodes,
                      size_type num_nodes_in_heap,
                      size_type first_node_in_mirror_heap,
                      size_type pos_hc,
                      size_type largest_node_in_biheap_hc) {
  bool is_odd_and_middle = (total_num_nodes % 2 == 1) && (pos_hc == total_num_nodes / 2);
  while ((pos_hc < first_node_in_mirror_heap) //While the node is NOT in the max heap
       || is_odd_and_middle) { //Or if total_num_nodes is odd and we've begun at the middle child.
    if (pos_hc > largest_node_in_biheap_hc) //If the node is not in the biheap.
      return ;

    auto left_child     = GetLeftChildInBiheap(pos_hc);
    auto right_child    = left_child + 1;
    bool is_left_child_valid  = (left_child <= largest_node_in_biheap_hc) && //Is the node in the biheap?
                                (left_child  < num_nodes_in_heap);            //Is the node in the min heap?
    bool is_right_child_valid = (right_child <= largest_node_in_biheap_hc) &&
                                (right_child < num_nodes_in_heap);
    if (!is_left_child_valid)
      return ;

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
  SiftUpMaxHeapMC(first, total_num_nodes, FLIP_COORDINATE(pos_hc), FLIP_COORDINATE(largest_node_in_biheap_hc));
  return ;
}

template<class RAI>
void SiftFromMaxToMinSimple(RAI first, size_type total_num_nodes,
                      size_type num_nodes_in_heap,
                      size_type first_node_in_mirror_heap,
                      size_type pos_mc,
                      size_type smallest_node_in_biheap_hc) {
  auto pos_hc     = FLIP_COORDINATE(pos_mc);
  bool is_odd_and_middle = (total_num_nodes % 2 == 1) && (pos_mc == total_num_nodes / 2);
  while ((pos_mc < first_node_in_mirror_heap) //While the node is NOT in the min heap
         || is_odd_and_middle) { //Or if total_num_nodes is odd and we've begun at the middle child.
    if (pos_hc < smallest_node_in_biheap_hc) //If the node is not in the biheap
      return ;

    auto left_child_mc  = GetLeftChildInBiheap(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP_COORDINATE(left_child_mc); //= pos_hc - pos_mc - 1;
    auto right_child_hc = left_child_hc - 1; //= FLIP_COORDINATE(right_child_mc)
    auto pos_it   = first + pos_hc;
    auto left_it  = first + left_child_hc;
    auto right_it = first + right_child_hc;

    bool is_left_child_valid  = (left_child_hc >= smallest_node_in_biheap_hc) && //Is the node in the biheap?
                                (left_child_mc < num_nodes_in_heap);             //Is the node in the max heap?
    bool is_right_child_valid = (right_child_hc >= smallest_node_in_biheap_hc) &&
                                (right_child_mc < num_nodes_in_heap);
    if (!is_left_child_valid)
      return ;

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
  SiftUpMinHeapHC(first, total_num_nodes, pos_hc, smallest_node_in_biheap_hc);
  return ;
}

template<class RAI>
void BiHeapifySimpleSinglePass(RAI first, size_type total_num_nodes,
                      size_type biheap_lower_bound_node_hc = 0,
                      size_type biheap_upper_bound_node_hc = 0,
                      size_type node_to_start_biheapification_at = static_cast<size_type>(-1)) {
  if (biheap_lower_bound_node_hc == 0 && biheap_upper_bound_node_hc == 0)
    biheap_upper_bound_node_hc = total_num_nodes - 1;
  //If it's small enough that it's easiest to just sort everything.
  if(biheap_upper_bound_node_hc - biheap_lower_bound_node_hc < 11) {
    std::sort(first + biheap_lower_bound_node_hc, first + (biheap_upper_bound_node_hc + 1));
    return ;
  }
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  auto first_node_in_mirror_heap  = total_num_nodes - num_nodes_in_heap - (total_num_nodes % 2);

  auto num_nodes_to_biheapify = biheap_upper_bound_node_hc - biheap_lower_bound_node_hc + 1;
  if (node_to_start_biheapification_at < biheap_lower_bound_node_hc
   || node_to_start_biheapification_at > biheap_upper_bound_node_hc
   || node_to_start_biheapification_at == static_cast<size_type>(-1)) {
    node_to_start_biheapification_at = biheap_lower_bound_node_hc
        + (num_nodes_to_biheapify / 2) + (num_nodes_to_biheapify % 2);
  }
  size_type smallest_node_in_biheap_hc = node_to_start_biheapification_at;
  size_type largest_node_in_biheap_hc  = node_to_start_biheapification_at - 1;

  if (total_num_nodes % 2 == 1)
    largest_node_in_biheap_hc++;

  while(smallest_node_in_biheap_hc > 0) {
    if (smallest_node_in_biheap_hc > biheap_lower_bound_node_hc) {
    --smallest_node_in_biheap_hc;
    SiftFromMinToMaxSimple<RAI>(first, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap,
                     smallest_node_in_biheap_hc, largest_node_in_biheap_hc);
    }
    if (largest_node_in_biheap_hc < biheap_upper_bound_node_hc) {
      ++largest_node_in_biheap_hc;
      SiftFromMaxToMinSimple<RAI>(first, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap,
                       FLIP_COORDINATE(largest_node_in_biheap_hc), (smallest_node_in_biheap_hc));
    }

  }
  return ;
}

template<class RAI>
void BiHeapifySimple(RAI first, size_type total_num_nodes) {
  do {
    BiHeapifySimpleSinglePass(first, total_num_nodes);
  } while(!IsBiheap(first, total_num_nodes, 0, total_num_nodes - 1, false));
  return ;
}

#undef FLIP_COORDINATE

#endif /* BIHEAPIFY_SIMPLE_H_ */
