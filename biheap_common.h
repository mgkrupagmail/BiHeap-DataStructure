/*
 * biheap_common.h
 *
 *  Created on: Jun 23, 2017
 *      Author: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 * This file contains function definitions that are used by other biheapify
 *  files, such as biheapify_even.h, biheapify_odd.h, and biheapify_simple.h.
 */

#ifndef BIHEAP_COMMON_H_
#define BIHEAP_COMMON_H_


#include <iostream>

namespace {
typedef std::size_t size_type;
}

#include "biheap_ostream.h"


#define FLIP_COORDINATE(a) (total_num_nodes - 1 - (a))

static inline size_type GetLeftChildInBiheap(size_type node) { return (2 * node) + 1; }
static inline size_type GetRightChildInBiheap(size_type node) {return 2 * (node + 1); }
static inline size_type GetParentInBiheapNotRoot(size_type node) { return static_cast<size_type> ((node - 1)/2); }
static inline size_type GetParentInBiheap(size_type node) { if(node == 0) return -1; else return GetParentInBiheapNotRoot(node); }
static inline size_type GetParentInBiheapZero(size_type node) { if(node == 0) return 0; else return GetParentInBiheapNotRoot(node); }

/*
 * This function is the most complicated part of implementing a biheapify
 *  algorithm and the value that it produces is a quantity fundamental to any
 *  biheap.
 * The following formulas were obtained by finding a repeating pattern in the
 *  number of nodes that max up a biheap's (min or max) heap, which is dependent
 *  on total_num_nodes.
 * Here we must separate num_nodes_in_heap (dependent on total_num_nodes) into
 *  two subsequences, num_nodes_in_heap_even and num_nodes_in_heap_odd, based on
 *  whether total_num_nodes is even or odd.
 * If num_nodes_in_heap_even[i] == #extended heap nodes in a biheap of size 2 * i,
 *  then for all i >= 0
 *  num_nodes_in_heap_even[i] == 4*((i-(i%3))/3) + num_nodes_in_heap_even[(i%3)]
 *                            == 4*((i-(i%3))/3) + (i%3) + 1 - (((i+2)%3)/2)
 *  holds and
 *  (*) num_nodes_in_heap[i + 3] == num_nodes_in_heap[i] + 4 for all i >= 0,
 *  where num_nodes_in_heap[0] == 0, num_nodes_in_heap[1] == 2,
 *        num_nodes_in_heap[2] == 3, (num_nodes_in_heap[3] == 4, ...)
 *
 * If num_nodes_in_heap_odd[i] == #extended heap nodes in a biheap of size
 * 2 * i - 1, (note the minus one)
 *  then for all i >= 1
 *   num_nodes_in_heap_odd[i] == 4*((i-(i%3))/3) + num_nodes_in_heap_odd[(i%3)]
 *                            == 4*((i-(i%3))/3) + (i%3)
 * Equivalently, if num_nodes_in_heap_odd_alt[j] == #extended heap nodes in
 *  a biheap of size 2 * j + 1, then for all j >= 0
 *  num_nodes_in_heap_odd_alt_alt[i] == 4*(((j-1)-((j-1)%3))/3) + num_nodes_in_heap_odd_alt[((j-1)%3)]
 *                                   == 4*(((j-1)-((j-1)%3))/3) + ((j-1)%3)
 * Also,
 *  num_nodes_in_heap_odd[i + 3] == num_nodes_in_heap_odd[i] + 4 for all i >= 0,
 *  where num_nodes_in_heap_odd[0] == 0, num_nodes_in_heap_odd[1] == 1,
 *        num_nodes_in_heap_odd[2] == 2, (num_nodes_in_heap_odd[3] == 4, ...)
 */
static inline size_type GetNumNodesInHeapContainedInBiheap(size_type total_num_nodes) {
  auto i = total_num_nodes / 2;
  size_type num_nodes_in_heap;
  if (total_num_nodes <= 2)
    num_nodes_in_heap = total_num_nodes;
  else if (total_num_nodes % 2 == 0)
    num_nodes_in_heap = 4*((i-(i%3))/3) + (i%3) + 1 - (((i+2)%3)/2);
  else {
    i--;
    num_nodes_in_heap = 4*((i-(i%3))/3) + (i%3);
  }
  return num_nodes_in_heap;
}

#define ISBIHEAP_OSTREAM_DEFAULT std::cout
//Calling Biheap(first, total_num_nodes, 0, 0) will change biheap_end_node_hc
// to equal total_num_nodes - 1.
template<class RAI>
static bool IsBiheap(RAI first, size_type total_num_nodes, size_type biheap_start_node_hc = 0,
                     size_type biheap_end_node_hc = 0, bool should_output_failure_message_to_cout = true,
                     std::ostream &ostrm = ISBIHEAP_OSTREAM_DEFAULT) {
  if (biheap_start_node_hc > biheap_end_node_hc) {
    ostrm << "WARNING In IsBiheap():" << biheap_start_node_hc
          << " = biheap_start_node_hc > biheap_end_node_hc = "
          << biheap_end_node_hc << ". Swapping values and continuing." << std::endl;
    std::swap(biheap_start_node_hc, biheap_end_node_hc);
  }
  if (biheap_start_node_hc == 0 && biheap_end_node_hc == 0)
    biheap_end_node_hc = total_num_nodes - 1;
  if (biheap_start_node_hc == biheap_end_node_hc)
    return true;

  auto num_nodes_in_biheap = biheap_end_node_hc - biheap_start_node_hc + 1;
  auto first_in_biheap = first + biheap_start_node_hc;
  if(num_nodes_in_biheap <= 1)
    return true;
  else if (num_nodes_in_biheap == 2)
    return *(first_in_biheap) <= *(first_in_biheap + 1);
  else if (num_nodes_in_biheap == 3)
    return *(first_in_biheap) <= *(first_in_biheap + 1)
        && *(first_in_biheap + 1) <= *(first_in_biheap + 2);

  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  size_type left_child;
  for (auto i = biheap_start_node_hc; GetLeftChildInBiheap(i) < num_nodes_in_heap
      && GetLeftChildInBiheap(i) <= biheap_end_node_hc; i++) {
    left_child = GetLeftChildInBiheap(i);
    //Check that the nodes first, ..., first + num_nodes_in_heap - 1 form a min heap.
    if (*(first + i) > *(first + left_child)) {
      if (should_output_failure_message_to_cout)
        ostrm << "Failed to be a MIN heap at pos_hc = " << i << " due to left_child= "
              << (left_child) << " \t*(first + pos_hc) = " << *(first + i)
              << " \t*(first + left_child) = " << *(first + left_child)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;
      return false;
    }
    //Check that the nodes FLIP_COORDINATE(first), ..., FLIP_COORDINATE(first + num_nodes_in_heap - 1) form a max heap.
    if (*(first + FLIP_COORDINATE(i)) < *(first + FLIP_COORDINATE(left_child))){
      if (should_output_failure_message_to_cout)
        ostrm << "Failed to be a MAX heap at pos_hc = " << (FLIP_COORDINATE(left_child))
              << " (pos_mc = " << i << ") due to left_child= " << (left_child)
              << " \t*(first + pos_hc) = " << *(first + FLIP_COORDINATE(i))
              << " \t*(first + left_child) = " << *(first + left_child)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;
      return false;
    }

    //If the right child is not in the heap or the biheap.
    if (left_child + 1 >= num_nodes_in_heap || left_child + 1 >= biheap_end_node_hc)
      break;

    if (*(first + i) > *(first + left_child + 1)) {
      if (should_output_failure_message_to_cout)
        ostrm << "Failed to be a MIN heap at pos_hc = " << i << " due to max heap right_child = "
              << (left_child + 1) << " \t*(first + pos_hc) = " << *(first + i)
              << " \t*(first + right_child) = " << *(first + left_child + 1)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;
      return false;
    }
    if (*(first + FLIP_COORDINATE(i)) < *(first + FLIP_COORDINATE(left_child + 1))){
      if (should_output_failure_message_to_cout)
        ostrm << "Failed to be a MAX heap at pos_hc = " << (FLIP_COORDINATE(left_child + 1))
              << " (pos_mc = " << i << ") due to max heap right_child = "
              << (left_child + 1) << " \t*(first + pos_hc) = "
              << *(first + FLIP_COORDINATE(i)) << " \t*(first + right_child) = "
              << *(first + left_child + 1) << " \ttotal_num_nodes = "
              << total_num_nodes << std::endl;
      return false;
    }
  }
  return true;
}

template<class RAI>
inline void SiftUpMaxHeapMC(RAI first, size_type total_num_nodes, size_type pos_mc, size_type smallest_node_in_biheap_mc) {
  size_type parent;
  if (pos_mc == 0 || (parent = GetParentInBiheapNotRoot(pos_mc)) < smallest_node_in_biheap_mc)
    return ;
  auto pos_it = first + FLIP_COORDINATE(pos_mc);
  do {
    auto parent_it = first + FLIP_COORDINATE(parent);
    if (*pos_it > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 && (parent = GetParentInBiheapNotRoot(pos_mc)) >= smallest_node_in_biheap_mc);
  return ;
}

template<class RAI>
inline void SiftUpMaxHeapHC(RAI first, size_type total_num_nodes, size_type pos_hc, size_type smallest_node_in_biheap_mc) {
  SiftUpMaxHeapMC(first, total_num_nodes, FLIP_COORDINATE(pos_hc), smallest_node_in_biheap_mc);
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapHC(RAI first, size_type total_num_nodes, size_type pos_hc, size_type smallest_node_in_biheap_hc) {
  size_type parent;
  if (pos_hc == 0 || (parent = GetParentInBiheapNotRoot(pos_hc)) < smallest_node_in_biheap_hc)
    return ;
  auto pos_it = first + pos_hc;
  do {
    auto parent_it = first + parent;
    if (*pos_it < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 && (parent = GetParentInBiheapNotRoot(pos_hc)) >= smallest_node_in_biheap_hc);
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapMC(RAI first, size_type total_num_nodes, size_type pos_mc, size_type smallest_node_in_biheap_hc) {
  SiftUpMinHeapHC(first, total_num_nodes, FLIP_COORDINATE(pos_mc), smallest_node_in_biheap_hc);
}

#undef ISBIHEAP_OSTREAM_DEFAULT
#endif /* BIHEAP_COMMON_H_ */
