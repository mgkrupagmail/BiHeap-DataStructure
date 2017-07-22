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
#include <cassert>

namespace {
typedef std::size_t size_type;
}

#include "biheap_ostream.h"

#define FLIP_COORDINATE(a) (total_num_nodes - 1 - (a))

inline size_type GetLeftChildInBiheap(size_type node) { return (2 * node) + 1; }
inline size_type GetRightChildInBiheap(size_type node) {return 2 * (node + 1); }
inline size_type GetParentInBiheapNotRoot(size_type node) {
  return static_cast<size_type> ((node - 1)/2);
}
inline size_type GetParentInBiheap(size_type node) {
  if(node == 0) return -1; else return GetParentInBiheapNotRoot(node);
}
inline size_type GetParentInBiheapZero(size_type node) {
  if(node == 0) return 0; else return GetParentInBiheapNotRoot(node);
}

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
 * If num_nodes_in_heap_odd[i] == #heap nodes in a biheap of size
 * 2 * i - 1, (note the minus one)
 *  then for all i >= 1
 *   num_nodes_in_heap_odd[i] == 4*((i-(i%3))/3) + num_nodes_in_heap_odd[(i%3)]
 * Equivalently, if total_num_nodes is odd then
 *   num_nodes_in_heap_odd[i] == 4 * (total_num_nodes / 6)
 *                               + (n < 4) ? (n + 1) / 2 : 4,
 *   where n = total_num_nodes % 6.
 * Also,
 *  num_nodes_in_heap_odd[i + 3] == num_nodes_in_heap_odd[i] + 4 for all i >= 0,
 *  where num_nodes_in_heap_odd[0] == 0, num_nodes_in_heap_odd[1] == 1,
 *        num_nodes_in_heap_odd[2] == 2, (num_nodes_in_heap_odd[3] == 4, ...)
 */
inline size_type GetNumNodesInHeapContainedInBiheap(size_type total_num_nodes) {
  auto i = total_num_nodes / 2;
  size_type num_nodes_in_heap;
  if (total_num_nodes <= 2)
    num_nodes_in_heap = total_num_nodes;
  else if (total_num_nodes % 2 == 0)
    num_nodes_in_heap = 4*((i-(i%3))/3) + (i%3) + 1 - (((i+2)%3)/2);
  else {
    auto n = total_num_nodes % 6;
    num_nodes_in_heap = 4 * (total_num_nodes / 6) + (n + 1) / 2 + (n == 5);
    //i.e. 4 * (total_num_nodes / 6) and if n == 1 then add 1, if n == 3 then
    // add 2 and if n == 5 then add 4, where these are the only possible
    // values of n. i.e.:
    //num_nodes_in_heap = 4 * (total_num_nodes / 6);
    //num_nodes_in_heap += (n < 4) ? (n + 1) / 2 : 4;
  }
  return num_nodes_in_heap;
}

/* Checks whether or not the first total_um_nodes given given the iterator
 *  first define a biheap.
 */
template<class RAI>
bool IsBiheap(RAI first, size_type total_num_nodes) {
  if (total_num_nodes <= 3) {
    if(total_num_nodes <= 1)
      return true;
    else if (total_num_nodes == 2)
      return *first <= *(first + 1);
    else if (total_num_nodes == 3)
      return *first <= *(first + 1) && *(first + 1) <= *(first + 2);
  }
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  //Check that the nodes first, ..., first + num_nodes_in_heap - 1 form
  // a min heap with min at first. This is half of the biheap condition.
  {
    auto i = 0;
    for (size_type right_child; (right_child = GetRightChildInBiheap(i))
                                                     < num_nodes_in_heap; i++) {
      assert(right_child == GetRightChildInBiheap(i));
      auto parent_val = *(first + i);
      //Check that the parent and left child satisfy the min heap condition.
      if (parent_val > *(first + (right_child - 1)))
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (parent_val > *(first + right_child))
        return false;
    }
    //If the min heap's last element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = GetLeftChildInBiheap(i)) < num_nodes_in_heap
          && *(first + i) > *(first + left_child))
        return false;
    }
  }
  //Check that the nodes FLIP_COORDINATE(first), ...,
  // FLIP_COORDINATE(first + num_nodes_in_heap - 1) form a max heap with max
  // at first + total_num_nodes - 1.
  {
    auto i = 0;
    for (size_type right_child; (right_child = GetRightChildInBiheap(i))
                                                     < num_nodes_in_heap; i++) {
      assert(GetRightChildInBiheap(i) < total_num_nodes);
      auto parent_val = *(first + FLIP_COORDINATE(i));
      auto mirror_left_child_hc = FLIP_COORDINATE(right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_val < *(first + mirror_left_child_hc))
        return false;

      //right_child_hc = FLIP_COORDINATE(right_child_mc)
      // = FLIP_COORDINATE(left_child + 1)
      // = total_num_nodes - 1 - (left_child + 1) = mirror_left_child_hc - 1.
      //assert(FLIP_COORDINATE(left_child + 1) == mirror_left_child_hc - 1);
      //Check that the parent and right child satisfy the max heap condition.
      if (parent_val < *(first + (mirror_left_child_hc - 1)))
        return false;
    }
    //If the max heap's last element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = GetLeftChildInBiheap(i)) < num_nodes_in_heap
          && *(first + FLIP_COORDINATE(i))
           < *(first + FLIP_COORDINATE(left_child)))
        return false;
    }
  }
  return true;
}

//Helper function for IsBiheap().
template<class RAI>
std::string GetIsBiHeapNotMinHeapErrorMessageHC(RAI first,
                                                size_type total_num_nodes,
                                                size_type parent_hc,
                                                bool is_left_child) {
  std::stringstream strm;
  std::string child_name = is_left_child ? std::string("left_child")
                                         : std::string("right_child");
  size_type child_index = GetLeftChildInBiheap(parent_hc) + !is_left_child;
  strm << "Failed to be a MIN  heap at parent_hc = " << parent_hc
       << " due to " << child_name << " = " << child_index
       << " \t*parent = " << *(first + parent_hc)
       << " \t*" << child_name << " = " << *(first + child_index)
       << " \ttotal_num_nodes = " << total_num_nodes << std::endl;
  return strm.str();
}

//Helper function for IsBiheap().
template<class RAI>
std::string GetIsBiHeapNotMaxHeapErrorMessageMC(RAI first,
                                                size_type total_num_nodes,
                                                size_type parent_mc,
                                                bool is_left_child) {
  std::stringstream strm;
  auto parent_hc = FLIP_COORDINATE(parent_mc);
  std::string child_name = is_left_child ? std::string("left_child")
                                         : std::string("right_child");
  size_type child_index_mc = GetLeftChildInBiheap(parent_mc) + !is_left_child;
  size_type child_index_hc = FLIP_COORDINATE(child_index_mc);
  strm << "Failed to be a MAX heap at parent_hc = " << parent_hc
       << " (parent_mc = " << parent_mc
       << ") due to " << child_name << "_hc = " << child_index_hc
       << " ("        << child_name << "_mc = " << child_index_mc
       << ") \t*parent = " << *(first + parent_hc)
       << " \t*" << child_name << " = " << *(first + child_index_hc)
       << " \ttotal_num_nodes = " << total_num_nodes << std::endl;
  return strm.str();
}

std::string GetIsBiheapStartAndEndWarningMessage(
    size_type biheap_start_node_hc = 0, size_type biheap_end_node_hc = 0) {
  std::stringstream strm;
  strm << "WARNING In IsBiheap():" << biheap_start_node_hc
       << " = biheap_start_node_hc > biheap_end_node_hc = "
       << biheap_end_node_hc << ". Swapping values and continuing."
       << std::endl;
  return strm.str();
}

#define ISBIHEAP_OSTREAM_DEFAULT std::cout
/* Given the iterator first, this checks whether or not the biheap condition is
 *  satisfied by all elements in
 *  first + [biheap_start_node_hc, biheap_end_node_hc] where this check is
 *  applied ONLY if BOTH the parent and child elements belong to this interval.
 * should_output_failure_messages - Setting this to true will cause an
 *   explanatory message to be outputted describing the elements that violate
 *   the biheap condition and how they violate it.
 * ostrm - Only used if should_output_failure_messages is true.
 * return_on_failure - Set this to false if you want to find ALL elements at
 *   which the biheap condition fails. Setting this to false is useful for
 *   debugging.
 * Returns true if it is a biheap and false otherwise.
 */
template<class RAI>
bool IsBiheap(RAI first, size_type total_num_nodes,
           size_type biheap_start_node_hc, size_type biheap_end_node_hc = 0,
           bool should_output_failure_messages = false,
           std::ostream &ostrm = ISBIHEAP_OSTREAM_DEFAULT,
           bool return_on_failure = true) {
  if (biheap_start_node_hc > biheap_end_node_hc) {
    ostrm << GetIsBiheapStartAndEndWarningMessage(biheap_start_node_hc,
                                                  biheap_end_node_hc);
    std::swap(biheap_start_node_hc, biheap_end_node_hc);
  }
  if (biheap_start_node_hc == 0 && biheap_end_node_hc == 0)
    biheap_end_node_hc = total_num_nodes - 1;

  {
    auto num_nodes_in_biheap = biheap_end_node_hc - biheap_start_node_hc + 1;
    if(total_num_nodes <= 1 || num_nodes_in_biheap <= 1)
      return true;
    auto first_in_biheap = first + biheap_start_node_hc;
    if (total_num_nodes == 2) //which implies that num_nodes_in_biheap == 2.
      return *(first_in_biheap) <= *(first_in_biheap + 1);
    //At this point, total_num_nodes > 2 and num_nodes_in_biheap > 1.
  }
  bool is_heap = true;
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  size_type left_child;
  for (auto i = biheap_start_node_hc; GetLeftChildInBiheap(i) <num_nodes_in_heap
      && GetLeftChildInBiheap(i) <= biheap_end_node_hc; i++) {
    left_child = GetLeftChildInBiheap(i);
    auto parent_val = *(first + i);
    auto mirror_parent_val = *(first + FLIP_COORDINATE(i));
    auto mirror_left_child_hc = FLIP_COORDINATE(left_child);
    //Check that the nodes first, ..., first + num_nodes_in_heap - 1 form
    // a min heap.
    if (parent_val > *(first + left_child)) {
      if (should_output_failure_messages)
        ostrm << GetIsBiHeapNotMinHeapErrorMessageHC(first, total_num_nodes, i,
                                                     true)
              << std::endl;
        /*ostrm << "Failed to be a MIN heap at parent_hc = " << i
              << " due to left_child = " << left_child
              << " \t*parent = " << parent_val
              << " \t*left_child = " << *(first + left_child)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;*/
      if (return_on_failure)
        return false;
      is_heap = false;
    }
    //Check that the nodes FLIP_COORDINATE(first), ...,
    // FLIP_COORDINATE(first + num_nodes_in_heap - 1) form a max heap.
    if (mirror_parent_val < *(first + mirror_left_child_hc)) {
      if (should_output_failure_messages)
        ostrm << GetIsBiHeapNotMaxHeapErrorMessageMC(first, total_num_nodes, i,
                                                     true)
              << std::endl;
        /*ostrm << "Failed to be a MAX heap at parent_hc = "
              <<  FLIP_COORDINATE(i)
              << " (parent_mc = " << i
              << ") due to max heap left_child_hc = " << mirror_left_child_hc
              << " (left_child_mc = " << left_child
              << ") \t*parent = " << mirror_parent_val
              << " \t*left_child = " << *(first + mirror_left_child_hc)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;*/
      if (return_on_failure)
        return false;
      is_heap = false;
    }

    //If the right child is not in the heap or the biheap.
    if (left_child + 1 >= num_nodes_in_heap
        || left_child + 1 >= biheap_end_node_hc)
      break;

    if (parent_val > *(first + (left_child + 1))) {
      if (should_output_failure_messages)
        ostrm << GetIsBiHeapNotMinHeapErrorMessageHC(first, total_num_nodes, i,
                                                     false)
              << std::endl;
        /*ostrm << "Failed to be a MIN heap at parent_hc = " << i
              << " due to min heap right_child = " << (left_child + 1)
              << " \t*parent = " << parent_val
              << " \t*right_child = " << *(first + (left_child + 1))
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;*/
      if (return_on_failure)
        return false;
      is_heap = false;
    }
    //right_child_hc = FLIP_COORDINATE(right_child_mc)
    // = FLIP_COORDINATE(left_child + 1)
    // = total_num_nodes - 1 - (left_child + 1) = mirror_left_child_hc - 1.
    //assert(FLIP_COORDINATE(left_child + 1) == mirror_left_child_hc - 1);
    if (mirror_parent_val < *(first + (mirror_left_child_hc - 1))) {
      if (should_output_failure_messages)
        ostrm << GetIsBiHeapNotMaxHeapErrorMessageMC(first, total_num_nodes, i,
                                                     false)
              << std::endl;
        /*ostrm << "Failed to be a MAX heap at parent_hc = "
              << (FLIP_COORDINATE(left_child + 1)) << " (parent_mc = " << i
              << ") due to max heap right_child_hc = "
              << (mirror_left_child_hc - 1)
              << " (right_child_mc = " << (left_child + 1)
              << ") \t*parent = " << mirror_parent_val
              << " \t*right_child = " << *(first + mirror_left_child_hc - 1)
              << " \ttotal_num_nodes = " << total_num_nodes << std::endl;*/
      if (return_on_failure)
        return false;
      is_heap = false;
    }
  }
  return is_heap;
}

template<class RAI>
inline void SiftUpMaxHeapMC(RAI first, size_type total_num_nodes,
                      size_type pos_mc, size_type smallest_node_in_biheap_mc) {
  size_type parent;
  if (pos_mc == 0 ||
      (parent = GetParentInBiheapNotRoot(pos_mc)) < smallest_node_in_biheap_mc)
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
  } while (pos_mc > 0 &&
     (parent = GetParentInBiheapNotRoot(pos_mc)) >= smallest_node_in_biheap_mc);
  return ;
}

template<class RAI>
inline void SiftUpMaxHeapHC(RAI first, size_type total_num_nodes,
                       size_type pos_hc, size_type smallest_node_in_biheap_mc) {
  SiftUpMaxHeapMC(first, total_num_nodes, FLIP_COORDINATE(pos_hc),
                  smallest_node_in_biheap_mc);
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapHC(RAI first, size_type total_num_nodes,
                       size_type pos_hc, size_type smallest_node_in_biheap_hc) {
  size_type parent;
  if (pos_hc == 0 ||
      (parent = GetParentInBiheapNotRoot(pos_hc)) < smallest_node_in_biheap_hc)
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
  } while (pos_hc > 0 &&
     (parent = GetParentInBiheapNotRoot(pos_hc)) >= smallest_node_in_biheap_hc);
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI>
inline void SiftUpMinHeapMC(RAI first, size_type total_num_nodes,
                       size_type pos_mc, size_type smallest_node_in_biheap_hc) {
  SiftUpMinHeapHC(first, total_num_nodes, FLIP_COORDINATE(pos_mc),
                  smallest_node_in_biheap_hc);
}

#undef ISBIHEAP_OSTREAM_DEFAULT
#endif /* BIHEAP_COMMON_H_ */
