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

#include <algorithm>
#include <cassert>
#include <iostream>

#define FLIP(a) (total_num_nodes - 1 - (a))

template<typename size_type = std::size_t>
inline size_type LeftChild(size_type node) { return (2 * node) + 1; }

template<typename size_type = std::size_t>
inline size_type RightChild(size_type node) {return 2 * (node + 1); }

//Note that even if size_type is unsigned, Parent(0) == 0.
template<typename size_type = std::size_t>
inline size_type Parent(size_type node) {
  return static_cast<size_type>((static_cast<long long>(node - 1))/2);
}

//Unlike Parent(), this assumes that node > 0.
//If size_type is unsigned, then ParentNotRoot(0) != 0.
template<typename size_type = std::size_t>
inline size_type ParentNotRoot(size_type node) {
  return static_cast<size_type>((node - 1)/2);
}

template<typename size_type = std::size_t>
inline size_type HeapSize(size_type total_num_nodes) {
  //assert((total_num_nodes - static_cast<size_type>(total_num_nodes / 3)) == (total_num_nodes - (total_num_nodes - (total_num_nodes % 3)) / 3));
  return total_num_nodes - static_cast<size_type>(total_num_nodes / 3);
}

/* Checks whether or not the first total_um_nodes given given the iterator
 *  first define a biheap.
 */
template<class RAI, typename size_type = std::size_t>
bool IsBiheap(RAI first, size_type total_num_nodes) {
  if (total_num_nodes <= 3) {
    if(total_num_nodes <= 1)
      return true;
    else if (total_num_nodes == 2)
      return *first <= *(first + 1);
    else if (total_num_nodes == 3)
      return *first <= *(first + 1) && *(first + 1) <= *(first + 2);
  }
  size_type num_nodes_in_heap = HeapSize<size_type>(total_num_nodes);
  //Check that the nodes first, ..., first + num_nodes_in_heap - 1 form
  // a min heap with min at first. This is half of the biheap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < num_nodes_in_heap; i++) {
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
      if ((left_child = LeftChild<size_type>(i)) < num_nodes_in_heap
          && *(first + i) > *(first + left_child))
        return false;
    }
  }
  //Check that the nodes FLIP(first), ...,
  // FLIP(first + num_nodes_in_heap - 1) form a max heap with max
  // at first + total_num_nodes - 1.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < num_nodes_in_heap; i++) {
      auto parent_val = *(first + FLIP(i));
      size_type mirror_left_child_hc = FLIP(right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_val < *(first + mirror_left_child_hc))
        return false;

      //right_child_hc = FLIP(right_child_mc)
      // = FLIPleft_child + 1)
      // = total_num_nodes - 1 - (left_child + 1) = mirror_left_child_hc - 1.
      //assert(FLIP(left_child + 1) == mirror_left_child_hc - 1);
      //Check that the parent and right child satisfy the max heap condition.
      if (parent_val < *(first + (mirror_left_child_hc - 1)))
        return false;
    }
    //If the max heap's last element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < num_nodes_in_heap
          && *(first + FLIP(i))
           < *(first + FLIP(left_child)))
        return false;
    }
  }
  return true;
}



template<class RAI, typename size_type = std::size_t>
inline void SiftUpMaxHeapMC(RAI first, size_type total_num_nodes,
                      size_type pos_mc, size_type smallest_node_in_biheap_mc) {
  size_type parent;
  if (pos_mc == 0 ||
      (parent = Parent<size_type>(pos_mc)) < smallest_node_in_biheap_mc)
    return ;
  auto pos_it = first + FLIP(pos_mc);
  do {
    auto parent_it = first + FLIP(parent);
    if (*pos_it > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
     (parent = Parent<size_type>(pos_mc)) >= smallest_node_in_biheap_mc);
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void SiftUpMaxHeapHC(RAI first, size_type total_num_nodes,
                       size_type pos_hc, size_type smallest_node_in_biheap_mc) {
  SiftUpMaxHeapMC<RAI, size_type>(first, total_num_nodes, FLIP(pos_hc),
                       smallest_node_in_biheap_mc);
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename size_type = std::size_t>
inline void SiftUpMinHeapHC(RAI first, size_type pos_hc,
                            size_type smallest_node_in_biheap_hc) {
  size_type parent;
  if (pos_hc == 0 ||
      (parent = Parent<size_type>(pos_hc)) < smallest_node_in_biheap_hc)
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
     (parent = Parent<size_type>(pos_hc)) >= smallest_node_in_biheap_hc);
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename size_type = std::size_t>
inline void SiftUpMinHeapMC(RAI first, size_type total_num_nodes,
                       size_type pos_mc, size_type smallest_node_in_biheap_hc) {
  SiftUpMinHeapHC<RAI, size_type>(first, total_num_nodes, FLIP(pos_mc),
                       smallest_node_in_biheap_hc);
}

#undef ISBIHEAP_OSTREAM_DEFAULT
#undef FLIP
#endif /* BIHEAP_COMMON_H_ */
