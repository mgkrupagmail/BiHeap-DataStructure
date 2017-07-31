/*
 * biheapify_even.h
 *
 *  Created on: Jun 22, 2017
 *      Author: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 */

#ifndef BIHEAPIFY_EVEN_H_
#define BIHEAPIFY_EVEN_H_

#include <algorithm>

#include "biheap_common.h"

/* Note that FlipCo(coord1) >= coord2 iff coord1 <= FlipCo(coord2)
 * (ditto for >, <=, and <).
 */
#define FLIP_COORDINATE(a) (total_num_nodes - 1 - (a))

namespace {
/* Assumes that total_num_nodes is even, that the node pos_mc belongs to the
 *  min heap, and that pos_hc <= largest_node_in_biheap_hc.
 */
template<class RAI>
inline void SiftFromMinToMaxEven(RAI first, size_type total_num_nodes,
                          size_type num_nodes_in_heap,
                          size_type first_node_in_mirror_heap,
                          size_type pos_hc,
                          size_type largest_node_in_biheap_hc) {
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child     = GetLeftChildInBiheap(pos_hc);
    auto right_child    = left_child + 1;
    auto left_it  = first + left_child;
    auto right_it = first + right_child;
    auto pos_it   = first + pos_hc;

    RAI smaller_it;
    if (*right_it < *left_it) {
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
  SiftUpMaxHeapHC<RAI>(first, total_num_nodes, pos_hc,
                       FLIP_COORDINATE(largest_node_in_biheap_hc));
  return ;
}

/* Assumes that total_num_nodes is even, that the node pos_mc belongs to the
 *  max heap, and that FLIP_COORDINATE(pos_mc) >= smallest_node_in_biheap_hc.
 */
template<class RAI>
inline void SiftFromMaxToMinEven(RAI first, size_type total_num_nodes,
                          size_type num_nodes_in_heap,
                          size_type first_node_in_mirror_heap,
                          size_type pos_mc,
                          size_type smallest_node_in_biheap_hc) {
  auto pos_hc = FLIP_COORDINATE(pos_mc);
  //While node pos_mc is NOT in the min heap.
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = GetLeftChildInBiheap(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= GetRightChildInBiheap(pos_mc);
    auto left_child_hc  = FLIP_COORDINATE(left_child_mc);
    auto right_child_hc = left_child_hc - 1; //= FLIP_COORDINATE(right_child_mc)
    auto pos_it   = first + pos_hc;
    auto left_it  = first + left_child_hc;
    auto right_it = first + right_child_hc;

    RAI larger_it;
    if (*right_it > *left_it) {
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
  SiftUpMinHeapHC<RAI>(first, pos_hc, smallest_node_in_biheap_hc);
  return ;
}

} //End anonymous namespace

/* This will BiHeapify the first total_num_nodes iterated by *first, ...,
 *  *(first + (total_num_nodes - 1)).
 * Assumes that total_num_nodes is even.
 * Experimentation (as opposed to proof) indicates that BiHeapifyEven() always
 *  succeeds in creating a biheap. Please contact the author if a
 *  counter-example is found in which case this algorithm can be made correct
 *  by mimicking the definitions of BiHeapifyOdd() and BiHeapifyOdd()
 *  found in biheapify_odd.h
 *
 * Remark:
 * (1) This algorithm has complexity O(total_num_nodes). To see why, recall the
 *  argument showing that the heapify operation has O(n) complexity (e.g. as
 *  found on pp. 115 - 116 of "The Algorithm Design Manual" 2nd edition); this
 *  argument generalizes to prove that this algorithm also runs in O(n) times.
 * The key observation is that the biheap is constructed by adding one node
 *  at a time with this node alternating between a node in a min heap and a
 *  node in the max heap. The complexity of this biheapification is easily
 *  seen to be twice the complexity of the above mentioned heapification
 *  operation plus a constant.
 * (2) It appears that only one call to BiHeapifyEven() is needed in order to
 *  obtain a biheap.
 */
template<class RAI>
void BiHeapifyEven(RAI first, size_type total_num_nodes) {
  if (total_num_nodes < 2)
    return ;
  auto num_nodes_in_heap = GetNumNodesInHeapContainedInBiheap(total_num_nodes);
  auto first_node_in_mirror_heap  = total_num_nodes - num_nodes_in_heap;

  size_type smallest_node_in_biheap_hc = total_num_nodes / 2;
  size_type largest_node_in_biheap_hc  = smallest_node_in_biheap_hc - 1;

  while (smallest_node_in_biheap_hc > 0) {
    smallest_node_in_biheap_hc--;
    SiftFromMinToMaxEven<RAI>(first, total_num_nodes, num_nodes_in_heap,
                         first_node_in_mirror_heap, smallest_node_in_biheap_hc,
                         largest_node_in_biheap_hc);
    largest_node_in_biheap_hc++;
    //Note that FLIP_COORDINATE(largest_node_in_biheap_hc) ==
    //                                               smallest_node_in_biheap_hc.
    SiftFromMaxToMinEven<RAI>(first, total_num_nodes, num_nodes_in_heap,
        first_node_in_mirror_heap, smallest_node_in_biheap_hc,
        smallest_node_in_biheap_hc);
  }
  return ;
}

/*This will BiHeapify all nodes in [0, total_num_nodes).
 * Assumes that total_num_nodes is even.
 */
template<class RAI>
inline void BiHeapifySafeEven(RAI first, size_type total_num_nodes) {
  if (total_num_nodes < 2)
    return ;

  /* Remarks:
   * (1) The following do {} while(); loop will stop due to the same reasoning
   *  that can be found in Remark (1) in BiHeapifyOdd()'s definition.
   *
   * (2) Experimentation shows that the only one call to
   *  BiHeapifyEven() is ever needed to form a biheap.
   */
  BiHeapifyEven(first, total_num_nodes);

  /* The while loop below should be redundant and is included just in case
   *  BiHeapifyEven() didn't successfully form a biheap. This is
   *  recommended if performance is important since IsBiheap() is an
   *  O(total_num_nodes) operation.
   */
  while (!IsBiheap(first, total_num_nodes)) {
    //The line below can be removed without affecting the correctness of the
    // algorithm.
    GetBiHeapifyFailureMessage(first, total_num_nodes);
    BiHeapifyEven(first, total_num_nodes);
  }
  return ;
}

/* Similar to BiHeapifyOdd(), except that it biheapifies an odd collection of
 *  elements using a pair of calls to BiHeapifyEven() in lieu of the
 *  single call to BiHeapifyOdd() that's done in BiHeapifyOdd().
 * It is not clear why this function empirically succeeds in biheapifying
 *  an odd number of elements.
 */
/* REMARK:
 * (1) If you were to replace the two calls to BiHeapifyEven() with
 *      two calls to BiHeapifyOdd(), then there is input data that the
 *      resulting function will NEVER succeed in biheapifying
 *      no matter how many times it repeatedly calls these two functions.
 *     This starts to happen when total_num_nodes >= 11.
 *     It is not clear why BiHeapifyUsingBiHeapifyEven() empirically always
 *      appears to produce a biheap while this fails to be true if
 *      BiHeapifyEven() is replaced with BiHeapifyOdd().
 */
template<class RAI>
void BiHeapifyUsingBiHeapifyEven(RAI first, size_type total_num_nodes) {
  if (total_num_nodes < 3) {
    if (total_num_nodes == 2 && *first > *(first + 1))
      std::iter_swap(first, first + 1);
    return ;
  }

  auto counter = 0u;
  do {
    BiHeapifyEven(first + 1, total_num_nodes - 1);
    BiHeapifyEven(first, total_num_nodes - 1);
    if (counter != 0 && counter % 3 == 0) {
      std::cerr << "BiHeapifyUsingBiHeapifyEven() could not biheapify after "
                << (counter + 1) << " call(s)." << std::endl;
      assert(false);
    }
    counter++;
  } while(!IsBiheap(first, total_num_nodes));
  return ;
}

#undef FLIP_COORDINATE

#endif /* BIHEAPIFY_EVEN_H_ */
