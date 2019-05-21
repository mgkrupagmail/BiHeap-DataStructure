/*
 * biqueue.h
 *
 *  Created on: Dec 2, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 *
 *  A BiQueue is a Double Ended Priority Queue based around
 *   the idea of the BiHeap and Fused BiHeap data structures.
 *
 *  It uses Fused BiHeaps to implement a double ended priority
 *   queue. It has amortized O(log N) insertions and
 *   amortized O(log N) deletions.
 *  Given a collection of elements, this double ended priority
 *   queue can be formed using O(N) swaps.
 *
 * Internally, BiQueue stores all objects in a std::vector<ValueType>
 *  object and when it needs more capacity will call std::vector's
 *  resize() method.
 *
 *  Example of using a BiQueue object:

#include "biqueue_standalone.h"

void BiQueueExample() {
  //Define the object.
  biqueue::BiQueue<int> biq;
  std::vector<int> vector = { 3, 6, 5, 4, 3, 2, 1, 2 };
  biqueue::BiQueue<int> biq2(vector.begin(), vector.end());

  //Insert an element into the object
  std::cout << "biq.insert(0) \t\t\t";
  biq.insert(0);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.insert(2) \t\t\t";
  biq.insert(2);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.insert(1) \t\t\t";
  biq.insert(1);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.insert(3) \t\t\t";
  biq.insert(3);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.insert(-1) \t\t\t";
  biq.insert(-1);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.popmax() \t\t\t";
  biq.popmax();
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  std::cout << "biq.popmin() \t\t\t";
  biq.popmin();
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;

  bool should_pop_max = true;
  std::cout << "biq.PopMinOrMax(should_pop_max)\t";
  //if should_pop_max == true then pop the max, otherwise pop the min.
  biq.PopMinOrMax(should_pop_max);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << "  \t2nd smallest = " << biq.second_smallest()
              << " \t2nd largest = " << biq.second_largest() << "  \tmax: " << biq.max() << " \tbiq.size() = " << biq.size()
              << " \tInternal biheap size = N = " << biq.biheap_size() << std::endl;
  return ;
}

 */

#ifndef BIQUEUE_STANDALONE_H_
#define BIQUEUE_STANDALONE_H_

#include <algorithm>
#include <cassert>
#include <vector>

namespace biqueue {

template<typename SizeType = std::size_t>
inline SizeType LeftChild(SizeType node) { return (2 * node) + 1; }

template<typename SizeType = std::size_t>
inline SizeType RightChild(SizeType node) {return 2 * (node + 1); }

//Note that even if SizeType is unsigned, Parent(0) == 0.
template<typename SizeType = std::size_t>
inline SizeType Parent(SizeType node) {
  return static_cast<SizeType>((static_cast<long long>(node) - 1) / 2);
}

//Unlike Parent(), this assumes that node > 0 (hence "NotRoot").
//Note that if SizeType is unsigned, then ParentNotRoot(0) != 0,
// which is the reason for assuming that node > 0.
//If node > 0 then the static_cast<long long> found in Parent()
// can be avoided.
template<typename SizeType = std::size_t>
inline SizeType ParentNotRoot(SizeType node) {
  return (node - 1) / 2;
}

template<typename SizeType = std::size_t>
inline SizeType HeapSize(SizeType N) {
  return N - static_cast<SizeType>(N / 3);
}

template<typename SizeType = std::size_t>
inline SizeType IndexOfLastMinHeapNodeGivenHeapSize(SizeType heap_size) {
  return heap_size - 1;
}

//Assumes that pos_mc is a node in the max heap.
//The node with max heap coordinate last_node_in_biheap_mc is a
// node in the BiHeap constructed so far such that if v is any
// node whose max heap coordinate is < last_node_in_biheap_mc,
// then v does NOT belong to the BiHeap constructed so far.
template<class RAI, typename SizeType = std::size_t>
inline void SiftUpMaxHeapMC(RAI V, SizeType N, SizeType N_minus1, SizeType pos_mc,
                            SizeType last_node_in_biheap_mc) {
  SizeType parent_mc;
  if (pos_mc == 0 ||
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) < last_node_in_biheap_mc)
    return ;
  auto pos_it    = V + (N_minus1 - pos_mc);
  auto pos_value = *pos_it;
  do {
    auto parent_it = V + (N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
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

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename SizeType = std::size_t>
inline void SiftUpMinHeapHC(RAI V, SizeType pos_hc,
                            SizeType first_node_in_biheap_hc) {
  SizeType parent_hc;
  if (pos_hc == 0 ||
      (parent_hc = ParentNotRoot<SizeType>(pos_hc)) < first_node_in_biheap_hc)
    return ;
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto parent_it = V + parent_hc;
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

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifySiftFromMinToMax(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc) {
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;

    //assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    decltype(pos_value) smaller_value = *left_it, right_value;
    RAI smaller_it;
    if (is_right_child_valid && (right_value = *right_it) < smaller_value) {
      smaller_value = right_value;
      smaller_it    = right_it;
      pos_hc        = right_child_hc;
    } else {
      smaller_it    = left_it;
      pos_hc        = left_child_hc;
    }
    if (smaller_value < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  } while (pos_hc < first_in_node) ;
  //Start sifting up the max heap.
  //At this point pos is an In node.
  auto pos_mc                 = N_minus1 - pos_hc;
  auto last_node_in_biheap_mc = N_minus1 - last_node_in_biheap_hc;
  SizeType parent_mc          = ParentNotRoot<SizeType>(pos_mc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_mc == 0 comparison since N > 8 implies this pos_mc > 0.
  if (pos_mc == 0 || parent_mc < last_node_in_biheap_mc)
    return ;
  do {
    auto parent_it = V + (N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
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

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifySiftFromMaxToMin(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc) {
  auto pos_hc    = N_minus1 - pos_mc;
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = N_minus1 - left_child_mc;//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;

    //assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
    bool is_right_child_valid = right_child_mc < heap_size;
    decltype(pos_value) larger_value = *left_it, right_value;
    RAI larger_it;
    if (is_right_child_valid && larger_value < (right_value = *right_it)) {
      larger_value = right_value;
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else {
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < larger_value){
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    }
    else
      return ;
  } while (pos_mc < first_in_node) ;
  //Start sifting up the min heap.
  //At this point pos is an In node.
  SizeType parent_hc = ParentNotRoot<SizeType>(pos_hc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_hc == 0 comparison since N > 8 implies this pos_hc > 0.
  if (pos_hc == 0 || parent_hc < first_node_in_biheap_hc)
    return ;
  do {
    auto parent_it = V + parent_hc;
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

/* This will BiHeapify all nodes in [0, N).
 * Assumes that N is odd.
 */
/*
 * Remark:
 *  This algorithm has complexity O(N). To see why, recall the
 *   argument showing that the heapify operation has O(n) complexity (e.g. as
 *   found on pp. 115 - 116 of "The Algorithm Design Manual" 2nd edition); this
 *   argument generalizes to prove that this algorithm also runs in O(n) times.
 *  The key observation is that the biheap is constructed by adding one node
 *   at a time with this node alternating between a node in a min heap and a
 *   node in the max heap. The complexity of this biheapification is easily
 *   seen to be twice the complexity of the above mentioned heapification
 *   operation plus a constant.
 */
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapify(RAI V, SizeType N) {
  if(N < 2)
    return ;
  SizeType heap_size     = HeapSize(N);
  SizeType first_in_node = N - heap_size;
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  SizeType last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                     - (N % 3 == 2);
  SizeType N_minus1                = N - 1;
  SizeType first_node_in_biheap_hc = N_minus1 - last_node_in_biheap_hc;
  //NOTE: It's not necessary to pass N_minus1, heap_size, or first_in_node to
  // these functions (since they can be computed from N) but we do so since
  // these variables won't actually be placed on the stack each time due to these
  // functions being both inlined and templates. Any half decent optimizer would
  // optimize away allocating stack space and copy these values.
  while (first_node_in_biheap_hc > 0) {
    BiHeapifySiftFromMinToMax<RAI, SizeType>(V, N, N_minus1, heap_size, first_in_node,
       --first_node_in_biheap_hc, last_node_in_biheap_hc);
    BiHeapifySiftFromMaxToMin<RAI, SizeType>(V, N, N_minus1, heap_size, first_in_node,
        N_minus1 - (++last_node_in_biheap_hc), first_node_in_biheap_hc);
  }
  return ;
}


/*
 * ================== START: Definition of lambda version of BiHeapify ====================
 */

//Assumes that pos_mc is a node in the max heap.
//The node with max heap coordinate last_node_in_biheap_mc is a
// node in the BiHeap constructed so far such that if v is any
// node whose max heap coordinate is < last_node_in_biheap_mc,
// then v does NOT belong to the BiHeap constructed so far.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void SiftUpMaxHeapMC(RAI V, SizeType N, SizeType N_minus1, SizeType pos_mc,
                            SizeType last_node_in_biheap_mc, LambdaType lambda) {
  SizeType parent_mc;
  if (pos_mc == 0 ||
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) < last_node_in_biheap_mc)
    return ;
  auto pos_it    = V + lambda(N, N_minus1 - pos_mc);
  auto pos_value = *pos_it;
  do {
    auto parent_it = V + lambda(N, N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
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

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void SiftUpMinHeapHC(RAI V, SizeType N, SizeType N_minus1, SizeType pos_hc,
                            SizeType first_node_in_biheap_hc, LambdaType lambda) {
  SizeType parent_hc;
  if (pos_hc == 0 ||
      (parent_hc = ParentNotRoot<SizeType>(pos_hc)) < first_node_in_biheap_hc)
    return ;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
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

//Assumes that
// (1) the node pos belongs to the min heap,
// (2) pos_hc <= last_node_in_biheap_hc,
// (3) pos is NOT an In node.
// (4) N > 2
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifySiftFromMinToMax(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  do {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);

    //assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    decltype(pos_value) smaller_value = *left_it, right_value;
    RAI smaller_it;
    if (is_right_child_valid && (right_value = *right_it) < smaller_value) {
      smaller_value = right_value;
      smaller_it    = right_it;
      pos_hc        = right_child_hc;
    } else {
      smaller_it    = left_it;
      pos_hc        = left_child_hc;
    }
    if (smaller_value < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  } while (pos_hc < first_in_node) ;
  //Start sifting up the max heap.
  //At this point pos is an In node.
  auto pos_mc                 = N_minus1 - pos_hc;
  auto last_node_in_biheap_mc = N_minus1 - last_node_in_biheap_hc;
  SizeType parent_mc          = ParentNotRoot<SizeType>(pos_mc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_mc == 0 comparison since N > 8 implies this pos_mc > 0.
  if (pos_mc == 0 || parent_mc < last_node_in_biheap_mc)
    return ;
  do {
    auto parent_it = V + lambda(N, N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
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

//Assumes that
// (1) the node pos belongs to the max heap,
// (2) FLIP(pos_mc) >= first_node_in_biheap_hc, and
// (3) pos is NOT an In node.
// (4) N > 2
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifySiftFromMaxToMin(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc    = N_minus1 - pos_mc;
  auto pos_it    = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  do {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = N_minus1 - left_child_mc;//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);

    //assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
    bool is_right_child_valid = right_child_mc < heap_size;
    decltype(pos_value) larger_value = *left_it, right_value;
    RAI larger_it;
    if (is_right_child_valid && larger_value < (right_value = *right_it)) {
      larger_value = right_value;
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else {
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < larger_value){
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    }
    else
      return ;
  } while (pos_mc < first_in_node) ;
  //Start sifting up the min heap.
  //At this point pos is an In node.
  SizeType parent_hc = ParentNotRoot<SizeType>(pos_hc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_hc == 0 comparison since N > 8 implies this pos_hc > 0.
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

/* This will BiHeapify all nodes in [0, N).
 * Assumes that N is odd.
 */
/*
 * Remark:
 *  This algorithm has complexity O(N). To see why, recall the
 *   argument showing that the heapify operation has O(n) complexity (e.g. as
 *   found on pp. 115 - 116 of "The Algorithm Design Manual" 2nd edition); this
 *   argument generalizes to prove that this algorithm also runs in O(n) times.
 *  The key observation is that the biheap is constructed by adding one node
 *   at a time with this node alternating between a node in a min heap and a
 *   node in the max heap. The complexity of this biheapification is easily
 *   seen to be twice the complexity of the above mentioned heapification
 *   operation plus a constant.
 */
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapify(RAI V, SizeType N, LambdaType lambda) {
  if (N % 3 == 2) { //Then sort the values of the double arrow.
    SizeType pmin_double_arrow_end_hc = (N - 2) / 3;
    auto pmin_double_arrow_end_it = V + lambda(N, pmin_double_arrow_end_hc);
    auto pmax_double_arrow_end_it = V + lambda(N, 2 * pmin_double_arrow_end_hc + 1);
    if (*pmax_double_arrow_end_it < *pmin_double_arrow_end_it)
      std::iter_swap(pmin_double_arrow_end_it, pmax_double_arrow_end_it);
  }
  if(N < 3)
    return ;
  SizeType heap_size     = HeapSize(N);
  SizeType first_in_node = N - heap_size;
  SizeType N_minus1      = N - 1;
  //Ignore all In arrows.
  SizeType last_node_in_biheap_hc  = heap_size - 1;
  SizeType first_node_in_biheap_hc = N - heap_size; //= N_minus1 - last_node_in_biheap_hc;
  while (first_node_in_biheap_hc > 0) {
    BiHeapifySiftFromMinToMax<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, --first_node_in_biheap_hc, //=pos_hc
        last_node_in_biheap_hc, lambda);
    (void)++last_node_in_biheap_hc;
    BiHeapifySiftFromMaxToMin<RAI, SizeType, LambdaType>(V, N, N_minus1,
        heap_size, first_in_node, first_node_in_biheap_hc, //= pos_mc
        first_node_in_biheap_hc, lambda);
  }
  return ;
}

/*
 * ================== END: Definition of lambda version of BiHeapify ====================
 */


template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftUpMaxHeapUnboundedMC(RAI V, SizeType N, SizeType N_minus1,
                            SizeType pos_mc, LambdaType lambda) {
  SizeType parent;
  auto pos_it = V + lambda(N, N_minus1 - pos_mc);
  auto pos_value = *pos_it;
  while (pos_mc > 0) {
    parent = ParentNotRoot<SizeType>(pos_mc);
    auto parent_it = V + lambda(N, N_minus1 - parent);
    if (*parent_it < pos_value) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename SizeType = std::size_t, class LambdaType>
inline void SiftUpMinHeapUnboundedHC(RAI V, SizeType N, SizeType N_minus1,
                            SizeType pos_hc, LambdaType lambda) {
  SizeType parent;
  auto pos_it = V + lambda(N, pos_hc);
  auto pos_value = *pos_it;
  while (pos_hc > 0) {
    parent = ParentNotRoot<SizeType>(pos_hc);
    auto parent_it = V + lambda(N, parent);
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  }
  return ;
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
    if (*parent_it < pos_value) {
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
    if (*parent_it < pos_value) {
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
      if (*V_1 < *V_0)
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
  if (*(V + lambda(N, (N - 1) - parent_of_pure_min_heap_double_arrow_node_hc))
      < *(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc)))
    std::iter_swap(V + lambda(N, parent_of_pure_min_heap_double_arrow_node_hc),
        V + lambda(N, (N - 1) - parent_of_pure_min_heap_double_arrow_node_hc));
  return ;
}

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define an Fused BiHeap.
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
      if (*(V + lambda(N, (right_child - 1))) < parent_value)
        return false;
      //Check that the parent and right child satisfy the min heap condition.
      if (*(V + lambda(N, right_child)) < parent_value) {
        if (!(is_N_mod_3_equal_to_2 && right_child == pmin_node_of_double_arrow))
          return false;
      }
    }
    //If the min heap's last non-In element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < first_in_node_hc
          && *(V + lambda(N, left_child)) < *(V + lambda(N, i)))
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
 * ================== START: Definition of FusedBiHeapifySift ====================
 */


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
      if (*maxh_parent_of_pmax_double_arrow_end_it < pos_value) {
        std::iter_swap(pos_it, maxh_parent_of_pmax_double_arrow_end_it);
        SiftUpMaxHeapUnboundedMC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      } else if (pos_value < *minh_parent_of_pmax_double_arrow_end_it) {
        std::iter_swap(pos_it, minh_parent_of_pmax_double_arrow_end_it);
        SiftUpMinHeapUnboundedHC<RAI, SizeType, LambdaType>(V, N, N_minus1, minh_parent_of_pmin_double_arrow_end_hc, lambda);
      }
      return ; //Return should be inside this if statement.
    } else if (pos_hc == pmax_double_arrow_end_hc && is_pmin_double_arrow_end_in_F) { //If !is_pmin_double_arrow_end_in_F then you can just sift up.
      if (*maxh_parent_of_pmax_double_arrow_end_it < pos_value) {
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
    if (*maxh_parent_of_pos_it < pos_value) {
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
    if (*maxh_parent_of_pos_it < pos_value) {
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



template<class ValueType, typename SizeType = std::size_t>
class BiQueue {
public:
  std::vector<ValueType> vec_;
  SizeType num_elements_; //The number of elements currently in the Fused BiHeap.
  SizeType N_;            //The size of the BiHeap that induced this Fused BiHeap.
                          //This is the maximum number of elements that the Fused BiHeap
                          // can hold before it needs to be resized to allow for
                          // the insertion of another element.
                          //N should always be even and non-zero.
  //SizeType N_minus_1;   //Frequently used value.
  //SizeType half_of_N;   //Frequently used value.
  SizeType F_first_hc_, F_last_hc_;
  //Invariants when num_elements > 0:
  // (1) vec_ holds a fused BiHeap (which could possibly also be a BiHeap).
  // (2) num_elements_ <= N_.
  // (3) num_elements_ <= vec_.size() or else num_elements <= 1.
  // (4) if lambda = get_index_lambda(), then the value of the
  //     node whose min heap coordinate is pos_hc is located
  //     at vec_[lambda(N_, i)] in memory.
  // (5) The graph is a BiHeap if and only if N_ == num_elements_.
  // (6) If F_first_hc_ < F_last_hc_ then N_ > 3 and
  //     the graph is a fused BiHeap graph on N_ nodes
  //     with each of the (F_last_hc_ + 1 - F_first_hc_) nodes
  //     F_first_hc_, F_first_hc_ + 1, ...., F_last_hc_
  //     fused.
  // (7) If F_last_hc_ < F_first_hc_ then it holds a BiHeap.
  // (8) Either F_last_hc_ < F_first_hc_
  //     or else Flip(F_last_hc_) + 1 >= F_first_hc_ >= Flip(F_last_hc_)
  //     (b) if num_elements_ == 1 then F_last_hc_ == F_first_hc_.
  //     (c) if F_last_hc_ == F_first_hc_ then either
  //         (i) num_elements_ == 1, in which case N_ == 1 or N_ == 2
  //             and the graph is a BiHeap (if N_ == 1) or else an
  //             Fused BiHeap (if N_ == 2) with F_last_hc_ = 1, or else
  //         (i) num_elements_ > 2, in which case the graph is a
  //             fused BiHeap graph fused at node F_first_hc_.
  // (9) N_ >= 2 where if num_elements_ <= 1 then N_ == 2.
  // (10) If vec_[i] does NOT store the value of any node then it will store 0
  //      and otherwise, it will store a positive number.
  //      - This is so that if the algorithm touches a node that is not
  //        in the fused BiHeap then this will be indicated by having
  //        a zero value where there should be a positive value.
  // (11) If num_elements_ > 1 then vec_[0] stores the minimum and
  //      vec_[1] stores the maximum. If num_elements_ == 1 then vec_[0]
  //      stores both the minimum and the maximum.
  // (12) For any index i, vec_[i] stores the value of a node in the
  //      fused BiHeap if and only if i < num_elements_.

  BiQueue() : num_elements_(0), N_(2), F_first_hc_(0), F_last_hc_(1) {
  }

  template<class Iterator>
  BiQueue(Iterator start, Iterator one_past_end) {
    num_elements_ = std::distance(start, one_past_end);
    if (num_elements_ <= 2) {
      N_            = 2;
      if (num_elements_ == 0) {
        F_first_hc_ = 0;
        F_last_hc_  = 1;
        return ;
      }
      vec_.resize(num_elements_);
      if (num_elements_ == 1) {
        F_first_hc_ = 1;
        F_last_hc_  = 1;
        vec_[0] = *start;
      } else {
        F_first_hc_ = 1;
        F_last_hc_  = 0;
        ValueType value_0 = *start;
        start++;
        ValueType value_1 = *start;
        if (value_1 < value_0) {
          vec_[0] = value_1;
          vec_[1] = value_0;
        }
        else {
          vec_[0] = value_0;
          vec_[1] = value_1;
        }
      }
      return ;
    }
    SizeType num_elements_evened_down = num_elements_ - (num_elements_ % 2);
    SizeType new_N = parent_heap_size(num_elements_evened_down);
    N_             = new_N;
    reserve(new_N);
    for (SizeType i = 0 ; start != one_past_end && i < num_elements_; i++, start++)
      vec_[i] = *start;
    //Unnecessary code, useful for checking correctness and
    // development of these algorithms but limits ValueType to numeric objects.
    //for (SizeType i = num_elements_; i < vec_.size(); i++)
    //  vec_[i] = static_cast<ValueType>(0); //Fill with 0 each value that doesn't store the value of a node in the fused BiHeap.
    F_first_hc_ = (num_elements_ + 1) / 2;
    F_last_hc_  = (N_ - 1) - (num_elements_ / 2);
    call_fused_biheapify();
    return ;
  }

  //Copy constructor
  BiQueue(const BiQueue<ValueType> &biq) {
    vec_ = biq.vec_;
    num_elements_ = biq.num_elements_;
    N_ = biq.N_;
    F_first_hc_ = biq.F_first_hc_;
    F_last_hc_ = biq.F_last_hc_;
  }

  //Assumes that vec_.size() is sufficiently large to store
  // the new parent BiHeap.
  inline void fused_biheapify() {
    N_ = parent_heap_size(num_elements_);
    if (num_elements_ <= 0)
      return ;
    F_first_hc_ = (num_elements_ + 1) / 2;        //Note that if num_elements == 1 then F_first_hc_ == 1, as desired.
    F_last_hc_  = (N_ - 1) - (num_elements_ / 2); //Note that if num_elements == 1 then F_last_hc_  == 1, as desired.
    if (num_elements_ > 1)
      call_fused_biheapify();
    return ;
  }

  inline SizeType biheap_size() const {
    return N_;
  }

  //Makes the num_elements_ into a BiHeap on num_elements_ nodes.
  inline void biheapify() {
    if (num_elements_ > 1) {
      call_biheapify(num_elements_);
    }
    return ;
  }

  inline void call_fused_biheapify() {
    call_fused_biheapify(N_);
    return ;
  }

  inline void call_fused_biheapify(SizeType new_N) {
    auto lambda = get_index_lambda(new_N);
    FusedBiHeapify<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(),
                                                                new_N, F_first_hc_, F_last_hc_, lambda);
    return ;
  }

  inline void call_fused_biheapify_sift(SizeType pos_hc) {
    auto lambda = get_index_lambda();
    FusedBiHeapifySift<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(),
                                                               N_, pos_hc, F_first_hc_, F_last_hc_, lambda);
    return ;
  }

  inline void call_biheapify() {
    call_biheapify(N_);
    return ;
  }

  inline void call_biheapify(SizeType new_N) {
    auto lambda = get_index_lambda(new_N);
    FusedBiHeapify<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(), new_N, new_N, 0, lambda);
    //Or equivalently:
    //BiHeapify<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(), new_N, lambda);
    return ;
  }

  inline SizeType capacity() const {
    return vec_.size();
  }

  inline void clear() noexcept {
     num_elements_ = 0;
     return ;
  }

  inline std::vector<ValueType>& data() {
    return vec_;
  }

  inline std::vector<ValueType>& data() const {
    return vec_;
  }

  inline bool empty() const {
    return num_elements_ <= 0;
  }

  inline void expand_parent_biheap() {
    SizeType new_N = parent_heap_size(N_);
    assert(N_ >= 2 && N_ % 2 == 0 && new_N > 2 && new_N % 2 == 0);
    assert(num_elements_ % 2 == 0 && N_ < new_N && num_elements_ <= N_ && num_elements_ + 1 < new_N);
    reserve(new_N);
    N_          = new_N;
    F_first_hc_ = (num_elements_ + 1) / 2;
    F_last_hc_  = (new_N - 1) - (num_elements_ / 2);
    /* //Unnecessary code, useful for checking correctness and
       // development but limits ValueType to numeric objects.
    auto index_lambda = get_index_lambda();
    //Fill the nodes that are not to be touched with 0's.
    for (SizeType i_hc = F_first_hc_; i_hc <= F_last_hc_; i_hc++) {
      SizeType vec_index_of_node_i = index_lambda(N_, i_hc);
      vec_[vec_index_of_node_i] = static_cast<ValueType>(0);
    }
    */
    call_fused_biheapify();
    return ;
  }

  inline auto get_index_lambda() const {
    return get_index_lambda(N_);
  }

  inline auto get_index_lambda(SizeType new_N) const {
    auto twice_new_N_minus1 = 2 * new_N - 1; // = 2 * (new_N - 1) + 1
    auto half_new_N = new_N / 2;
    return [twice_new_N_minus1, half_new_N](SizeType N_local, SizeType pos_hc) -> SizeType {
      SizeType twice_pos_hc = 2 * pos_hc;
      if (pos_hc < half_new_N)
        return twice_pos_hc;
      else
        return twice_new_N_minus1 - twice_pos_hc;// = 2 * (new_N_minus1 - pos_hc) + 1 = 2 * pos_mc + 1
    };
  }

  //To see that insert() has amortized O(log(N_)) complexity, note
  // that after this structure is first created, if one were to
  // continue calling insert(), then it would FusedBiHeapifySift()
  // (an O(log N) operation) approximately (i.e. plus or minus 2) N_ / 3
  // times before it would have to call FusedBiHeapify(), where since
  // FusedBiHeapify() is an O(N) operation, there is some constant C
  // such that FusedBiHeapify() performs no more than C * N_ operations.
  //If each call to FusedBiHeapifySift() performs at most D log N
  // operations, then at most D*((N_ / 3) + 2)*log(N_) + C N_ operations
  // will have been performed. Dividing by N_, shows that the amortized
  // complexity is O(log(N_)).
  //Note that since num_elements_ always satisfies
  // (2 * N_) / 3 - 4 <= num_elements_ <= N_,
  // insert() also has amortized O(log(num_elements_)) complexity.
  inline void insert(ValueType value) {
    if (num_elements_ <= 1) {
      assert(N_ == 2);
      reserve(num_elements_ + 1);
      F_first_hc_ = 1;
      F_last_hc_  = 1;
      if (num_elements_ == 0) {
        vec_[0] = value;
      } else if (num_elements_ == 1) {
        if (vec_[0] < value)
          vec_[1] = value;
        else {
          vec_[1] = vec_[0];
          vec_[0] = value;
        }
      }
      num_elements_++;
      return ;
    }
    if (num_elements_ == N_) {//If we need to make room for the new value
      expand_parent_biheap();
      assert((N_ == 2 || N_ % 3 != 2) && num_elements_ < N_ && F_first_hc_ < F_last_hc_);
    }
    assert(num_elements_ < N_ && F_first_hc_ <= F_last_hc_);
    SizeType placement_node_hc;
    //If there are as many nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
    if (num_elements_ % 2 == 0) { //Then "un-fuse" node F_first_hc_ to place the value there.
      placement_node_hc = F_first_hc_;
      assert(F_first_hc_ <= (N_ - 1) - F_last_hc_);
      F_first_hc_++;
    } else { //Then "un-fuse" node F_last_hc_ to place the value there.
      placement_node_hc = F_last_hc_;
      assert(F_first_hc_ > (N_ - 1) - F_last_hc_);
      F_last_hc_--;
    }
    vec_[num_elements_] = value;
    num_elements_++;
    call_fused_biheapify_sift(placement_node_hc);
    return ;
  }

  //Assumes that num_elements_ > 0.
  inline ValueType max() const {
    return vec_[num_elements_ > 1];
    /* The above is short for:
    if (num_elements_ == 1)
      return vec_[0];
    return vec_[1];*/
  }

  //Assumes that num_elements_ > 0.
  inline ValueType min() const {
    return vec_[0];
  }

  //Assumes that even_N is even and positive.
  static inline SizeType parent_heap_size_even(SizeType even_N) {
    SizeType half_N             = even_N / 2;
    SizeType three_times_half_N = 3 * half_N;
    if (half_N % 2 == 1)
      return three_times_half_N + 1;
    else
      return three_times_half_N;
  }

  //Assumes that even_N is even and positive.
  static inline SizeType parent_heap_size(SizeType N) {
    if (N <= 1)
      return 2;
    else
      return parent_heap_size_even(N - (N % 2));
  }

  //Assumes that num_elements_ > 0.
  //pop_index should be either 0 or 1 for the min or max, respectively.
  //The argument that this function has amortized O(log(N_)) complexity
  // is analogous to the argument used to show that insert() also
  // has amortized O(log(N_)) complexity.
  void PopMinOrMax(SizeType pop_index) {
    assert(pop_index == 0 || pop_index == 1);
    if (num_elements_ == 1) {
      num_elements_ = 0;
      F_first_hc_   = 0;
      F_last_hc_    = 1;
      N_            = 2;
      return ;
    } else if (num_elements_ == 2) {
      num_elements_ = 1;
      if (pop_index == 0)
        std::iter_swap(vec_.begin(), vec_.begin() + 1); //Swap the min and max.
      F_first_hc_   = 1;
      F_last_hc_    = 1;
      N_            = 2;
      return ;
    }
    SizeType heap_size = HeapSize<SizeType>(N_);
    SizeType first_in_node = N_ - heap_size;
    if (F_first_hc_ == first_in_node) { //If we can not remove any more In nodes from the fused BiHeap.
      assert(num_elements_ % 2 == 0);
      N_ = num_elements_;
      call_biheapify(); //BiHeapify it.
      F_first_hc_  = ((N_ + 1) / 2);
      F_last_hc_   = F_first_hc_;
    } else {
      if (num_elements_ % 2 == 1) { //If there are more nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
        assert(F_first_hc_ > (N_ - 1) - F_last_hc_);
        F_first_hc_--;
      } else {
        assert(F_first_hc_ <= (N_ - 1) - F_last_hc_);
        F_last_hc_++;
      }
    }
    num_elements_--;
    SizeType pop_node_hc = 0;
    if (pop_index == 1)     //If we're to pop the max node
      pop_node_hc = N_ - 1; //then this node has min heap coordinate N_ - 1.
    std::iter_swap(vec_.begin() + pop_index, vec_.begin() + num_elements_);
    call_fused_biheapify_sift(pop_node_hc); //Sift the element into place.
    return ;
  }

  //Assumes that num_elements_ > 0.
  inline void popmax() {
    PopMinOrMax(1);
    return ;
  }

  //Assumes that num_elements_ > 0.
  inline void popmin() {
    PopMinOrMax(0);
    return ;
  }

  //Expands the size of the container to be at least new_expanded_vec_size.
  //Does not call BiHeapify or Fused BiHeapify.
  inline void reserve(SizeType new_expanded_vec_size) {
    if (capacity() < new_expanded_vec_size)
      vec_.resize(new_expanded_vec_size);
    return ;
  }

  //Resizes the container to size hold num_elements_ elements and
  // calls call_biheapify() if necessary.
  inline void shrink_to_fit(bool should_call_biheapify = true) {
    if (vec_.size() > num_elements_) {
      vec_.resize(num_elements_);
      N_ = num_elements_;
      if (should_call_biheapify)
        call_biheapify();
    }
    return ;
  }

  inline SizeType size() const {
    return num_elements_;
  }

  //Returns the second largest element in the BiQueue.
  //This is an O(const) function.
  //Assumes that num_elements_ >= 1.
  //Although note a standard feature of double ended queues,
  // it is included because of how easily this element can be found.
  //If there is only 1 element in the BiQueue then it returns that element.
  inline ValueType second_largest() const {
    if (num_elements_ <= 5) {
      if (num_elements_ <= 2) {
        return vec_[0];
      } else if (num_elements_ == 3) {
        return vec_[2];
      } else if (num_elements_ == 4) {
        return (vec_[2] < vec_[3]) ? vec_[3] : vec_[2];
      } else { //Else num_elements_ == 5
        return (vec_[3] < vec_[4]) ? vec_[4] : vec_[3];
      }
    } else {
      return (vec_[3] < vec_[5]) ? vec_[5] : vec_[3];
    }
  }

  //Returns the second smallest element in the BiQueue.
  //This is an O(const) function.
  //Assumes that num_elements_ >= 1.
  //Although note a standard feature of double ended queues,
  // it is included because of how easily this element can be found.
  //If there is only 1 element in the BiQueue then it returns that element.
  inline ValueType second_smallest() const {
    if (num_elements_ <= 4) {
      /* The following if statement if short for:
      if (num_elements_ == 1) {
        return vec_[0];
      } else if (num_elements_ == 2) {
        return vec_[1];
      } else if (num_elements_ == 3) {
        return vec_[2];
      }*/
      if (num_elements_ <= 3) {
        return vec_[num_elements_ - 1];
      } else { //Else num_elements_ == 4
        return (vec_[2] < vec_[3]) ? vec_[2] : vec_[3];
      }
    } else {
      return (vec_[2] < vec_[4]) ? vec_[2] : vec_[4];
    }
  }
};

} //End namespace biqueue


#endif /* BIQUEUE_STANDALONE_H_ */
