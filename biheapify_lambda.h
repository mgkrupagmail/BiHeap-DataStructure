/*
 * biheapify_lambda.h
 *
 * This is a generalization of biheapify.h where BiHeapify()
 *  uses a lambda to access objects.
 * Note that RAI = Random Access Iterator
 *
 *  Created on: Nov 14, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_LAMBDA_H_
#define BIHEAPIFY_LAMBDA_H_

#include <algorithm>

//#define FLIP(a) ((N - 1) - (a))

/*
 * ================== END: Definition of lambda version of BiHeapify ====================
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

/*
 * ================== START: Definition of BiHeapifyJumpMiddle and derived functions ====================
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
inline void BiHeapifyJumpMiddle(RAI V, SizeType distance, SizeType num_nodes_to_biheapify) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (distance == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    BiHeapify<RAI, SizeType>(V, distance);
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
    BiHeapify<RAI, SizeType>(V, num_nodes_to_biheapify, index_lambda);
  } else {
    auto index_lambda = [half_of_num_nodes, distance_minus_num_nodes]
                   (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return i;
      else
        return distance_minus_num_nodes + i; // = (distance - 1) - [(local_N - 1) - i]
    };
    BiHeapify<RAI, SizeType>(V, num_nodes_to_biheapify, index_lambda);
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifyJumpMiddleInwards(RAI V, SizeType original_N, SizeType N) {
  while (N > 1) {
    BiHeapifyJumpMiddle<RAI, SizeType>(V, original_N, N);
    if (N <= 3)
      break ;
    SizeType heap_size     = HeapSize(N);
    SizeType first_in_node = N - heap_size;
    SizeType num_in_nodes  = heap_size - first_in_node;
    //assert(N - num_in_nodes ==  2 * first_in_node);
    N                       = num_in_nodes;
    V                      += first_in_node;
    original_N             -= 2 * first_in_node;
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyJumpMiddle(RAI V, SizeType distance, SizeType num_nodes_to_biheapify, LambdaType lambda) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (distance == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    BiHeapify<RAI, SizeType, LambdaType>(V, distance, lambda);
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
    BiHeapify<RAI, SizeType, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  } else {
    auto index_lambda =
        [half_of_num_nodes, distance_minus_num_nodes, num_nodes_to_biheapify, lambda]
         (SizeType local_N, SizeType i) -> SizeType {
      if (i < half_of_num_nodes)
        return lambda(num_nodes_to_biheapify, i);
      else
        return lambda(num_nodes_to_biheapify, distance_minus_num_nodes + i); // = (distance - 1) - [(local_N - 1) - i]
    };
    BiHeapify<RAI, SizeType, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  }
  return ;
}

template<class RAI, typename SizeType = std::size_t, typename LambdaType>
inline void BiHeapifyJumpMiddleInwards(RAI V, SizeType original_N, SizeType N, LambdaType lambda) {
  while (N > 1) {
    BiHeapifyJumpMiddle<RAI, SizeType, LambdaType>(V, original_N, N, lambda);
    if (N <= 3)
      break ;
    SizeType heap_size     = HeapSize(N);
    SizeType first_in_node = N - heap_size;
    SizeType num_in_nodes  = heap_size - first_in_node;
    //assert(N - num_in_nodes ==  2 * first_in_node);
    N                       = num_in_nodes;
    V                      += first_in_node;
    original_N             -= 2 * first_in_node;
  }
  return ;
}

/*
 * ================== END: Definition of BiHeapifyJumpMiddle and derived functions ====================
 */

/*
 * ================== START: Check if nodes form a BiHeap ====================
 */

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define a BiHeap.
 */
template<class RAI, typename SizeType = std::size_t, typename LambdaType>
bool IsBiHeap(RAI V, SizeType N, LambdaType lambda) {
  if (N <= 1)
    return true;
  SizeType heap_size = HeapSize<SizeType>(N);
  //Check that the nodes V, ..., V + heap_size - 1 form
  // a min heap with min at V. This is half of the biheap condition.
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + lambda(N, i));
      //Check that the parent and left child satisfy the min heap condition.
      if (*(V + lambda(N, right_child - 1)) < parent_value)
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (*(V + lambda(N, right_child)) < parent_value)
        return false;
    }
    //If the min heap's last element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < heap_size
          && *(V + lambda(N, left_child)) < *(V + lambda(N, i)))
        return false;
    }
  }
  //Check that the nodes FLIP(V), ...,
  // FLIP(V + heap_size - 1) form a max heap with max
  // at V + N - 1.
  SizeType N_minus1 = N - 1;
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + lambda(N, N_minus1 - i));
      SizeType mirror_left_child_hc = N_minus1 - (right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, mirror_left_child_hc)))
        return false;

      //right_child_hc = FLIP(right_child_mc) = FLIP(left_child + 1)
      // = N - 1 - (left_child + 1) = mirror_left_child_hc - 1.
      //Check that the parent and right child satisfy the max heap condition.
      if (parent_value < *(V + lambda(N, mirror_left_child_hc - 1)))
        return false;
    }
    //If the max heap's last element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < heap_size
          && *(V + lambda(N, N_minus1 - i)) < *(V + lambda(N, N_minus1 - left_child)))
        return false;
    }
  }
  return true;
}

/*
 * ================== END: Check if nodes form a BiHeap ====================
 */

//#undef FLIP

#endif /* BIHEAPIFY_LAMBDA_H_ */
