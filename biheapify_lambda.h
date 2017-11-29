/*
 * biheapify_lambda.h
 *
 *  Created on: Nov 14, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_LAMBDA_H_
#define BIHEAPIFY_LAMBDA_H_

#include <algorithm>
#include <cassert>

#define FLIP(a) ((N - 1) - (a))

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define a BiHeap.
 */
template<class RAI, typename size_type = std::size_t, typename LambdaType>
bool IsBiHeap(RAI V, size_type N, LambdaType lambda) {
  if (N <= 1)
    return true;
  size_type heap_size = HeapSize<size_type>(N);
  //Check that the nodes V, ..., V + heap_size - 1 form
  // a min heap with min at V. This is half of the biheap condition.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + lambda(N, i));
      //Check that the parent and left child satisfy the min heap condition.
      if (parent_value > *(V + lambda(N, right_child - 1)))
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (parent_value > *(V + lambda(N, right_child)))
        return false;
    }
    //If the min heap's last element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < heap_size
          && *(V + lambda(N, i)) > *(V + lambda(N, left_child)))
        return false;
    }
  }
  //Check that the nodes FLIP(V), ...,
  // FLIP(V + heap_size - 1) form a max heap with max
  // at V + N - 1.
  {
    size_type i = 0;
    for (size_type right_child; (right_child = RightChild<size_type>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + lambda(N, FLIP(i)));
      size_type mirror_left_child_hc = FLIP(right_child - 1);
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
      size_type left_child;
      if ((left_child = LeftChild<size_type>(i)) < heap_size
          && *(V + lambda(N, FLIP(i))) < *(V + lambda(N, FLIP(left_child))))
        return false;
    }
  }
  return true;
}

//Assumes that pos_mc is a node in the max heap.
//The node with max heap coordinate last_node_in_biheap_mc is a
// node in the BiHeap constructed so far such that if v is any
// node whose max heap coordinate is < last_node_in_biheap_mc,
// then v does NOT belong to the BiHeap constructed so far.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void SiftUpMaxHeapMC(RAI V, size_type N, size_type pos_mc,
                            size_type last_node_in_biheap_mc, LambdaType lambda) {
  size_type parent;
  if (pos_mc == 0 ||
      (parent = ParentNotRoot<size_type>(pos_mc)) < last_node_in_biheap_mc)
    return ;
  auto pos_it = V + lambda(N, FLIP(pos_mc));
  do {
    auto parent_it = V + lambda(N, FLIP(parent));
    if (*pos_it > *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
     (parent = ParentNotRoot<size_type>(pos_mc)) >= last_node_in_biheap_mc);
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void SiftUpMaxHeapHC(RAI V, size_type N, size_type pos_hc,
                            size_type last_node_in_biheap_mc, LambdaType lambda) {
  SiftUpMaxHeapMC<RAI, size_type, LambdaType>(V, N, FLIP(pos_hc), last_node_in_biheap_mc, lambda);
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void SiftUpMinHeapHC(RAI V, size_type N, size_type pos_hc,
                            size_type first_node_in_biheap_hc, LambdaType lambda) {
  size_type parent;
  if (pos_hc == 0 ||
      (parent = ParentNotRoot<size_type>(pos_hc)) < first_node_in_biheap_hc)
    return ;
  auto pos_it = V + lambda(N, pos_hc);
  do {
    auto parent_it = V + lambda(N, parent);
    if (*pos_it < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 &&
     (parent = ParentNotRoot<size_type>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifySiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_in_node,
                                      size_type pos_hc,
                                      size_type last_node_in_biheap_hc,
                                      LambdaType lambda) {
  while (pos_hc < first_in_node) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    auto pos_it         = V + lambda(N, pos_hc);

    //assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    RAI smaller_it;
    if (is_right_child_valid && *right_it < *left_it) {
      smaller_it = right_it;
      pos_hc     = right_child_hc;
    } else {
      smaller_it = left_it;
      pos_hc     = left_child_hc;
    }
    if (*pos_it > *smaller_it)
      std::iter_swap(pos_it, smaller_it);
    else
      return ;
  }
  SiftUpMaxHeapMC<RAI, size_type, LambdaType>(V, N, FLIP(pos_hc),
                                  FLIP(last_node_in_biheap_hc), lambda);
  return ;
}

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifySiftFromMaxToMin(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_in_node,
                                      size_type pos_mc,
                                      size_type first_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_in_node) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = V + lambda(N, pos_hc);
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);

    //assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
    bool is_right_child_valid = right_child_mc < heap_size;
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
  SiftUpMinHeapHC<RAI, size_type, LambdaType>(V, N, pos_hc, first_node_in_biheap_hc, lambda);
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
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapify(RAI V, size_type N, LambdaType lambda) {
  if(N < 2)
    return ;
  size_type heap_size     = HeapSize(N);
  size_type first_in_node = N - heap_size;
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  size_type last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                         - (N % 3 == 2);
  size_type first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    BiHeapifySiftFromMinToMax<RAI, size_type, LambdaType>(V, N, heap_size,
       first_in_node, --first_node_in_biheap_hc, last_node_in_biheap_hc,
       lambda);
    BiHeapifySiftFromMaxToMin<RAI, size_type, LambdaType>(V, N, heap_size,
       first_in_node, FLIP(++last_node_in_biheap_hc), first_node_in_biheap_hc,
       lambda);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapify(RAI V, RAI one_past_last, LambdaType lambda) {
  BiHeapify<RAI, size_type, LambdaType>(V, std::distance(V, one_past_last), lambda);
}

/*
 * ================== START: Definition of BiHeapifyJumpMiddle and derived functions ====================
 */


/* In short, given a list of N nodes: V, ..., V + (N - 1)
 *  this takes the first h := num_nodes_to_biheapify / 2 nodes
 *  and the last h nodes and makes them into a BiHeap.
 * Specifically, if num_nodes_to_biheapify is even then it makes the nodes
 *  V, ..., V + (h-1), V + Flip(h-1), ..., V + (N - 1)
 *  into a BiHeap while if num_nodes_to_biheapify is odd then it makes
 *  V, ..., V + (h-1), V + (N / 2), V + Flip(h-1), ..., V + (N - 1)
 *  into a BiHeap (with minimum at V and maximum at V + (N - 1)).
 * If num_nodes_to_biheapify is odd then it is assumed that N is odd
 *  (since otherwise there is no natural way to distinguish a middle node)
 *  and the middle node V + (N / 2) will become the
 *  middle node this BiHeap of size num_nodes_to_biheapify.
 * Assumes that num_nodes_to_biheapify <= N.
 */
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyJumpMiddle(RAI V, size_type N, size_type num_nodes_to_biheapify) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (N == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    BiHeapify<RAI, size_type>(V, N);
    return ;
  }
  size_type half_of_num_nodes = num_nodes_to_biheapify / 2;
  size_type N_minus_num_nodes = N - num_nodes_to_biheapify;
  if (num_nodes_to_biheapify % 2 == 1 && N % 2 == 1) {
    size_type half_of_N = N / 2;
    //Note that local_N will equal num_nodes_to_biheapify in the subsequent call to BiHeapify.
    auto lambda = [N, half_of_N, half_of_num_nodes, N_minus_num_nodes]
         (size_type local_N, size_type i) -> size_type {
      if (i < half_of_num_nodes)
        return i;
      else if (i != half_of_num_nodes) //So i > half_of_num_nodes
        return N_minus_num_nodes + i; // = (N - 1) - [(local_N - 1) - i]
      else //i == half_of_num_nodes
        return half_of_N;
    };
    BiHeapify<RAI, size_type>(V, num_nodes_to_biheapify, lambda);
  } else {
    auto lambda = [N, half_of_num_nodes, N_minus_num_nodes]
                   (size_type local_N, size_type i) -> size_type {
      if (i < half_of_num_nodes)
        return i;
      else
        return N_minus_num_nodes + i; // = (N - 1) - [(local_N - 1) - i]
    };
    BiHeapify<RAI, size_type>(V, num_nodes_to_biheapify, lambda);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyJumpMiddleInwards(RAI V, size_type original_N, size_type N) {
  while (N > 1) {
    BiHeapifyJumpMiddle<RAI, size_type>(V, original_N, N);
    if (N <= 3)
      break ;
    size_type heap_size     = HeapSize(N);
    size_type first_in_node = N - heap_size;
    size_type num_in_nodes  = heap_size - first_in_node;
    //assert(N - num_in_nodes ==  2 * first_in_node);
    N                       = num_in_nodes;
    V                      += first_in_node;
    original_N             -= 2 * first_in_node;
  }
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyJumpMiddle(RAI V, size_type N, size_type num_nodes_to_biheapify, LambdaType lambda) {
  if (num_nodes_to_biheapify < 2)
    return ;
  if (N == num_nodes_to_biheapify) { //Then there are no middle nodes to jump.
    BiHeapify<RAI, size_type, LambdaType>(V, N, lambda);
    return ;
  }
  size_type half_of_num_nodes = num_nodes_to_biheapify / 2;
  size_type N_minus_num_nodes = N - num_nodes_to_biheapify;
  if (num_nodes_to_biheapify % 2 == 1 && N % 2 == 1) {
    size_type lambda_of_half_of_N = lambda(num_nodes_to_biheapify, N / 2);
    //Note that local_N will equal num_nodes_to_biheapify in the subsequent call to BiHeapify.
    auto index_lambda =
        [N, half_of_num_nodes, N_minus_num_nodes, num_nodes_to_biheapify, lambda_of_half_of_N, lambda]
         (size_type local_N, size_type i) -> size_type {
      if (i < half_of_num_nodes)
        return lambda(num_nodes_to_biheapify, i);
      else if (i != half_of_num_nodes) //So i > half_of_num_nodes
        return lambda(num_nodes_to_biheapify, N_minus_num_nodes + i); // = (N - 1) - [(local_N - 1) - i]
      else //i == half_of_num_nodes
        return lambda_of_half_of_N;
    };
    BiHeapify<RAI, size_type, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  } else {
    auto index_lambda =
        [N, half_of_num_nodes, N_minus_num_nodes, num_nodes_to_biheapify, lambda]
         (size_type local_N, size_type i) -> size_type {
      if (i < half_of_num_nodes)
        return lambda(num_nodes_to_biheapify, i);
      else
        return lambda(num_nodes_to_biheapify, N_minus_num_nodes + i); // = (N - 1) - [(local_N - 1) - i]
    };
    BiHeapify<RAI, size_type, decltype(index_lambda)>(V, num_nodes_to_biheapify, index_lambda);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifyJumpMiddleInwards(RAI V, size_type original_N, size_type N, LambdaType lambda) {
  while (N > 1) {
    BiHeapifyJumpMiddle<RAI, size_type, LambdaType>(V, original_N, N, lambda);
    if (N <= 3)
      break ;
    size_type heap_size     = HeapSize(N);
    size_type first_in_node = N - heap_size;
    size_type num_in_nodes  = heap_size - first_in_node;
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

#undef FLIP

#endif /* BIHEAPIFY_LAMBDA_H_ */
