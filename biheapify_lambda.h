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
 *  V define a biheap.
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
      (parent = Parent<size_type>(pos_mc)) < last_node_in_biheap_mc)
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
     (parent = Parent<size_type>(pos_mc)) >= last_node_in_biheap_mc);
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
      (parent = Parent<size_type>(pos_hc)) < first_node_in_biheap_hc)
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
     (parent = Parent<size_type>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapifySiftFromMinToMax(RAI V, size_type N,
                                      size_type heap_size,
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_hc,
                                      size_type last_node_in_biheap_hc,
                                      LambdaType lambda) {
  while (pos_hc < first_node_in_mirror_heap) {
    auto left_child_hc  = LeftChild<size_type>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);
    auto pos_it         = V + lambda(N, pos_hc);

    assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
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
                                      size_type first_node_in_mirror_heap,
                                      size_type pos_mc,
                                      size_type first_node_in_biheap_hc,
                                      LambdaType lambda) {
  auto pos_hc = FLIP(pos_mc);
  while (pos_mc < first_node_in_mirror_heap) {
    auto left_child_mc  = LeftChild<size_type>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = FLIP(left_child_mc);//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto pos_it         = V + lambda(N, pos_hc);
    auto left_it        = V + lambda(N, left_child_hc);
    auto right_it       = V + lambda(N, right_child_hc);

    assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
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
  size_type heap_size                  = HeapSize(N);
  size_type first_node_in_mirror_heap  = N - heap_size;
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  size_type last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                         - (N % 3 == 2);
  size_type first_node_in_biheap_hc = FLIP(last_node_in_biheap_hc);
  while (first_node_in_biheap_hc > 0) {
    BiHeapifySiftFromMinToMax<RAI, size_type, LambdaType>(V, N, heap_size,
        first_node_in_mirror_heap,
        --first_node_in_biheap_hc, last_node_in_biheap_hc, lambda);
    BiHeapifySiftFromMaxToMin<RAI, size_type, LambdaType>(V, N, heap_size,
        first_node_in_mirror_heap,
        FLIP(++last_node_in_biheap_hc), first_node_in_biheap_hc, lambda);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t, typename LambdaType>
inline void BiHeapify(RAI V, RAI one_past_last, LambdaType lambda) {
  BiHeapify<RAI, size_type, LambdaType>(V, std::distance(V, one_past_last), lambda);
}

#undef FLIP



#endif /* BIHEAPIFY_LAMBDA_H_ */
