/*
 * biheapify_inwards_and_outwards.h
 *
 * The O(N) function BiHeapifyInwardsAndOutwards() appears to, based on
 *  empirical testing, always emplace the median(s) into the middle of
 *  the list while simultaneously arranging the elements so that all
 *  elements before (resp. after) the middle element are less than or
 *  equal to (resp. greater than or equal to) the median value(s).
 *
 *  Created on: Nov 28, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_INWARDS_AND_OUTWARDS_H_
#define BIHEAPIFY_INWARDS_AND_OUTWARDS_H_

#include "biheapify.h"
#include "biheapify_lambda.h"
#include "almost_biheapify.h"
#include "biheapify_flip_ordered.h"

#include <algorithm>
#include <cassert>

//Note that this is an O(N) operation since the single recursive call
// is passed a value, next_num_nodes, for its value of N that satisfies
// next_num_nodes <= 2 * N / 3 and all of the BiHeapify calls are O(N)
// operations.
template<class RAI, typename size_type = std::size_t>
inline size_type BiHeapifyOutwardsFromInNodesSingleRecursion(const RAI V, const size_type original_N, size_type N) {
  if (N < 2)
    return 0;

  BiHeapifyJumpMiddleInwards<RAI, size_type>(V, original_N, N);
  if (N < 4)
    return 1;
  const size_type heap_size     = HeapSize<size_type>(N);
  const size_type first_in_node = N - heap_size;
  size_type start_hc = first_in_node;
  size_type one_past_end_hc = N / 2;
  if (N % 2 == 1) {
    assert(N == original_N);
    one_past_end_hc++; //Then incorporate the middle node.
  }
  BiHeapify<RAI, size_type>(V,                                  one_past_end_hc);
  BiHeapify<RAI, size_type>(V + (original_N - one_past_end_hc), one_past_end_hc);
  if (*(V + (one_past_end_hc - 1)) > *(V + (original_N - one_past_end_hc))) {
    std::iter_swap(V + (one_past_end_hc - 1), V + (original_N - one_past_end_hc));
  }
  const size_type one_past_mid_point = (one_past_end_hc + start_hc) / 2 + 1;
  BiHeapify<RAI, size_type>(V,                                     one_past_mid_point);
  BiHeapify<RAI, size_type>(V + (original_N - one_past_mid_point), one_past_mid_point);
  if (*(V + (one_past_mid_point - 1)) > *(V + (original_N - one_past_mid_point))) {
    std::iter_swap(V + (one_past_mid_point - 1), V + (original_N - one_past_mid_point));
  }

  //Ignore all of B_N's In nodes, which leaves 2 * first_in_nodes for the subsequent
  // recursive call to work on.
  const size_type next_num_nodes    = 2 * first_in_node;
  const size_type next_start_offset = BiHeapifyOutwardsFromInNodesSingleRecursion<RAI, size_type>(V, original_N, next_num_nodes);

  const RAI new_V                        = V + next_start_offset;
  const size_type new_N                  = N - 2 * next_start_offset;
  const size_type first_in_node_of_new_N = new_N - HeapSize<size_type>(new_N);
  const size_type return_value           = next_start_offset + first_in_node_of_new_N;
  if (new_N == 1)
    return return_value;
  const size_type new_original_N         = original_N - 2 * next_start_offset;



  assert(4 * next_start_offset + 2 >= 3 * first_in_node && 2 * next_start_offset >= first_in_node); //But replacing " + 2" with " + 1" leads to false assertions.
  assert(next_start_offset <= first_in_node && next_start_offset < first_in_node + (N % 3));
  assert(next_start_offset < first_in_node || N == 4);
  //NOTE: next_start_offset < first_in_node if and only if N != 4.
  assert(7 * next_start_offset <= 2 * N);
  {
    assert((4 * next_start_offset + 4 >= N));
    //Thus 4 * (next_start_offset + 1) >= N  so next_start_offset >= N / 4.0 - 1
    assert((next_start_offset + 1 >= N / 4));
    //Thus next_start_offset >= FLOOR( N / 4 ) - 1
    if (N != 8) {
      assert((next_start_offset >= N / 4));
      //Thus if N != 8 then next_start_offset >= FLOOR( N / 4 )
    }
  }
  //NOTE: 21 * (first_in_node - next_start_offset) approaches N from above as N gets larger.
  if (first_in_node - next_start_offset != 0) // This happens if and only if N = 4
    assert(N < 21 * (first_in_node - next_start_offset));
  else
    assert(N == 4);
  //Thus, for N != 4, N / 21.0 < first_in_node - next_start_offset

  {
    auto value = 21 * (first_in_node - next_start_offset) - N;
    if (N != 4) {
      assert(1 <= value && value <= 24); //Note that value can equal any integer between (and including) 1 and 24.
      //IMPORTANT: So 1 <= 21 * (first_in_node - next_start_offset) - N <= 24
      //           so 0 <= (first_in_node - next_start_offset) - (N + 1) / 21 <= 23 / 21
      //           so distance(first_in_node - next_start_offset, (N + 1) / 21) <= 23 / 21
      //           so (N + 1) / 21 <= (first_in_node - next_start_offset) <= (N + 24) / 21
      //THUS first_in_node - next_start_offset is approximately = (N + 1) / 21
      //THUS next_start_offset is approx. = first_in_node - (N + 1) / 21, which is approx. = 2 * N / 7
      //THUS new_N = N - 2 * next_start_offset is approx. N - 4 * N / 7 = 3 * N / 7

      //The following establishes a bijection between (value % 3) and (N % 3):
      if (value % 3 == 0) assert(N % 3 == 0);
      if (value % 3 == 1) assert(N % 3 == 2);
      if (value % 3 == 2) assert(N % 3 == 1);
      assert((2 * value) % 3 == N % 3);
      assert((2 * N) % 3 == value % 3);
    }
  }

  BiHeapifyJumpMiddleInwards<RAI, size_type>(new_V, new_original_N, new_N);

  return return_value;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyOutwardsFromInNodesSingleRecursion(RAI V, size_type N) {
  if (N < 2)
    return ;
  BiHeapifyOutwardsFromInNodesSingleRecursion<RAI, size_type>(V, N, N);
  return ;
}


//Note that this IS an O(N) operation. This is despite the fact
// that each call to this function calls itself twice.
//See "BiHeaps and Pivot Selection.pdf" for details.
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsAndOutwardsSingle(RAI V, size_type N) {
  if (N < 8) {
    if (N > 1)
      std::sort(V, V + N);
    return;
  }
  size_type heap_size     = HeapSize(N);
  size_type first_in_node = N - heap_size;
  for (size_type i = 0; i < 2; i++) {
    BiHeapifyOutwardsFromInNodesSingleRecursion<RAI, size_type>(V, N);
    AlmostBiHeapify<RAI, size_type>(V, N);
    BiHeapify<RAI, size_type>(V, N);
    BiHeapifyInwardsAndOutwardsSingle<RAI, size_type>(V + first_in_node, heap_size - first_in_node); //Recursive call.
  }
  return ;
}

//This function is not currently used.
//Empirical testing indicates that is an O(N) function and I strongly suspect
// the the proof that BiHeapifyInwardsAndOutwardsSingle() is O(N) generalizes
// to prove that this function is O(N).
//Returns the offset from the first node, V, that the parent call should use in
// its call to BiHeapify.
template<class RAI, typename size_type = std::size_t>
inline size_type BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion(const RAI V, const size_type original_N, size_type N) {
  if (N < 2)
    return 0;
  const size_type heap_size     = HeapSize<size_type>(N);
  const size_type first_in_node = N - heap_size;

  BiHeapifyJumpMiddleInwards<RAI, size_type>(V, original_N, N);
  if (N < 4)
    return 1;
  size_type start_hc = first_in_node;
  size_type one_past_end_hc = N / 2;
  if (N % 2 == 1) {
    assert(N == original_N);
    one_past_end_hc++; //Then incorporate the middle node.
  }
  if (N == 5) {
    std::sort(V, V + 5);
    return 1;
  }
  BiHeapify<RAI, size_type>(V,                                  one_past_end_hc);
  BiHeapify<RAI, size_type>(V + (original_N - one_past_end_hc), one_past_end_hc);
  if (*(V + (one_past_end_hc - 1)) > *(V + (original_N - one_past_end_hc))) {
    std::iter_swap(V + (one_past_end_hc - 1), V + (original_N - one_past_end_hc));
  }
  const size_type one_past_mid_point = (one_past_end_hc + start_hc) / 2 + 1;
  BiHeapify<RAI, size_type>(V,                                     one_past_mid_point);
  BiHeapify<RAI, size_type>(V + (original_N - one_past_mid_point), one_past_mid_point);
  if (*(V + (one_past_mid_point - 1)) > *(V + (original_N - one_past_mid_point))) {
    std::iter_swap(V + (one_past_mid_point - 1), V + (original_N - one_past_mid_point));
  }

  //Ignore all of B_N's In nodes, which leaves 2 * first_in_nodes to work on.
  const size_type next_num_nodes    = 2 * first_in_node;
  const size_type next_start_offset = BiHeapifyOutwardsFromInNodesSingleRecursion<RAI, size_type>(V, original_N, next_num_nodes);
  //NOTE: next_start_offset is approx. = first_in_node - (N + 1) / 21, which is approx. = 2 * N / 7

  const size_type new_N                  = N - 2 * next_start_offset;
  const size_type first_in_node_of_new_N = new_N - HeapSize<size_type>(new_N);
  const size_type return_value     = next_start_offset + first_in_node_of_new_N;
  if (new_N == 1)
    return return_value;
  const RAI new_V                        = V + next_start_offset;
  const size_type new_original_N         = original_N - 2 * next_start_offset;

  {
    if (first_in_node - next_start_offset != 0) // This happens if and only if N = 4
      assert(N < 21 * (first_in_node - next_start_offset));
    else
      assert(N == 4);
    const auto value = 21 * (first_in_node - next_start_offset) - N;
    if (N != 4)
      assert(1 <= value && value <= 24);
    //Thus next_start_offset is approx. = first_in_node - (N + 1) / 21, which is approx. = 2 * N / 7
    assert(first_in_node >= next_start_offset && 2 * (next_start_offset) >= first_in_node && 2 * (next_start_offset) <= (3 * N) / 2);
  }

  BiHeapifyJumpMiddleInwards<RAI, size_type>(new_V, new_original_N, new_N);
  //Note that if we were to replace the above call with a recursive call
  // to this same function, then this function would not be O(N).
  //So that this function will be O(N), we instead perform the
  // last recursive call on the In nodes' In nodes, which is
  // what we we compute starting with computing first_in_node_of_new_N.

  //Nodes V + next_start_offset, ..., V + (first_in_node - 1), .., and their flips
  size_type N_new_start_to_first_in = 2 * (first_in_node - next_start_offset); //
  if (first_in_node > next_start_offset + 1) { //If there are at least two nodes
    if (first_in_node - next_start_offset != 0) // This happens if and only if N = 4
      assert(N < 21 * (first_in_node - next_start_offset));
    else
      assert(N == 4);
    BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion<RAI, size_type>(new_V, new_original_N, N_new_start_to_first_in);
    BiHeapifyJumpMiddleInwards<RAI, size_type>(new_V, new_original_N, new_N);
  }

  //size_type first_in_node_of_new_N  = new_N - HeapSize<size_type>(new_N);
  //NOTE: new_N is approx. = 3 * N / 7 so HeapSize<size_type>(new_N) is approx. = N / 7
  size_type final_rec_offset       = next_start_offset + first_in_node_of_new_N;
  if (2 * final_rec_offset >= N)
    return return_value;
  size_type final_rec_N            = N - 2 * final_rec_offset;
  size_type final_rec_original_N   = original_N - 2 * final_rec_offset;
  RAI final_rec_V                  = V + final_rec_offset;
  if (final_rec_N > 1) {
    BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion<RAI, size_type>(final_rec_V, final_rec_original_N, final_rec_N);
    //The following inequality will guarantee that this function is O(N).
    assert(21 * (next_num_nodes + N_new_start_to_first_in + final_rec_N) < 20 * (N + 4));
  }
  return return_value;
}

//This function is not currently used.
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion(RAI V, size_type N) {
  if (N < 2)
    return ;
  BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion<RAI, size_type>(V, N, N);
  return ;
}

//This function is not currently used.
//Note that this is an O(N) operation if the same is true of
// BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion.
// See "BiHeaps and Pivot Selection.pdf" for details.
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsAndOutwardsDouble(RAI V, size_type N) {
  if (N < 8) {
    if (N > 1)
      std::sort(V, V + N);
    return;
  }
  size_type heap_size     = HeapSize(N);
  size_type first_in_node = N - heap_size;
  for (size_type i = 0; i < 2; i++) {
    BiHeapifyOutwardsFromInNodesAndInwardsDoubleRecursion<RAI, size_type>(V, N);
    AlmostBiHeapify<RAI, size_type>(V, N);
    BiHeapify<RAI, size_type>(V, N);
    BiHeapifyInwardsAndOutwardsDouble<RAI, size_type>(V + first_in_node, heap_size - first_in_node);
  }
  return ;
}

template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyInwardsAndOutwards(RAI V, size_type N) {
  BiHeapifyInwardsAndOutwardsSingle<RAI, size_type>(V, N); //Note that this is an O(N) operation.
  //BiHeapifyInwardsAndOutwardsDouble<RAI, size_type>(V, N);
  return ;
}


#endif /* BIHEAPIFY_INWARDS_AND_OUTWARDS_H_ */
