/*
 * biqueue_with_verification.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIQUEUE_WITH_VERIFICATION_H_
#define BIQUEUE_WITH_VERIFICATION_H_

#include <algorithm>
#include <vector>

#include "biheapify.h"
#include "almost_biheapify.h"
#include "biheap_sift.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <random>


template<class ValueType, typename SizeType = std::size_t>
class BiQueueWithVerification {

public:
  std::vector<ValueType> vec_;
  SizeType num_elements_; //The number of elements currently in the almost BiHeap.
  SizeType N_;            //The size of the BiHeap that induced this almost BiHeap.
                          //This is the maximum number of elements that the almost BiHeap
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
  //             almost BiHeap (if N_ == 2) with F_last_hc_ = 1, or else
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

  template<class Iterator>
  BiQueueWithVerification(Iterator start, Iterator one_past_end) {
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
        auto value_0 = *start;
        start++;
        auto value_1 = *start;
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
    N_          = new_N;
    expand_container(new_N);
    for (SizeType i = 0 ; start != one_past_end && i < num_elements_; i++, start++)
      vec_[i] = *start;
    for (SizeType i = num_elements_; i < vec_.size(); i++)
      vec_[i] = static_cast<ValueType>(0); //Fill with 0 each value that doesn't store the value of a node in the fused BiHeap.
    F_first_hc_ = (num_elements_ + 1) / 2;
    F_last_hc_  = (N_ - 1) - (num_elements_ / 2);
    call_almost_biheapify();
    return ;
  }

  inline void call_almost_biheapify() {
    call_almost_biheapify(N_);
    return ;
  }

  inline void call_almost_biheapify(SizeType new_N) {
    auto lambda = get_index_lambda(new_N);
    FusedBiHeapify<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(),
                                                                new_N, F_first_hc_, F_last_hc_, lambda);
    return ;
  }

  inline void call_almost_biheapify_sift(SizeType pos_hc) {
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
    BiHeapify<typename std::vector<ValueType>::iterator, SizeType, decltype(lambda)>(vec_.begin(), new_N, lambda);
    return ;
  }

  inline SizeType capacity() const {
    return vec_.size();
  }

  inline void expand_parent_biheap() {
    SizeType new_N = parent_heap_size(N_);
    assert(N_ >= 2 && N_ % 2 == 0 && new_N > 2 && new_N % 2 == 0);
    assert(num_elements_ % 2 == 0 && N_ < new_N && num_elements_ <= N_ && num_elements_ + 1 < new_N);
    expand_container(new_N);
    N_          = new_N;
    F_first_hc_ = (num_elements_ + 1) / 2;
    F_last_hc_  = (new_N - 1) - (num_elements_ / 2);
    auto index_lambda = get_index_lambda();
    //Fill the nodes that are not to be touched with 0's.
    for (SizeType i_hc = F_first_hc_; i_hc <= F_last_hc_; i_hc++) {
      SizeType vec_index_of_node_i = index_lambda(N_, i_hc);
      vec_[vec_index_of_node_i] = static_cast<ValueType>(0);
    }
    call_almost_biheapify();
    return ;
  }

  inline void expand_container(SizeType new_expanded_vec_size) {
    if (capacity() < new_expanded_vec_size)
      vec_.resize(new_expanded_vec_size);
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
      expand_container(num_elements_ + 1);
      F_first_hc_ = 1;
      F_last_hc_  = 1;
      if (num_elements_ == 0) {
        vec_[0] = value;
      } else if (num_elements_ == 1) {
        if (value > vec_[0])
          vec_[1] = value;
        else {
          vec_[1] = vec_[0];
          vec_[0] = value;
        }
      }
      num_elements_++;
      return ;
    }
    if (num_elements_ == N_) {//If we need to make more room.
      expand_parent_biheap();
      assert((N_ == 2 || N_ % 3 != 2) && num_elements_ < N_ && F_first_hc_ < F_last_hc_);
    }
    assert(num_elements_ < N_ && F_first_hc_ <= F_last_hc_);
    SizeType placement_index, placement_node_hc;
    //If there are as many nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
    if (F_first_hc_ <= (N_ - 1) - F_last_hc_) { //Then "un-fuse" node F_first_hc_ to place the value there.
      placement_node_hc = F_first_hc_;
      F_first_hc_++;
    } else { //Then "un-fuse" node F_last_hc_ to place the value there.
      placement_node_hc = F_last_hc_;
      F_last_hc_--;
    }
    {
      auto index_lambda = get_index_lambda();
      placement_index = index_lambda(N_, placement_node_hc);
      vec_[placement_index] = value;
    }
    num_elements_++;
    call_almost_biheapify_sift(placement_node_hc);
    return ;
  }

  //Assumes that vec_.size() is sufficiently large to store
  // the new almost BiHeap.
  inline void almost_biheapify() {
    N_ = parent_heap_size(num_elements_);
    if (num_elements_ <= 0)
      return ;
    F_first_hc_ = (num_elements_ + 1) / 2;        //Note that if num_elements == 1 then F_first_hc_ == 1, as desired.
    F_last_hc_  = (N_ - 1) - (num_elements_ / 2); //Note that if num_elements == 1 then F_last_hc_  == 1, as desired.
    if (num_elements_ > 1)
      call_almost_biheapify();
    return ;
  }

  //Assumes that vec_.size() is sufficiently large to store
  // the new almost BiHeap.
  inline void biheapify() {
    if (num_elements_ > 1) {
      call_biheapify();
    }
    return ;
  }

  inline SizeType biheap_size() const {
    return N_;
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
    SizeType swap_index, swap_node_hc;
    if (F_first_hc_ == first_in_node) { //If we can not remove any more In nodes from the fused BiHeap.
      assert(num_elements_ % 2 == 0);
      N_ = num_elements_;
      call_biheapify(); //BiHeapify it.
      F_first_hc_  = ((N_ + 1) / 2);
      F_last_hc_   = F_first_hc_;
      swap_node_hc = F_first_hc_;
    } else {
      if (F_first_hc_ > (N_ - 1) - F_last_hc_) { //If there are more nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
        F_first_hc_--;
        swap_node_hc = F_first_hc_;
      } else {
        F_last_hc_++;
        swap_node_hc = F_last_hc_;
      }
    }
    num_elements_--;
    SizeType pop_node_hc = 0;
    if (pop_index == 1)     //If we're to pop the max node
      pop_node_hc = N_ - 1; //then this node has min heap coordinate N_ - 1.
    auto index_lambda = get_index_lambda();
    //Check that node pop_node_hc is stored in the vector vec at index pop_index_.
    assert(pop_index == index_lambda(N_, pop_node_hc));
    swap_index = index_lambda(N_, swap_node_hc);
    std::iter_swap(vec_.begin() + pop_index, vec_.begin() + swap_index);
    //Indicate that this node is no longer in the almost BiHeap by setting it to 0.
    vec_[swap_index] = static_cast<ValueType>(0);
    //Sift the element into place.
    call_almost_biheapify_sift(pop_node_hc);
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

  inline void shrink_to_fit() const {
    if (size() > N_)
      vec_.resize(size());
    return ;
  }

  inline SizeType size() const {
    return num_elements_;
  }

  //This method is for verifying the functionality of the algorithm.
  inline auto GetMaxBySearchingThroughVector() const {
    auto max_value = vec_[0];
    for (SizeType i = 0; i < num_elements_; i++) {
      if (vec_[i] > max_value)
        max_value = vec_[i];
    }
    return max_value;
  }

  //This method is for verifying the functionality of the algorithm.
  inline auto GetMinBySearchingThroughVector() const {
    auto min_value = vec_[0];
    for (SizeType i = 0; i < num_elements_; i++) {
      if (vec_[i] < min_value)
        min_value = vec_[i];
    }
    return min_value;
  }

  //This method is for verifying the functionality of the algorithm.
  //If num_elements > 1 then after popping the max, what should
  // the next maximum be? This function gives the answer.
  inline auto GetNextMaxBySearchingThroughVector() const {
    auto max_value = GetMaxBySearchingThroughVector();
    auto second_largest_value = static_cast<ValueType>(0);
    SizeType num_of_times_max_value_is_repated = 0;
    for (SizeType i = 0; i < num_elements_; i++) {
      if (vec_[i] == max_value)
        num_of_times_max_value_is_repated++;
      else if (vec_[i] > second_largest_value)
        second_largest_value = vec_[i];
    }
    if (num_of_times_max_value_is_repated > 1) //If the maximum value is repeated at least twice.
      return max_value; //Then the next maximum value should be the same as the current maximum value.
    else
      return second_largest_value;
  }

  //This method is for verifying the functionality of the algorithm.
  //If num_elements > 1 then after popping the min, what should
  // the next minimum be? This function gives the answer.
  inline auto GetNextMinBySearchingThroughVector() const {
    auto min_value = GetMinBySearchingThroughVector();
    auto second_smallest_value = std::numeric_limits<ValueType>::max();
    SizeType num_of_times_min_value_is_repated = 0;
    for (SizeType i = 0; i < num_elements_; i++) {
      if (vec_[i] == min_value)
        num_of_times_min_value_is_repated++;
      else if (vec_[i] < second_smallest_value)
        second_smallest_value = vec_[i];
    }
    if (num_of_times_min_value_is_repated > 1) //If the minimum value is repeated at least twice.
      return min_value; //Then the next minimum value should be the same as the current minimum value.
    else
      return second_smallest_value;
  }

  inline bool AreMinAndMaxCorrect() const {
    if (num_elements_ == 0)
      return true; //Then there's nothing to check.
    auto min_by_searching = GetMinBySearchingThroughVector();
    auto max_by_searching = GetMaxBySearchingThroughVector();
    //If the condition inside assert() is false, then execution
    // immediately terminates in error.
    assert(min() == min_by_searching && min() != 0 && min_by_searching != 0); //Should never be 0 since only positive values are inserted.
    assert(max() == max_by_searching && max() != 0 && max_by_searching != 0);
    return true;
  }
};

#include <numeric>
#include <random>


template<class T>
T GetRandomPositiveIntForBiQueue(T last_value = 1) {
  std::random_device rnd_device;
  std::mt19937 generator(rnd_device());
  //Using last_value + 10 makes repeated values more likely,
  // which are important special cases to check.
  std::uniform_int_distribution<T> dist(1, last_value + 10);
  return dist(generator);
}


//This function will randomly add and removes values (by popping the min or the max)
// to and from a double ended BiHeap queue. After check each operation it check the
// vector to make sure that it has the correct min and max.
void PerformSingleBiQueueVerificationTest(int num_random_iterations = 64) {
  std::vector<int> vec({5, 4, 3, 2, 1});
  BiQueueWithVerification<int> biq(vec.begin(), vec.end());
  biq.AreMinAndMaxCorrect();
  int next_value = vec.size() + 1;
  for (int i = 0; i < num_random_iterations; i++) {
    if (biq.size() == 0 || GetRandomPositiveIntForBiQueue(next_value) % 2 == 0) {
      auto value = (GetRandomPositiveIntForBiQueue(next_value++)) + 1;
      biq.insert(value);
      assert(biq.biheap_size() == 2 || biq.biheap_size() % 3 != 2); //Make sure we don't have to deal with a double arrow.
      biq.AreMinAndMaxCorrect();
    }
    if (biq.size() > 0 && GetRandomPositiveIntForBiQueue(next_value) % 5 == 0) {
      const auto original_size = biq.size(); //The number of elements in the BiQueue.
      //Before changing anything, what should be next minimum value be?
      const auto next_min = biq.GetNextMinBySearchingThroughVector();
      biq.popmin();
      assert(biq.size() == original_size - 1);
      assert(biq.biheap_size() == 2 || biq.biheap_size() % 3 != 2); //Make sure we don't have to deal with a double arrow.
      if (biq.size() > 0)
        biq.AreMinAndMaxCorrect();
      if (original_size > 1) //If there were at least two elements before the pop.
        assert(next_min == biq.min());
    }
    if (biq.size() > 0 && GetRandomPositiveIntForBiQueue(next_value) % 5 == 0) {
      const auto original_size = biq.size(); //The number of elements in the BiQueue.
      //Before changing anything, what should be next maximum value be?
      const auto next_max = biq.GetNextMaxBySearchingThroughVector();
      biq.popmax();
      assert(biq.size() == original_size - 1);
      assert(biq.biheap_size() == 2 || biq.biheap_size() % 3 != 2); //Make sure we don't have to deal with a double arrow.
      if (biq.size() > 0)
        biq.AreMinAndMaxCorrect();
      if (original_size > 1) //If there were at least two elements before the pop.
        assert(next_max == biq.max());
    }
  }
  return ;
}

//If num_random_iterations == 0 then the current test number will
// be passed as input to PerformSingleBiQueueVerificationTest().
void BiQueueVerificationTests(bool verbose = true, int number_of_test = 30000, int num_random_iterations = 0) {
  if (verbose) {
    std::cout << "Started PerformSingleBiQueueVerificationTest(). Will perform " << number_of_test << "tests.\n";
    std::cout.flush();
  }
  for (int i = 0; i < number_of_test; i++) {
    int num_random_its = num_random_iterations;
    if (num_random_iterations == 0)
      num_random_its = i;
    if (verbose) {
      std::cout << "Started PerformSingleBiQueueVerificationTest(num_random_iterations = "
                << num_random_its << ") \t";
      std::cout.flush();
    }
    PerformSingleBiQueueVerificationTest(num_random_its);
    if (verbose) {
      std::cout << "Finished PerformSingleBiQueueVerificationTest(num_random_iterations = "
                << num_random_its << ")\n";
      std::cout.flush();
    }
  }
  return ;
}


#endif /* BIQUEUE_WITH_VERIFICATION_H_ */
