/*
 * biqueue.h
 *
 *  Created on: Dec 5, 2017
 *
 *  A BiQueue is a Double Ended Queue based around the idea of
 *   the BiHeap and almost BiHeap data structures.
 *
 *  It uses almost BiHeaps to implement a double ended priority
 *   queue. It has amortized O(log N) insertions and
 *   amortized O(log N) deletions.
 *  Given a collection of elements, this double ended priority
 *   queue can be formed using O(N) swaps.
 *
 *  Example of using a BiQueue object:

void BiQueueExample() {
  //Define the object.
  BiQueue<int> biq;
  std::vector<int> vector = { 1, 2, 3, 4, 5, 6 };
  BiQueue<int> biq2(vector.begin(), vector.end());

  //Insert an element into the object
  std::cout << "biq.insert(0) \t\t";
  biq.insert(0)
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  std::cout << "biq.insert(2) \t\t";
  biq.insert(2)
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  std::cout << "biq.insert(1) \t\t";
  biq.insert(1)
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  std::cout << "biq.insert(3) \t\t";
  biq.insert(3)
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  std::cout << "biq.popmax() \t\t";
  biq.popmax()
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  std::cout << "biq.popmin() \t\t";
  biq.popmin()
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;

  bool should_pop_max = true;
  std::cout << "biq.PopMinOrMax(should_pop_max) \t\t";
  //if should_pop_max == true then pop the max, otherwise pop the min.
  biq.PopMinOrMax(should_pop_max);
  if (!biq.empty())
    std::cout << "Min: " << biq.min() << " \tmax: " << biq.max() << std::endl;
  return 0;
}

 *
 *  Created on: Dec 2, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIQUEUE_H_
#define BIQUEUE_H_


#include <algorithm>
#include <cassert>
#include <vector>

#include "biheapify.h"
#include "almost_biheapify.h"
#include "biheap_sift.h"

template<class ValueType, typename SizeType = std::size_t>
class BiQueue {
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
    //Unnecessary code, useful for checking correctness.
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
    auto index_lambda = get_index_lambda();
    /* //Unnecessary code, useful for checking correctness.
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
  //Does not call BiHeapify or almost BiHeapify.
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
};

#endif /* BIQUEUE_H_ */
