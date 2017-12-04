/*
 * biqueue_demonstration.h
 *
 *  A BiQueue is a Double Ended Queue based around the idea of
 *   the BiHeap and almost BiHeap data structures.
 *  This file contains a simple BiQueue implementation specifically
 *   designed to help users understand how the push and pop
 *   algorithm works. It was not meant to be used in any application
 *   so it was not designed with efficiency in mind.
 *
 *  It uses almost BiHeaps to implement a double ended priority
 *   queue. It has amortized O(log N) insertions and
 *   amortized O(log N) deletions.
 *  Given a collection of elements, this double ended priority
 *   queue can be formed using O(N) swaps.
 *
 *  Created on: Dec 2, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIQUEUE_DEMONSTRATION_H_
#define BIQUEUE_DEMONSTRATION_H_

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#include "biheapify.h"
#include "almost_biheapify.h"
#include "biheap_sift.h"

#include "biheap_ostream.h"

template<class ValueType, typename SizeType = std::size_t>
class BiQueueDemonstrationClass {
public:
  bool verbose_ = true;
  std::vector<ValueType> vec_;
  SizeType num_elements_; //The number of elements currently in the almost BiHeap.
  //Note that num_elements_ always satisfies (2 * N_) / 3 - 4 <= num_elements_ <= N_.
  SizeType N_;            //The size of the BiHeap that induced this almost BiHeap.
                          //This is the maximum number of elements that the almost BiHeap
                          // can hold before it needs to be resized to allow for
                          // the insertion of another element.
                          //N should always be even and non-zero.
  //SizeType N_minus_1;   //Frequently used value.
  //SizeType half_of_N;   //Frequently used value.
  SizeType F_first_hc_, F_last_hc_;


  template<class Iterator>
  BiQueueDemonstrationClass(Iterator start, Iterator one_past_end, bool verbose = true) : verbose_(verbose) {
    Initialize<Iterator>(start, one_past_end);
    return ;
  }


  template<class Iterator>
  void Initialize(Iterator start, Iterator one_past_end) {
    num_elements_ = std::distance(start, one_past_end);
    std::vector<ValueType> vec;
    if (num_elements_ <= 2) { //Take care of this special case.
      if (num_elements_ == 2)
        vec.resize(4);
      else
        vec.resize(2);
      N_ = vec.size();  //BiHeap should always have even non-zero size.
      if (num_elements_ == 1) {
        vec[0] = *start;
      } else if (num_elements_ == 2) {
        F_first_hc_ = 1;
        F_last_hc_  = 2;
        vec[0] = *start;
        auto value_0 = *start;
        start++;
        auto value_1 = *start;
        if (value_0 <= value_1) {
          vec[0] = value_0;
          vec[N_ - 1] = value_1;
        }
        else {
          vec[N_ - 1] = value_0;
          vec[0] = value_1;
        }
      }
      vec_ = std::move(vec);
      return ;
    }
    SizeType num_elements_evened_down = num_elements_ - (num_elements_ % 2);
    N_                                = ParentHeapSize(num_elements_evened_down);
    vec.resize(N_);
    SizeType half_of_num_elements     = (num_elements_ + 1) / 2;
    F_first_hc_                       = half_of_num_elements;
    SizeType num_elements_remaining   = num_elements_ - half_of_num_elements;
    F_last_hc_                        = (N_ - 1) - (num_elements_ / 2);
    for (SizeType i = 0 ; start != one_past_end && i < F_first_hc_; i++, start++)
      vec[i] = *start;
    for ( ; start != one_past_end; num_elements_remaining--, start++)
      vec[N_ - num_elements_remaining] = *start;
    //Zero out all fused nodes to make their location clear.
    for (SizeType i = F_first_hc_; i <= F_last_hc_; i++)
      vec[i] = 0;
    //Form an almost BiHeap.
    AlmostBiHeapify<typename std::vector<ValueType>::iterator, SizeType>(vec.begin(), N_, F_first_hc_, F_last_hc_);
    vec_ = std::move(vec);
    return ;
  }

  template<class Iterator>
  void Initialize(Iterator start, SizeType new_N) {
    Initialize<Iterator>(start, start + new_N);
    return ;
  }

  //Assumes that even_N is even and positive.
  static SizeType ParentHeapSize(SizeType even_N) {
    SizeType half_N             = even_N / 2;
    SizeType three_times_half_N = 3 * half_N;
    if (half_N % 2 == 1)
      return three_times_half_N + 1;
    else
      return three_times_half_N;
  }

  //To see that insert() has amortized O(log(N_)) complexity, note
  // that after this structure is first created, if one were to
  // continue calling insert(), then it would AlmostBiHeapifySift()
  // (an O(log N) operation) approximately (i.e. plus or minus 2) N_ / 3
  // times before it would have to call AlmostBiHeapify(), where since
  // AlmostBiHeapify() is an O(N) operation, there is some constant C
  // such that AlmostBiHeapify() performs no more than C * N_ operations.
  //If each call to AlmostBiHeapifySift() performs at most D log N
  // operations, then at most D*((N_ / 3) + 2)*log(N_) + C N_ operations
  // will have been performed. Dividing by N_, shows that the amortized
  // complexity is O(log(N_)).
  //Note that since num_elements_ always satisfies
  // (2 * N_) / 3 - 4 <= num_elements_ <= N_,
  // insert() also has amortized O(log(num_elements_)) complexity.
  void insert(ValueType value) {
    if (num_elements_ <= 1) {
      if (num_elements_ == 0) {
        vec_[0] = value;
        N_ = 2;
        F_first_hc_ = 1;
        F_last_hc_ = 1;
        num_elements_++;
      } else if (num_elements_ == 1) {
        if (vec_.size() < 4)
          vec_.resize(4);
        N_ = 4;
        if (value > vec_[0])
          vec_[3] = value;
        else {
          vec_[3] = vec_[0];
          vec_[0] = value;
        }
        F_first_hc_ = 1;
        F_last_hc_  = 2;
        num_elements_++;
      }
      return ;
    }
    if (num_elements_ == N_) {//If we need to make more room.
      Initialize<typename std::vector<ValueType>::iterator>(vec_.begin(), vec_.begin() + num_elements_);
      assert(N_ == 2 || N_ % 3 != 2);
      if (verbose_) {
        std::cout << "\n******************ENLARGED BIHEAP******************** The New Almost BiHeap is:\n";
        PrintAlmostBiHeap();
        PrintVariableInformation();
      }
    }
    SizeType placement_index;
    if (F_first_hc_ <= (N_ - 1) - F_last_hc_) { //If there are as many nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
      placement_index = F_first_hc_;
      F_first_hc_++;
    } else {
      placement_index = F_last_hc_;
      F_last_hc_--;
    }
    vec_[placement_index] = value;
    if (verbose_) {
      std::cout << "Placed the value " << value << " at node " << placement_index
                << ", which has just been inserted into the this data structure.\n";
      PrintAlmostBiHeap();
      std::cout << "We still need to sift this element into its correct location using AlmostBiHeapifySift()." << std::endl;
    }
    num_elements_++;
    AlmostBiHeapifySift<typename std::vector<ValueType>::iterator, SizeType>(vec_.begin(), N_, placement_index, F_first_hc_, F_last_hc_);
    return ;
  }

  //Assumes that num_elements_ > 0.
  //pop_index should be either 0 or N_ - 1.
  //The argument that this function has amortized O(log(N_)) complexity
  // is analogous to the argument used to show that insert() also
  // has amortized O(log(N_)) complexity.
  void PopMinOrMax(SizeType pop_index) {
    if (num_elements_ == 1) {
      num_elements_ = 0;
      F_first_hc_   = 0;
      F_last_hc_    = 1;
      N_            = 2;
      return ;
    } else if (num_elements_ == 2) {
      num_elements_ = 1;
      if (pop_index == 0)
        std::iter_swap(vec_.begin(), vec_.begin() + (N_ - 1));
      F_first_hc_   = 0;
      F_last_hc_    = 0;
      N_            = 2;
      return ;
    }
    SizeType heap_size = HeapSize<SizeType>(N_);
    SizeType first_in_node = N_ - heap_size;
    SizeType swap_index;
    if (F_first_hc_ == first_in_node) { //If we can not remove any more In nodes.
      //Shift the pure max heap to the left so as to be contiguous with
      // the pure min heap in memory.
      SizeType i = first_in_node;
      SizeType j = heap_size; //one past last In node.
      while (j < N_) {
        vec_[i] = vec_[j];
        i++;
        j++;
      }
      N_ = i;
      //BiHeapify it.
      BiHeapify<typename std::vector<ValueType>::iterator, SizeType>(vec_.begin(), N_);
      F_first_hc_ = ((N_ + 1) / 2);
      F_last_hc_  = F_first_hc_;
      swap_index  = F_first_hc_;
      pop_index   = N_ - 1;
      if (verbose_) {
        std::cout << "\n******************SHRUNK ALMOST BIHEAP********************"
                  << " The New Almost BiHeap is:\n";
        PrintAlmostBiHeap();
        PrintVariableInformation();
      }
    } else {
      if (F_first_hc_ > (N_ - 1) - F_last_hc_) { //If there are more nodes to the left of F_first_hc_ as there are to the right of F_last_hc_.
        F_first_hc_--;
        swap_index = F_first_hc_;
      } else {
        F_last_hc_++;
        swap_index = F_last_hc_;
      }
    }
    num_elements_--;
    if (verbose_) {
      std::cout << "Swapping values with min heap coordinates "
                << swap_index << " and " << pop_index
                << " and then applying AlmostBiHeapifySift() to sift nodes "
                << pop_index << ".\n"
                << "The node that has been removed from the Almost BiHeap is now indicated by a 0."
                << std::endl;
    }
    std::iter_swap(vec_.begin() + pop_index, vec_.begin() + swap_index);
    //Indicate that this node is no longer in the almost BiHeap by setting it to 0.
    vec_[swap_index] = static_cast<ValueType>(0);
    //Sift the element into place.
    AlmostBiHeapifySift<typename std::vector<ValueType>::iterator, SizeType>(vec_.begin(), N_, pop_index, F_first_hc_, F_last_hc_);
    return ;
  }

  //Assumes that num_elements_ > 0.
  void PopMin() {
    PopMinOrMax(0);
    return ;
  }

  //Assumes that num_elements_ > 0.
  void PopMax() {
    PopMinOrMax(N_ - 1);
    return ;
  }

  ValueType Min() const {
    return vec_[0];
  }

  ValueType Max() const {
    if (num_elements_ == 1)
      return vec_[0];
    return vec_[N_ - 1];
  }

  SizeType size() {
    return num_elements_;
  }

  void PrintAlmostBiHeap() {
    PrintBiHeap(vec_.begin(), N_, false);
  }

  void PrintVariableInformation() {
    SizeType heap_size     = HeapSize<SizeType>(N_);
    SizeType first_in_node = N_ - heap_size;
    std::cout << "vec_.size() = "         << vec_.size()
              << " \tN_ = "               << N_
              << " \tnum_elements_ = "    << num_elements_
              << " \theap_size = "        << heap_size
              << " \tfirst_in_node = "    << first_in_node
              << " \tF_first_hc_ = "      << F_first_hc_
              << " \tF_last_hc_ = "       << F_last_hc_
              << " \tFlip(F_last_hc_) = " << ((N_ - 1) - F_last_hc_)
              << std::endl;
  }
};

template<class T>
T GetRandomPositiveInt(T last_value = 1) {
  std::random_device rnd_device;
  std::mt19937 generator(rnd_device());
  std::uniform_int_distribution<T> dist(1, 99);
  return dist(generator);
}

void DEBQDemonstration(bool verbose = true, int num_random_iterations = 64) {
  std::vector<int> vec({5, 4, 3, 2, 1});
  BiQueueDemonstrationClass<int> biq(vec.begin(), vec.end(), verbose);
  int next_value = 6;
  if (verbose) {
    std::cout << "NOTE: The top half is the pure min heap and the bottom half is the pure max heap.\n";
    std::cout << "The node at the very top is node 0 (the minimum) and the node at the very bottom is "
              << "the maximum.\n";
    std::cout << "These are just the values of the nodes of the almost BiHeap graphs.\n"
              << "To view the edges of the BiHeap graph, use the functions found in biheap_tikz_graph.h \n"
              << " to generate LaTeX code that you can then compile to produce an actual image.\n";
    std::cout << "The value of each element actually in the doubled ended BiHeap priority queue is POSITIVE.\n"
              << "The 0's represent unused nodes that are not part of the current almost BiHeap.\n";
    std::cout << "The coordinate of each node is represented by its min heap coordinate (which is\n"
              << " the same thing as its location/index in the vector containing this data)\n";
    std::cout << "The initial values of the vector (not a BiHeap) are:" << '\n';
    PrintBiHeap(vec.begin(), vec.size(), false);
    std::cout << "-------------------------------------------------------------------------------------" << std::endl;
    std::cout << "After construction of DEBQ object:" << '\n';
    biq.PrintAlmostBiHeap();
    std::cout << "\n\n_____________________________________________________________________________________\n\n\n";
  }
  for (int i = 0; i < num_random_iterations; i++) {
    if (biq.size() == 0 || GetRandomPositiveInt(next_value) % 2 == 0) {
      auto value = (GetRandomPositiveInt(next_value++) % 99) + 1;
      if (verbose) {
        std::cout << "->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->" << std::endl;
        std::cout << "Values before call: ";
        biq.PrintVariableInformation();
        std::cout << "Calling insert(" << value << "). ";
      }
      biq.insert(value);
      if (verbose) {
        std::cout << "After call to insert(" << value << "):\n";
        biq.PrintAlmostBiHeap();
        std::cout << "Values after call: ";
        biq.PrintVariableInformation();
        std::cout << "->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->" << std::endl;
      }
    }
    if (biq.size() > 0 && GetRandomPositiveInt(next_value) % 5 == 0) {
      if (verbose) {
        std::cout << "______________________________________________________________________________________"<< std::endl;
        std::cout << "Values before call: ";
        biq.PrintVariableInformation();
        std::cout << "Calling PopMin() where the Min Value was " << biq.Min() << ". \n";
      }
      biq.PopMin();
      if (verbose) {
        std::cout << "After call to PopMin(): \n";
        biq.PrintAlmostBiHeap();
        std::cout << "Values after call: ";
        biq.PrintVariableInformation();
        std::cout << "______________________________________________________________________________________"<< std::endl;
      }
    }
    if (biq.size() > 0 && GetRandomPositiveInt(next_value) % 5 == 0) {
      if (verbose) {
        std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
        std::cout << "Values before call: ";
        biq.PrintVariableInformation();
        std::cout << "Calling PopMax() where the Max Value was " << biq.Max() << ". \n";
      }
      biq.PopMax();
      if (verbose) {
        std::cout << "After call to PopMax(): \n";
        biq.PrintAlmostBiHeap();
        std::cout << "Values after call: ";
        biq.PrintVariableInformation();
        std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
      }
    }
  }
  if (verbose)
    std::cout << "Done" << std::endl;
  return ;
}


#endif /* BIQUEUE_DEMONSTRATION_H_ */
