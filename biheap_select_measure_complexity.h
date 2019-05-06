/*
 * biheap_select_measure_complexity.h
 *
 *  Created on: May 2, 2019
 *      Author: Matthew Gregory Krupa
 *      Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAP_SELECT_MEASURE_COMPLEXITY_H_
#define BIHEAP_SELECT_MEASURE_COMPLEXITY_H_

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <string>
#include <vector>

//#include "biheap_common.h"
#include "biheap_ostream.h"
#include "biheapify.h"
#include "biheap_select.h"

std::size_t num_swaps_this_test = 0;

//Emplaces the element V + desired index.
// That is, it rearranges V, V + 1, ..., V+ (N - 1) so that
// (1) *(V + i) <= *(V + desired_index) for all 0 <= i < desired_index, and
// (2) *(V + i) >= *(V + desired_index) for all desired_index < i < N.
template<class RAI, typename size_type = std::size_t>
void BiHeapSelectMeasureComplexity(RAI V, size_type N, size_type desired_index) {
  assert(desired_index < N || !(std::cout << "desired_index = " << desired_index << " \tN = " << N << std::endl));
  while(N > 0) {
    if (N <= BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_THIS) {
      //Note that since N <= BIHEAP_SELECT_SORT_IF_LENGTH_IS_LESS_THAN_THIS, a constant, this sort operation is effectively O(constant).
      num_swaps_this_test += N * N; //Use N * N to approximate this constant since the sort algorithm is O(N log N).
      std::sort<RAI>(V, V + N);
      break ;
    }
    assert(desired_index >= 0);
    assert(desired_index < N || !(std::cout << "desired_index = " << desired_index << " \tN = " << N << std::endl));

    num_swaps_this_test += (7.0l / 3.0l) * N; //It suffices to add N here since it is known that BiHeapify (called next)
                              // is O(N) and that a call to BiHeapify(V, N) performs at most 7N/3 swaps.
    BiHeapify<RAI, size_type>(V, N);

    std::pair<size_type, size_type> next_range_pair = BiHeapSelectGetRangeToFindMedianOf<size_type>(N, desired_index);
    size_type next_range_start = next_range_pair.first;
    size_type next_range_end = next_range_pair.second;
    assert(next_range_end >= next_range_start);
    size_type next_range_length = (next_range_end + 1) - next_range_start;
    size_type middle_of_next_range = next_range_start + (next_range_length / 2);
    //std::cout << "next_range_start = " << next_range_start << " \tnext_range_end = " << next_range_end << " \tnext_range_length = " << next_range_length << " \tmiddle_of_next_range = " << middle_of_next_range << std::endl;
    assert(next_range_start <= middle_of_next_range);
    assert(middle_of_next_range <= next_range_end);

    //For the following call, we do NOT add middle_of_next_range - next_range_start to num_swaps_this_test since the addition
    // to num_swaps_this_test of the appropriate quantity will be taken care of within the subsequence call.
    //Find the median of the range.
    BiHeapSelectMeasureComplexity<RAI, size_type>(V + next_range_start, next_range_length, middle_of_next_range - next_range_start);

    num_swaps_this_test += N; //We add N here again since DutchPartition()is an O(N) algorithm where
                              // the call DutchPartition(V, N, middle_of_next_range) performs at most N swaps.

    //Dutch flag partition all N elements about the median of the range that was found by the above call to BiHeapSelect().
    std::pair<size_type, size_type> pivot_range = DutchPartition<RAI, size_type>(V, N, middle_of_next_range);

    auto a = pivot_range.first;
    auto b = pivot_range.second;
    // a is the index that is 1 + (the index of last element that is < the pivot value), and
    // b is the index of the first element that is > the pivot value.
    if (desired_index < a) {
      N = a;
      assert(desired_index < N);
    } else if (desired_index >= b) {
      V = V + b;
      N -= b;
      desired_index -= b;
      assert(desired_index < N);
    } else { //a <= diesred_index < b so we've emplaced the desired index.
      break ;
    }
  }
  return ;
}











//***********************************************
// MEASURING COMPLEXITY
//***********************************************


template<class ValueType, typename size_type = std::size_t>
struct MeasureComplexityOfBiHeapifySelectMediansOnGivenSize {
  size_type N_;
  std::vector<ValueType> vec_;//, vec_original_;
  static long double largest_ratio_of_num_swaps_over_N_so_far_;
  static size_type N_of_largest_ratio_of_num_swaps_over_N_so_far_;
  bool verbose_ = true;
  ValueType lower_, upper_;
  size_type left_index_;
  size_type right_index_;

  size_type total_num_test_runs_ = 0;
  long double largest_ratio_of_num_swaps_over_N_for_this_N_ = 0.0;

  MeasureComplexityOfBiHeapifySelectMediansOnGivenSize(size_type N) : N_(N) {
    vec_ = std::vector<ValueType>(N_);
    left_index_ = N_ / 2;
    right_index_ = left_index_;
    if (N_ % 2 == 0)
      left_index_--;
    lower_ = std::numeric_limits<ValueType>::min();
    upper_ = std::numeric_limits<ValueType>::max();
    //vec_original_ = vec_;
  }

  template<class Iterator, typename T>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, T a, T b) {
    std::random_device rnd_device;
    std::mt19937 generator(rnd_device());
    std::uniform_int_distribution<T> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, float a,
                             float b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<float> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, double a,
                             double b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<double> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, long double a,
                             long double b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<long double>  dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  void FillVecWithRandomNumbers() {
    FillWithRandomNumbers(vec_.begin(), vec_.end(), lower_, upper_);
    //vec_original_ = vec_;
  }

  //Note that is FillVecWithRandomNumbers() is used then the vast majority of the time
  // all of vec_'s values will be distinct from each other.
  //So to investigate what happens when there are frequently repeated values,
  // the following function may be used in place of FillVecWithRandomNumbers().
  void FillVecWithRandomNumbersSizeScaledInterval(long double scale = 10.0l) {
    if (scale <= 0.0l)
      scale = 1.0l;
    ValueType max;
    if (scale == 1.0l)
      max = static_cast<ValueType>(N_); //Done so as to minimize the conversion error that is introduced.
    else if (scale < 1.0l || (scale > 1.0l && static_cast<long double>(std::numeric_limits<ValueType>::max() / scale) < static_cast<long double>(N_)))
      max = static_cast<ValueType>(scale * N_);
    else
      max = std::numeric_limits<ValueType>::max();
    FillWithRandomNumbers(vec_.begin(), vec_.end(), static_cast<ValueType>(0), max);
    //vec_original_ = vec_;
  }

  size_type GetNumberOfElementsGreaterThanOrEqualToValue(const ValueType value) {
    size_type counter = 0;
    auto first = vec_.begin();
    for (size_type i = 0; i < N_; i++) {
      if (*(first + i) >= value)
        counter++;
    }
    return counter;
  }

  size_type GetNumberOfElementsLessThanOrEqualToValue(const ValueType value) {
    size_type counter = 0;
    auto first = vec_.begin();
    for (size_type i = 0; i < N_; i++) {
      if (*(first + i) <= value)
        counter++;
    }
    return counter;
  }

  void Analyze() {
    long double ratio_of_num_swaps_over_N = static_cast<long double>(num_swaps_this_test) / static_cast<long double>(N_);

    if (ratio_of_num_swaps_over_N > largest_ratio_of_num_swaps_over_N_for_this_N_){
      largest_ratio_of_num_swaps_over_N_for_this_N_ = ratio_of_num_swaps_over_N;
    }
    if (ratio_of_num_swaps_over_N > largest_ratio_of_num_swaps_over_N_so_far_) {
      largest_ratio_of_num_swaps_over_N_so_far_ = ratio_of_num_swaps_over_N;
      N_of_largest_ratio_of_num_swaps_over_N_so_far_ = N_;
      std::cout << "New largest_ratio_of_num_swaps_over_N_so_far_ = " << largest_ratio_of_num_swaps_over_N_so_far_ << " \t for N = " << N_of_largest_ratio_of_num_swaps_over_N_so_far_ << std::endl;
    }
  }

  void RunBiHeapSelect() {
    BiHeapSelectMeasureComplexity(vec_.begin(), N_, N_ / 2);
    total_num_test_runs_++;
    return ;
  }

  void PerformTest(long double scale = 0.0l) {
    if (scale <= 0.0l)
      FillVecWithRandomNumbers();
    else
      FillVecWithRandomNumbersSizeScaledInterval(scale);
    num_swaps_this_test = 0;
    RunBiHeapSelect();
    Analyze();
    return ;
  }

  std::string GetDescriptionOfResult() {
    std::stringstream sstrm;
    sstrm << "N = " << N_ << " \t(N+1)/2 = " << ((N_ + 1) / 2) << " \tNumber of tests = " << total_num_test_runs_;
    sstrm << " \tmax {(num swaps/N} = " << largest_ratio_of_num_swaps_over_N_for_this_N_;
    sstrm << " \tMax Ever {(num swaps/N} = " << largest_ratio_of_num_swaps_over_N_so_far_ << " \t occurred for N = " << N_of_largest_ratio_of_num_swaps_over_N_so_far_;

    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double MeasureComplexityOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::largest_ratio_of_num_swaps_over_N_so_far_ = 0.0l;
template<class ValueType, typename size_type>
size_type MeasureComplexityOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::N_of_largest_ratio_of_num_swaps_over_N_so_far_ = 1;

template<class ValueType = int, typename size_type = std::size_t>
void MeasureBiHeapSelectMediansComplexity(size_type num_test = (1u << 17), size_type N_start = 20, size_type N_one_past_end = (1u << 17)) {
  if (N_start >= N_one_past_end) {
    N_one_past_end = N_start + 1;
  }
  for (size_type N = N_start; N < N_one_past_end; N += N / 3) {
    MeasureComplexityOfBiHeapifySelectMediansOnGivenSize<ValueType> test(N);
    test.verbose_ = false;
    for (std::size_t i = 0; i < num_test; i++) {
      test.PerformTest();
    }
    //The two key quantities to look at are:
    //smallest min(# elements <= pivot, # elements >= pivot), and
    //average of min(# elements <= pivot, # elements >= pivot)
    std::cout << test.GetDescriptionOfResult() << std::endl;
  }
}





template<class ValueType, typename size_type = std::size_t>
struct MeasureComplexityOfBiHeapifySelectOnGivenSize {
  size_type N_, pivot_position_;
  std::vector<ValueType> vec_;//, vec_original_;
  std::vector<long double> permitted_pivot_positions_;
  long double pivot_average_ = 0.0l;
  static long double largest_ratio_of_num_swaps_over_N_so_far_;
  static size_type N_of_largest_ratio_of_num_swaps_over_N_so_far_;
  bool verbose_ = true;
  ValueType lower_, upper_;
  size_type left_index_;
  size_type right_index_;

  size_type total_num_test_runs_ = 0;
  long double largest_ratio_of_num_swaps_over_N_for_this_N_ = 0.0;

  MeasureComplexityOfBiHeapifySelectOnGivenSize(size_type N, std::vector<long double> permitted_pivot_positions = std::vector<long double>()) : N_(N),
      permitted_pivot_positions_(permitted_pivot_positions) {
    vec_ = std::vector<ValueType>(N_);
    left_index_ = N_ / 2;
    right_index_ = left_index_;
    if (N_ % 2 == 0)
      left_index_--;
    lower_ = std::numeric_limits<ValueType>::min();
    upper_ = std::numeric_limits<ValueType>::max();
    //vec_original_ = vec_;
  }

  template<class Iterator, typename T>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, T a, T b) {
    std::random_device rnd_device;
    std::mt19937 generator(rnd_device());
    std::uniform_int_distribution<T> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, float a,
                             float b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<float> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, double a,
                             double b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<double> dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  template<class Iterator>
  void FillWithRandomNumbers(Iterator start, Iterator one_past_end, long double a,
                             long double b) {
    std::random_device rnd_device;
    std::mt19937_64 generator(rnd_device());
    std::uniform_real_distribution<long double>  dist(a, b);
    for (auto it = start; it != one_past_end; it++)
      *it = dist(generator);
    return ;
  }

  void FillVecWithRandomNumbers() {
    FillWithRandomNumbers(vec_.begin(), vec_.end(), lower_, upper_);
    //vec_original_ = vec_;
  }

  //Note that is FillVecWithRandomNumbers() is used then the vast majority of the time
  // all of vec_'s values will be distinct from each other.
  //So to investigate what happens when there are frequently repeated values,
  // the following function may be used in place of FillVecWithRandomNumbers().
  void FillVecWithRandomNumbersSizeScaledInterval(long double scale = 10.0l) {
    if (scale <= 0.0l)
      scale = 1.0l;
    ValueType max;
    if (scale == 1.0l)
      max = static_cast<ValueType>(N_); //Done so as to minimize the conversion error that is introduced.
    else if (scale < 1.0l || (scale > 1.0l && static_cast<long double>(std::numeric_limits<ValueType>::max() / scale) < static_cast<long double>(N_)))
      max = static_cast<ValueType>(scale * N_);
    else
      max = std::numeric_limits<ValueType>::max();
    FillWithRandomNumbers(vec_.begin(), vec_.end(), static_cast<ValueType>(0), max);
    //vec_original_ = vec_;
  }

  size_type GetNumberOfElementsGreaterThanOrEqualToValue(const ValueType value) {
    size_type counter = 0;
    auto first = vec_.begin();
    for (size_type i = 0; i < N_; i++) {
      if (*(first + i) >= value)
        counter++;
    }
    return counter;
  }

  size_type GetNumberOfElementsLessThanOrEqualToValue(const ValueType value) {
    size_type counter = 0;
    auto first = vec_.begin();
    for (size_type i = 0; i < N_; i++) {
      if (*(first + i) <= value)
        counter++;
    }
    return counter;
  }

  void Analyze() {
    long double ratio_of_num_swaps_over_N = static_cast<long double>(num_swaps_this_test) / static_cast<long double>(N_);

    if (ratio_of_num_swaps_over_N > largest_ratio_of_num_swaps_over_N_for_this_N_){
      largest_ratio_of_num_swaps_over_N_for_this_N_ = ratio_of_num_swaps_over_N;
    }
    if (ratio_of_num_swaps_over_N > largest_ratio_of_num_swaps_over_N_so_far_) {
      largest_ratio_of_num_swaps_over_N_so_far_ = ratio_of_num_swaps_over_N;
      N_of_largest_ratio_of_num_swaps_over_N_so_far_ = N_;
      std::cout << "New largest_ratio_of_num_swaps_over_N_so_far_ = " << largest_ratio_of_num_swaps_over_N_so_far_ << " \t for N = " << N_of_largest_ratio_of_num_swaps_over_N_so_far_ << std::endl;
    }
  }

  void RunBiHeapSelect() {
    BiHeapSelectMeasureComplexity(vec_.begin(), N_, N_ / 2);
    total_num_test_runs_++;
    return ;
  }

  void PerformTest(long double scale = 0.0l) {
    if (scale <= 0.0l)
      FillVecWithRandomNumbers();
    else
      FillVecWithRandomNumbersSizeScaledInterval(scale);
    //Select the pivot position
    if (permitted_pivot_positions_.size() == 0) { //If the pivot position is to be selected at random
      if (vec_.size() >= 1) {
        std::random_device rnd_device;
        std::mt19937 generator(rnd_device());
        std::uniform_int_distribution<size_type> dist(0, vec_.size() - 1);
        pivot_position_ = dist(generator);
      } else {
        pivot_position_ = 0;
        assert(N_ == 0);
      }
    } else { //else the pivot position is to be selected from the list permitted_pivot_positions_,
             // where each element of this vector should be between 0.0 and 1.0
      std::random_device rnd_device;
      std::mt19937 generator(rnd_device());
      std::uniform_int_distribution<size_type> dist(0, permitted_pivot_positions_.size() - 1);
      auto pivot_position_fraction = permitted_pivot_positions_[dist(generator)];
      assert(pivot_position_fraction >= 0.0l && pivot_position_fraction <= 1.0l);
      pivot_position_ = pivot_position_fraction * N_;
      if (pivot_position_ >= N_) {
        pivot_position_ = N_ - 1;
      } else if (pivot_position_ < 0) {
        pivot_position_ = 0;
      }
    }
    assert(pivot_position_ >= 0 && pivot_position_ < N_);
    //Update the average values of all pivots
    pivot_average_ = (pivot_average_ * static_cast<long double>(total_num_test_runs_) + pivot_position_) / static_cast<long double>(total_num_test_runs_ + 1);
    num_swaps_this_test = 0;
    RunBiHeapSelect();
    Analyze();
    return ;
  }

  std::string GetDescriptionOfResult() {
    std::stringstream sstrm;
    sstrm << "N = " << N_
          << " \t(N-1)/2 = " << (static_cast<long double>(N_ - 1) / 2.0l)
          << " \tpivot_average_ = " << pivot_average_
          << " \tNum tests = " << total_num_test_runs_;
    sstrm << " \tMax {(num swaps/N} = " << largest_ratio_of_num_swaps_over_N_for_this_N_;
    sstrm << " \tMax Ever {(num swaps/N} = " << largest_ratio_of_num_swaps_over_N_so_far_
          << " occurred for N = " << N_of_largest_ratio_of_num_swaps_over_N_so_far_;

    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double MeasureComplexityOfBiHeapifySelectOnGivenSize<ValueType, size_type>::largest_ratio_of_num_swaps_over_N_so_far_ = 0.0l;
template<class ValueType, typename size_type>
size_type MeasureComplexityOfBiHeapifySelectOnGivenSize<ValueType, size_type>::N_of_largest_ratio_of_num_swaps_over_N_so_far_ = 1;

template<class ValueType = int, typename size_type = std::size_t>
void MeasureBiHeapSelectComplexity(size_type num_tests = (1u << 17), size_type N_start = 20, size_type N_one_past_end = (1u << 17)) {
  if (N_start >= N_one_past_end) {
    N_one_past_end = N_start + 1;
  }
  for (size_type N = N_start; N < N_one_past_end; N += N / 3) {
    MeasureComplexityOfBiHeapifySelectOnGivenSize<ValueType> test(N);
    test.verbose_ = false;
    for (std::size_t i = 0; i < num_tests; i++) {
      test.PerformTest();
    }
    std::cout << test.GetDescriptionOfResult() << std::endl;
  }
}

#endif /* BIHEAP_SELECT_MEASURE_COMPLEXITY_H_ */
