/*
 * biheapify_inwards_pivot_testing.h
 *
 *  Created on: Nov 3, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAPIFY_INWARDS_PIVOT_TESTING_H_
#define BIHEAPIFY_INWARDS_PIVOT_TESTING_H_

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <sstream>
#include <vector>

#include "biheapify_inwards.h"

template<class ValueType, typename size_type = std::size_t>
struct TestBiHeapifyInwardsOnGivenSize {
  size_type N_;
  std::vector<ValueType> vec_;//, vec_original_;
  bool verbose_ = true;
  ValueType lower_, upper_;
  size_type left_index_;
  size_type right_index_;

  size_type total_num_elements_less_than_or_equal_to_pivot_ = 0;
  size_type total_num_elements_greater_than_or_equal_to_pivot_ = 0;
  size_type total_min_num_elements_ge_and_le_pivot_ = 0;
  size_type total_num_test_runs_ = 0;
  size_type smallest_total_num_elements_less_than_or_equal_to_pivot_ = 0;
  size_type smallest_total_num_elements_greater_than_or_equal_to_pivot_ = 0;
  static long double smallest_ever_num_elements_less_than_or_equal_ratio_;
  static long double smallest_ever_num_elements_greater_than_or_equal_ratio_;
  static long double smallest_ever_of_smaller_of_num_le_and_ge_ratios_;
  static size_type N_for_smallest_ever_of_smaller_of_num_le_and_ge_;
  static size_type smaller_of_num_le_and_ge_for_smallest_ever_;
  long double smallest_ever_num_elements_less_than_or_equal_ratio_this_size_ = 1.0l;
  long double smallest_num_elements_greater_than_or_equal_ratio_this_size_ = 1.0l;
  long double smallest_of_smaller_of_num_le_and_ge_ratios_this_size_ = 1.0;

  TestBiHeapifyInwardsOnGivenSize(size_type N) : N_(N) {
    vec_ = std::vector<ValueType>(N_);
    //vec_original_ = vec_;
    left_index_ = N_ / 2;
    right_index_ = left_index_;
    smallest_total_num_elements_less_than_or_equal_to_pivot_   = N_;
    smallest_total_num_elements_greater_than_or_equal_to_pivot_ = N_;
    if (N_ % 2 == 0)
      left_index_--;
    lower_ = std::numeric_limits<ValueType>::min();
    upper_ = std::numeric_limits<ValueType>::max();
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
    auto first = vec_.begin();
    auto left_middle_value = *(first + left_index_);
    auto right_middle_value = *(first + right_index_);

    size_type num_elements_less_than_or_equal    = GetNumberOfElementsLessThanOrEqualToValue(left_middle_value);
    size_type num_elements_greater_than_or_equal = GetNumberOfElementsGreaterThanOrEqualToValue(right_middle_value);
    size_type min_num_elements_ge_and_le_pivot   = num_elements_less_than_or_equal < num_elements_greater_than_or_equal ? num_elements_less_than_or_equal : num_elements_greater_than_or_equal;
    total_num_elements_less_than_or_equal_to_pivot_    += num_elements_less_than_or_equal;
    total_num_elements_greater_than_or_equal_to_pivot_ += num_elements_greater_than_or_equal;
    total_min_num_elements_ge_and_le_pivot_            += min_num_elements_ge_and_le_pivot;

    long double num_elements_less_than_or_equal_ratio    = static_cast<long double>(num_elements_less_than_or_equal) / static_cast<long double>(N_);
    long double num_elements_greater_than_or_equal_ratio = static_cast<long double>(num_elements_greater_than_or_equal) / static_cast<long double>(N_);
    long double smaller_of_num_le_and_ge_ratios = num_elements_less_than_or_equal_ratio < num_elements_greater_than_or_equal_ratio ? num_elements_less_than_or_equal_ratio : num_elements_greater_than_or_equal_ratio;

    if (num_elements_greater_than_or_equal < smallest_total_num_elements_greater_than_or_equal_to_pivot_)
      smallest_total_num_elements_greater_than_or_equal_to_pivot_ = num_elements_greater_than_or_equal;
    if (num_elements_less_than_or_equal < smallest_total_num_elements_less_than_or_equal_to_pivot_)
      smallest_total_num_elements_less_than_or_equal_to_pivot_ = num_elements_less_than_or_equal;

    if (num_elements_less_than_or_equal_ratio < smallest_ever_num_elements_less_than_or_equal_ratio_)
      smallest_ever_num_elements_less_than_or_equal_ratio_ = num_elements_less_than_or_equal_ratio;
    if (num_elements_greater_than_or_equal_ratio < smallest_ever_num_elements_greater_than_or_equal_ratio_)
      smallest_ever_num_elements_greater_than_or_equal_ratio_ = num_elements_greater_than_or_equal_ratio;

    if (smaller_of_num_le_and_ge_ratios < smallest_ever_of_smaller_of_num_le_and_ge_ratios_) {
      if (verbose_)
        std::cout << "smallest min(# elements <= pivot, # elements >= pivot) ever seen = "
                  << std::setprecision(30) << smaller_of_num_le_and_ge_ratios << std::endl;
      smallest_ever_of_smaller_of_num_le_and_ge_ratios_ = smaller_of_num_le_and_ge_ratios;
      N_for_smallest_ever_of_smaller_of_num_le_and_ge_ = N_;
      smaller_of_num_le_and_ge_for_smallest_ever_ = min_num_elements_ge_and_le_pivot;
    }

    if (num_elements_less_than_or_equal_ratio < smallest_ever_num_elements_less_than_or_equal_ratio_this_size_)
      smallest_ever_num_elements_less_than_or_equal_ratio_this_size_ = num_elements_less_than_or_equal_ratio;
    if (num_elements_greater_than_or_equal_ratio < smallest_num_elements_greater_than_or_equal_ratio_this_size_)
      smallest_num_elements_greater_than_or_equal_ratio_this_size_ = num_elements_greater_than_or_equal_ratio;

    if (smaller_of_num_le_and_ge_ratios < smallest_of_smaller_of_num_le_and_ge_ratios_this_size_) {
      if (verbose_)
        std::cout << "smallest min(# elements <= pivot, # elements >= pivot) for this size = "
                  << std::setprecision(30) << smaller_of_num_le_and_ge_ratios << std::endl;
      smallest_of_smaller_of_num_le_and_ge_ratios_this_size_ = smaller_of_num_le_and_ge_ratios;
    }
  }

  void RunBiHeapifyInwards() {
    BiHeapifyInwards(vec_.begin(), N_);
    total_num_test_runs_++;
  }

  void RunBiHeapifyInwardsNicerMath() {
    BiHeapifyInwardsNicerMath(vec_.begin(), N_);
    total_num_test_runs_++;
  }

  void PerformTest(long double scale = 0.0l) {
    if (scale <= 0.0l)
      FillVecWithRandomNumbers();
    else
      FillVecWithRandomNumbersSizeScaledInterval(scale);
    RunBiHeapifyInwards();
    Analyze();
  }

  void PerformTestWithBiHeapifyInwardsNicerMath(long double scale = 0.0l) {
    if (scale <= 0.0l)
      FillVecWithRandomNumbers();
    else
      FillVecWithRandomNumbersSizeScaledInterval(scale);
    RunBiHeapifyInwardsNicerMath();
    Analyze();
  }

  //If total_num_runs_to_print !=0 then it will output that the number of tests
  // performed is total_num_runs_to_print rather than total_num_test_runs_.
  std::string GetStartOfDescriptionOfResult(size_type total_num_runs_to_print = 0) {
    std::stringstream sstrm;
    size_type total_num_runs = total_num_test_runs_;
    if (total_num_runs_to_print != 0)
      total_num_runs = total_num_runs_to_print;
    sstrm << "N = " << N_ << " \tNumber of tests performed = " << total_num_runs << '\n';
    sstrm << "For this size N = " << N_ << '\n';
    sstrm << "smallest min(# elements <= pivot, # elements >= pivot) is ";
    return sstrm.str();
  }

  std::string GetDescriptionOfResult(bool include_start_of_description = true,
                                     bool output_smallest_num_less_than_info = true,
                                     bool output_smallest_num_greater_than_info = true) {
    std::stringstream sstrm;
    if (include_start_of_description)
      sstrm << GetStartOfDescriptionOfResult();
    size_type min_of_num_le_and_num_ge = smallest_total_num_elements_less_than_or_equal_to_pivot_;
    if (smallest_total_num_elements_greater_than_or_equal_to_pivot_ < min_of_num_le_and_num_ge)
      min_of_num_le_and_num_ge = smallest_total_num_elements_greater_than_or_equal_to_pivot_;

    long double min_of_num_le_and_num_ge_ratio = static_cast<long double>(min_of_num_le_and_num_ge) / static_cast<long double>(N_);
    sstrm << min_of_num_le_and_num_ge;
    sstrm << " \t(" << (100.0l * min_of_num_le_and_num_ge_ratio) << "% of N)" << '\n';
    if (output_smallest_num_less_than_info) {
      long double smallest_total_num_elements_less_than_or_equal_to_pivot_ratio = static_cast<long double>(smallest_total_num_elements_less_than_or_equal_to_pivot_) / static_cast<long double>(N_);
      sstrm << "smallest # of elements <= pivot is                        " << smallest_total_num_elements_less_than_or_equal_to_pivot_;
      sstrm << " \t(" << (100.0l * smallest_total_num_elements_less_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    }
    if (output_smallest_num_greater_than_info) {
      long double smallest_total_num_elements_greater_than_or_equal_to_pivot_ratio = static_cast<long double>(smallest_total_num_elements_greater_than_or_equal_to_pivot_) / static_cast<long double>(N_);
      sstrm << "smallest # of elements >= pivot is                        " << smallest_total_num_elements_greater_than_or_equal_to_pivot_;
      sstrm << " \t(" << (100.0l * smallest_total_num_elements_greater_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    }

    long double ave_total_min_num_elements_ge_and_le_pivot_ = static_cast<long double>(total_min_num_elements_ge_and_le_pivot_) / static_cast<long double>(total_num_test_runs_);
    long double ave_total_min_num_elements_ge_and_le_pivot_ratio = ave_total_min_num_elements_ge_and_le_pivot_ / static_cast<long double>(N_);
    sstrm << "average of min(# elements <= pivot, # elements >= pivot) is " << ave_total_min_num_elements_ge_and_le_pivot_;
    sstrm << " \t(" << (100.0l * ave_total_min_num_elements_ge_and_le_pivot_ratio) << "% of N)" << '\n';
    if (output_smallest_num_less_than_info) {
      long double ave_num_elements_less_than_or_equal_to_pivot_ = static_cast<long double>(total_num_elements_less_than_or_equal_to_pivot_) / static_cast<long double>(total_num_test_runs_);
      long double ave_num_elements_less_than_or_equal_to_pivot_ratio = ave_num_elements_less_than_or_equal_to_pivot_ / static_cast<long double>(N_);
      sstrm << "average  # of elements <= pivot is                          " << ave_num_elements_less_than_or_equal_to_pivot_;
      sstrm << " \t(" << (100.0l * ave_num_elements_less_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    }
    if (output_smallest_num_greater_than_info) {
      long double ave_num_elements_greater_than_or_equal_to_pivot_ = static_cast<long double>(total_num_elements_greater_than_or_equal_to_pivot_) / static_cast<long double>(total_num_test_runs_);
      long double ave_num_elements_greater_than_or_equal_to_pivot_ratio = ave_num_elements_greater_than_or_equal_to_pivot_ / static_cast<long double>(N_);
      sstrm << "average  # of elements >= pivot is                          " << ave_num_elements_greater_than_or_equal_to_pivot_;
      sstrm << " \t(" << (100.0l * ave_num_elements_greater_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    }

    sstrm << "Over all sizes N tested so far:\n";
    sstrm << "smallest ever min(# elements <= pivot, # elements >= pivot)/N is " << (smallest_ever_of_smaller_of_num_le_and_ge_ratios_)
          << " = " << smaller_of_num_le_and_ge_for_smallest_ever_ << " / " << N_for_smallest_ever_of_smaller_of_num_le_and_ge_ << '\n';
    if (output_smallest_num_less_than_info) {
      sstrm << "smallest ever (# of elements <= pivot)/N is " << (smallest_ever_num_elements_less_than_or_equal_ratio_) << '\n';
    }
    if (output_smallest_num_greater_than_info) {
      sstrm << "smallest ever (# of elements >= pivot)/N is " << (smallest_ever_num_elements_greater_than_or_equal_ratio_) << '\n';
    }
    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_ever_num_elements_less_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_ever_num_elements_greater_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_ever_of_smaller_of_num_le_and_ge_ratios_ = 1.0l;
template<class ValueType, typename size_type>
size_type TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::N_for_smallest_ever_of_smaller_of_num_le_and_ge_ = 1;
template<class ValueType, typename size_type>
size_type TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smaller_of_num_le_and_ge_for_smallest_ever_ = 1;

//To only increase N linearly set N_exponential_growth_divisor equal to 0.
//If N < N_switch_to_exponential_growth_threshold then the next value of N
// that will be tested will be N + N_linear_increment_size, otherwise
// (assuming that N_exponential_growth_divisor != 0)
// the next value of N that will be tested will be
// N + (N / N_exponential_growth_divisor).
template<class ValueType = int, typename size_type = std::size_t>
void MeasureBiHeapifyInwardsPivotProperties(size_type num_tests = (1u << 20),
                                            size_type N_start = 101,
                                            size_type N_one_past_end = (1u << 17),
                                            size_type N_exponential_growth_divisor = 11,
                                            size_type N_linear_increment_size = 1,
                                            size_type N_switch_to_exponential_growth_threshold = (1u << 9)) {
  assert(N_linear_increment_size > 0);
  for (size_type N = N_start; N < N_one_past_end; ) {
    TestBiHeapifyInwardsOnGivenSize<ValueType> test(N);
    test.verbose_ = false;
    std::cout << test.GetStartOfDescriptionOfResult(num_tests);
    std::cout.flush();
    for (std::size_t i = 0; i < num_tests; i++) {
      test.PerformTestWithBiHeapifyInwardsNicerMath();
      //test.PerformTest();
    }
    //The two key quantities to look at are:
    //smallest min(# elements <= pivot, # elements >= pivot), and
    //average of min(# elements <= pivot, # elements >= pivot)
    std::cout << test.GetDescriptionOfResult(false) << std::endl;
    if (N_exponential_growth_divisor == 0
        || static_cast<size_type>(N / N_exponential_growth_divisor) == 0
        || N < N_switch_to_exponential_growth_threshold)
      N += N_linear_increment_size; //Grow N linearly
    else
      N += N / N_exponential_growth_divisor; //Grow N exponentially
  }
}

#endif /* BIHEAPIFY_INWARDS_PIVOT_TESTING_H_ */
