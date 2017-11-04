/*
 * biheapify_inwards_pivot_testing.h
 *
 *  Created on: Nov 3, 2017
 *      Author: Matthew Gregory Krupa
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
  std::vector<ValueType> vec_, vec_original_;
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
  static long double smallest_num_elements_less_than_or_equal_ratio_;
  static long double smallest_num_elements_greater_than_or_equal_ratio_;
  static long double smallest_of_smaller_of_num_le_and_ge_ratios_;
  long double smallest_num_elements_less_than_or_equal_ratio_this_size_ = 1.0l;
  long double smallest_num_elements_greater_than_or_equal_ratio_this_size_ = 1.0l;
  long double smallest_of_smaller_of_num_le_and_ge_ratios_this_size_ = 1.0;

  TestBiHeapifyInwardsOnGivenSize(size_type N) : N_(N) {
    vec_ = std::vector<ValueType>(N_);
    vec_original_ = vec_;
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
    vec_original_ = vec_;
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

  void RunBiHeapifyInwards() {
    BiHeapifyInwards(vec_.begin(), N_);
    total_num_test_runs_++;
  }

  void Analyze() {
    auto first = vec_.begin();
    auto left_middle_value = *(first + left_index_);
    auto right_middle_value = *(first + right_index_);

    size_type num_elements_less_than_or_equal    = GetNumberOfElementsLessThanOrEqualToValue(left_middle_value);
    size_type num_elements_greater_than_or_equal = GetNumberOfElementsGreaterThanOrEqualToValue(right_middle_value);
    size_type min_num_elements_ge_and_le_pivot_  = num_elements_less_than_or_equal < num_elements_greater_than_or_equal ? num_elements_less_than_or_equal : num_elements_greater_than_or_equal;
    total_num_elements_less_than_or_equal_to_pivot_    += num_elements_less_than_or_equal;
    total_num_elements_greater_than_or_equal_to_pivot_ += num_elements_greater_than_or_equal;
    total_min_num_elements_ge_and_le_pivot_            += min_num_elements_ge_and_le_pivot_;

    long double num_elements_less_than_or_equal_ratio    = static_cast<long double>(num_elements_less_than_or_equal) / static_cast<long double>(N_);
    long double num_elements_greater_than_or_equal_ratio = static_cast<long double>(num_elements_greater_than_or_equal) / static_cast<long double>(N_);
    long double smaller_of_num_le_and_ge_ratios = num_elements_less_than_or_equal_ratio < num_elements_greater_than_or_equal_ratio ? num_elements_less_than_or_equal_ratio : num_elements_greater_than_or_equal_ratio;

    if (num_elements_greater_than_or_equal < smallest_total_num_elements_greater_than_or_equal_to_pivot_)
      smallest_total_num_elements_greater_than_or_equal_to_pivot_ = num_elements_greater_than_or_equal;
    if (num_elements_less_than_or_equal < smallest_total_num_elements_less_than_or_equal_to_pivot_)
      smallest_total_num_elements_less_than_or_equal_to_pivot_ = num_elements_less_than_or_equal;

    if (num_elements_less_than_or_equal_ratio < smallest_num_elements_less_than_or_equal_ratio_)
      smallest_num_elements_less_than_or_equal_ratio_ = num_elements_less_than_or_equal_ratio;
    if (num_elements_greater_than_or_equal_ratio < smallest_num_elements_greater_than_or_equal_ratio_)
      smallest_num_elements_greater_than_or_equal_ratio_ = num_elements_greater_than_or_equal_ratio;

    if (smaller_of_num_le_and_ge_ratios < smallest_of_smaller_of_num_le_and_ge_ratios_) {
      if (verbose_)
        std::cout << "smallest min(# elements <= pivot, # elements >= pivot) ever seen = "
                  << std::setprecision(30) << smaller_of_num_le_and_ge_ratios << std::endl;
      smallest_of_smaller_of_num_le_and_ge_ratios_ = smaller_of_num_le_and_ge_ratios;
    }

    if (num_elements_less_than_or_equal_ratio < smallest_num_elements_less_than_or_equal_ratio_this_size_)
      smallest_num_elements_less_than_or_equal_ratio_this_size_ = num_elements_less_than_or_equal_ratio;
    if (num_elements_greater_than_or_equal_ratio < smallest_num_elements_greater_than_or_equal_ratio_this_size_)
      smallest_num_elements_greater_than_or_equal_ratio_this_size_ = num_elements_greater_than_or_equal_ratio;

    if (smaller_of_num_le_and_ge_ratios < smallest_of_smaller_of_num_le_and_ge_ratios_this_size_) {
      if (verbose_)
        std::cout << "smallest min(# elements <= pivot, # elements >= pivot) for this size = "
                  << std::setprecision(30) << smaller_of_num_le_and_ge_ratios << std::endl;
      smallest_of_smaller_of_num_le_and_ge_ratios_this_size_ = smaller_of_num_le_and_ge_ratios;
    }
  }

  void PerformTest() {
    FillVecWithRandomNumbers();
    RunBiHeapifyInwards();
    Analyze();
  }

  std::string GetDescriptionOfResult() {
    std::stringstream sstrm;
    sstrm << "N = " << N_ << " \tNumber of tests performed = " << total_num_test_runs_ << '\n';
    sstrm << "For this size N = " << N_ << '\n';
    size_type min_of_num_le_and_num_ge = smallest_total_num_elements_less_than_or_equal_to_pivot_;
    if (smallest_total_num_elements_greater_than_or_equal_to_pivot_ < min_of_num_le_and_num_ge)
      min_of_num_le_and_num_ge = smallest_total_num_elements_greater_than_or_equal_to_pivot_;

    long double min_of_num_le_and_num_ge_ratio = static_cast<long double>(min_of_num_le_and_num_ge) / static_cast<long double>(N_);
    sstrm << "smallest min(# elements <= pivot, # elements >= pivot) is " << min_of_num_le_and_num_ge;
    sstrm << " \t(" << (100.0l * min_of_num_le_and_num_ge_ratio) << "% of N)" << '\n';
    long double smallest_total_num_elements_less_than_or_equal_to_pivot_ratio = static_cast<long double>(smallest_total_num_elements_less_than_or_equal_to_pivot_) / static_cast<long double>(N_);
    sstrm << "smallest # of elements <= pivot is                        " << smallest_total_num_elements_less_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * smallest_total_num_elements_less_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    long double smallest_total_num_elements_greater_than_or_equal_to_pivot_ratio = static_cast<long double>(smallest_total_num_elements_greater_than_or_equal_to_pivot_) / static_cast<long double>(N_);
    sstrm << "smallest # of elements >= pivot is                        " << smallest_total_num_elements_greater_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * smallest_total_num_elements_greater_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';

    long double ave_total_min_num_elements_ge_and_le_pivot_ = static_cast<long double>(total_min_num_elements_ge_and_le_pivot_) / static_cast<long double>(total_num_test_runs_);
    long double ave_total_min_num_elements_ge_and_le_pivot_ratio = ave_total_min_num_elements_ge_and_le_pivot_ / static_cast<long double>(N_);
    sstrm << "average of min(# elements <= pivot, # elements >= pivot) is " << ave_total_min_num_elements_ge_and_le_pivot_;
    sstrm << " \t(" << (100.0l * ave_total_min_num_elements_ge_and_le_pivot_ratio) << "% of N)" << '\n';
    long double ave_num_elements_less_than_or_equal_to_pivot_ = static_cast<long double>(total_num_elements_less_than_or_equal_to_pivot_) / static_cast<long double>(total_num_test_runs_);
    long double ave_num_elements_less_than_or_equal_to_pivot_ratio = ave_num_elements_less_than_or_equal_to_pivot_ / static_cast<long double>(N_);
    sstrm << "average  # of elements <= pivot is                          " << smallest_total_num_elements_less_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * ave_num_elements_less_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    long double ave_num_elements_greater_than_or_equal_to_pivot_ = static_cast<long double>(total_num_elements_greater_than_or_equal_to_pivot_) / static_cast<long double>(total_num_test_runs_);
    long double ave_num_elements_greater_than_or_equal_to_pivot_ratio = ave_num_elements_greater_than_or_equal_to_pivot_ / static_cast<long double>(N_);
    sstrm << "average  # of elements >= pivot is                          " << smallest_total_num_elements_greater_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * ave_num_elements_greater_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';

    //sstrm << '\n';
    sstrm << "Over all sizes N tested so far:\n";
    sstrm << "smallest ever min(# elements <= pivot, # elements >= pivot)/N is " << (smallest_of_smaller_of_num_le_and_ge_ratios_) << '\n';
    sstrm << "smallest ever (# of elements <= pivot)/N is " << (smallest_num_elements_less_than_or_equal_ratio_) << '\n';
    sstrm << "smallest ever (# of elements >= pivot)/N is " << (smallest_num_elements_greater_than_or_equal_ratio_) << '\n';
    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_num_elements_less_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_num_elements_greater_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestBiHeapifyInwardsOnGivenSize<ValueType, size_type>::smallest_of_smaller_of_num_le_and_ge_ratios_ = 1.0l;

void MeasureBiHeapifyInwardsPivotProperties() {
  for (std::size_t N = 110; N < (1u << 17); N += N / 10) {
    TestBiHeapifyInwardsOnGivenSize<int> test(N);
    test.verbose_ = false;
    for (std::size_t i = 0; i < (1u << 15); i++) {
      test.PerformTest();
    }
    //The two key quantities to look at are:
    //smallest min(# elements <= pivot, # elements >= pivot), and
    //average of min(# elements <= pivot, # elements >= pivot)
    std::cout << test.GetDescriptionOfResult() << std::endl;
  }
}

#endif /* BIHEAPIFY_INWARDS_PIVOT_TESTING_H_ */
