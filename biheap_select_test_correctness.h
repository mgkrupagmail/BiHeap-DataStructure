/*
 * biheap_select_test_correctness.h
 *
 *  Created on: May 6, 2019
 *      Author: Matthew Gregory Krupa
 *      Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAP_SELECT_TEST_CORRECTNESS_H_
#define BIHEAP_SELECT_TEST_CORRECTNESS_H_

#include <vector>

#include "biheap_select.h"


//***********************************************
// TESTING FOR CORRECTNESS
//***********************************************


template<class ValueType, typename size_type = std::size_t>
struct TestCorrectnessOfBiHeapifySelectMediansOnGivenSize {
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

  TestCorrectnessOfBiHeapifySelectMediansOnGivenSize(size_type N) : N_(N) {
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

  void RunBiHeapSelect() {
    BiHeapSelect(vec_.begin(), N_, N_ / 2);
    if (N_ % 2 == 0) {
      BiHeapify(vec_.begin(), N_ / 2); //Emplace the maximum of V+0, ..., V + ((N_/2) - 1)
    }
    total_num_test_runs_++;
    return ;
  }

  void PerformTest(long double scale = 0.0l) {
    if (scale <= 0.0l)
      FillVecWithRandomNumbers();
    else
      FillVecWithRandomNumbersSizeScaledInterval(scale);
    RunBiHeapSelect();
    Analyze();
    return ;
  }

  std::string GetDescriptionOfResult() {
    std::stringstream sstrm;
    sstrm << "N = " << N_ << " \t(N+1)/2 = " << ((N_ + 1) / 2) << " \tNumber of tests performed = " << total_num_test_runs_ << '\n';
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
    sstrm << "average  # of elements <= pivot is                          " << ave_num_elements_less_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * ave_num_elements_less_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';
    long double ave_num_elements_greater_than_or_equal_to_pivot_ = static_cast<long double>(total_num_elements_greater_than_or_equal_to_pivot_) / static_cast<long double>(total_num_test_runs_);
    long double ave_num_elements_greater_than_or_equal_to_pivot_ratio = ave_num_elements_greater_than_or_equal_to_pivot_ / static_cast<long double>(N_);
    sstrm << "average  # of elements >= pivot is                          " << ave_num_elements_greater_than_or_equal_to_pivot_;
    sstrm << " \t(" << (100.0l * ave_num_elements_greater_than_or_equal_to_pivot_ratio) << "% of N)" << '\n';

    sstrm << "Over all sizes N tested so far:\n";
    sstrm << "smallest ever min(# elements <= pivot, # elements >= pivot)/N is " << (smallest_ever_of_smaller_of_num_le_and_ge_ratios_)
          << " = " << smaller_of_num_le_and_ge_for_smallest_ever_ << " / " << N_for_smallest_ever_of_smaller_of_num_le_and_ge_ << '\n';
    sstrm << "smallest ever (# of elements <= pivot)/N is " << (smallest_ever_num_elements_less_than_or_equal_ratio_) << '\n';
    sstrm << "smallest ever (# of elements >= pivot)/N is " << (smallest_ever_num_elements_greater_than_or_equal_ratio_) << '\n';
    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::smallest_ever_num_elements_less_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::smallest_ever_num_elements_greater_than_or_equal_ratio_ = 1.0l;
template<class ValueType, typename size_type>
long double TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::smallest_ever_of_smaller_of_num_le_and_ge_ratios_ = 1.0l;
template<class ValueType, typename size_type>
size_type TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::N_for_smallest_ever_of_smaller_of_num_le_and_ge_ = 1;
template<class ValueType, typename size_type>
size_type TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType, size_type>::smaller_of_num_le_and_ge_for_smallest_ever_ = 1;

#define DEFAULT_num_test (1u << 19)
#define DEFAULT_N_start 20
#define DEFAULT_N_one_past_end (1u << 23)
#define DEFAULT_permitted_pivot_positions std::vector<long double>()

template<class ValueType = int, typename size_type = std::size_t>
void TestCorrectnessOfBiHeapSelectMedians(size_type num_test = DEFAULT_num_test,
    size_type N_start = DEFAULT_N_start,
    size_type N_one_past_end = DEFAULT_N_one_past_end) {
  if (num_test <= 0)
    num_test = DEFAULT_num_test;
  if (N_start <= 0)
    N_start = DEFAULT_N_start;
  if (N_one_past_end <= 0)
    N_one_past_end = DEFAULT_N_one_past_end;

  for (size_type N = N_start; N < N_one_past_end; N += N / 11) {
    TestCorrectnessOfBiHeapifySelectMediansOnGivenSize<ValueType> test(N);
    test.verbose_ = false;
    for (std::size_t i = 0; i < num_test; i++) {
      test.PerformTest();
    }
    //The two key quantities to look at are:
    //smallest min(# elements <= pivot, # elements >= pivot), and
    //average of min(# elements <= pivot, # elements >= pivot)
    std::cout << test.GetDescriptionOfResult();
    std::cout.flush();
  }
}








/////////////////////////////////////////////////////////////////////////
// This test randomly selects a position in a random vector and emplaces
//  the value at that position.
/////////////////////////////////////////////////////////////////////////

template<class ValueType, typename size_type = std::size_t>
struct TestCorrectnessOfBiHeapSelectOnGivenSize {
  size_type N_, total_num_test_runs_, total_num_test_successes_, pivot_position_;
  std::vector<long double> permitted_pivot_positions_;
  long double pivot_average_ = 0.0l;
  static size_type N_of_smallest_ever_success_percentage_;
  long double success_percentage_;
  static long double smallest_ever_success_percentage_;
  std::vector<ValueType> vec_;//, vec_original_;
  bool verbose_ = true;
  ValueType lower_, upper_;

  TestCorrectnessOfBiHeapSelectOnGivenSize(size_type N, std::vector<long double> permitted_pivot_positions = std::vector<long double>()) : N_(N),
      permitted_pivot_positions_(permitted_pivot_positions) {
    vec_ = std::vector<ValueType>(N_);
    total_num_test_successes_ = 0;
    total_num_test_runs_ = 0;
    success_percentage_ = 1.0l;
    if (N_of_smallest_ever_success_percentage_ == 0) {//If this is the first time this class has been instantiated.
      N_of_smallest_ever_success_percentage_ = N_;
    }
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
    auto pivot_value = vec_[pivot_position_];

    size_type num_elements_less_than_or_equal    = GetNumberOfElementsLessThanOrEqualToValue(pivot_value);
    size_type num_elements_greater_than_or_equal = GetNumberOfElementsGreaterThanOrEqualToValue(pivot_value);
    if (num_elements_less_than_or_equal >= pivot_position_ + 1 && num_elements_greater_than_or_equal >= (N_ - pivot_position_)) {
      total_num_test_successes_++;
    }
    if (total_num_test_successes_ == total_num_test_runs_) {
      success_percentage_ = 1.0l;
    } else {
      success_percentage_ = static_cast<long double>(total_num_test_successes_) / static_cast<long double>(total_num_test_runs_);
      assert(success_percentage_ < 1.0l); //Make sure that we don't indicate any false successes due to rounding error.
    }
    if (success_percentage_ < smallest_ever_success_percentage_) {
      smallest_ever_success_percentage_ = success_percentage_;
      N_of_smallest_ever_success_percentage_ = N_;
    }

    return ;
  }

  void RunBiHeapSelect() {
    BiHeapSelect(vec_.begin(), N_, pivot_position_);
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
    RunBiHeapSelect();
    Analyze();
    return ;
  }

  std::string GetDescriptionOfResult() {
    std::stringstream sstrm;
    sstrm << "N = " << N_
          << " \t(N-1)/2 = " << (static_cast<long double>(N_ - 1) / 2.0l)
          << " \tpivot_average_ = " << pivot_average_
          << " \tNum successes/Num tests = " << total_num_test_successes_
          << " / " << total_num_test_runs_
          << " \t= " << success_percentage_
          << " \tMin ever % successes = " << smallest_ever_success_percentage_
          << " \t occurred when N = " << N_of_smallest_ever_success_percentage_
          << '\n';
    return sstrm.str();
  }
};

template<class ValueType, typename size_type>
long double TestCorrectnessOfBiHeapSelectOnGivenSize<ValueType, size_type>::smallest_ever_success_percentage_ = 1.0l;
template<class ValueType, typename size_type>
size_type TestCorrectnessOfBiHeapSelectOnGivenSize<ValueType, size_type>::N_of_smallest_ever_success_percentage_ = 0;

template<class ValueType = int, typename size_type = std::size_t>
void TestCorrectnessOfBiHeapSelect(size_type num_test = DEFAULT_num_test,
    size_type N_start = DEFAULT_N_start,
    size_type N_one_past_end = DEFAULT_N_one_past_end,
    std::vector<long double> permitted_pivot_positions = DEFAULT_permitted_pivot_positions) {
  if (num_test <= 0)
    num_test = DEFAULT_num_test;
  if (N_start <= 0)
    N_start = DEFAULT_N_start;
  if (N_one_past_end <= 0)
    N_one_past_end = DEFAULT_N_one_past_end;

  for (size_type N = N_start; N < N_one_past_end; N += N / 11) {
    TestCorrectnessOfBiHeapSelectOnGivenSize<ValueType> test(N, permitted_pivot_positions);
    test.verbose_ = false;
    for (std::size_t i = 0; i < num_test; i++) {
      test.PerformTest();
    }
    std::cout << test.GetDescriptionOfResult();
    std::cout.flush();
  }
}



#endif /* BIHEAP_SELECT_TEST_CORRECTNESS_H_ */
