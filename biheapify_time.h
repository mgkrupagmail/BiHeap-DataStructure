/*
 * biheapify_time.h
 *
 *  Created on: Jun 27, 2017
 *      Author: Matthew Gregory Krupa
 *
 *  The TimeBiHeapifies() functions times the various BiHeapify functions.
 *  An example call can be found in main().
 */

#ifndef BIHEAPIFY_TIME_H_
#define BIHEAPIFY_TIME_H_

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "biheap_common.h"
#include "biheapify.h"
#include "biheapify_simple.h"
#include "random_helpers.h"

template<class T>
std::chrono::nanoseconds TimeBiHeapifySafeGivenVec(T &vec, const T &vec_original,
                         std::size_t total_num_nodes, int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  {
  BiHeapifySafe(vec.begin(), total_num_nodes);
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter <
                                   num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapifySafe(vec.begin(), total_num_nodes);
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T>
std::chrono::nanoseconds TimeBiHeapifyGivenVec(T &vec,
                             const T &vec_original, std::size_t total_num_nodes,
                             int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  {
    if (total_num_nodes % 2 == 0){
      BiHeapifyEven(vec.begin(), total_num_nodes);
    } else {
      BiHeapifyOdd(vec.begin(), total_num_nodes);
  }
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter <
                                   num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapify(vec.begin(), total_num_nodes);
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T>
std::chrono::nanoseconds::rep TimeBiHeapifySafeOnGivenSize(
                          std::size_t total_num_nodes, int num_vecs_to_try = 1,
                          int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                  std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    vec_original = vec;
    total += TimeBiHeapifySafeGivenVec(vec, vec_original, total_num_nodes,
                                   num_repititions_per_vec);
  }
  return total.count();
}

template<class T>
std::chrono::nanoseconds::rep TimeBiHeapifyOnGivenSize(
                           std::size_t total_num_nodes, int num_vecs_to_try = 1,
                           int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                  std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    vec_original = vec;
    total += TimeBiHeapifyGivenVec(vec, vec_original, total_num_nodes,
                                             num_repititions_per_vec);
  }
  return total.count();
}

std::string GetTimeString(std::chrono::nanoseconds::rep time_rep, long long divisor) {
  std::stringstream strm;
  strm << time_rep/divisor << " ns  = "
       << (time_rep/divisor)/1000 << " mus  = "
       << (time_rep/divisor)/1000000 << " ms ";
  return strm.str();
}

template<class T> void TimeBiHeapifies(long start_total_num_nodes,
                                       long end_total_num_nodes,
                                       long num_vecs_to_try = 1,
                                       long num_repititions_per_vec = 1,
                                       long increment_size = 1,
                                       long long divisor = 0) {
  if (divisor == 0)
    divisor = num_vecs_to_try * num_repititions_per_vec;
  else if (divisor < 0)
    divisor = num_repititions_per_vec;
  auto num_sizes_tested = (end_total_num_nodes + 1 - start_total_num_nodes)
                                                               / increment_size;
  std::vector<std::chrono::nanoseconds::rep> biheapify_times(num_sizes_tested,
                                           std::chrono::nanoseconds{0}.count());
  for (auto i = start_total_num_nodes; i <= end_total_num_nodes; i +=
                                                              increment_size) {
    std::cout << "total_num_nodes = \t" << i << '\n';

    std::cout << "BiHeapify() ave       = ";
    std::cout.flush();
    auto biheapify_single_pass = TimeBiHeapifyOnGivenSize<T>(i,
                                      num_vecs_to_try, num_repititions_per_vec);
    std::cout << GetTimeString(biheapify_single_pass, divisor);

    std::cout << "  \tBiHeapifySafe() ave             = ";
    std::cout.flush();
    auto biheapify = TimeBiHeapifySafeOnGivenSize<T>(i, num_vecs_to_try,
                                                 num_repititions_per_vec);
    biheapify_times[(i - start_total_num_nodes) / increment_size] = biheapify;
    std::cout << GetTimeString(biheapify, divisor) << std::endl;
  }

  std::cout << "\nTime differences between consecutive biheap sizes:\n";
  for (auto i = 1u; i < biheapify_times.size(); i++) {
    auto diff = static_cast<std::size_t>(biheapify_times[i]
                                         - biheapify_times[i - 1]);
    std::cout << diff << '\n';
  }
  std::cout .flush();
}

#endif /* BIHEAPIFY_TIME_H_ */
