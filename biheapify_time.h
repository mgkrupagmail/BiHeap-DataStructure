/*
 * biheapify_time.h
 *
 *  Created on: Jun 27, 2017
 *      Author: Matthew Gregory Krupa
 *
 *  The TimeBiHeapifies() functions times the various BiHeapify functions. An
 *  example call can be found in main().
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

template<class T> std::chrono::nanoseconds TimeBiHeapifyGivenVec(T &vec, const T &vec_original, std::size_t total_num_nodes, int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  {
  BiHeapify(vec.begin(), total_num_nodes);
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter < num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapify(vec.begin(), total_num_nodes);
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T> std::chrono::nanoseconds TimeBiHeapifySinglePassGivenVec(T &vec, const T &vec_original, std::size_t total_num_nodes, int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  {
  BiHeapifyOddSinglePass(vec.begin(), total_num_nodes);
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter < num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (total_num_nodes % 2 == 0){
      BiHeapifyEvenSinglePass(vec.begin(), total_num_nodes);
    } else {
      BiHeapifyOddSinglePass(vec.begin(), total_num_nodes);
    }
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T> std::chrono::nanoseconds TimeBiHeapifySimpleGivenVec(T &vec, const T &vec_original, std::size_t total_num_nodes, int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  {
  BiHeapifySimple(vec.begin(), total_num_nodes);
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter < num_repititions_per_vec; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapifySimple(vec.begin(), total_num_nodes);
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return total;
}

template<class T> std::chrono::nanoseconds TimeBiHeapifySimpleSinglePassGivenVec(T &vec, const T &vec_original, std::size_t total_num_nodes, int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  {
  BiHeapifySimpleSinglePass(vec.begin(), total_num_nodes);
  vec = vec_original;
  }
  for (int num_repititions_counter = 0; num_repititions_counter < num_repititions_per_vec; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapifySimpleSinglePass(vec.begin(), total_num_nodes);
    total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return total;
}

template<class T> std::chrono::nanoseconds::rep TimeBiHeapifyOnGivenSize(std::size_t total_num_nodes,
                                                                         int num_vecs_to_try = 1,
                                                                         int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillVectorWithRandomNumbers<T>(vec, std::numeric_limits<T>::min(), std::numeric_limits<T>::max()); //(vec, 0, total_num_nodes);
    vec_original = vec;
    total += TimeBiHeapifyGivenVec(vec, vec_original, total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

template<class T> std::chrono::nanoseconds::rep TimeBiHeapifySinglePassOnGivenSize(std::size_t total_num_nodes,
                                                                         int num_vecs_to_try = 1,
                                                                         int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillVectorWithRandomNumbers<T>(vec, std::numeric_limits<T>::min(), std::numeric_limits<T>::max()); //(vec, 0, total_num_nodes);
    vec_original = vec;
    total += TimeBiHeapifySinglePassGivenVec(vec, vec_original, total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

template<class T> std::chrono::nanoseconds::rep TimeBiHeapifySimpleOnGivenSize(std::size_t total_num_nodes,
                                                                         int num_vecs_to_try = 1,
                                                                         int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillVectorWithRandomNumbers<T>(vec, std::numeric_limits<T>::min(), std::numeric_limits<T>::max()); //(vec, 0, total_num_nodes);
    vec_original = vec;
    total += TimeBiHeapifySimpleGivenVec(vec, vec_original, total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

template<class T> std::chrono::nanoseconds::rep TimeBiHeapifySimpleSinglePassOnGivenSize(std::size_t total_num_nodes,
                                                                         int num_vecs_to_try = 1,
                                                                         int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillVectorWithRandomNumbers<T>(vec, std::numeric_limits<T>::min(), std::numeric_limits<T>::max()); //(vec, 0, total_num_nodes);
    vec_original = vec;
    total += TimeBiHeapifySimpleSinglePassGivenVec(vec, vec_original, total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

template<class T> void TimeBiHeapifies(int start_total_num_nodes, int end_total_num_nodes,
                                       int num_vecs_to_try = 1, int num_repititions_per_vec = 1,
                                       int increment_size = 1) {
  long double num_repititions_per_total_num_nodes = static_cast<long double>(num_repititions_per_vec);//= static_cast<long double>(num_vecs_to_try * num_repititions_per_vec);
  auto num_sizes_tested = (end_total_num_nodes + 1 - start_total_num_nodes) / increment_size;
  std::vector<std::chrono::nanoseconds::rep> biheapify_times(num_sizes_tested, std::chrono::nanoseconds{0}.count());
  for (int i = start_total_num_nodes; i <= end_total_num_nodes; i += increment_size) {
    std::cout << "total_num_nodes = \t" << i << std::endl;
    auto biheapify = TimeBiHeapifyOnGivenSize<T>(i, num_vecs_to_try, num_repititions_per_vec);
    biheapify_times[(i - start_total_num_nodes) / increment_size] = biheapify;
    std::cout << "BiHeapify() ave       = " << static_cast<long long>(biheapify/num_repititions_per_total_num_nodes) << " ns ";
    std::cout << " = " << static_cast<long long>((biheapify/num_repititions_per_total_num_nodes)/1000) << " mus ";
    std::cout << " = " << static_cast<long long>((biheapify/num_repititions_per_total_num_nodes)/1000000) << " ms ";
    std::cout.flush();
    auto biheapify_single_pass = TimeBiHeapifySinglePassOnGivenSize<T>(i, num_vecs_to_try, num_repititions_per_vec);
    std::cout << "  \tBiHeapifySinglePass() ave       = " << static_cast<long long>(biheapify_single_pass/num_repititions_per_total_num_nodes) << " ns ";
    std::cout << " = " << static_cast<long long>((biheapify_single_pass/num_repititions_per_total_num_nodes)/1000) << " mus ";
    std::cout << " = " << static_cast<long long>((biheapify_single_pass/num_repititions_per_total_num_nodes)/1000000) << " ms ";
    std::cout << std::endl;
    auto biheapify_simple = TimeBiHeapifySimpleOnGivenSize<T>(i, num_vecs_to_try, num_repititions_per_vec);
    std::cout << "BiHeapifySimple() ave = " << static_cast<long long>(biheapify_simple/num_repititions_per_total_num_nodes)  << " ns ";
    std::cout << " = " << static_cast<long long>((biheapify_simple/num_repititions_per_total_num_nodes)/1000)  << " mus ";
    std::cout << " = " << static_cast<long long>((biheapify_simple/num_repititions_per_total_num_nodes)/1000000)  << " ms ";
    std::cout.flush();
    auto biheapify_simple_single_pass = TimeBiHeapifySimpleSinglePassOnGivenSize<T>(i, num_vecs_to_try, num_repititions_per_vec);
    std::cout << "  \tBiHeapifySimpleSinglePass() ave = " << static_cast<long long>(biheapify_simple_single_pass/num_repititions_per_total_num_nodes) << " ns ";
    std::cout << " = " << static_cast<long long>((biheapify_simple_single_pass/num_repititions_per_total_num_nodes)/1000) << " mus ";
    std::cout << " = " << static_cast<long long>((biheapify_simple_single_pass/num_repititions_per_total_num_nodes)/1000000) << " ms.";
    std::cout << std::endl;
  }

  std::cout << "\nTime differences between consecutive biheap sizes" << std::endl;
  for (auto i = 1u; i < biheapify_times.size(); i++) {
    std::cout << static_cast<long double>(biheapify_times[i] - biheapify_times[i - 1]) << "\n";
  }
  std::cout .flush();
}



#endif /* BIHEAPIFY_TIME_H_ */
