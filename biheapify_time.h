/*
 * biheapify_time.h
 *
 *  Created on: Jun 27, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Copyright Matthew Gregory Krupa
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

#include "../../CppBiHeapMedian/src/biheapify_simple.h"
#include "biheapify.h"
#include "biheapify_single_pass_success_rate.h"
#include "random_helpers.h"
/* Similar to BiHeapifyOdd(), except that it biheapifies an odd collection of
 *  elements using a pair of calls to BiHeapifyEven() in lieu of the
 *  single call to BiHeapifyOdd() that's done in BiHeapifyOdd().
 * It is not clear why this function empirically succeeds in biheapifying
 *  an odd number of elements.
 */
/* REMARK:
 * (1) If you were to replace the two calls to BiHeapifyEven() with
 *      two calls to BiHeapifyOdd(), then there is input data that the
 *      resulting function will NEVER succeed in biheapifying
 *      no matter how many times it repeatedly calls these two functions.
 *     This starts to happen when total_num_nodes >= 11.
 *     It is not clear why BiHeapifyUsingBiHeapifyEven() empirically always
 *      appears to produce a biheap while this fails to be true if
 *      BiHeapifyEven() is replaced with BiHeapifyOdd().
 */
template<class RAI, typename size_type = std::size_t>
inline void BiHeapifyUsingBiHeapifyEven(RAI first, size_type total_num_nodes) {
  if (total_num_nodes < 3) {
    if (total_num_nodes == 2 && *first > *(first + 1))
      std::iter_swap(first, first + 1);
    return ;
  }

  auto counter = 0u;
  do {
    BiHeapify<RAI, size_type>(first + 1, total_num_nodes - 1);
    BiHeapify<RAI, size_type>(first, total_num_nodes - 1);
    if (counter != 0 && counter % 3 == 0) {
      std::cerr << "BiHeapifyUsingBiHeapifyEven() could not biheapify after "
                << (counter + 1) << " call(s)." << std::endl;
      assert(false);
    }
    counter++;
  } while(!IsBiHeap<RAI, size_type>(first, total_num_nodes));
  return ;
}

template<class T, typename size_type = std::size_t>
std::chrono::nanoseconds TimeBiHeapifyGivenVec(T &vec,
                             const T &vec_original, std::size_t total_num_nodes,
                             int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  //Start at -1 in order to load the instructions into the cache
  // (this does affect timing).
  for (int num_repititions_counter = -1; num_repititions_counter <
                                   num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapify<typename T::iterator, size_type>(vec.begin(), total_num_nodes);
    if (num_repititions_counter >= 0)
      total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T, typename size_type = std::size_t>
std::chrono::nanoseconds TimeBiHeapifyUsingBiHeapifyEvenGivenVec(T &vec,
                             const T &vec_original, std::size_t total_num_nodes,
                             int num_repititions = 1) {
  std::chrono::nanoseconds total{0};
  for (int num_repititions_counter = -1; num_repititions_counter <
                                   num_repititions; num_repititions_counter++) {
    auto start_time = std::chrono::high_resolution_clock::now();
    BiHeapifyUsingBiHeapifyEven<typename T::iterator, size_type>(vec.begin(), total_num_nodes);
    if (num_repititions_counter >= 0)
      total += std::chrono::high_resolution_clock::now() - start_time;
    vec = vec_original;
  }
  return static_cast<std::chrono::nanoseconds>(total);
}

template<class T, typename size_type = std::size_t>
std::chrono::nanoseconds::rep TimeBiHeapifyOnGivenSize(
                           std::size_t total_num_nodes, int num_vecs_to_try = 1,
                           int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                  std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    vec_original = vec;
    total += TimeBiHeapifyGivenVec<std::vector<T>, size_type>(vec, vec_original,
                                      total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

template<class T, typename size_type = std::size_t>
std::chrono::nanoseconds::rep TimeBiHeapifyUsingBiHeapifyEvenOnGivenSize(
                           std::size_t total_num_nodes, int num_vecs_to_try = 1,
                           int num_repititions_per_vec = 1) {
  std::chrono::nanoseconds total{0};
  std::vector<T> vec(total_num_nodes), vec_original(total_num_nodes);
  for (int i = 0; i < num_vecs_to_try; i++) {
    randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                  std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    vec_original = vec;
    total += TimeBiHeapifyUsingBiHeapifyEvenGivenVec<std::vector<T>, size_type>(vec,
                            vec_original, total_num_nodes, num_repititions_per_vec);
  }
  return total.count();
}

std::string GetTimeString(std::chrono::nanoseconds::rep time_rep, long long divisor) {
  std::stringstream strm;
  strm << time_rep/divisor << " ns  = " << (time_rep/divisor)/1000 << " mus  = "
       << (time_rep/divisor)/1000000 << " ms ";
  return strm.str();
}

template<class T, typename size_type = std::size_t>
void TimeBiHeapifies(long start_total_num_nodes,
                     long end_total_num_nodes,
                     long num_vecs_to_try = 1,
                     long num_repititions_per_vec = 1,
                     long increment_size = 1,
                     long long divisor = 0,
                     bool should_time_BiHeapify = true,
                     bool should_time_BiHeapifySafe = true,
                     bool should_time_BiHeapifyUsingEven = true) {
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
    std::cout << "total_num_nodes = " << i << " \t";

    if (should_time_BiHeapify) {
      std::cout << "BiHeapify() ave = ";
      std::cout.flush();
      auto biheapify = TimeBiHeapifyOnGivenSize<T, size_type>(i, num_vecs_to_try,
                                                   num_repititions_per_vec);
      biheapify_times[(i - start_total_num_nodes) / increment_size] = biheapify;
      std::cout << GetTimeString(biheapify, divisor);
    }

    if (should_time_BiHeapifyUsingEven) {
      if (should_time_BiHeapify || should_time_BiHeapifySafe)
        std::cout << "  \t";
      std::cout << "BiHeapifyUsingBiHeapifyEven() ave = ";
      std::cout.flush();
      auto biheapify_using_even = TimeBiHeapifyUsingBiHeapifyEvenOnGivenSize<T, size_type>(
                                   i, num_vecs_to_try, num_repititions_per_vec);

      if (!should_time_BiHeapify && !should_time_BiHeapifySafe)
        biheapify_times[(i - start_total_num_nodes) / increment_size] =
                                                           biheapify_using_even;
      std::cout << GetTimeString(biheapify_using_even, divisor);
    }
    std::cout << std::endl;
  }

  std::cout << "\nTime differences between consecutive biheap sizes:\n";
  for (auto i = 1u; i < biheapify_times.size(); i++) {
    auto diff = static_cast<size_type>(biheapify_times[i]
                                         - biheapify_times[i - 1]);
    std::cout << diff << '\n';
  }
  std::cout .flush();
}

#endif /* BIHEAPIFY_TIME_H_ */
