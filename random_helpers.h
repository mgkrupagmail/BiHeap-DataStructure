/*
 * random_helpers.h
 *
 *  Created on: May 5, 2017
 *      Author: Matthew Gregory Krupa
 */

#ifndef RANDOM_HELPERS_H_
#define RANDOM_HELPERS_H_

#include <iostream>
#include <random>
#include <vector>

namespace randomhelpers {

//This general case is designed to work with any type T for which
// std::is_integral<A>::value == true.
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
void FillWithRandomNumbers(Iterator start, Iterator one_past_end, float a, float b) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<float> dist(a, b);
  for (auto it = start; it != one_past_end; it++)
    *it = dist(generator);
  return ;
}

template<class Iterator>
void FillWithRandomNumbers(Iterator start, Iterator one_past_end, double a, double b) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<double> dist(a, b);
  for (auto it = start; it != one_past_end; it++)
    *it = dist(generator);
  return ;
}

template<class Iterator>
void FillWithRandomNumbers(Iterator start, Iterator one_past_end, long double a, long double b) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<long double>  dist(a, b);
  for (auto it = start; it != one_past_end; it++)
    *it = dist(generator);
  return ;
}

} //END NAMESPACE randomhelpers

#endif /* RANDOM_HELPERS_H_ */
