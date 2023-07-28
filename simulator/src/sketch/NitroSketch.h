/**
 * @file NitroSketch
 * @author 
 * @brief Implementation of NitroSketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>

namespace OmniSketch::Sketch {
/**
 * @brief NitroSketch
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 * @tparam T        type of the counter
 */

template <int32_t key_len, typename T, 
          typename hash_t = Hash::AwareHash> 
class NitroSketch : public SketchBase<key_len, T> {
public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  NitroSketch(int depth, int width);
  /**
   * @brief Release the pointer
   *
   */
  ~NitroSketch();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T value);

  void alwaysLineRateUpdate(const FlowKey<key_len> &flowkey, T value);

  void alwaysCorrectUpdate(const FlowKey<key_len> &flowkey, T value);

  T query(const FlowKey<key_len> &flowkey) const;

  void adjustUpdateProb(double traffic_rate);

  void clear();
  
  size_t size() const;

private:
  int depth_;
  int width_;

  T **array_;
  hash_t *hash_fns_;
  double *square_sum_; // maintain the square sum of each hash table
                       // used in isLineRateUpdate()

  int32_t next_bucket_; // index of the next update packet
  int32_t next_packet_; // number of skipped packets

  double update_prob_; // sampling probability
  std::default_random_engine generator;

  bool line_rate_enable_; // enable always line rate
  double switch_thresh_;  // switch to always line rate update

  static double update_probs[8];

  void getNextUpdate(double prob); // update next_bbucket_ and next_packet_
  bool isLineRateUpdate();         // judge whether enable line rate update

  void __do_update(const FlowKey<key_len> &flowkey, T value, double prob);
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
double NitroSketch<key_len, T, hash_t>::update_probs[8] = {
    1.0, 1.0 / 2, 1.0 / 4, 1.0 / 8, 1.0 / 16, 1.0 / 32, 1.0 / 64, 1.0 / 128};

template <int32_t key_len, typename T, typename hash_t>
NitroSketch<key_len, T, hash_t>::NitroSketch(int depth, int width)
    : depth_(depth), width_(Util::NextPrime(width)) {

  switch_thresh_ = (1.0 + std::sqrt(11.0 / width_)) * width_ * width_;

  line_rate_enable_ = false; // always correct update by default
  update_prob_ = 1.0;

  next_packet_ = 1;
  next_bucket_ = 0;

  hash_fns_ = new hash_t[depth_ * 2];

  square_sum_ = new double[depth_]();
  array_ = new T *[depth_];
  array_[0] = new T[depth_ * width_]();
  for (int i = 1; i < depth_; i++) {
    array_[i] = array_[i - 1] + width_;
  }
}

template <int32_t key_len, typename T, typename hash_t>
NitroSketch<key_len, T, hash_t>::~NitroSketch() {
  delete[] hash_fns_;
  delete[] array_[0];
  delete[] array_;
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::alwaysLineRateUpdate(
    const FlowKey<key_len> &flowkey, T value) {
  __do_update(flowkey, value, update_prob_);
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::alwaysCorrectUpdate(
    const FlowKey<key_len> &flowkey, T value) {
  if (isLineRateUpdate()) {
    __do_update(flowkey, value, update_prob_);
  } else {
    __do_update(flowkey, value, 1.0);
  }
}

template <int32_t key_len, typename T, typename hash_t>
T NitroSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{
  T median;
  T values[depth_];
  for (int i = 0; i < depth_; i++) {
    int index = hash_fns_[i](flowkey) % width_;
    values[i] = array_[i][index] *
                (2 * static_cast<int>(hash_fns_[depth_ + i](flowkey) & 1) - 1);
  }
  std::sort(values, values + depth_);
  if (depth_ & 1) {
    median = std::abs(values[depth_ / 2]);
  } else {
    median = std::abs((values[depth_ / 2 - 1] + values[depth_ / 2]) / 2);
  }
  return median;
}

template <int32_t key_len, typename T, typename hash_t>
std::size_t NitroSketch<key_len, T, hash_t>::size() const{
  return sizeof(NitroSketch<key_len, T, hash_t>) + depth_ * 2 * sizeof(hash_t) +
         depth_ * sizeof(T *) + depth_ * width_ * sizeof(T);
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::clear() {
  std::fill(array_[0], array_[0] + depth_ * width_);
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::__do_update(const FlowKey<key_len> &flowkey,
                                         T value, double prob) {
  next_packet_--; // skip packets
  if (next_packet_ == 0) {
    int i;
    for (;;) {
      i = next_bucket_;
      int index = hash_fns_[i](flowkey) % width_;

      double delta =
          1.0 * value / prob *
          (2 * static_cast<int>(hash_fns_[depth_ + i](flowkey) & 1) - 1);

      square_sum_[i] += (2.0 * array_[i][index] + delta) * delta;
      array_[i][index] += static_cast<T>(delta);

      getNextUpdate(prob);
      if (next_packet_ > 0) {
        break;
      }
    }
  }
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::getNextUpdate(double prob) {
  int sample = 1;
  if (prob < 1.0) {
    std::geometric_distribution<int> dist(prob);
    sample = 1 + dist(generator);
  }
  next_bucket_ = next_bucket_ + sample;
  next_packet_ = ((int)(next_bucket_ / depth_));
  next_bucket_ %= depth_;
}

template <int32_t key_len, typename T, typename hash_t>
bool NitroSketch<key_len, T, hash_t>::isLineRateUpdate() {
  if (line_rate_enable_) {
    return true;
  } else {
    double median;
    double values[depth_];
    for (int i = 0; i < depth_; i++) {
      values[i] = square_sum_[i];
    }
    std::sort(values, values + depth_);
    if (depth_ & 1) {
      median = values[depth_ / 2];
    } else {
      median = (values[depth_ / 2 - 1] + values[depth_ / 2]) / 2;
    }
    if (median >= switch_thresh_) {
      std::cout << "line rate update enable\n";
      line_rate_enable_ = true;
    }
    return line_rate_enable_;
  }
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::adjustUpdateProb(double traffic_rate) {
  int log_rate = static_cast<int>(std::log2(traffic_rate));
  int update_index = std::max(0, std::min(log_rate, 7));
  update_prob_ = update_probs[update_index];
}

template <int32_t key_len, typename T, typename hash_t>
void NitroSketch<key_len, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T value) {
  __do_update(flowkey, value, update_prob_);
}

} // namespace OmniSketch::Sketch
