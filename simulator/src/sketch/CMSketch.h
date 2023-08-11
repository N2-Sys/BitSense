/**
 * @file CMSketch.h
 * @author dromniscience (dr_abc@pku.edu.cn)
 * @brief Implementation of Count Min Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>
#include <iostream>

namespace OmniSketch::Sketch {
/**
 * @brief Count Min Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class CMSketch : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;
  hash_t *hash_fns;
  T **counter;

  CMSketch(const CMSketch &) = delete;
  CMSketch(CMSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  CMSketch(int32_t depth_, int32_t width_);
  /**
   * @brief Release the pointer
   *
   */
  ~CMSketch();
  /**
   * @brief Update a flowkey with certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Query a flowkey
   *
   */
  T query(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const override;
  /**
   * @brief Reset the sketch
   *
   */
  void clear();
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
CMSketch<key_len, T, hash_t>::CMSketch(int32_t depth_, int32_t width_)
    : depth(depth_), width(Util::NextPrime(width_)) {

  hash_fns = new hash_t[depth];
  // Allocate continuous memory
  counter = new T *[depth];
  counter[0] = new T[depth * width](); // Init with zero
  for (int32_t i = 1; i < depth; ++i) {
    counter[i] = counter[i - 1] + width;
  }
}

template <int32_t key_len, typename T, typename hash_t>
CMSketch<key_len, T, hash_t>::~CMSketch() {
  delete[] hash_fns;
  delete[] counter[0];
  delete[] counter;
}

template <int32_t key_len, typename T, typename hash_t>
void CMSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
  for (int32_t i = 0; i < depth; ++i) {
    int32_t index = hash_fns[i](flowkey) % width;
    counter[i][index] += val;
  }
}

template <int32_t key_len, typename T, typename hash_t>
T CMSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  static bool cnt_distrib = true;
  if (cnt_distrib) {
    std::vector<double> distrib(32);
    for (int32_t i = 0; i < depth; ++i)
      for (int32_t j = 0; j < width; ++j) {
        for (int32_t k = 0; k < 32; ++k) {
          if (counter[i][j] >= (1 << k))
            distrib[k] += 1.0;
          else
            break;
        }
      }
    for (int32_t k = 0; k < 32; ++k) {
      std::cout << distrib[k] / depth / width << " ";
      if (distrib[k] == 0.0) {
        std::cout << std::endl;
        break;
      }
    }
    cnt_distrib = false;
  }
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < depth; ++i) {
    int32_t index = hash_fns[i](flowkey) % width;
    min_val = std::min(min_val, counter[i][index]);
  }
  return min_val;
}

template <int32_t key_len, typename T, typename hash_t>
size_t CMSketch<key_len, T, hash_t>::size() const {
  return sizeof(*this)                // instance
         + sizeof(hash_t) * depth     // hashing class
         + sizeof(T) * depth * width; // counter
}

template <int32_t key_len, typename T, typename hash_t>
void CMSketch<key_len, T, hash_t>::clear() {
  std::fill(counter[0], counter[0] + depth * width, 0);
}

} // namespace OmniSketch::Sketch
