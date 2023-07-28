/**
 * @file UnivMon.h
 * @author dromniscience XierLabber (you@domain.com)
 * @brief Implementation of UnivMon
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once
#include <iostream>

#include <common/hash.h>
#include <common/sketch.h>
#include <sketch/CountSketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief UnivMon
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class UnivMon : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;
  int32_t logn;
  hash_t *hash_fns;
  CountSketch<key_len, T, hash_t> **sketch;
  Data::Estimation<key_len> *flows;

  UnivMon(const UnivMon &) = delete;
  UnivMon(UnivMon &&) = delete;

public:
  /**
   * @brief Construct by specifying depth, width and $\log n$, where $n$ is the
   * number of flows to insert.
   *
   */
  UnivMon(int32_t depth_, int32_t width_, int32_t log_n);
  /**
   * @brief Release the pointer
   *
   */
  ~UnivMon();
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
  int32_t getHash(const FlowKey<key_len> &flowkey, int32_t layer) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
int32_t UnivMon<key_len, T, hash_t>::getHash(const FlowKey<key_len> &flowkey, int32_t layer) const{
  if(layer == 0){
    return 1;
  }
  return (hash_fns[layer - 1](flowkey) & 1);
}

template <int32_t key_len, typename T, typename hash_t>
UnivMon<key_len, T, hash_t>::UnivMon(int32_t depth_, int32_t width_,
                                     int32_t log_n)
    : depth(depth_), width(Util::NextPrime(width_)), logn(log_n) {

  hash_fns = new hash_t[logn - 1];
  // Allocate continuous memory
  sketch = new CountSketch<key_len, T, hash_t> *[logn];
  for (int32_t i = 0; i < logn; ++i){
    sketch[i] = new CountSketch<key_len, T, hash_t>(depth, width_);
    width_ = std::max(1, width_ / 2);
  }
  flows = new Data::Estimation<key_len>[logn];
}

template <int32_t key_len, typename T, typename hash_t>
UnivMon<key_len, T, hash_t>::~UnivMon() {
  if (hash_fns)
    delete[] hash_fns;
  if (sketch) {
    for (int32_t i = 0; i < logn; ++i)
      delete sketch[i];
    delete[] sketch;
  }
  if (flows)
    delete[] flows;
}

template <int32_t key_len, typename T, typename hash_t>
void UnivMon<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                         T val) {
  for (int32_t i = 0; i < logn; ++i) {
    if (getHash(flowkey, i)) {
      sketch[i]->update(flowkey, val);
    } else
      break;
  }
}

template <int32_t key_len, typename T, typename hash_t>
T UnivMon<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  int level;
  for(level = 0; level < logn; level++){
    if(!getHash(flowkey, level)){
      break;
    }
  }
  level--;
  T ret = sketch[level]->query(flowkey);
  for(int i = level - 1; i >= 0; i--){
    ret = 2 * ret - sketch[i]->query(flowkey);
  }
  return ret;
}

/*
template <int32_t key_len, typename T, typename hash_t>
T UnivMon<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  static bool cnt_distrib = true;
  int32_t total = 0;
  if (cnt_distrib) {
    std::vector<double> distrib(32);
    for (int32_t t = 0; t < logn; ++t) {
      total += sketch[t]->getDepth() * sketch[t]->getWidth();
      for (int32_t i = 0; i < sketch[t]->getDepth(); ++i)
        for (int32_t j = 0; j < sketch[t]->getWidth(); ++j) {
          for (int32_t k = 0; k < 32; ++k) {
            if (std::abs(sketch[t]->getCnt(i,j)) >= (1 << k))
              distrib[k] += 1.0;
            else
              break;
          }
        }
    }
    for (int32_t k = 0; k < 32; ++k) {
      std::cout << distrib[k] / total << " ";
      if (distrib[k] == 0.0) {
        std::cout << std::endl;
        break;
      }
    }
    cnt_distrib = false;
  }
  return 0;
}
*/

template <int32_t key_len, typename T, typename hash_t>
size_t UnivMon<key_len, T, hash_t>::size() const {
  size_t total = sizeof(*this); // instance
  for (int32_t i = 0; i < logn; ++i)
    total += sketch[i]->size(); // L2 HH
  total += sizeof(hash_t) * (logn - 1);
  return total;
}

template <int32_t key_len, typename T, typename hash_t>
void UnivMon<key_len, T, hash_t>::clear() {
  // std::fill(counter[0], counter[0] + depth * width, 0);
}

} // namespace OmniSketch::Sketch


