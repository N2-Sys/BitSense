/**
 * @file BSCountSketch.h
 * @author dromniscience (you@domain.com)
 * @brief Implementation of Count Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/bitsense.h>
#include <common/sketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief Count Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class BSCountSketch : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;
  hash_t *hash_fns;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  BitSense<no_layer, T, hash_t> *ch;

  BSCountSketch(const BSCountSketch &) = delete;
  BSCountSketch(BSCountSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth and width
   *
   */
  BSCountSketch(int32_t depth_, int32_t width_, double cnt_no_ratio,
                const std::vector<size_t> &width_cnt,
                const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~BSCountSketch();
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
  int32_t getDepth() const;
  int32_t getWidth() const;
  int32_t getIndex(int32_t i, int32_t j) const;
  T getCnt(int32_t i, int32_t j);
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t BSCountSketch<key_len, no_layer, T, hash_t>::getDepth() const {
  return depth;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t BSCountSketch<key_len, no_layer, T, hash_t>::getWidth() const {
  return width;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t BSCountSketch<key_len, no_layer, T, hash_t>::getIndex(int32_t i,
                                                              int32_t j) const {
  return i * width + j;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSCountSketch<key_len, no_layer, T, hash_t>::getCnt(int32_t i, int32_t j) {
  return ch->getCnt(getIndex(i, j));
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSCountSketch<key_len, no_layer, T, hash_t>::BSCountSketch(
    int32_t depth_, int32_t width_, double cnt_no_ratio,
    const std::vector<size_t> &width_cnt, const std::vector<size_t> &no_hash)
    : depth(depth_), width(Util::NextPrime(width_)), ch(nullptr),
      width_cnt(width_cnt), no_hash(no_hash) {

  // The first depth hash functions: CM
  // The last depth hash function: signed bit
  hash_fns = new hash_t[depth * 2];
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(depth) * width);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // BS
  ch = new BitSense<no_layer, T, hash_t>(no_cnt, width_cnt, no_hash,
                                                 true);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSCountSketch<key_len, no_layer, T, hash_t>::~BSCountSketch() {
  delete[] hash_fns;
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSCountSketch<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  for (int i = 0; i < depth; ++i) {
    int idx = hash_fns[i](flowkey) % width;
    ch->updateCnt(
        getIndex(i, idx),
        val * (static_cast<int>(hash_fns[depth + i](flowkey) & 1) * 2 - 1));
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSCountSketch<key_len, no_layer, T, hash_t>::query(
    const FlowKey<key_len> &flowkey) const {
  T values[depth];
  for (int i = 0; i < depth; ++i) {
    int idx = hash_fns[i](flowkey) % width;
    values[i] = ch->getCnt(getIndex(i, idx)) *
                (static_cast<int>(hash_fns[depth + i](flowkey) & 1) * 2 - 1);
  }
  std::sort(values, values + depth);
  if (!(depth & 1)) { // even
    return std::abs((values[depth / 2 - 1] + values[depth / 2]) / 2);
  } else { // odd
    return std::abs(values[depth / 2]);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSCountSketch<key_len, no_layer, T, hash_t>::size() const {
  ch->print_rate("");
  return sizeof(*this)                // instance
         + sizeof(hash_t) * depth * 2 // hashing class
         + ch->size();                // counter
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSCountSketch<key_len, no_layer, T, hash_t>::clear() {
  ch->clear();
}

} // namespace OmniSketch::Sketch
