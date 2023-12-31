/**
 * @file BSCMSketch.h
 * @author dromniscience (dr_abc@pku.edu.cn)
 * @brief Implementation of Count Min Sketch with BitSense
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
 * @brief Count Min Sketch with BS
 *
 * @tparam key_len  length of flowkey
 * @tparam no_layer layer of BS
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class BSCMSketch : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t width;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  hash_t *hash_fns;
  BitSense<no_layer, T, hash_t> *ch;

  BSCMSketch(const BSCMSketch &) = delete;
  BSCMSketch(BSCMSketch &&) = delete;

public:
  /**
   * @brief Construct by specifying depth, width and ch parameters
   *
   * @param depth   depth of CM
   * @param width   width of CM
   * @param cnt_no_width  ratio of the number of counters in two adjacent layers
   * (should be in (0, 1))
   * @param width_cnt   Width of counters on each layer
   * @param no_hash     #hash between adjacent layers
   *
   */
  BSCMSketch(int32_t depth, int32_t width, double cnt_no_ratio,
             const std::vector<size_t> &width_cnt,
             const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~BSCMSketch();
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSCMSketch<key_len, no_layer, T, hash_t>::BSCMSketch(
    int32_t depth, int32_t width, double cnt_no_ratio,
    const std::vector<size_t> &width_cnt, const std::vector<size_t> &no_hash)
    : depth(depth), width(Util::NextPrime(width)), ch(nullptr),
      width_cnt(width_cnt), no_hash(no_hash) {

  hash_fns = new hash_t[this->depth];
  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in BS should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->depth) * this->width);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // BS
  ch = new BitSense<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, true);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSCMSketch<key_len, no_layer, T, hash_t>::~BSCMSketch() {
  delete[] hash_fns;
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSCMSketch<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  for (int32_t i = 0; i < depth; ++i) {
    int32_t index = hash_fns[i](flowkey) % width;
    ch->updateCnt(i * width + index, val);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSCMSketch<key_len, no_layer, T, hash_t>::query(
    const FlowKey<key_len> &flowkey) const {
  T min_val = std::numeric_limits<T>::max();
  for (int32_t i = 0; i < depth; ++i) {
    int32_t index = hash_fns[i](flowkey) % width;
    min_val = std::min(min_val, ch->getCnt(i * width + index));
  }
  return min_val;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSCMSketch<key_len, no_layer, T, hash_t>::size() const {
  // ch->print_rate("COUNT MIN BS");
  return sizeof(*this)            // instance
         + depth * sizeof(hash_t) // hashing class
         + ch->size();            // ch
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSCMSketch<key_len, no_layer, T, hash_t>::clear() {
  ch->clear();
}

} // namespace OmniSketch::Sketch
