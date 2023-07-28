/**
 * @file BSUnivMon.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of BSUnivMon
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once
#include <iostream>

#include <common/hash.h>
#include <common/bitsense.h>
#include <common/sketch.h>
#include <sketch/CountSketch.h>

namespace OmniSketch::Sketch {
/**
 * @brief BSUnivMon
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
class BSUnivMon : public SketchBase<key_len, T> {
private:
  int32_t depth;
  int32_t total_length;
  int32_t *width;
  int32_t *width_idx;
  int32_t logn;
  hash_t *global_hash_fns;
  hash_t **CS_hash_fns;
  hash_t **CS_update_hash_fns;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;

  BitSense<no_layer, T, hash_t> *ch;

  Data::Estimation<key_len> *flows;

  BSUnivMon(const BSUnivMon &) = delete;
  BSUnivMon(BSUnivMon &&) = delete;

public:
  /**
   * @brief Construct by specifying depth, width and $\log n$, where $n$ is the
   * number of flows to insert.
   *
   */
  BSUnivMon(int32_t depth_, int32_t width_, int32_t log_n, double cnt_no_ratio,
            const std::vector<size_t> &width_cnt,
            const std::vector<size_t> &no_hash);
  /**
   * @brief Release the pointer
   *
   */
  ~BSUnivMon();
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
  int32_t getGlobalHash(const FlowKey<key_len> &flowkey, int32_t layer) const;
  int32_t getHashIdx(int32_t sketch_idx, int32_t hash_idx) const;
  int32_t getBSIdx(int32_t sketch_id, int32_t r, int32_t c) const;
  void updateSketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey, T val);
  T querySketch(int32_t sketch_idx, const FlowKey<key_len> &flowkey) const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSUnivMon<key_len, no_layer, T, hash_t>::querySketch(
    int32_t sketch_idx, const FlowKey<key_len> &flowkey) const {
  T values[depth];
  for (int i = 0; i < depth; ++i) {
    int32_t idx = CS_hash_fns[sketch_idx][i](flowkey) % width[sketch_idx];
    int32_t chIdx = getBSIdx(sketch_idx, i, idx);
    values[i] =
        ch->getCnt(chIdx) *
        (static_cast<int>(CS_update_hash_fns[sketch_idx][i](flowkey) & 1) * 2 -
         1);
  }
  std::sort(values, values + depth);
  if (!(depth & 1)) { // even
    return std::abs((values[depth / 2 - 1] + values[depth / 2]) / 2);
  } else { // odd
    return std::abs(values[depth / 2]);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSUnivMon<key_len, no_layer, T, hash_t>::updateSketch(
    int32_t sketch_idx, const FlowKey<key_len> &flowkey, T val) {
  for (int i = 0; i < depth; i++) {
    int32_t idx = CS_hash_fns[sketch_idx][i](flowkey) % width[sketch_idx];
    T update_val =
        val *
        (static_cast<int>(CS_update_hash_fns[sketch_idx][i](flowkey) & 1) * 2 -
         1);
    ch->updateCnt(getBSIdx(sketch_idx, i, idx), update_val);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t BSUnivMon<key_len, no_layer, T, hash_t>::getBSIdx(int32_t sketch_id,
                                                          int32_t r,
                                                          int32_t c) const {
  return (width_idx[sketch_id] + c) + r * total_length;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int32_t BSUnivMon<key_len, no_layer, T, hash_t>::getGlobalHash(
    const FlowKey<key_len> &flowkey, int32_t layer) const {
  if (layer == 0) {
    return 1;
  }
  return (global_hash_fns[layer - 1](flowkey) & 1);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSUnivMon<key_len, no_layer, T, hash_t>::BSUnivMon(
    int32_t depth_, int32_t width_, int32_t log_n, double cnt_no_ratio,
    const std::vector<size_t> &width_cnt, const std::vector<size_t> &no_hash)
    : depth(depth_), logn(log_n), width_cnt(width_cnt), no_hash(no_hash) {

  global_hash_fns = new hash_t[logn - 1];
  width_ = (Util::NextPrime(width_));

  CS_hash_fns = new hash_t *[logn];
  CS_update_hash_fns = new hash_t *[logn];

  int32_t running_idx = 0;
  width = new int32_t[logn];
  width_idx = new int32_t[logn];

  for (int i = 0; i < logn; i++) {
    width[i] = Util::NextPrime(width_);
    width_idx[i] = running_idx;
    running_idx += width[i];
    width_ = std::max(1, width_ / 2);
    CS_hash_fns[i] = new hash_t[depth];
    CS_update_hash_fns[i] = new hash_t[depth];
  }

  total_length = running_idx;

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in BS should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->depth) * running_idx);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // BS
  ch = new BitSense<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, true);

  flows = new Data::Estimation<key_len>[logn];
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSUnivMon<key_len, no_layer, T, hash_t>::~BSUnivMon() {
  if (global_hash_fns)
    delete[] global_hash_fns;
  for (int i = 0; i < depth; i++) {
    delete[] CS_hash_fns[i];
    delete[] CS_update_hash_fns[i];
  }
  delete[] CS_hash_fns;
  delete[] CS_update_hash_fns;
  delete[] width;
  delete[] width_idx;
  delete[] ch;
  if (flows)
    delete[] flows;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSUnivMon<key_len, no_layer, T, hash_t>::update(
    const FlowKey<key_len> &flowkey, T val) {
  for (int32_t i = 0; i < logn; ++i) {
    if (getGlobalHash(flowkey, i)) {
      updateSketch(i, flowkey, val);
    } else
      break;
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSUnivMon<key_len, no_layer, T, hash_t>::query(
    const FlowKey<key_len> &flowkey) const {
  static bool cnt_distrib = true;
  if (cnt_distrib) {
    std::vector<double> distrib(32);
    for (int32_t i = 0; i < no_cnt[0]; ++i) {
      T val = ch->getOriginalCnt(i);
      for (int32_t k = 0; k < 32; ++k) {
        if (std::abs(val) >= (1 << k))
          distrib[k] += 1.0;
        else
          break;
      }
    }
    for (int32_t k = 0; k < 32; ++k) {
      std::cout << distrib[k] / no_cnt[0] << " ";
      if (distrib[k] == 0.0) {
        std::cout << std::endl;
        break;
      }
    }
    cnt_distrib = false;
  }
  int level;
  for (level = 0; level < logn; level++) {
    if (!getGlobalHash(flowkey, level)) {
      break;
    }
  }
  level--;
  T ret = querySketch(level, flowkey);
  for (int i = level - 1; i >= 0; i--) {
    ret = 2 * ret - querySketch(i, flowkey);
  }
  return ret;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSUnivMon<key_len, no_layer, T, hash_t>::size() const {
  ch->print_rate("UnivMon BS");
  return sizeof(*this) + sizeof(hash_t) * (logn + 2 * logn * depth - 1) +
         sizeof(int32_t) * 2 * logn + ch->size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSUnivMon<key_len, no_layer, T, hash_t>::clear() {
  // std::fill(counter[0], counter[0] + depth * width, 0);
}

} // namespace OmniSketch::Sketch
