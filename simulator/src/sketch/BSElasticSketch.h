/**
 * @file BSElasticSketch.h
 * @author XierLabber (you@domain.com)
 * @brief Implementation of BitSense-optimized Elastic Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

#include <sketch/BSCMSketch.h>

#define JUDGE_IF_SWAP(min_val, guard_val)                                      \
  ((guard_val) > ((min_val) << 3)) // delta = 8

namespace OmniSketch::Sketch {
/**
 * @brief BS Elastic Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
class BSElasticSketch : public SketchBase<key_len, T> {
private:
  struct Entry {
    FlowKey<key_len> flowkey_;
    bool flag_;
    Entry() : flowkey_(), flag_(false) {}
    bool isEmpty() {
      FlowKey<key_len> emptyKey{};
      return flowkey_ == emptyKey;
    }
  };
  // heavy part
  Entry **buckets_;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  int32_t heavy_cm_r, heavy_cm_w;

  BitSense<no_layer, T, hash_t> *ch;

  hash_t *hash_fns;

  int32_t num_buckets_;
  int32_t num_per_bucket_; // each bucket has num_per_bucket_ entries

  hash_t hash_h_;
  // light part
  BSCMSketch<key_len, no_layer, T, hash_t> cm_;

public:
  BSElasticSketch(int32_t num_buckets, int32_t num_per_bucket, int32_t l_depth,
                int32_t l_width, double cnt_no_ratio,
                const std::vector<size_t> &width_cnt,
                const std::vector<size_t> &no_hash, 
                const int32_t heavy_cm_r,
                const int32_t heavy_cm_w,
                double cm_cnt_no_ratio, 
                const std::vector<size_t> &cm_width_cnt, 
                const std::vector<size_t> &cm_no_hash);
  ~BSElasticSketch();

  int heavypartInsert(const FlowKey<key_len> &flowkey, T val,
                      FlowKey<key_len> &swap_key, T &swap_val);

  void lightpartInsert(const FlowKey<key_len> &flowkey, T val);
  void update(const FlowKey<key_len> &flowkey, T val);
  T heavypartQuery(const FlowKey<key_len> &flowkey, bool &flag) const;
  T lightpartQuery(const FlowKey<key_len> &flowkey) const;
  T query(const FlowKey<key_len> &flowkey) const;
  size_t size() const;
  void clear();
  int32_t getBSIdx(int32_t i, int32_t j) const{
    return i * num_per_bucket_ + j;
  }
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSElasticSketch<key_len, no_layer, T, hash_t>::BSElasticSketch(int32_t num_buckets,
                                                 int32_t num_per_bucket,
                                                 int32_t l_depth,
                                                 int32_t l_width, double cnt_no_ratio,
                                                 const std::vector<size_t> &width_cnt,
                                                 const std::vector<size_t> &no_hash, 
                                                 const int32_t heavy_cm_r,
                                                 const int32_t heavy_cm_w,
                                                 double cm_cnt_no_ratio, 
                                                 const std::vector<size_t> &cm_width_cnt, 
                                                 const std::vector<size_t> &cm_no_hash)
    : num_buckets_(Util::NextPrime(num_buckets)),
      width_cnt(width_cnt), no_hash(no_hash), 
      num_per_bucket_(num_per_bucket), cm_(l_depth, l_width, cm_cnt_no_ratio, cm_width_cnt, cm_no_hash) {
  buckets_ = new Entry *[num_buckets_];
  buckets_[0] = new Entry[num_buckets_ * num_per_bucket_]();
  for (int i = 1; i < num_buckets_; ++i) {
    buckets_[i] = buckets_[i - 1] + num_per_bucket_;
  }

  // check ratio
  if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in BS should be in (0, 1), but got " +
                            std::to_string(cnt_no_ratio) + " instead.");
  }
  // prepare no_cnt
  no_cnt.push_back(static_cast<size_t>(this->num_buckets_) * this->num_per_bucket_);
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = no_cnt.back();
    no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
  }
  // BS
  ch = new BitSense<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                 this->no_hash, false, true, 
                                                 heavy_cm_r, heavy_cm_w);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSElasticSketch<key_len, no_layer, T, hash_t>::size() const {
  ch->print_rate("HEAVY PART");
  return sizeof(*this) 
         + (key_len + 0.125) * num_buckets_ * num_per_bucket_ 
         + ch->size()
         + cm_.size();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSElasticSketch<key_len, no_layer, T, hash_t>::clear() {
  ch->clear();
  cm_.clear();
  std::fill(buckets_[0], buckets_[0] + num_buckets_ * num_per_bucket_, 0);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSElasticSketch<key_len, no_layer, T, hash_t>::~BSElasticSketch() {
  delete[] buckets_[0];
  delete[] buckets_;
  delete[] ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
int BSElasticSketch<key_len, no_layer, T, hash_t>::heavypartInsert(
    const FlowKey<key_len> &flowkey, T val, FlowKey<key_len> &swap_key,
    T &swap_val) {

  int32_t index = hash_h_(flowkey) % num_buckets_;
  int32_t matched = -1;
  int32_t empty = -1;
  int32_t min_counter = 0;
  T min_counter_val = ch->getEstCnt(getBSIdx(index, 0));
  for (int32_t i = 0; i < num_per_bucket_ - 1; ++i) {
    int32_t chIdx = getBSIdx(index, i);
    if (buckets_[index][i].flowkey_ == flowkey) { // flowkey is in the bucket
      matched = i;
      ch->updateCnt(chIdx, val);
      return 0;
    }
    if (buckets_[index][i].isEmpty()) {
      buckets_[index][i].flowkey_ = flowkey;
      ch->updateCnt(chIdx, val);
      return 0;
    }
    if (min_counter_val > ch->getEstCnt(chIdx)) {
      min_counter = i;
      min_counter_val = ch->getEstCnt(chIdx);
    }
  }
  /* update guard val and comparison */
  int32_t chIdx = getBSIdx(index, num_per_bucket_ - 1);
  T guard_val = ch->getEstCnt(chIdx);
  guard_val += 1;

  if (!JUDGE_IF_SWAP(min_counter_val, guard_val)) {
    ch->resetCnt(chIdx, guard_val);
    return 2;
  } else {
    int32_t minBSIdx = getBSIdx(index, min_counter);
    swap_key = buckets_[index][min_counter].flowkey_;
    swap_val = ch->getEstCnt(minBSIdx);
    ch->resetCnt(chIdx, 0);
    buckets_[index][min_counter].flowkey_ = flowkey;
    ch->resetCnt(minBSIdx, val);
    buckets_[index][min_counter].flag_ = true;
    return 1;
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSElasticSketch<key_len, no_layer, T, hash_t>::lightpartInsert(
    const FlowKey<key_len> &flowkey, T val) {
  cm_.update(flowkey, val);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSElasticSketch<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {

  FlowKey<key_len> swap_key;
  T swap_val;
  int result = heavypartInsert(flowkey, val, swap_key, swap_val);
  switch (result) {
  case 0:
    return;
  case 1: {
    lightpartInsert(swap_key, swap_val);
    return;
  }
  case 2: {
    lightpartInsert(flowkey, val);
    return;
  }
  default:
    printf("error return value !\n");
    exit(1);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSElasticSketch<key_len, no_layer, T, hash_t>::heavypartQuery(
    const FlowKey<key_len> &flowkey, bool &flag) const {
  int index = hash_h_(flowkey) % num_buckets_;
  for (int i = 0; i < num_per_bucket_ - 1; ++i) {
    int32_t chIdx = getBSIdx(index, i);
    if (buckets_[index][i].flowkey_ == flowkey) {
      flag = buckets_[index][i].flag_;
      return ch->getCnt(chIdx);
    }
  }
  return 0;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSElasticSketch<key_len, no_layer, T, hash_t>::lightpartQuery(
    const FlowKey<key_len> &flowkey) const {
  return cm_.query(flowkey);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSElasticSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const {
  
  bool flag = false;
  T heavy_result = heavypartQuery(flowkey, flag);
  T light_result = 0;
  if (heavy_result == 0 || flag == true) {
    light_result = lightpartQuery(flowkey);
  }
  return heavy_result + light_result;
}

} // namespace OmniSketch::Sketch
