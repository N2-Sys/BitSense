/**
 * @file BSFlowRadar.h
 * @author XierLabber (you@domain.com)
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/bitsense.h>
#include <common/sketch.h>
#include <sketch/BloomFilter.h>

namespace OmniSketch::Sketch {
/**
 * @brief Flow Radar
 *
 * @tparam key_len  length of flowkey
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class BSFlowRadar : public SketchBase<key_len, T> {
private:
  FlowKey<key_len>* flowXOR;

  std::vector<size_t> flow_no_cnt;
  std::vector<size_t> flow_width_cnt;
  std::vector<size_t> flow_no_hash;
  BitSense<no_layer, T, hash_t>* flow_ch;

  std::vector<size_t> packet_no_cnt;
  std::vector<size_t> packet_width_cnt;
  std::vector<size_t> packet_no_hash;  
  BitSense<no_layer, T, hash_t>* packet_ch;

  const int32_t num_bitmap;
  const int32_t num_bit_hash;
  const int32_t num_count_table;
  const int32_t num_count_hash;
  int32_t num_flows;

  hash_t *hash_fns;
  BloomFilter<key_len, hash_t> *flow_filter;

  BSFlowRadar(const BSFlowRadar &) = delete;
  BSFlowRadar(BSFlowRadar &&) = delete;

  struct CountTableEntry {
    FlowKey<key_len> flowXOR;
    T flow_count;
    T packet_count;
    CountTableEntry() : flowXOR(), flow_count(0), packet_count(0) {}
  };

public:
  /**
   * @brief Construct a new Flow Radar object
   *
   * @param flow_filter_size Number of bits in flow filter (a Bloom Filter)
   * @param flow_filter_hash Number of hash functions in flow filter
   * @param count_table_size Number of elements in count table
   * @param count_table_hash Number of hash functions in count table
   */
  BSFlowRadar(int32_t flow_filter_size, int32_t flow_filter_hash,
            int32_t count_table_size, int32_t count_table_hash, 
             double flow_cnt_no_ratio,
             const std::vector<size_t> &flow_width_cnt,
             const std::vector<size_t> &flow_no_hash, 
             double packet_cnt_no_ratio,
             const std::vector<size_t> &packet_width_cnt,
             const std::vector<size_t> &packet_no_hash);
  /**
   * @brief Destructor
   *
   */
  ~BSFlowRadar();
  /**
   * @brief Update a flowkey with a certain value
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Decode flowkey and its value
   *
   */
  Data::Estimation<key_len, T> decode() override;
  /**
   * @brief Reset the sketch
   *
   */
  void clear();
  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const override;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSFlowRadar<key_len, no_layer, T, hash_t>::BSFlowRadar(
             int32_t flow_filter_size, int32_t flow_filter_hash,
             int32_t count_table_size, int32_t count_table_hash, 
             double flow_cnt_no_ratio,
             const std::vector<size_t> &flow_width_cnt,
             const std::vector<size_t> &flow_no_hash, 
             double packet_cnt_no_ratio,
             const std::vector<size_t> &packet_width_cnt,
             const std::vector<size_t> &packet_no_hash)
    : num_bitmap(Util::NextPrime(flow_filter_size)),
      num_bit_hash(flow_filter_hash),
      num_count_table(Util::NextPrime(count_table_size)),
      num_count_hash(count_table_hash), num_flows(0),
      flow_width_cnt(flow_width_cnt), flow_no_cnt(), 
      packet_width_cnt(packet_width_cnt), packet_no_cnt() {
  // check ratio
  if (flow_cnt_no_ratio <= 0.0 || flow_cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in BS should be in (0, 1), but got " +
                            std::to_string(flow_cnt_no_ratio) + " instead.");
  }
  if (packet_cnt_no_ratio <= 0.0 || packet_cnt_no_ratio >= 1.0) {
    throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                            "layers in BS should be in (0, 1), but got " +
                            std::to_string(packet_cnt_no_ratio) + " instead.");
  }
  hash_fns = new hash_t[num_count_hash];
  // flow filter
  flow_filter = new BloomFilter<key_len, hash_t>(num_bitmap, num_bit_hash);
  // count table
  flowXOR = new FlowKey<key_len>[num_count_table];
  // prepare no_cnt
  flow_no_cnt.push_back(static_cast<size_t>(num_count_table));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = flow_no_cnt.back();
    flow_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * flow_cnt_no_ratio)));
  }
  packet_no_cnt.push_back(static_cast<size_t>(num_count_table));
  for (int32_t i = 1; i < no_layer; ++i) {
    size_t last_layer = packet_no_cnt.back();
    packet_no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * packet_cnt_no_ratio)));
  }
  // BS
  flow_ch = new BitSense<no_layer, T, hash_t>(flow_no_cnt, 
                                                      flow_width_cnt,
                                                      flow_no_hash);
  packet_ch = new BitSense<no_layer, T, hash_t>(packet_no_cnt, 
                                                        packet_width_cnt,
                                                        packet_no_hash);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSFlowRadar<key_len, no_layer, T, hash_t>::~BSFlowRadar() {
  delete[] hash_fns;
  delete flow_filter;
  delete[] flowXOR;
  delete[] flow_ch;
  delete[] packet_ch;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSFlowRadar<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                           T val) {
  bool exist = flow_filter->lookup(flowkey);
  // a new flow
  if (!exist) {
    flow_filter->insert(flowkey);
    num_flows++;
  }

  for (int32_t i = 0; i < num_count_hash; i++) {
    int32_t index = hash_fns[i](flowkey) % num_count_table;
    // a new flow
    if (!exist) {
      flow_ch->updateCnt(index, 1);
      flowXOR[index] ^= flowkey;
    }
    // increment packet count
    packet_ch->updateCnt(index, val);
  }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
Data::Estimation<key_len, T> BSFlowRadar<key_len, no_layer, T, hash_t>::decode() {

  CountTableEntry* count_table = new CountTableEntry[num_count_table];
  for(int i = 0; i < num_count_table; i++)
  {
    count_table[i].flowXOR = flowXOR[i];
    count_table[i].flow_count = flow_ch->getCnt(i);
    count_table[i].packet_count = packet_ch->getCnt(i);
  }

  // an optimized implementation
  class CompareFlowCount {
  public:
    bool operator()(CountTableEntry *ptr1, CountTableEntry *ptr2) const {
      if (ptr1->flow_count == ptr2->flow_count) {
        return ptr1 < ptr2;
      } else
        return ptr1->flow_count < ptr2->flow_count;
    }
  };
  std::set<CountTableEntry *, CompareFlowCount> set;
  for (int32_t i = 0; i < num_count_table; ++i) {
    set.insert(count_table + i);
  }

  Data::Estimation<key_len, T> est;
  while (!set.empty()) {
    int32_t index = *set.begin() - count_table;
    T value = count_table[index].flow_count;
    // no decodable flow count
    if (value > 1) {
      break;
    }

    set.erase(set.begin());
    // ignore vacant counts
    if (value == 0)
      continue;

    FlowKey<key_len> flowkey = count_table[index].flowXOR;
    T size = count_table[index].packet_count;
    for (int i = 0; i < num_count_hash; ++i) {
      int l = hash_fns[i](flowkey) % num_count_table;
      set.erase(count_table + l);
      count_table[l].flow_count--;
      count_table[l].packet_count -= size;
      count_table[l].flowXOR ^= flowkey;
      set.insert(count_table + l);
    }
    est[flowkey] = size;
  }

  delete[] count_table;

  return est;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSFlowRadar<key_len, no_layer, T, hash_t>::size() const {
  flow_ch->print_rate("FLOW");
  packet_ch->print_rate("PACKET");
  return sizeof(*this)                                 // instance
         + num_count_hash * sizeof(hash_t)             // hashing class
         + num_count_table * key_len
         + flow_ch->size()
         + packet_ch->size()                           // ch
         + flow_filter->size();                        // flow filter
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSFlowRadar<key_len, no_layer, T, hash_t>::clear() {
  // reset flow counter
  num_flows = 0;
  // reset flow filter
  flow_filter->clear();
  // reset count table
  delete[] flowXOR;
  flowXOR = new FlowKey<key_len>[num_count_table];
  flow_ch->clear();
  packet_ch->clear();
}

} // namespace OmniSketch::Sketch