/**
 * @file BSPRSketch.h
 * @author XierLabber (you@domain.com)
 * @brief BSPR-Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/hash.h>
#include <common/bitsense.h>
#include <common/sketch.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include <map>

#define BYTE(n) ((n) >> 3)
#define BIT(n) ((n)&7)
#define FILTER_LENGTH(n) ((n + 7) >> 3)

#define TEST_DECODE_TIME

namespace OmniSketch::Sketch {
/**
 * @brief BSPR Sketch
 *
 * @tparam key_len  length of flowkey
 * @tparam no_layer numbers of layers of BS
 * @tparam T        type of the counter
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
class BSPRSketch : public SketchBase<key_len, T> {
private:
  int32_t counter_length;
  int32_t counter_hash_num;
  hash_t* counter_hash_func;

  std::vector<size_t> no_cnt;
  std::vector<size_t> width_cnt;
  std::vector<size_t> no_hash;
  int32_t ch_cm_r, ch_cm_w;
  BitSense<no_layer, T ,hash_t> *ch;

  int32_t filter_length;
  int32_t filter_hash_num;
  hash_t* filter_hash_func;
  uint8_t* filter;

  T phi;

  std::vector<FlowKey<key_len>> recorded_keys;
  std::map<FlowKey<key_len>, T> recovered_value;

  /**
   * @brief Set a bit
   *
   */
  void setBit(int32_t pos) { filter[BYTE(pos)] |= (1 << BIT(pos)); }
  /**
   * @brief Fetch a bit
   *
   * @return `true` if it is `1`; `false` otherwise.
   */
  bool getBit(int32_t pos) const { return (filter[BYTE(pos)] >> BIT(pos)) & 1; }
  /**
   * @brief Fetch a bit, and set that bit afer fetching
   *
   * @return `true` if it is `1`; `false` otherwise.
   */
  bool SetGetBit(int32_t pos) {
    uint8_t ans = filter[BYTE(pos)];
    filter[BYTE(pos)] = ans | (1 << BIT(pos));
    return (ans >> BIT(pos)) & 1;
  }

public:
  /**
   * @brief Construct by specifying counter_length, 
   * counter_hash_num, filter_length and filter_hash_num
   *
   */
  BSPRSketch(int32_t counter_length, int32_t counter_hash_num, 
           int32_t filter_length, int32_t filter_hash_num, 
           double cnt_no_ratio,
           const std::vector<size_t> &width_cnt,
           const std::vector<size_t> &no_hash, 
           const int32_t ch_cm_r, const int32_t ch_cm_w,
           T phi = 10);
  /**
   * @brief Release the pointer
   *
   */
  ~BSPRSketch();
  /**
   * @brief Update a flowkey with certain value
   *        val is a useless parameter here, as we only count
   *        the number of times flowkey occurs
   *
   */
  void update(const FlowKey<key_len> &flowkey, T val) override;
  /**
   * @brief Recover the values of every recorded keys
   *
   */
  void recover();
  /**
   * @brief Get the size of the sketch
   *
   */
  size_t size() const override;
  /**
   * @brief Query a flowkey
   *
   */
  T query(const FlowKey<key_len> &flowkey) const override;
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
BSPRSketch<key_len, no_layer, T, hash_t>::BSPRSketch(int32_t counter_length, 
    int32_t counter_hash_num, int32_t filter_length, 
    int32_t filter_hash_num, 
    double cnt_no_ratio,
    const std::vector<size_t> &width_cnt,
    const std::vector<size_t> &no_hash, 
    const int32_t ch_cm_r, const int32_t ch_cm_w, T phi): counter_length(Util::NextPrime(counter_length)),
    counter_hash_num(counter_hash_num), filter_length(Util::NextPrime(filter_length)),
    filter_hash_num(filter_hash_num), phi(phi),
    width_cnt(width_cnt), no_hash(no_hash) {  
    
    counter_hash_func = new hash_t[this->counter_hash_num];
    filter_hash_func = new hash_t[this->filter_hash_num];

    // check ratio
    if (cnt_no_ratio <= 0.0 || cnt_no_ratio >= 1.0) {
      throw std::out_of_range("Out of Range: Ratio of #counters of adjacent "
                              "layers in BS should be in (0, 1), but got " +
                              std::to_string(cnt_no_ratio) + " instead.");
    }
    // prepare no_cnt
    no_cnt.push_back(static_cast<size_t>(this->counter_length));
    for (int32_t i = 1; i < no_layer; ++i) {
      size_t last_layer = no_cnt.back();
      no_cnt.push_back(Util::NextPrime(std::ceil(last_layer * cnt_no_ratio)));
    }
    // BS
    ch = new BitSense<no_layer, T, hash_t>(no_cnt, this->width_cnt,
                                                   this->no_hash, false, true, 
                                                   ch_cm_r, ch_cm_w);

    filter = new uint8_t[FILTER_LENGTH(this->filter_length)]();

    recorded_keys.clear();
    recovered_value.clear();
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
BSPRSketch<key_len, no_layer, T, hash_t>::~BSPRSketch(){
    delete[] counter_hash_func;
    delete[] filter_hash_func;
    delete[] ch;
    delete[] filter;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSPRSketch<key_len, no_layer, T, hash_t>::update(const FlowKey<key_len> &flowkey,
                                          T val) {
    bool is_first_item = false;
    bool is_new_item = false;

    for(int i = 0; i < counter_hash_num; i++){
        size_t idx = counter_hash_func[i](flowkey) % counter_length;
        if(ch->getEstCnt(idx) <= phi){
            is_first_item = true;
        }
        ch->updateCnt(idx, val);
    }

    if(is_first_item){
        for(int i = 0; i < filter_hash_num; i++){
            size_t idx = filter_hash_func[i](flowkey) % filter_length;
            if(!getBit(idx)){
                is_new_item = true;
                setBit(idx);
            }
        }
    }

    if(is_new_item){
        recorded_keys.push_back(flowkey);
    }
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSPRSketch<key_len, no_layer, T, hash_t>::recover(){

#ifdef TEST_DECODE_TIME
  auto MY_TIMER = std::chrono::microseconds::zero();                              \
  auto MY_TICK = std::chrono::steady_clock::now();                                \
  auto MY_TOCK = std::chrono::steady_clock::now();
#endif

  int32_t key_num = recorded_keys.size();

  // solver
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
      solver_sparse;

  std::vector<double> recorded_b(counter_length);
  for(int i = 0; i < counter_length; i++)
  {
    recorded_b[i] = (double)ch->getCnt(i);
  }
  Eigen::VectorXd X(key_num),
      b = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(recorded_b.data(),
                                                        recorded_b.size());
  Eigen::SparseMatrix<double> A(counter_length, key_num);
  std::vector<Eigen::Triplet<double>> tripletlist;
  for(int i = 0; i < key_num; i++)
  {
    FlowKey<key_len> the_key = recorded_keys[i];
    for(int j = 0; j < counter_hash_num; j++)
    {
        size_t idx = counter_hash_func[j](the_key) % counter_length;
        tripletlist.push_back(Eigen::Triplet<double>(idx, i, 1.0));
    }
  }
  // duplicates are summed up, see
  // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0
  A.setFromTriplets(tripletlist.begin(), tripletlist.end());
  A.makeCompressed();
  solver_sparse.compute(A);
  X = solver_sparse.solve(b);
  recovered_value.clear();
  for (size_t i = 0; i < key_num; i++) {
    if(X[i] > 0)
    {
      recovered_value[recorded_keys[i]] = (static_cast<T>(X[i] + 0.5));
    }
    else
    {
      recovered_value[recorded_keys[i]] = 0;
    }
  }

#ifdef TEST_DECODE_TIME
  MY_TOCK = std::chrono::steady_clock::now();
  MY_TIMER = std::chrono::duration_cast<std::chrono::microseconds>(MY_TOCK - MY_TICK);
  printf("\nDECODE COST %jdms FOR ALL\n", static_cast<intmax_t>(MY_TIMER.count()));
#endif

}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
size_t BSPRSketch<key_len, no_layer, T, hash_t>::size() const{
    ch->print_rate("");
    return ch->size()
           + (filter_length >> 3)
           + (counter_hash_num + filter_hash_num) * sizeof(hash_t)
           + sizeof(*this);
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
T BSPRSketch<key_len, no_layer, T, hash_t>::query(const FlowKey<key_len> &flowkey) const{
    if(recovered_value.empty()){
        const_cast<BSPRSketch<key_len, no_layer, T, hash_t>* >(this)->recover();
    }
    auto cit = recovered_value.find(flowkey);
    if (cit != recovered_value.end()) {
        return cit->second;
    }
    return 0;
}

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void BSPRSketch<key_len, no_layer, T, hash_t>::clear(){
    ch->clear();
    std::fill(filter, filter + FILTER_LENGTH(filter_length), 0);
    recorded_keys.clear();
    recovered_value.clear();
}

} // namespace OmniSketch::Sketch