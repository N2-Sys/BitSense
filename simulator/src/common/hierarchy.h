/**
 * @file hierarchy.h
 * @author dromniscience (you@domain.com)
 * @brief Counter Hierarchy
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include "hash.h"
#include "utils.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include <map>

#define RECORD_ACCESS_TIME
#define TEST_DECODE_TIME
// #define SKIP_HASH
// #define USE_CLP

namespace OmniSketch::Sketch {
/**
 * @brief Use the counter hierarchy to better save space while preserving
 * accuracy!
 *
 * @tparam no_layer   Number of layers in CH
 * @tparam T          Counter Type of CH
 * @tparam hash_t     Hashing classes used internally
 *
 * @note
 * - In CH, counters are serialized, so it is the user's job to convert the
 * index of a possibly multi-dimensional array into a unique serial number.
 * Since CH uses 0-based array internally, this serial number is chosen to
 * have type `size_t`.
 * - CH uses **lazy update policy**. That is, only when you try to get a counter
 * value in CH will the updates genuinely be propagated to the higher layers
 * and then decoded.
 * - Note that CH cannot bring better accuracy. Ideally, if the first layer in
 * CH never overflows, your sketch achieves the same accuracy as does without
 * it.
 *
 * @warning
 * - To apply CH, make sure values of the original counters are always
 * non-negative during the whole process, though negative update is OK.
 * (i.e., only *cash register case* and the *non-negative case*
 * in *turnstile model* are allowed) Otherwise, errors of decoding can be
 * unbounded.
 * - It is highly recommended that on each layer the number of counters is
 * prime.
 */
template <int32_t no_layer, typename T, typename hash_t = Hash::AwareHash>
class CounterHierarchy {
private:
  using CarryOver = std::map<std::size_t, T>;
  /**
   * @brief Number of counters on each layer, from low to high.
   *
   */
  const std::vector<size_t> no_cnt;
  /**
   * @brief Width of counters on each layer, from low to high.
   *
   */
  const std::vector<size_t> width_cnt;
  /**
   * @brief Number of hash function used on each layer, from low to
   * high (except for the last layer)
   *
   */
  const std::vector<size_t> no_hash;
#ifndef SKIP_HASH
  /**
   * @brief array of hashing classes
   *
   */
  std::vector<hash_t> *hash_fns;
#endif
  /**
   * @brief counters in CH
   *
   */
  std::vector<Util::DynamicIntX<T>> *cnt_array;
  /**
   * @brief Status bits
   *
   */
  boost::dynamic_bitset<uint8_t> *status_bits;
  /**
   * @brief Original counters
   *
   */
  std::vector<T> original_cnt;
  /**
   * @brief get decoded counter
   * @details `double` will round to `T` after decoding each layer. The reason
   * why `double` here is to facilitate NZE decoding.
   */
  std::vector<double> decoded_cnt;
#ifndef RECORD_ACCESS_TIME
  /**
   * @brief For lazy update policy
   *
   */
  CarryOver lazy_update;
#else
  /**
   * @brief Times of Accessing DynamicIntX
   *
   */
  int32_t access_time;
  /**
   * @brief Total Update Times
   *
   */
  int32_t total_update_time;

  double est_ARE;
  int32_t est_ERROR_TIME;
  int32_t est_TIMES;

#endif
  /**
   * @brief A time-saving optimization
   * @details If no upper layer counter is set, there is no need to decode.
   *
   */
  bool need_to_decode;
  bool use_negative_counters;
  bool have_decoded;

  bool use_cm_sketch;

  size_t cm_row, cm_width;

  std::vector<Util::DynamicIntX<T>> *cm_sketch;
  std::vector<hash_t> cm_hash;

private:
  /**
   * @brief Update a layer (aggregation)
   *
   * @details An overflow error would be thrown if there is an overflow at the
   * last layer. For the other layers, since a counter may first oveflow and
   * then be substracted to withdraw any carry over, the number of overflows of
   * counter whose status bit is set is not assumed to be 1 at least.
   *
   * @param layer   the current layer
   * @param updates updates to be propagated to the current layer
   * @return updates to be propagated to the next layer
   */
  void updateSegment(const int32_t layer, const size_t index, const T val);
  /**
   * @brief Update a layer (aggregation)
   *
   * @details An overflow error would be thrown if there is an overflow at the
   * last layer. For the other layers, since a counter may first oveflow and
   * then be substracted to withdraw any carry over, the number of overflows of
   * counter whose status bit is set is not assumed to be 1 at least.
   *
   * @param layer   the current layer
   * @param updates updates to be propagated to the current layer
   * @return updates to be propagated to the next layer
   */
  [[nodiscard]] CarryOver updateLayer(const int32_t layer, CarryOver &&updates);
  /**
   * @brief Decode a layer
   *
   * @param layer   the current layer to decode
   * @param higher  decoded results of the higher layer
   * @return results of the current layer
   */
  [[nodiscard]] std::vector<double>
  decodeLayer(const int32_t layer, std::vector<double> &&higher) const;

public:
  /**
   * @brief Construct by specifying detailed architectural parameters
   *
   * @param no_cnt      number of counters on each layer, from low to high
   * @param width_cnt   width of counters on each layer, from low to high
   * @param no_hash     number of hash functions used on each layer, from low
   * to high (except for the last layer)
   *
   * @details The meaning of the three parameters stipulates the following
   * requirements:
   * - size of `no_cnt` should equal `no_layer`.
   * - size of `width_cnt` should equal `no_layer`.
   * - size of `no_hash` should equal `no_layer - 1`.
   *
   * If any of these is violated, an exception would be thrown. Other
   * conditions that trigger an exception:
   * - Items in these vectors contains a 0.
   * - `no_layer <= 0`
   * - Sum of `width_cnt` exceeds `sizeof(T) * 8`. This constraint is imposed to
   * guarantee proper shifting of counters when decoding.
   */
  CounterHierarchy(const std::vector<size_t> &no_cnt,
                   const std::vector<size_t> &width_cnt,
                   const std::vector<size_t> &no_hash,
                   bool use_negative_counters = false,
                   bool use_cm_sketch = false, size_t cm_row = -1,
                   size_t cm_width = -1);
  /**
   * @brief Destructor
   *
   */
  ~CounterHierarchy();
  /**
   * @brief Update a counter
   *
   * @param index Serialized index of a counter. It is the user's job to get the
   * index serialized in advance.
   * @param val   Value to be updated
   *
   * @note
   * - Use the lazy update policy.
   * - An out-of-range exception would be thrown if `index` is out of range.
   */
  void updateCnt(size_t index, T val);
  /**
   * @brief Get the value of counters in CH
   *
   * @details An overflow exception would be thrown if there is an overflow at
   * the last layer. An out-of-range exception would be thrown if the index is
   * out of range.
   *
   * @param index Serialized index of a counter. It is the user's job to get the
   * index serialized in advance.
   */
  T getCnt(size_t index);
  /**
   * @brief Get the original value of counters.
   *
   * @details I.e., the value of the counter without CH.
   *
   * @param index Serialized index of a counter. It is the user's job to get
   * the index serialized in advance.
   */
  T getOriginalCnt(size_t index) const;
  /**
   * @brief Size of CH.
   *
   */
  size_t size() const;
  /**
   * @brief Size of counters without CH.
   *
   */
  size_t originalSize() const;
  /**
   * @brief Reset CH.
   *
   */
  void clear();
  /**
   * @brief Get the value of status bits in CH
   *
   */
  uint8_t getStatus(int32_t idx) const;
  /**
   * @brief Output Overflow rates and correct rates
   *
   */
  void print_rate(const char *name);
  void highest_bit_add(int32_t val);
  /**
   * @brief return -1 if we think counter[layer][idx] < 0, return 1
   *        otherwise
   *
   */
#ifdef SKIP_HASH
  int32_t my_hash(int32_t layer, int32_t hash_id, int32_t counter_id) const;
#endif
  int guess_sign(int32_t negative_weight, int32_t positive_weight, int32_t idx,
                 int32_t layer = 0);
  T get_current_cnt(size_t idx);
  T getEstCnt(int32_t idx, int32_t layer = 0);
  void resetCnt(size_t index, T val);
  T getTotalCnt(int32_t idx, int32_t layer = 0);
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

#ifdef SKIP_HASH
template <int32_t no_layer, typename T, typename hash_t>
int32_t
CounterHierarchy<no_layer, T, hash_t>::my_hash(int32_t layer, int32_t hash_id,
                                               int32_t counter_id) const {
  int32_t offset = (6 * hash_id + 6);
  int32_t mask = (1 << offset) - 1;
  int32_t tmp = counter_id & mask;
  counter_id = (counter_id >> offset) ^ (tmp << (32 - offset));
  return counter_id % no_cnt[layer + 1];
}
#endif

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::updateSegment(const int32_t layer,
                                                          const size_t index,
                                                          const T val) {
  access_time++;
  T overflow = cnt_array[layer][index] + val;
  if (overflow) {
    need_to_decode = true;
    // mark status bits
    status_bits[layer][index] = true;
    if (use_cm_sketch && layer == 0) {
      for (size_t i = 0; i < cm_row; ++i) {
        size_t ind = cm_hash[i](index) % cm_width;
        cm_sketch[i][ind].setVal(cm_sketch[i][ind].getVal() + overflow);
      }
    }
    if (layer == no_layer - 1) { // last layer
      throw std::overflow_error(
          "Counter overflow at the last layer in CH, overflow by " +
          std::to_string(overflow) + ".");
    } else { // hash to upper-layer counters
      for (size_t i = 0; i < no_hash[layer]; i++) {
#ifndef SKIP_HASH
        std::size_t new_index = hash_fns[layer][i](index) % no_cnt[layer + 1];
#else
        std::size_t new_index = my_hash(layer, i, index);
#endif
        updateSegment(layer + 1, new_index, overflow);
      }
    }
  }
}

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::highest_bit_add(int32_t val) {
  int32_t tmp = (val >= 0) ? 1 : -1;
  for (int t = 0; t < no_layer; t++) {
    for (int i = 0; i < no_cnt[t]; i++) {
      cnt_array[t][i] + (1 << (width_cnt[t] - 1));
    }
  }
}

template <int32_t no_layer, typename T, typename hash_t>
typename CounterHierarchy<no_layer, T, hash_t>::CarryOver
CounterHierarchy<no_layer, T, hash_t>::updateLayer(const int32_t layer,
                                                   CarryOver &&updates) {
  // A time-saving optimization
  if (layer > 0 && !updates.empty()) {
    need_to_decode = true;
  }

  CarryOver ret; // aggregate all updates on the current layer
  for (const auto &kv : updates) {
    T overflow = cnt_array[layer][kv.first] + kv.second;
    if (overflow) {
      // mark status bits
      status_bits[layer][kv.first] = true;
      if (layer == no_layer - 1) { // last layer
        throw std::overflow_error(
            "Counter overflow at the last layer in CH, overflow by " +
            std::to_string(overflow) + ".");
      } else { // hash to upper-layer counters
        for (size_t i = 0; i < no_hash[layer]; i++) {
#ifndef SKIP_HASH
          std::size_t index = hash_fns[layer][i](kv.first) % no_cnt[layer + 1];
#else
          std::size_t index = my_hash(layer, i, kv.first);
#endif
          ret[index] += overflow;
        }
      }
    }
  }
  return ret;
}

template <int32_t no_layer, typename T, typename hash_t>
std::vector<double> CounterHierarchy<no_layer, T, hash_t>::decodeLayer(
    const int32_t layer, std::vector<double> &&higher) const {
  // make sure the size match
  if (higher.size() != no_cnt[layer + 1]) {
    throw std::length_error("Size Error: Expect a vector of size " +
                            std::to_string(no_cnt[layer + 1]) +
                            ", but got one of size " +
                            std::to_string(higher.size()) + " instead.");
  }

#ifdef USE_CLP
  // solver
  OsiClpSolverInterface model;
  int numcols = no_cnt[layer];
  int numrows = no_cnt[layer + 1];
  double *obj = new double[numcols + 100];
  double *rowlb = new double[numrows + 100];
  double *rowub = new double[numrows + 100];
  double *collb = new double[numcols + 100];
  double *colub = new double[numcols + 100];

  int *my_start = new int[numcols + 1 + 100];
  int *my_index = new int[numcols * no_hash[layer] + 1 + 100];
  double *values = new double[numcols * no_hash[layer] + 100];

  for (int i = 0; i < numcols; i++) {
    obj[i] = -1.0;
    collb[i] = -__DBL_MAX__;
    colub[i] = __DBL_MAX__;
  }
  for (int i = 0; i < numrows; i++) {
    rowlb[i] = higher[i];
    rowub[i] = higher[i];
  }

  int elem_num = 0;
  for (int i = 0; i < numcols; i++) {
    my_start[i] = elem_num;
    if (!status_bits[layer][i]) {
      continue;
    }

    std::map<int, double> tmp;
    for (int j = 0; j < no_hash[layer]; j++) {
#ifndef SKIP_HASH
      int k = hash_fns[layer][j](i) % no_cnt[layer + 1];
#else
      int k = my_hash(layer, j, i);
#endif
      tmp[k] += 1.0;
    }
    for (auto item : tmp) {
      my_index[elem_num] = item.first;
      values[elem_num++] = item.second;
    }
  }
  my_start[numcols] = elem_num;

  model.loadProblem(numcols, numrows, my_start, my_index, values, collb, colub,
                    obj, rowlb, rowub);
  for (int i = 0; i < numcols; i++) {
    model.setInteger(i); // Sets xi to integer
  }
  model.setObjSense(-1.0); // Maximise

  CbcModel solver(model);
  solver.branchAndBound();
  const double *ans = solver.getColSolution();

  std::vector<double> ret(ans, ans + numcols);
  for (int i = 0; i < no_cnt[layer]; ++i) {
    ret[i] = status_bits[layer][i]
                 ? ((static_cast<int>(ret[i])) << width_cnt[layer])
                 : 0;
    ret[i] += cnt_array[layer][i].getVal();
  }

  delete[] obj;
  delete[] rowlb;
  delete[] rowub;
  delete[] collb;
  delete[] colub;
  delete[] my_start;
  delete[] my_index;
  delete[] values;
  printf("REACH HERE12!\n");

  return ret;

#else
  // solver
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
      solver_sparse;

  Eigen::VectorXd X(no_cnt[layer]),
      b = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(higher.data(),
                                                        higher.size());
  Eigen::SparseMatrix<double> A(no_cnt[layer + 1], no_cnt[layer]);
  std::vector<Eigen::Triplet<double>> tripletlist;

  for (size_t i = 0; i < no_cnt[layer]; i++) {
    if (!status_bits[layer][i])
      continue; // no overflow
    // hash to higher-layer counter
    for (size_t j = 0; j < no_hash[layer]; ++j) {
#ifndef SKIP_HASH
      size_t k = hash_fns[layer][j](i) % no_cnt[layer + 1];
#else
      size_t k = my_hash(layer, j, i);
#endif
      tripletlist.push_back(Eigen::Triplet<double>(k, i, 1.0));
    }
  }
  // duplicates are summed up, see
  // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0
  A.setFromTriplets(tripletlist.begin(), tripletlist.end());
  A.makeCompressed();
  solver_sparse.compute(A);
  X = solver_sparse.solve(b);

  std::vector<double> ret(no_cnt[layer]);
  for (size_t i = 0; i < no_cnt[layer]; ++i) {
    if (X[i] > 0) {
      ret[i] = status_bits[layer][i]
                   ? static_cast<double>(static_cast<T>(X[i] + 0.5)
                                         << width_cnt[layer])
                   : 0.0;
    } else {
      ret[i] = status_bits[layer][i]
                   ? static_cast<double>(static_cast<T>(X[i] - 0.5)
                                         << width_cnt[layer])
                   : 0.0;
    }
    ret[i] += cnt_array[layer][i].getVal();
  }
  return ret;
#endif
}

template <int32_t no_layer, typename T, typename hash_t>
CounterHierarchy<no_layer, T, hash_t>::CounterHierarchy(
    const std::vector<size_t> &no_cnt, const std::vector<size_t> &width_cnt,
    const std::vector<size_t> &no_hash, bool use_negative_counters,
    bool use_cm_sketch_, size_t cm_row_, size_t cm_width_)
    : no_cnt(no_cnt), width_cnt(width_cnt), no_hash(no_hash),
      use_negative_counters(use_negative_counters), need_to_decode(false),
      have_decoded(false), use_cm_sketch(use_cm_sketch_), cm_row(cm_row_),
      cm_width(use_cm_sketch_ ? Util::NextPrime(cm_width_) : -1),
      cm_sketch(NULL), est_TIMES(0), est_ARE(0), est_ERROR_TIME(0) {
  // validity check
  if (use_cm_sketch_) {
    if (cm_row_ <= 0) {
      throw std::invalid_argument(
          "Invalid Template Argument: `cm_row` must >= 1, got " +
          std::to_string(cm_row_) + ".");
    }
    if (cm_width_ <= 0) {
      throw std::invalid_argument(
          "Invalid Template Argument: `cm_width` must >= 1, got " +
          std::to_string(cm_width_) + ".");
    }
  }
  if (no_layer < 1) {
    throw std::invalid_argument(
        "Invalid Template Argument: `no_layer` must > 1, got " +
        std::to_string(no_layer) + ".");
  }
  if (no_cnt.size() != no_layer) {
    throw std::invalid_argument(
        "Invalid Argument: `no_cnt` should be of size " +
        std::to_string(no_layer) + ", but got size " +
        std::to_string(no_cnt.size()) + ".");
  }
  if (width_cnt.size() != no_layer) {
    throw std::invalid_argument(
        "Invalid Argument: `width_cnt` should be of size " +
        std::to_string(no_layer) + ", but got size " +
        std::to_string(width_cnt.size()) + ".");
  }
  if (no_hash.size() != no_layer - 1) {
    throw std::invalid_argument(
        "Invalid Argument: `no_hash` should be of size " +
        std::to_string(no_layer - 1) + ", but got size " +
        std::to_string(no_hash.size()) + ".");
  }
  for (auto i : no_cnt) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `no_cnt`.");
    }
  }
  for (auto i : width_cnt) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `width_cnt`.");
    }
  }
  for (auto i : no_hash) {
    if (i == 0) {
      throw std::invalid_argument(
          "Invalid Argument: There is a zero in `no_hash`.");
    }
  }
  size_t length = 0;
  for (auto i : width_cnt) {
    size_t tmp = length + i;
    if (tmp < length || tmp > sizeof(T) * 8) {
      throw std::invalid_argument(
          "Invalid Argument: Aggregate length of `width_cnt` is too large.");
    }
    length = tmp;
  }

  // allocate in heap
#ifndef SKIP_HASH
  hash_fns = new std::vector<hash_t>[no_layer - 1];
  for (int32_t i = 0; i < no_layer - 1; ++i) {
    hash_fns[i] = std::vector<hash_t>(no_hash[i]);
  }
#endif
  cnt_array = new std::vector<Util::DynamicIntX<T>>[no_layer];
  for (int32_t i = 0; i < no_layer; ++i) {
    cnt_array[i] = std::vector<Util::DynamicIntX<T>>(no_cnt[i], {width_cnt[i]});
  }
  status_bits = new boost::dynamic_bitset<uint8_t>[no_layer];
  for (int32_t i = 0; i < no_layer; ++i) {
    status_bits[i].resize(no_cnt[i], false);
  }
  // original counters, value initialized
  original_cnt.resize(no_cnt[0]);
  for (int i = 0; i < no_cnt[0]; i++) {
    original_cnt[i] = 0;
  }
  // decoded counters, value initialized
  decoded_cnt.resize(no_cnt[0]);
  // set counters to make them record negtive values
  if (use_negative_counters) {
    highest_bit_add(1);
  }
#ifdef RECORD_ACCESS_TIME
  access_time = 0;
  total_update_time = 0;
#endif
  if (use_cm_sketch) {
    int32_t length = 0;
    for (int i = 1; i < no_layer; i++) {
      length += width_cnt[i];
    }
    cm_sketch = new std::vector<Util::DynamicIntX<T>>[cm_row];
    for (int32_t i = 0; i < cm_row; ++i)
      cm_sketch[i].resize(cm_width, length);
    if (use_negative_counters) {
      int offset = (1 << (length - 1));
      for (int i = 0; i < cm_row; i++) {
        for (int j = 0; j < cm_width; j++) {
          cm_sketch[i][j] + offset;
        }
      }
    }
    cm_hash.resize(cm_row);
  }
}

template <int32_t no_layer, typename T, typename hash_t>
CounterHierarchy<no_layer, T, hash_t>::~CounterHierarchy() {
#ifndef SKIP_HASH
  if (hash_fns)
    delete[] hash_fns;
#endif
  if (cnt_array)
    delete[] cnt_array;
  if (status_bits)
    delete[] status_bits;
  if (use_cm_sketch) {
    delete[] cm_sketch;
  }
}

template <int32_t no_layer, typename T, typename hash_t>
uint8_t CounterHierarchy<no_layer, T, hash_t>::getStatus(int32_t idx) const {
  return status_bits[0].test(idx);
}

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::updateCnt(size_t index, T val) {
  if (index >= no_cnt[0]) {
    throw std::out_of_range("Index Out of Range: Should be in [0, " +
                            std::to_string(no_cnt[0] - 1) + "], but got " +
                            std::to_string(index) + " instead.");
  }
  have_decoded = false;
#ifndef RECORD_ACCESS_TIME
  // lazy update policy
  lazy_update[index] += val;
#else
  updateSegment(0, index, val);
  total_update_time++;
#endif
  // original counters
  original_cnt[index] += val;
}

template <int32_t no_layer, typename T, typename hash_t>
T CounterHierarchy<no_layer, T, hash_t>::getCnt(size_t index) {
  if (index >= no_cnt[0]) {
    throw std::out_of_range("Index Out of Range: Should be in [0, " +
                            std::to_string(no_cnt[0] - 1) + "], but got " +
                            std::to_string(index) + " instead.");
  }

#ifndef RECORD_ACCESS_TIME
  // lazy update
  if (!lazy_update.empty()) {
    for (int32_t i = 0; i < no_layer; i++) {
      lazy_update = updateLayer(i, std::move(lazy_update)); // throw exception
    }
    lazy_update.clear();
  }
#endif

  // A time-saving optimization
  if (have_decoded || !need_to_decode) {
    if (!have_decoded)
      return cnt_array[0][index].getVal();
    else
      return static_cast<T>(decoded_cnt[index]);
  } else { // decode
    printf("\nDECODER CALLED!\n");
#ifdef TEST_DECODE_TIME
    auto MY_TIMER = std::chrono::microseconds::zero();
    auto MY_TICK = std::chrono::steady_clock::now();
    auto MY_TOCK = std::chrono::steady_clock::now();
#endif
    decoded_cnt = std::vector<double>(no_cnt.back());
    for (size_t i = 0; i < no_cnt.back(); ++i) {
      if (!use_negative_counters) {
        decoded_cnt[i] =
            static_cast<double>(cnt_array[no_layer - 1][i].getVal());
      } else {
        decoded_cnt[i] =
            static_cast<double>(cnt_array[no_layer - 1][i].getVal()) -
            (1 << (width_cnt[no_layer - 1] - 1));
      }
    }
    for (int32_t i = no_layer - 2; i >= 0; i--) {
      decoded_cnt = decodeLayer(i, std::move(decoded_cnt));
      if (use_negative_counters) {
        for (int j = 0; j < no_cnt[i]; j++) {
          decoded_cnt[j] -= (1 << (width_cnt[i] - 1));
        }
      }
    }
#ifdef TEST_DECODE_TIME
    MY_TOCK = std::chrono::steady_clock::now();
    MY_TIMER = std::chrono::duration_cast<std::chrono::microseconds>(MY_TOCK -
                                                                     MY_TICK);
    printf("\nDECODE COST %lldms\n", static_cast<int64_t>(MY_TIMER.count()));
#endif
    printf("\nDECODER END!\n");
    have_decoded = true;
    return static_cast<T>(decoded_cnt[index]);
  }
  // always return

  /*
  // lazy update
  if (!lazy_update.empty()) {
    for (int32_t i = 0; i < no_layer; i++) {
      lazy_update = updateLayer(i, std::move(lazy_update)); // throw exception
    }
    lazy_update.clear();
    // A time-saving optimization
    if (!need_to_decode) {
      return cnt_array[0][index].getVal();
    } else { // decode
      printf("\nDECODER CALLED!\n");
      decoded_cnt = std::vector<double>(no_cnt.back());
      for (size_t i = 0; i < no_cnt.back(); ++i) {
        if(!use_negative_counters)
        {
          decoded_cnt[i] =
              static_cast<double>(cnt_array[no_layer - 1][i].getVal());
        }
        else
        {
          decoded_cnt[i] =
            static_cast<double>(cnt_array[no_layer - 1][i].getVal()) - (1 <<
  (width_cnt[no_layer - 1] - 1));
        }
      }
      for (int32_t i = no_layer - 2; i >= 0; i--) {
        decoded_cnt = decodeLayer(i, std::move(decoded_cnt));
        if(use_negative_counters)
        {
          for(int j = 0; j < no_cnt[i]; j++)
          {
            decoded_cnt[j] -= (1 << (width_cnt[i] - 1));
          }
        }
      }
      printf("\nDECODER END!\n");
      return static_cast<T>(decoded_cnt[index]);
    }
    // always return
  }
  // no lazy update
  if (!need_to_decode)
    return cnt_array[0][index].getVal();
  return static_cast<T>(decoded_cnt[index]);
  */
}

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::resetCnt(size_t index, T val) {
  T est = getEstCnt(index, 0);
  updateSegment(0, index, -est);
  total_update_time++;
  // sketch.erase(index);
  status_bits[0][index] = false;
  updateSegment(0, index, val);
  original_cnt[index] = val;
}

template <int32_t no_layer, typename T, typename hash_t>
T CounterHierarchy<no_layer, T, hash_t>::getOriginalCnt(size_t index) const {
  if (index >= no_cnt[0]) {
    throw std::out_of_range("Index Out of Range: Should be in [0, " +
                            std::to_string(no_cnt[0] - 1) + "], but got " +
                            std::to_string(index) + " instead.");
  }
  return original_cnt[index];
}

template <int32_t no_layer, typename T, typename hash_t>
size_t CounterHierarchy<no_layer, T, hash_t>::size() const {
  // counters + status bits
  size_t tot = 0; // first in bits
  for (int32_t i = 0; i < no_layer; ++i) {
    tot += no_cnt[i] * width_cnt[i];
    tot += no_cnt[i];
  }
  tot >>= 3; // to bytes
  // hash functions
#ifndef SKIP_HASH
  for (int32_t i = 0; i < no_layer - 1; ++i) {
    tot += sizeof(hash_t) * no_hash[i];
  }
#endif
  if (use_cm_sketch) {
    int32_t length = 0;
    for (int i = 1; i < no_layer; i++) {
      length += width_cnt[i];
    }
    tot += ((length * cm_width * cm_row) >> 3);
    tot += cm_row * sizeof(hash_t);
  }
  return tot;
}

template <int32_t no_layer, typename T, typename hash_t>
size_t CounterHierarchy<no_layer, T, hash_t>::originalSize() const {
  return sizeof(T) * no_cnt[0];
}

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::clear() {
  // reset counters
  for (int32_t i = 0; i < no_layer; ++i) {
    cnt_array[i] = std::vector<Util::DynamicIntX<T>>(no_cnt[i], width_cnt[i]);
  }
  // reset status bits
  for (int32_t i = 0; i < no_layer; ++i) {
    status_bits[i].reset();
  }
  // reset original counters
  original_cnt = std::vector<T>(no_cnt[0]);
  // // reset decoded counters
  // decoded_cnt = std::vector<double>(no_cnt[0]);
  // reset lazy_update
#ifndef RECORD_ACCESS_TIME
  lazy_update.clear();
#endif
  // reset tag
  need_to_decode = false;
  have_decoded = false;

  if (use_cm_sketch) {
    if (!use_negative_counters) {
      for (size_t i = 0; i < cm_row; ++i)
        for (size_t j = 0; j < cm_width; ++j)
          cm_sketch[i][j].setVal(0);
    } else {
      int length = 0;
      for (int i = 1; i < no_layer; i++) {
        length += width_cnt[i];
      }
      int offset = (1 << (length - 1));
      for (size_t i = 0; i < cm_row; ++i)
        for (size_t j = 0; j < cm_width; ++j)
          cm_sketch[i][j].setVal(offset);
    }
  }
}

template <int32_t no_layer, typename T, typename hash_t>
void CounterHierarchy<no_layer, T, hash_t>::print_rate(const char *name) {
  size_t length = no_cnt[0];
  int32_t overflow_num = 0;
  int32_t correct_num = 0;
  double are = 0;
  for (int32_t i = 0; i < length; i++) {
    if (getStatus(i) == true) {
      overflow_num++;
      are += 1.0 * std::abs(getCnt(i) - getOriginalCnt(i)) /
             ((getOriginalCnt(i) >> width_cnt[0]) << width_cnt[0]);
    }
    if (getCnt(i) == getOriginalCnt(i)) {
      correct_num++;
    }
  }
  if (use_cm_sketch) {
    printf("\n %s EST ERROR %d / %d = %lf TIMES! ARE: %lf\n", name,
           est_ERROR_TIME, est_TIMES, (double)est_ERROR_TIME / est_TIMES,
           est_ARE / est_TIMES);
  }
#ifdef RECORD_ACCESS_TIME
  printf("\n%s AVERAGE UPDATE RATE: %lf\n", name,
         (double)access_time / total_update_time);
#endif
  printf("\n%s OVERFLOW RATE: %lf\n\n", name, (double)overflow_num / length);
  printf("\n%s CORRECT NUM: %lf\n\n", name, (double)correct_num / length);
  printf("\n%s ARE: %lf\n\n", name, (double)are / length);
}

template <int32_t no_layer, typename T, typename hash_t>
int CounterHierarchy<no_layer, T, hash_t>::guess_sign(int32_t negative_weight,
                                                      int32_t positive_weight,
                                                      int32_t idx,
                                                      int32_t layer) {
  if (!use_negative_counters) {
    return 1;
  }

#ifndef RECORD_ACCESS_TIME
  // lazy update
  if (!lazy_update.empty()) {
    for (int32_t i = 0; i < no_layer; i++) {
      lazy_update = updateLayer(i, std::move(lazy_update)); // throw exception
    }
    lazy_update.clear();
  }
#endif

  // The last layer should not overflow

  if (status_bits[layer][idx]) {
    int result = 0;
    for (int i = 0; i < no_hash[layer]; i++) {
#ifndef SKIP_HASH
      result +=
          guess_sign(negative_weight, positive_weight,
                     hash_fns[layer][i](idx) % width_cnt[layer + 1], layer + 1);
#else
      result += guess_sign(negative_weight, positive_weight,
                           my_hash(layer, i, idx), layer + 1);
#endif
    }
    if (result < 0) {
      return -negative_weight;
    } else {
      return positive_weight;
    }
  } else {
    if (cnt_array[layer][idx].getVal() < (1 << (width_cnt[layer] - 1))) {
      return -negative_weight;
    }
    return positive_weight;
  }
}

template <int32_t no_layer, typename T, typename hash_t>
T CounterHierarchy<no_layer, T, hash_t>::get_current_cnt(size_t idx) {
  if (idx >= no_cnt[0]) {
    throw std::out_of_range("Index Out of Range: Should be in [0, " +
                            std::to_string(no_cnt[0] - 1) + "], but got " +
                            std::to_string(idx) + " instead.");
  }

#ifndef RECORD_ACCESS_TIME
  // lazy update
  if (!lazy_update.empty()) {
    for (int32_t i = 0; i < no_layer; i++) {
      lazy_update = updateLayer(i, std::move(lazy_update)); // throw exception
    }
    lazy_update.clear();
  }
#endif

  return cnt_array[0][idx].getVal();
}
template <int32_t no_layer, typename T, typename hash_t>
T CounterHierarchy<no_layer, T, hash_t>::getTotalCnt(int32_t idx,
                                                     int32_t layer) {
  T ans = 0;
  if (!status_bits[layer][idx]) {
    ans = cnt_array[layer][idx].getVal();
  } else {
    T result = (1 << (width_cnt[layer + 1]));
    for (int i = 0; i < no_hash[layer]; i++) {
#ifndef SKIP_HASH
      result = std::min(
          result,
          getTotalCnt(hash_fns[layer][i](idx) % no_cnt[layer + 1], layer + 1));
#else
      result = std::min(result, getTotalCnt(my_hash(layer, i, idx), layer + 1));
#endif
    }
    ans = ((result << (width_cnt[layer])) + cnt_array[layer][idx].getVal());
  }
  return ans;
}

template <int32_t no_layer, typename T, typename hash_t>
T CounterHierarchy<no_layer, T, hash_t>::getEstCnt(int32_t idx, int32_t layer) {
  if (!use_cm_sketch) {
    if (!status_bits[layer][idx]) {
      return cnt_array[layer][idx].getVal();
    } else {
      T result = (1 << (width_cnt[layer + 1]));
      for (int i = 0; i < no_hash[layer]; i++) {
#ifndef SKIP_HASH
        result = std::min(
            result,
            getEstCnt(hash_fns[layer][i](idx) % no_cnt[layer + 1], layer + 1));
#else
        result = std::min(result, getEstCnt(my_hash(layer, i, idx), layer + 1));
#endif
      }
      return ((result << (width_cnt[layer])) + cnt_array[layer][idx].getVal());
    }
  } else {
    // assert(layer == 0);
    T val = cnt_array[0][idx].getVal();
    if (status_bits[0][idx]) {
      if (!use_negative_counters) {
        T val_min = std::numeric_limits<T>::max();
        for (size_t i = 0; i < cm_row; ++i) {
          size_t ind = cm_hash[i](idx) % cm_width;
          val_min = std::min(val_min, cm_sketch[i][ind].getVal());
        }
        // val += (sketch.at(index) << width_cnt[0]);
        val += (val_min << width_cnt[0]);
      } else {
        T result[cm_row];
        for (size_t i = 0; i < cm_row; ++i) {
          size_t ind = cm_hash[i](idx) % cm_width;
          result[i] = cm_sketch[i][ind].getVal();
        }
        int length = 0;
        for (int i = 1; i < no_layer; i++) {
          length += width_cnt[i];
        }
        int offset = (1 << (length - 1));
        std::sort(result, result + cm_row);
        val += ((result[cm_row / 2] - offset) << width_cnt[0]);
      }
    }
    if (use_negative_counters) {
      val -= ((1 << (width_cnt[0] - 1)));
    } else {
      val = std::min(val, getTotalCnt(idx));
    }
    est_TIMES++;
    T ori = getOriginalCnt(idx);
    if (val != ori) {
      est_ERROR_TIME++;
      est_ARE += (ori == 0) ? 0 : std::abs((ori - val) / ori);
    }
    return val;
  }
}

} // namespace OmniSketch::Sketch