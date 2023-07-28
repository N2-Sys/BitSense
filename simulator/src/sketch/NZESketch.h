/**
 * @file NZESketch.h
 * @author XierLabber
 * @brief Implementation of NZE Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <common/hash.h>
#include <common/sketch.h>
#include <common/utils.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include <map>
#include "BloomFilter.h"

namespace OmniSketch::Sketch {
template <int32_t key_len, typename T, typename hash_t> 
class NZESketch : public SketchBase<key_len, T> {

class HashTableElem{
public:
  FlowKey<key_len> f;
  T c;
  T d;
  bool empty() { return (c == 0) && (d == 0); }
};

  HashTableElem* HT;         // hash table
  hash_t HTHashFunc;
  int32_t HTLength;

  std::map<FlowKey<key_len>, T> EvictResult;
  std::map<FlowKey<key_len>, T> DecodedValue;
  std::vector<FlowKey<key_len>> RecordedKeys;

  int32_t BFBitsNum;
  int32_t BFHashNum;
  BloomFilter<key_len, hash_t> BF;

  int32_t FSdepth, FSwidth;
  T** FS;
  hash_t* FSHashFunc;

  bool have_decoded;

public:
  NZESketch(int32_t HTLength, int32_t BFBitsNum, int32_t BFHashNum, 
      int32_t FSdepth, int32_t FSwidth);
  ~NZESketch();
  NZESketch(NZESketch &&) = delete;
  NZESketch &operator=(const NZESketch &) = delete;
  NZESketch &operator=(NZESketch &&) = delete;

  void update(const FlowKey<key_len> &flow_key, T val);
  void recover();
  T query(const FlowKey<key_len> &flow_key) const;

  void clear();
  size_t size() const;
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
NZESketch<key_len, T, hash_t>::NZESketch(
    int32_t HTLength, int32_t BFBitsNum, int32_t BFHashNum, 
    int32_t FSdepth, int32_t FSwidth) : HTLength(Util::NextPrime(HTLength)), 
    BFBitsNum(Util::NextPrime(BFBitsNum)), BFHashNum(BFHashNum), 
    FSdepth(FSdepth), FSwidth(Util::NextPrime(FSwidth)), BF(BFBitsNum, BFHashNum){
    
    HT = new HashTableElem[this->HTLength]();
    FS = new T* [this->FSdepth];
    FS[0] = new T[this->FSdepth * this->FSwidth]();
    for(int i = 1; i < this->FSdepth; i++){
        FS[i] = FS[i - 1] + this->FSwidth;
    }
    FSHashFunc = new hash_t[this->FSdepth];

    have_decoded = false;
}

template <int32_t key_len, typename T, typename hash_t>
NZESketch<key_len, T, hash_t>::~NZESketch(){
    delete[] HT;
    delete[] FS[0];
    delete[] FS;
    delete[] FSHashFunc;
    EvictResult.swap();
    DecodedValue.swap();
    RecordedKeys.swap();
    BF.~BloomFilter();
}

template <int32_t key_len, typename T, typename hash_t>
void NZESketch<key_len, T, hash_t>::clear(){
    std::fill(HT, HT + HTLength, {0, 0, 0});
    std::fill(FS[0], FS[0] + FSwidth * FSdepth, 0);
    BF.clear();
    EvictResult.clear();
    DecodedValue.clear();
    RecordedKeys.clear();
}

template <int32_t key_len, typename T, typename hash_t>
void NZESketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flow_key, T val){
  size_t idx = HTHashFunc(flow_key) % HTLength;
  HashTableElem* bucket = HT + idx;
  if(bucket->empty()){
    bucket->f = flow_key;
    bucket->c = val;
    bucket->d = 0;
  }
  else if(HT[idx].f == flow_key){
    bucket->c += val;
  }
  else{
    bucket->d = bucket->d + val;
    if(bucket->d >= bucket->c){
      EvictResult[bucket->f] += bucket->c;
      bucket->f = flow_key;
      bucket->c = val;
      bucket->d = 0;
    }
    else{
      for(int i = 0; i < FSdepth; i++){
        size_t idx = FSHashFunc[i](flow_key) % FSwidth;
        FS[i][idx] += val;
      }
      if(!BF.lookup(flow_key)){
        RecordedKeys.push_back(flow_key);
        BF.insert(flow_key);
      }
    }
  }
}

template <int32_t key_len, typename T, typename hash_t>
void NZESketch<key_len, T, hash_t>::recover(){
  std::sort(RecordedKeys.begin(), RecordedKeys.end());
  RecordedKeys.erase(std::unique(RecordedKeys.begin(), RecordedKeys.end()), RecordedKeys.end());

  uint32_t w = FSwidth;
  int M = FSdepth * w, N = RecordedKeys.size();

  Eigen::VectorXd X(N), b(M);
  Eigen::SparseMatrix<double> A(M, N);
  std::vector<Eigen::Triplet<double>> tripletlist;

  for(int i = 0, j; i < FSdepth; i++){
    for(int j = 0; j < RecordedKeys.size(); j++){
      int idx = i * w + FSHashFunc[i](RecordedKeys[j]) % w;
      tripletlist.push_back(Eigen::Triplet<double>(idx, j, 1));
    }
    for(int j = 0; j < w; j++){
      b(i * w + j) = FS[i][j];
    }
  }

  A.setFromTriplets(tripletlist.begin(), tripletlist.end());
  A.makeCompressed();
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> Solver_sparse;
  Solver_sparse.compute(A);
  X = Solver_sparse.solve(b);

  DecodedValue.clear();

  for(int i = 0; i < N; i++){
    T ans = static_cast<T>(X[i] + 0.5);
    DecodedValue[RecordedKeys[i]] = (ans <= 0)? 1 : ans;
  }

  have_decoded = true;
}

template <int32_t key_len, typename T, typename hash_t>
T NZESketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flow_key) const{
  if(!have_decoded){
    const_cast<NZESketch<key_len, T, hash_t>* >(this)->recover();
  }
  T ans = 0;
  HashTableElem* elem = HT + (HTHashFunc(flow_key) % HTLength);
  if(elem->f == flow_key){
    ans += elem->c;
  }

  auto cit = EvictResult.find(flow_key);
  if (cit != EvictResult.end()) {
      ans += cit->second;
  }

  auto cit_ = DecodedValue.find(flow_key);
  if (cit_ != DecodedValue.end()) {
      ans += cit_->second;
  }

  return ans;
}

template <int32_t key_len, typename T, typename hash_t>
size_t NZESketch<key_len, T, hash_t>::size() const{
    return sizeof(*this)
           + HTLength * (2 * sizeof(T) + sizeof(FlowKey<key_len>))
           + BF.size()
           + FSdepth * sizeof(hash_t)
           + FSdepth * FSwidth * sizeof(T);
}

} // namespace OmniSketch::Sketch