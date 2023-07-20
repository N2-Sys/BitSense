/**
 * @file StreamSummary.h
 * @author XierLabber (you@domain.com)
 * @brief Stream Summary
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <common/hash.h>
#include <common/sketch.h>
#include <stdint.h>

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
class hash_table_elem;

template <int32_t key_len, typename T, typename hash_t>
class bucket_list_elem;

template <int32_t key_len, typename T, typename hash_t>
class StreamSummary;

template <int32_t key_len, typename T, typename hash_t>
class hash_table_elem{
public:
    bucket_list_elem<key_len, T, hash_t>* parent;
    hash_table_elem<key_len, T, hash_t>* next;
    hash_table_elem<key_len, T, hash_t>* prev;
    hash_table_elem<key_len, T, hash_t>* next_hash_elem;
    FlowKey<key_len> first;
    T second;
    hash_table_elem(): second(-1), parent(NULL), next(NULL), next_hash_elem(NULL){};
};

template <int32_t key_len, typename T, typename hash_t>
class bucket_list_elem{
public:
    T value;
    bucket_list_elem<key_len, T, hash_t>* prev;
    bucket_list_elem<key_len, T, hash_t>* next;
    hash_table_elem<key_len, T, hash_t>* child;
    bucket_list_elem(): value(-1), prev(NULL), next(NULL), child(NULL){};
};

template <int32_t key_len, typename T, typename hash_t>
class StreamSummary{
private:
    hash_table_elem<key_len, T, hash_t>** hash_table;
    bucket_list_elem<key_len, T, hash_t>* bucket_list_head;
    bucket_list_elem<key_len, T, hash_t>* bucket_list_tail;
    int32_t num_threshold;
    int32_t hash_table_length;
    hash_t hash_func;
    int32_t size_;

    void create_new_bucket(bucket_list_elem<key_len, T, hash_t>* prev, hash_table_elem<key_len, T, hash_t>* h);
    void insert_bucket_list_with_hash_table(hash_table_elem<key_len, T, hash_t>* h);
    void delete_bucket(bucket_list_elem<key_len, T, hash_t>* b);
    void push_forward(hash_table_elem<key_len, T, hash_t>* h);

public:
    void init(int32_t num_threshold_, double hash_table_alpha);
    void destroy();
    void emplace(FlowKey<key_len> key, T val);
    void insert(FlowKey<key_len> key, T val);
    void erase(hash_table_elem<key_len, T, hash_t>* h);
    void increment(hash_table_elem<key_len, T, hash_t>* h, T delta);
    hash_table_elem<key_len, T, hash_t>* find(FlowKey<key_len> key);
    hash_table_elem<key_len, T, hash_t>* get_least_elem();
    Data::Estimation<key_len, T> getHeavyHitter(double val_threshold) const;

    int32_t size() const;
    size_t memory_size() const;
    int32_t get_hashtable_length() const { return hash_table_length;}
    T get_hashtable_val(int32_t idx) const { return hash_table[idx].second;}
    FlowKey<key_len> get_hashtable_key(int32_t idx) const { return hash_table[idx].first;}

};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::init(int32_t num_threshold_, double hash_table_alpha){
  num_threshold = num_threshold_;
  hash_table_length = Util::NextPrime((int32_t)(num_threshold * hash_table_alpha));
  hash_table = new hash_table_elem<key_len, T, hash_t>* [hash_table_length];
  for(int i = 0; i < hash_table_length; i++)
  {
    hash_table[i] = NULL;
  }
  bucket_list_head = new bucket_list_elem<key_len, T, hash_t>;
  bucket_list_tail = new bucket_list_elem<key_len, T, hash_t>;
  bucket_list_head->next = bucket_list_tail;
  bucket_list_tail->prev = bucket_list_head;
  bucket_list_tail->value = 0xfffffff;
  bucket_list_head->value = -1;
  size_ = 0;
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::destroy(){
  for(int i = 0; i < hash_table_length; i++)
  {
    hash_table_elem<key_len, T, hash_t>* ptr = hash_table[i];
    while(ptr != NULL)
    {
      hash_table_elem<key_len, T, hash_t>* tmp = ptr;
      ptr = ptr->next_hash_elem;
      delete tmp;
    }
  }
  delete[] hash_table;
  bucket_list_elem<key_len, T, hash_t>* ptr = bucket_list_head;
  while(ptr != bucket_list_tail)
  {
    ptr = ptr->next;
    delete ptr->prev;
  }
  delete bucket_list_tail;
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::insert(FlowKey<key_len> key, T val){
  hash_table_elem<key_len, T, hash_t>* h = find(key);
  if(h == NULL)
  {
    if(size_ < num_threshold)
    {
      emplace(key, val);
    }
    else
    {
      hash_table_elem<key_len, T, hash_t>* min_h = get_least_elem();
      if(min_h == NULL || min_h->second < val){
        erase(min_h);
        emplace(key, val);
      }
    }
  }
  else
  {
    increment(h, val - h->second);
  }
}

template <int32_t key_len, typename T, typename hash_t>
hash_table_elem<key_len, T, hash_t>* StreamSummary<key_len, T, hash_t>::find(
  FlowKey<key_len> key){
  int32_t index = hash_func(key) % hash_table_length;
  hash_table_elem<key_len, T, hash_t>* ptr = hash_table[index];
  while(ptr != NULL)
  {
    if(ptr->first == key)
    {
      return ptr;
    }
    ptr = ptr->next_hash_elem;
  }
  return NULL;
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::create_new_bucket(
  bucket_list_elem<key_len, T, hash_t>* prev, 
  hash_table_elem<key_len, T, hash_t>* h){

  #ifdef DEBUG
  assert(prev != NULL);
  assert(h != NULL);
  #endif

  bucket_list_elem<key_len, T, hash_t>* tmp = new bucket_list_elem<key_len, T, hash_t>;
  tmp->next = prev->next;
  tmp->prev = prev;
  prev->next->prev = tmp;
  prev->next = tmp;
  tmp->child = h;
  tmp->value = h->second;
  h->next = h;
  h->prev = h;
  h->parent = tmp;
  #ifdef DEBUG
  assert(h->parent != NULL);
  #endif
  /*
  if(tmp->value == 1){
    printf"VAL 1 CREATE! %p, h: %p\n",tmp, h);
  }
  */
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::insert_bucket_list_with_hash_table(
  hash_table_elem<key_len, T, hash_t>* h){
  #ifdef DEBUG
  assert(h != NULL);
  #endif

  bucket_list_elem<key_len, T, hash_t>* ptr = bucket_list_head;
  while(ptr->value < h->second)
  {
      ptr = ptr->next;
  }
  #ifdef DEBUG
  assert(ptr != NULL);
  #endif
  if(ptr->value == h->second)
  {
      h->parent = ptr;
      h->next = ptr->child->next;
      h->prev = ptr->child;
      h->next->prev = h;
      ptr->child->next = h;
  }
  else
  {
      create_new_bucket(ptr->prev, h);
  }
  #ifdef DEBUG
  assert(h->parent != NULL);
  assert(h->parent->child != NULL);
  assert(h->parent->child->parent != NULL);
  #endif
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::emplace(FlowKey<key_len> key, T val){
  int32_t index = hash_func(key) % hash_table_length;
  hash_table_elem<key_len, T, hash_t>* ptr = hash_table[index];
  if(ptr == NULL)
  {
    hash_table[index] = new hash_table_elem<key_len, T, hash_t>;
    ptr = hash_table[index];
  }
  else
  {
    while(ptr->next_hash_elem != NULL)
    {
      ptr = ptr->next_hash_elem;
    }
    ptr->next_hash_elem = new hash_table_elem<key_len, T, hash_t>;
    ptr = ptr->next_hash_elem;
  }
  ptr->first = key;
  ptr->second = val;
  insert_bucket_list_with_hash_table(ptr);
  size_++;
  #ifdef DEBUG
  assert(ptr->parent != NULL);
  assert(ptr->parent->child->parent != NULL);
  #endif
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::delete_bucket(
  bucket_list_elem<key_len, T, hash_t>* b){

  /*
  if(b->value == 1){
    // printf"delete bucket: %p\n",b);
  }
  */
  #ifdef DEBUG
  assert(b != NULL);
  assert(b->prev != NULL);
  assert(b->next != NULL);
  #endif
  // printf("B IS %p\n",b);
  // printf("B PREV IS %p\n",b->prev);
  // printf("B NEXT IS %p\n",b->next);
  b->next->prev = b->prev;
  b->prev->next = b->next;
  delete b;

}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::push_forward(
  hash_table_elem<key_len, T, hash_t>* h){
  
  if(h->parent->value == 1){
    // printf"push called, pushed: %p, parent: %p\n",h,h->parent);
  }

  bucket_list_elem<key_len, T, hash_t>* need_to_delete = NULL;
  // printf("PUSH: REACH HERE! h is %p\n",h);
  // printf("PUSH: REACH HERE! h parent is %p\n",h->parent);
  #ifdef DEBUG
  assert(h != NULL);
  assert(h->parent != NULL);
  assert(h->prev != NULL);
  assert(h->next != NULL);
  assert(h->parent->child != NULL);
  #endif

  if(h->parent->child == h){
    if(h->next != h){
      hash_table_elem<key_len, T, hash_t>* ptr = h->prev;
      ptr->next = h->next;
      h->next->prev = ptr;
      h->parent->child = ptr;
      #ifdef DEBUG
      assert(ptr->parent != NULL);
      assert(ptr->parent == h->parent);
      #endif
      /*
      if(h->parent->value == 1)
      {
        printf"%p NEW CHILD: %p\n", ptr->parent, ptr);
      }
      */
    }
    else{
      need_to_delete = h->parent;
      /*
      if(need_to_delete->value == 1){
        printf"VAL 1 BUT DELETE! %p, cur h: %p\n",need_to_delete,h);
      }
      */
    }
  }
  else{
    hash_table_elem<key_len, T, hash_t>* ptr = h->prev;
    h->next->prev = ptr;
    ptr->next = h->next;
  }

  #ifdef DEBUG
  if(need_to_delete != NULL){
    assert(h->parent->child->parent != NULL);
  }
  #endif

  // printf("PUSH: REACH HERE2!\n");
  // printf("H->PARENT->NEXT IS %p\n",h->parent->next);
  // printf("H->PARENT->NEXT->CHILD IS %p\n",h->parent->next->child);
  if(h->parent->next->value == h->second){
    bucket_list_elem<key_len, T, hash_t>* ptr = h->parent->next;
    #ifdef DEBUG
    assert(ptr != NULL);
    assert(ptr->child != NULL);
    assert(ptr->child->next != NULL);
    #endif
    h->parent = ptr;
    h->next = ptr->child->next;
    h->prev = ptr->child;
    h->next->prev = h;
    ptr->child->next = h;
  }
  else{
    // printf("HERE!?\n");
    create_new_bucket(h->parent, h);
  }
  if(need_to_delete != NULL){
    delete_bucket(need_to_delete);
  }

  #ifdef DEBUG
  assert(h->parent != NULL);
  #endif

  // printf("PUSH: REACH HERE3!\n");
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::erase(
  hash_table_elem<key_len, T, hash_t>* h){
  if(h->parent->value == 1){
    // printf"erase called! %p is erased, its parent is %p\n", h, h->parent);
  }
  if(h->parent == NULL)
  {
    // printf"H: %p, H->parent: %p,H->VAL: %d\n",h,h->parent,h->second);
  }
  #ifdef DEBUG
  assert(h->parent != NULL);
  #endif
  if(h == NULL)
  {
      return;
  }
  if(h->next == h)
  {
  // printf("REACH HERE! 0\n");
      delete_bucket(h->parent);
  }
  else
  {
  // printf("REACH HERE! 1\n");
      hash_table_elem<key_len, T, hash_t>* ptr = h->prev;
      h->next->prev = ptr;
      ptr->next = h->next;
      h->parent->child = h->next;
  }
  // printf("REACH HERE! 2\n");
  int32_t idx = hash_func(h->first) % hash_table_length;
  hash_table_elem<key_len, T, hash_t>* ptr = hash_table[idx];
  // printf("REACH HERE! 3\n");
  if(ptr == h)
  {
  // printf("REACH HERE! 4\n");
    hash_table[idx] = h->next_hash_elem;
  }
  else
  {
  // printf("REACH HERE! 5\n");
    while(ptr->next_hash_elem != h)
    {
      ptr = ptr->next_hash_elem;
    }
    ptr->next_hash_elem = h->next_hash_elem;
  }
  // printf("REACH HERE! 6\n");
  delete h;
  size_--;
  // printf("REACH HERE! 7\n");
}

template <int32_t key_len, typename T, typename hash_t>
hash_table_elem<key_len, T, hash_t>* 
StreamSummary<key_len, T, hash_t>::get_least_elem(){
  bucket_list_elem<key_len, T, hash_t>* ptr = bucket_list_head->next;
  if(ptr->value == 1){
    // printf"get called! parent: %p\n",ptr);
  }
  if(ptr->value == 0xfffffff)
  {
      return NULL;
  }
  #ifdef DEBUG
  assert(ptr->child != NULL);
  if(ptr->child->parent == NULL){
    // printf"err ptr: %p\nerr ptr val: %d, err child val: %d\nerror ptr->child: %p\n", ptr, ptr->value, ptr->child->second, ptr->child);
  }
  assert(ptr->child->parent != NULL);
  #endif
  return ptr->child;
}

template <int32_t key_len, typename T, typename hash_t>
void StreamSummary<key_len, T, hash_t>::increment(
  hash_table_elem<key_len, T, hash_t>* h, T delta){
  
  h->second += delta;
  push_forward(h);
}

template <int32_t key_len, typename T, typename hash_t>
int32_t StreamSummary<key_len, T, hash_t>::size() const{
  return size_;
}

template <int32_t key_len, typename T, typename hash_t>
Data::Estimation<key_len, T> 
StreamSummary<key_len, T, hash_t>::getHeavyHitter(double val_threshold) const{
  Data::Estimation<key_len, T> heavy_Hitter;
  int found = 0;
  for(int i=0; i < hash_table_length; i++)
  {
      hash_table_elem<key_len, T, hash_t>* ptr = hash_table[i];
      while(ptr != NULL)
      {
      T val = ptr->second;
      if(val >= val_threshold)
      {
          heavy_Hitter[ptr->first] = ptr->second;
          found++;
      }
      ptr = ptr->next_hash_elem;
      }
  }
  return heavy_Hitter;
}

template <int32_t key_len, typename T, typename hash_t>
size_t StreamSummary<key_len, T, hash_t>::memory_size() const{
  return sizeof(*this) + 
         hash_table_length * sizeof(hash_table_elem<key_len, T, hash_t>*) + 
         num_threshold * (4 * sizeof(void *) + sizeof(T) + sizeof(FlowKey<key_len>));
}

} // namespace OmniSketch::Sketch