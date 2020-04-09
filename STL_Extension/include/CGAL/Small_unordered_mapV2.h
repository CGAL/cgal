// Copyright (c) 2019 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_SMALL_UNORDERED_MAPV2_H
#define CGAL_SMALL_UNORDERED_MAPV2_H

#include <boost/unordered_map.hpp>
#include <array>
#include <iostream>
#include <bitset>

//#define CGAL_SMALL_UNORDERED_MAP_STATS
namespace CGAL {


template <typename K, typename T, typename H, unsigned int M, unsigned int Factor>
class Small_unordered_mapV2 {
#ifdef    CGAL_SMALL_UNORDERED_MAP_STATS
  std::array<int, 20> collisions = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif
  static constexpr int B = M * Factor ; // the number of bins
  int head = B;
  // mutable std::array<boost::int8_t, B> occupied {}; // all 0
  mutable std::bitset<B> occupied;
  std::array<std::pair<K, T>, B> data;
  const H hash = {};
  boost::unordered_map<K, T> * big = nullptr;
  std::size_t N = 0; // the number of stored elements
public:

  Small_unordered_mapV2()
  {}

#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
  ~Small_unordered_mapV2()
  {
    std::cout << "N = "<< N << std::endl;
    int total = 0;
    std::cout << "0 " << collisions[0] << std::endl;
    for (int i = 1; i < 20; i++) {
      total += collisions[i];
      if (collisions[i] != 0) {
        std::cout << i << " " << collisions[i] << std::endl;
      }
    }
    std::cout << "Total: " << total << " " << 100 * (double(total) / double(total + collisions[0])) << "%" << std::endl;
  }
#else
  ~Small_unordered_mapV2()
  {
    if(big != nullptr){
      delete big;
    }
  }
#endif

  /// Set only once for a key
  void set(const K& k, const T& t)
  {
    if (N < M) {
      unsigned int h = hash(k) % B;
      unsigned i = h;
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
      int collision = 0;
#endif
      do {
        if (occupied[i] == false) {
          occupied[i] = true;
          data[i].first = k;
          data[i].second = t;
          if(i < head){
            head = i;
          }
#ifdef  CGAL_SMALL_UNORDERED_MAP_STATS
          if (collision > 19) {
            std::cerr << collision << " collisions" << std::endl;
          }
          else {
            ++collisions[collision];
          }
#endif
          ++N;
          return;
        }
        i = (i + 1) % B;
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
        ++collision;
#endif
      } while (i != h);
      CGAL_error();
    }
    else if (N == M) {
      big = new boost::unordered_map<K, T>(); 
      for(int pos = head; pos < B; ++pos){
        if(occupied[pos]){
          big->insert(data[pos]);
        }
      }
      (*big)[k] = t;
    }
    else {
      (*big)[k] = t;
    }
    ++N;
  }

  
  void insert(const std::pair<K,T>& p)
  {
    set(p.first,p.second);
  }
  

  const T& at(const K& k) const
  {
    if (N <= M) {
      unsigned int h = hash(k) % B;
      unsigned int i = h;
      do {
        if ((occupied[i]) && (data[i].first == k)) {
          return data[i].second;
        }
        i = (i + 1) % B;
      } while (i != h);
      CGAL_error();
    }
    else {
      return big.at(k);
    }
  }

  struct const_iterator;

  const_iterator find(const K& k) const
  {
     if (N <= M) {
       unsigned int h = hash(k) % B;
      unsigned int i = h;
      do {
        if(! occupied[i]){
          // the element is not in the map
          return end();
        }
        if (data[i].first == k) {
          return const_iterator(*this,i);
        }
        i = (i + 1) % B;
      } while (i != h);
      CGAL_error();
    }
    else {
      return const_iterator(*this, big->find(k));
    }
  }

  
  std::size_t size() const
  {
    if(big != nullptr){
      return big->size();
    }
    return N;
  }

  /*
  void clear()
  {
    head = B;
    // occupied.fill(0);
    if(big != nullptr){
      delete big;
      big = nullptr;
    }
    N = 0;
  }
  */

  struct const_iterator {
    typedef std::pair<K, T>           value_type;
    typedef const std::pair<K, T>&    reference;
    typedef std::size_t               size_type;
    typedef std::ptrdiff_t            difference_type;
    typedef std::forward_iterator_tag iterator_category;

    const Small_unordered_mapV2& map;
    int pos;
    
    typedef typename boost::unordered_map<K, T>::const_iterator Bigit;
    Bigit bigit;
    bool big = false;

    const_iterator(const Small_unordered_mapV2& map)
      : map(map), pos(B), big(map.N > M)
    {
      if(big){
        bigit = map.big->end();
      }
    }

    const_iterator(const Small_unordered_mapV2& map, int pos)
      : map(map), pos(pos), big(map.N > M)
    {
      if(big){
        bigit = map.big->begin();
      }
    }

    const_iterator(const Small_unordered_mapV2& map, const Bigit& bigit)
      : map(map), pos(B), bigit(bigit), big(true)
    {}
    
    bool operator==(const const_iterator& other) const
    {
      if(big){
        return bigit == other.bigit;
      }else{
        return pos == other.pos;
      }
    }

    bool operator!=(const const_iterator& other) const
    {
      return !((*this) == other);
    }
    const_iterator operator++()
    {
      if (big) {
        ++bigit;
      } else {
        if(pos != B){
          do {
            ++pos;
          }while((pos !=B) && (! map.occupied[pos]));
        }
      }
      return *this;
    }

    value_type operator*()const
    {
      if (big) {
        return *bigit;
      }
      else {
        return map.data[pos];
      }
    }
  };


  const_iterator begin() const
  {
    return const_iterator(*this, head);
  }

  const_iterator end() const
  {
    return const_iterator(*this);
  }


  friend struct const_iterator;
};

} // namespace CGAL
#endif // CGAL_SMALL_UNORDERED_MAPV2_H
