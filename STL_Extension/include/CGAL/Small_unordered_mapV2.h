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


template <typename Mesh, typename K, typename T, int M>
class Small_unordered_mapV2 {

  const Mesh& mesh;
  // std::array<std::pair<K, T>, M> data;
  std::array<T, M> data;

  boost::unordered_map<K, T> * big = nullptr;
  int N = 0; // the number of stored elements
public:

  Small_unordered_mapV2(const Mesh& mesh)
    : mesh(mesh)
  {}


  ~Small_unordered_mapV2()
  {
    if(big != nullptr){
      delete big;
    }
  }


  /// Set only once for a key
  void set(const K& k, const T& t)
  {
    if (N < M) {
      //data[N].first = k;
      //data[N].second = t;
      data[N] = t;
      ++N;
    }else{
      if (N == M) {
        big = new boost::unordered_map<K, T>(); 
        for(int pos = 0; pos < M; ++pos){
          // big->insert(data[pos]);
          big->insert(std::make_pair(target(data[pos],mesh),data[pos]));
        }
        (*big)[k] = t;
      }else{
        (*big)[k] = t;
      }
      ++N;
    }
  }

  
  void insert(const std::pair<K,T>& p)
  {
    set(p.first,p.second);
  }
  

  
  struct const_iterator;

  const_iterator find(const K& k) const
  {
    if (N <= M) {
      for(int i =0; i < N; ++i){
        if (target(data[i],mesh) == k) {
          return const_iterator(*this,i);
        }
      }
      return end();
    }else{
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
      : map(map), pos(M), big(map.N > M)
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
      : map(map), pos(M), bigit(bigit), big(true)
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
      } else if(pos < map.N-1){
        ++pos;
      }else {
          CGAL_assertion((pos == map.N-1) || (pos = M));
          pos = M;
      }
      return *this;
    }

    value_type operator*()const
    {
      if (big) {
        return *bigit;
      }
      else {
        return std::make_pair(target(map.data[pos],map.mesh),map.data[pos]);
      }
    }
  };


  const_iterator begin() const
  {
    if(N==0){
      return const_iterator(*this);
    }
    return const_iterator(*this, 0);
  }

  const_iterator end() const
  {
    return const_iterator(*this);
  }


  friend struct const_iterator;
};

} // namespace CGAL
#endif // CGAL_SMALL_UNORDERED_MAPV2_H
