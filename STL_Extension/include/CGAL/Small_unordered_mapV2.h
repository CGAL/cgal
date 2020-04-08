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

//#define CGAL_SMALL_UNORDERED_MAP_STATS
namespace CGAL {


template <typename K, typename T, typename H, unsigned int M, unsigned int Factor>
class Small_unordered_mapV2 {
#ifdef    CGAL_SMALL_UNORDERED_MAP_STATS
  std::array<int, 20> collisions = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif
  static constexpr int B = M * Factor; // the number of bins
  int head = -2;
  mutable std::array<int, B>    occupied;
  std::array<int, B>            unfreelist;
  std::array<std::pair<K, T>, B> data;
  const H hash = {};
  boost::unordered_map<K, T> big;
  std::size_t N = 0; // the number of stored elements
public:

  Small_unordered_mapV2()
  {
    occupied.fill(-1);
  }

#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
  ~Small_unordered_mapV2()
  {
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
        if (occupied[i] == -1) {
          occupied[i] = 1;
          data[i].first = k;
          data[i].second = t;
          unfreelist[i] = head;
          head = i;
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
      int pos = head;
      while(pos != -2){
        big.insert(data[pos]);
        pos = unfreelist[pos];
      }
      big[k] = t;
    }
    else {
      big[k] = t;
    }
    ++N;
  }


  const T& get(const K& k) const
  {
    if (N < M) {
      unsigned int h = hash(k) % B;
      unsigned int i = h;
      do {
        if ((occupied[i] == 1) && (data[i].first == k)) {
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

  std::size_t size() const
  {
    return N;
  }

  void clear()
  {
    head = -2;
    occupied.fill(-1);
    big.clear();
    N = 0;
  }

  struct iterator {
    typedef std::pair<K, T>           value_type;
    typedef const std::pair<K, T>&    reference;
    typedef std::size_t               size_type;
    typedef std::ptrdiff_t            difference_type;
    typedef std::forward_iterator_tag iterator_category;

    const Small_unordered_mapV2& map;
    int pos;
    typename boost::unordered_map<K, T>::const_iterator bigit;
    bool big = false;

    iterator(const Small_unordered_mapV2& map)
      : map(map), pos(-2), bigit(map.big.end()), big(map.N > M)
    {}

    iterator(const Small_unordered_mapV2& map, int pos)
      : map(map), pos(pos), bigit(map.big.begin()), big(map.N > M)
    {}

    bool operator==(const iterator& other) const
    {
      if(big){
        return bigit == other.bigit;
      }else{
        return pos == other.pos;
      }
    }

    bool operator!=(const iterator& other) const
    {
      return !((*this) == other);
    }
    iterator operator++()
    {
      if (big) {
        ++bigit;
      }
      else {
        pos = map.unfreelist[pos];
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


  iterator begin() const
  {
    return iterator(*this, head);
  }

  iterator end() const
  {
    return iterator(*this);
  }

  void clear(const iterator it)
  {
    occupied[it.pos] = -1;
  }

  friend struct iterator;
};

} // namespace CGAL
#endif // CGAL_SMALL_UNORDERED_MAPV2_H
