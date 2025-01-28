<<<<<<< HEAD
// Copyright (c) 2019  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
=======
// Copyright (c) 2019 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
>>>>>>> cgal/master
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_SMALL_UNORDERED_MAP_H
#define CGAL_SMALL_UNORDERED_MAP_H

#include <array>
#include <iostream>

//#define CGAL_SMALL_UNORDERED_MAP_STATS
namespace CGAL {
<<<<<<< HEAD
  
template <typename K, typename T, typename H, unsigned int N>
class Small_unordered_map{
#ifdef    CGAL_SMALL_UNORDERED_MAP_STATS
  std::array<int,10> collisions = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#endif
  static const unsigned int M = N<<3;
  int head = -2;
  mutable std::array<int, M>           occupied;
  std::array<int, M>           unfreelist; 
  std::array<std::pair<K,T>, M> data;
  const H hash = {};

  
=======


template <typename K, typename T, typename H, unsigned int M>
class Small_unordered_map{
#ifdef    CGAL_SMALL_UNORDERED_MAP_STATS
  std::array<int,20> collisions = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#endif
  int head = -2;
  mutable std::array<int, M>    occupied;
  std::array<int, M>            unfreelist;
  std::array<std::pair<K,T>, M> data;
  const H hash = {};

>>>>>>> cgal/master
public:
  Small_unordered_map()
  {
    occupied.fill(-1);
  }

<<<<<<< HEAD
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS  
  ~Small_unordered_map()
  {
    for(int i = 0; i < 10; i++){
      std::cout << i << " " << collisions[i] << std::endl;
    }
  }
#endif
  
  /// Set only once for a key and not more than N 
=======
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
  ~Small_unordered_map()
  {
    int total = 0;
    std::cout << "0 " << collisions[0] << std::endl;
    for(int i = 1; i < 20; i++){
      total += collisions[i];
      if(collisions[i] != 0){
        std::cout << i << " " << collisions[i] << std::endl;
      }
    }
    std::cout << "Total: " << total << " " << 100 * (double(total) / double(total + collisions[0])) << "%" << std::endl;
  }
#endif

  /// Set only once for a key and not more than N
>>>>>>> cgal/master
  void set(const K& k, const T& t)
  {
    unsigned int h  = hash(k)%M;
    unsigned i = h;
<<<<<<< HEAD
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS    
    int collision = 0;
#endif     
=======
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
    int collision = 0;
#endif
>>>>>>> cgal/master
    do {
      if(occupied[i]== -1){
        occupied[i] = 1;
        data[i].first = k;
        data[i].second = t;
        unfreelist[i] = head;
        head = i;
<<<<<<< HEAD
#ifdef  CGAL_SMALL_UNORDERED_MAP_STATS        
        if(collision>9){
          std::cerr << collision << " collisions" << std::endl;
        }
        ++collisions[collision];
#endif        
        return;
      }
      i = (i+1)%M;
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS   
      ++collision;
#endif      
=======
#ifdef  CGAL_SMALL_UNORDERED_MAP_STATS
        if(collision>19){
          std::cerr << collision << " collisions" << std::endl;
        }else{
          ++collisions[collision];
        }
#endif
        return;
      }
      i = (i+1)%M;
#ifdef CGAL_SMALL_UNORDERED_MAP_STATS
      ++collision;
#endif
>>>>>>> cgal/master
    }while(i != h);
    CGAL_error();
  }

<<<<<<< HEAD
  
  const T& get(const K& k) const
  {
    unsigned int h  = hash(k)%M;
    int i = h;
=======
  // et only once as it is erased
  const T& get_and_erase(const K& k) const
  {
    unsigned int h  = hash(k)%M;
    unsigned int i = h;
>>>>>>> cgal/master
    do{
      if((occupied[i] == 1) && (data[i].first == k)){
        occupied[i] = -1;
        return data[i].second;
      }
      i = (i+1)%M;
    }while(i != h);
    CGAL_error();
  }

  void clear()
  {
    head = -2;
<<<<<<< HEAD
    //    occupied.fill(-1);
=======
    // without erase we would have to call occupied.fill(-1); which is costly
>>>>>>> cgal/master
  }

  struct iterator {
    const Small_unordered_map& map;
    int pos;

    iterator(const Small_unordered_map& map)
      : map(map),pos(-2)
    {}

    iterator(const Small_unordered_map& map, int pos)
      :  map(map), pos(pos)
    {}

    bool operator==(const iterator& other) const
    {
      return pos == other.pos;
    }

    bool operator!=(const iterator& other) const
    {
      return pos != other.pos;
    }
    iterator operator++()
    {
      pos = map.unfreelist[pos];
      return *this;
    }
<<<<<<< HEAD
    
=======

>>>>>>> cgal/master
    const std::pair<K,T>& operator*()const
    {
      return map.data[pos];
    }
  };

  iterator begin() const
  {
    return iterator(*this,head);
  }

  iterator end() const
  {
    return iterator(*this);
  }

  void clear(const iterator it)
  {
    occupied[it.pos] = -1;
  }
<<<<<<< HEAD
  
=======

>>>>>>> cgal/master
  friend struct iterator;
};

} // namespace CGAL
#endif // CGAL_SMALL_UNORDERED_MAP_H
