// Copyright (c) 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BLOCKLIST_H
#define CGAL_BLOCKLIST_H

#include <list>

#include <CGAL/array.h>

namespace CGAL {

  template < typename T, int N>
  class Block_list{
    typedef CGAL::cpp11::array<T,N> Block;
    int index;
    std::list<Block> blocks;
    typedef typename std::list<Block>::iterator List_iterator;
    typedef typename Block::iterator Block_iterator;

  public:    
    Block_list()
      : index(N)
    {}

    std::size_t size() const
    {
      return (index == N)? (blocks.size()*N) : ((blocks.size()-1)*N+index);
    }


    T* push_back(const T& t)
    {
      if(index==N){
        Block b;
        blocks.push_back(b);
        index = 0;
      }
      Block& b = blocks.back();
      b[index] = t;
      T* ptr = &b[index];
      ++index;
      return ptr;
    }
    


    void clear()
    {
      blocks.clear();
      index = N;
    }

    void print(std::ostream& os)
    {
      std::size_t s = blocks.size();
      std::size_t i = 0;
      for(std::list<Block>::iterator it = blocks.begin(); it!=  blocks.end(); ++it, ++i){
        Block& b = *it;
        os << "Block "<< i << std::endl;
        std::size_t n = (i+1 == s)? index : N;
        for(std::size_t j=0; j < n; j++){
          os << b[j] << std::endl;
        }
      }
    }
};

} // namespace CGAL

#endif 
