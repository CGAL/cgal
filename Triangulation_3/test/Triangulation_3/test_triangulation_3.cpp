// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Francois Rebufat

#include <CGAL/Triangulation_3.h>
#include <cassert>

bool del = false;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <iostream>
#include <CGAL/Unique_hash_map.h>
// Explicit instantiation of the whole class :
template class CGAL::Triangulation_3<K>;

template<typename Vb>
struct Custom_vertex_base : public Vb
{
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Custom_vertex_base<Vb2>           Other;
  };
  
  std::ostream& write_data(std::ostream& os,
                           const CGAL::Unique_hash_map<typename Vb::Vertex_handle, std::size_t > &
                           )const 
  {
    os << this->info() << std::endl;
    return os;
  }
  
  std::istream& read_data(std::istream& is,
                          std::vector<typename Vb::Cell_handle>&,
                          std::vector<typename Vb::Vertex_handle>&)
  {
    is >> this->info();
    return is;
  }
};

template<typename Cb>
struct Custom_cell_base : public Cb
{
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other       Cb2;
    typedef Custom_cell_base<Cb2>  Other;
  };

  Custom_cell_base():Cb(){}
  Custom_cell_base(typename Cb::Vertex_handle v0, typename Cb::Vertex_handle v1, typename Cb::Vertex_handle v2, typename Cb::Vertex_handle v3):Cb(v0, v1, v2, v3){}
  
  std::ostream& write_data(std::ostream& os,
                           const CGAL::Unique_hash_map<typename Cb::Vertex_handle, std::size_t > &
                           )const
  {
    os << this->info() << std::endl;
    return os;
  }
  
  std::istream& read_data(std::istream& is,
                          std::vector<typename Cb::Cell_handle>&,
                          std::vector<typename Cb::Vertex_handle>&)
  {
    is >> this->info();
    return is;
  }
};

typedef Custom_vertex_base<CGAL::Triangulation_vertex_base_with_info_3<int, K> > Vertex_base;
typedef Custom_cell_base<CGAL::Triangulation_cell_base_with_info_3<int, K> > Cell_base;

int main()
{
  typedef CGAL::Triangulation_3<K>                               Cls3;

  _test_cls_triangulation_3( Cls3() );

  // Test operator== between triangulations of different Tds types.

  typedef CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vertex_base, Cell_base > > Cls3_2;

  assert(Cls3() == Cls3_2());
  
  typedef Cls3_2::Vertex_handle                                      Vertex_handle;
  typedef Cls3_2::Point                                              Point;
  
  Cls3_2 T;
  Vertex_handle v0 = T.insert(Point(0,0,0));
  Vertex_handle v1 = T.insert(Point(1,0,0));
  Vertex_handle v2 = T.insert(Point(0,1,0));
  Vertex_handle v3 = T.insert(Point(0,0,1));
  Vertex_handle v4 = T.insert(Point(2,2,2));
  Vertex_handle v5 = T.insert(Point(-1,0,1));
  v0->info() = 0;
  v1->info() = 1;
  v2->info() = 2;
  v3->info() = 3;
  v4->info() = 4;
  v5->info() = 5;
  
  int i=0;
  for(auto it = T.cells_begin();
      it != T.cells_end(); ++it)
  {
    it->info() = i++;
  }
  
  
  std::ofstream ofs("triangulation_output");
  ofs << T;
  Cls3_2 T2;
  std::ifstream ifs("triangulation_output");
  if(!ifs)
    return 1;
  ifs >> T2;
  auto vit = T.vertices_begin();
  auto cit = T.cells_begin();
  auto vit2 = T2.vertices_begin();
  vit++; //first vertex is infinite
  vit2++; //first vertex is infinite
  for(;vit2 != T2.vertices_end(); ++vit2)
  {
    CGAL_assertion(vit2->info() == vit->info());
    ++vit;
  }

  for(auto it2 = T2.cells_begin();it2 != T2.cells_end(); ++it2)
  {
    CGAL_assertion(it2->info() == cit->info());
    ++cit;
  }
  
  return 0;
}
