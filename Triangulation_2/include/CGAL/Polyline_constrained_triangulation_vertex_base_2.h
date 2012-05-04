// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Philipp Moeller

#ifndef CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_VERTEX_BASE_2_H


namespace CGAL {

template<class Vb>
class Polyline_constraint_hierarchy_vertex_base 
  : public Vb
{
  typedef Vb                                            Base;
  typedef typename Base::Triangulation_data_structure   Tds;
public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef Triangulation_hierarchy_vertex_base_2<Vb2>         Other;
  };

  Polyline_constraint_hierarchy_vertex_base() 
    : Base(), fixed(false), id(-1), cost(-1.0) 
  {}
  
  bool fixed;
  int id;
  double cost;
};

} // CGAL

#endif /* CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_VERTEX_BASE_2_H */
