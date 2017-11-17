// Copyright (c) 2012  GeometryFactory Sarl (France)
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_COMBINATORIAL_MAP_FOR_COREFINEMENT_H
#define CGAL_INTERNAL_COMBINATORIAL_MAP_FOR_COREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Cell_attribute.h>
#include <set>

namespace CGAL{
  namespace internal_IOP{

template <class Polyhedron>
struct Volume_info{
  std::set<Polyhedron*> outside;
  std::set<Polyhedron*> inside;
  bool is_empty;
  Volume_info():is_empty(false){}
};

struct Volume_on_merge
{
  template <class Attribute>
  void operator() (Attribute& a1,const Attribute& a2) const
  {
    CGAL_assertion(!a1.info().is_empty && !a2.info().is_empty);
    std::copy(a2.info().outside.begin(),a2.info().outside.end(),std::inserter(a1.info().outside,a1.info().outside.begin()));
    std::copy(a2.info().inside.begin(),a2.info().inside.end(),std::inserter(a1.info().inside,a1.info().inside.begin()));
  }
};

#ifndef NDEBUG
struct Point_on_merge
{
  template <class Attribute>
  void operator() (Attribute& a1,const Attribute& a2) const
  {
    CGAL_assertion(a1.point()==a2.point() );
    CGAL_USE(a1); CGAL_USE(a2);
  }
};
#endif


template < class Refs, class T, class Point_,
           class Functor_on_merge_=CGAL::Null_functor,
           class Functor_on_split_=CGAL::Null_functor >
class My_cell_attribute_with_point :
  public CGAL::Cell_attribute_without_info<Refs, T, Functor_on_merge_, Functor_on_split_>
{
   Point_ mpoint;
public:
  typedef Point_            Point;
  typedef Functor_on_merge_ Functor_on_merge;
  typedef Functor_on_split_ Functor_on_split;

  My_cell_attribute_with_point(){}
  My_cell_attribute_with_point(const Point& apoint) : mpoint(apoint) {}
  Point& point()             { return mpoint; }
  const Point& point() const { return mpoint; }  
    
};

template <typename Traits_,class Polyhedron>
struct Item_with_points_and_volume_info
{
  static const unsigned int dimension = 3;
  static const unsigned int NB_MARKS = 32;

  template<class CMap>
  struct Dart_wrapper
  {
    typedef Traits_                    Traits;
    typedef typename Traits::FT        FT;
    typedef typename Traits::Point_3   Point;
    typedef typename Traits::Vector_3  Vector;
    #ifndef NDEBUG
    typedef My_cell_attribute_with_point<CMap,CGAL::Tag_true,Point,Point_on_merge>                 Vertex_attribute;
    #else
    typedef My_cell_attribute_with_point<CMap,CGAL::Tag_true,Point>                                Vertex_attribute;
    #endif
    typedef CGAL::Cell_attribute<CMap,Volume_info<Polyhedron>,CGAL::Tag_true,Volume_on_merge >     Volume_attribute;
    typedef CGAL::cpp11::tuple< Vertex_attribute,
                                void,
                                void,
                                Volume_attribute>    Attributes;
  };
};

} } //namespace CGAL::internal_IOP

#endif //CGAL_INTERNAL_COMBINATORIAL_MAP_FOR_COREFINEMENT_H
