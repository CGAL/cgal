// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_2_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
#define CGAL_PERIODIC_2_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H

#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Periodic_2_hyperbolic_triangulation_dummy.h>

namespace CGAL {
  
template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                         Triangulation_vertex_base_2<Gt> > >
class Periodic_2_Delaunay_hyperbolic_triangulation_2 : public Delaunay_triangulation_2<Gt,Tds>
{
public:
  typedef Periodic_2_Delaunay_hyperbolic_triangulation_2<Gt, Tds> Self;
  typedef Delaunay_triangulation_2<Gt,Tds> Base;
  
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Face_handle   Face_handle;
  typedef typename Base::Edge          Edge;
  
  typedef typename Base::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Base::Face_circulator Face_circulator;
  
  typedef Gt Geom_traits;
  typedef typename Geom_traits::FT            FT;
  typedef typename Geom_traits::Point_2       Point_2;
  typedef typename Geom_traits::Segment_2     Segment;
  
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::tds;
#endif
  
  Periodic_2_Delaunay_hyperbolic_triangulation_2(const Gt& gt = Gt())
    : Delaunay_triangulation_2<Gt,Tds>(gt) 
  {
    insert_dummy_points();
  }
  
  Periodic_2_Delaunay_hyperbolic_triangulation_2(
    const Periodic_2_Delaunay_hyperbolic_triangulation_2<Gt,Tds> &tr)
    : Delaunay_triangulation_2<Gt,Tds>(tr)
  { CGAL_triangulation_postcondition( this->is_valid() ); }
  
  /*vector<Vertex_handle>*/void insert_dummy_points();
  
  Object
  dual(const Finite_edges_iterator& ei) const {
    return this->dual(*ei);
  }
  
  Object
  dual(const Edge &e) const {
    // default implementation
    assert(false);
    return make_object(e);
  }
  
  // Clears the triangulation and initializes it again.
  void clear() {
    tds().clear();
    init_tds();
  }
  
private:
  // Initializes the triangulation data structure
  void init_tds() {
    this->_infinite_vertex = tds().insert_first();
  }
  
  void paste_together_opposite_sides(const std::vector<Vertex_handle>& on_vertex, const std::vector<Vertex_handle>& on_boundary);
};

} // namespace CGAL

#include <CGAL/Periodic_2_hyperbolic_triangulation_dummy.h>

#endif // CGAL_PERIODIC_2_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
