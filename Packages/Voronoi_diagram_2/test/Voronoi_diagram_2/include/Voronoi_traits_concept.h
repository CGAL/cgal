// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_TRAITS_CONCEPT_H
#define CGAL_VORONOI_TRAITS_CONCEPT_H 1

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>

CGAL_BEGIN_NAMESPACE

//=========================================================================

template<class DG>
class Voronoi_traits_concept
{
 public:
  typedef DG                               Delaunay_graph;

  typedef CGAL::Object    Object;
  typedef Tag_false       Has_nearest_site_2;
  typedef Tag_false       Has_get_conflicts;
  typedef Tag_false       Has_insert;
  typedef Tag_false       Has_remove;

  struct Edge_degeneracy_tester
  {
    typedef DG  Delaunay_graph;

    typedef typename DG::Edge                   Edge;
    typedef typename DG::Face_handle            Face_handle;
    typedef typename DG::Edge_circulator        Edge_circulator;
    typedef typename DG::All_edges_iterator     All_edges_iterator;
    typedef typename DG::Finite_edges_iterator  Finite_edges_iterator;

    typedef bool           result_type;
    typedef Arity_tag<2>   Arity;

    bool operator()(const Delaunay_graph&, const Edge&) const {
      return false;
    }

    bool operator()(const Delaunay_graph&, const Face_handle&, int) const {
      return false;
    }

    bool operator()(const Delaunay_graph&, const Edge_circulator&) const {
      return false;
    }

    bool operator()(const Delaunay_graph&, const All_edges_iterator&) const {
      return false;
    } 
  };

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DG>
  Face_degeneracy_tester;

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

  typedef typename DG::Point_2        Point_2;
  typedef typename DG::Site_2         Site_2;
  typedef typename DG::Vertex_handle  Vertex_handle;
  typedef typename DG::Face_handle    Face_handle;

  struct Get_site_2 {
    typedef Site_2         result_type;
    typedef Vertex_handle  value_type;
    typedef Arity_tag<1>   Arity;

    result_type operator()(const Vertex_handle& v) const {
      return Site_2();
    }
  };

  Get_site_2 get_site_2_object() const { return Get_site_2(); }

  struct Get_point_2 {
    typedef Point_2       result_type;
    typedef Face_handle   value_type;
    typedef Arity_tag<1>  Arity;

    Point_2 operator()(const Face_handle& f) const {
      return Point_2();
    }
  };

  Get_point_2 get_point_2_object() const { return Get_point_2(); }

 private:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_TRAITS_CONCEPT_H
