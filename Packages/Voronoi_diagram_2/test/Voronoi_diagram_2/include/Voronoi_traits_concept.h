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
#include <CGAL/Voronoi_diagram_2/Voronoi_traits_functors.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>

CGAL_BEGIN_NAMESPACE

//=========================================================================

template<class DG>
class Voronoi_traits_concept
{
private:
  typedef Voronoi_traits_concept<DG>       Self;

public:
  typedef DG                               Delaunay_graph;

  typedef CGAL::Object    Object;
  typedef Tag_false       Has_nearest_site_2;
  typedef Tag_false       Has_site_inserter;
  typedef Tag_false       Has_remove;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Identity_edge_degeneracy_tester<DG>
  Edge_degeneracy_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_degeneracy_tester<DG>
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

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<Site_2,DG,Tag_false>
  Access_site_2;

  Access_site_2 access_site_2_object() const { return Access_site_2(); }

  struct Construct_dual_point_2 {
    typedef Point_2       result_type;
    typedef typename DG::Face_handle   Face_handle;
    typedef Arity_tag<1>  Arity;

    Point_2 operator()(const Face_handle& f) const {
      return Point_2();
    }
  };

  Construct_dual_point_2 construct_dual_point_2_object() const {
    return Construct_dual_point_2();
  }

  void clear() {
    e_tester_.clear();
    f_tester_.clear();
  }

  void swap(Self& other) {
    e_tester_.swap(other.e_tester_);
    f_tester_.swap(other.f_tester_);
  }

  bool is_valid() const {
    return e_tester_.is_valid() && f_tester_.is_valid();
  }

private:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_TRAITS_CONCEPT_H
