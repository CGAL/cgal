// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_ADAPTATION_TRAITS_CONCEPT_H
#define CGAL_ADAPTATION_TRAITS_CONCEPT_H 1

#include <CGAL/tags.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>

namespace CGAL {

//=========================================================================

template<class DG>
class Adaptation_traits_concept
{
private:
  typedef Adaptation_traits_concept<DG>       Self;

public:
  typedef DG                               Delaunay_graph;

  typedef CGAL::Object    Object;
  typedef Tag_false       Has_nearest_site_2;
  typedef Tag_false       Has_site_inserter;
  typedef Tag_false       Has_remove;


  typedef typename DG::Point_2        Point_2;
  typedef typename DG::Site_2         Site_2;
  typedef typename DG::Vertex_handle  Delaunay_vertex_handle;
  typedef typename DG::Face_handle    Delaunay_face_handle;
  typedef typename DG::Edge           Delaunay_edge;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<Site_2,DG,Tag_false>
  Access_site_2;

  Access_site_2 access_site_2_object() const { return Access_site_2(); }

  struct Construct_Voronoi_point_2 {
    typedef Point_2       result_type;
    typedef typename DG::Face_handle   Face_handle;

    Point_2 operator()(const Face_handle& ) const {
      return Point_2();
    }
  };

  Construct_Voronoi_point_2 construct_Voronoi_point_2_object() const {
    return Construct_Voronoi_point_2();
  }

#if 0
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
#endif
};


} //namespace CGAL


#endif // CGAL_ADAPTATION_TRAITS_CONCEPT_H
