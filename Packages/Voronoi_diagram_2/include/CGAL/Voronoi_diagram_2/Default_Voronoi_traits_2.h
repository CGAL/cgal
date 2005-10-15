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

#ifndef CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
#define CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Cached_degeneracy_testers.h>
#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class DG>
struct Default_face_degeneracy_tester
{
  // tests whether a face has zero area
  typedef DG                                      Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle  Vertex_handle;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

  bool operator()(const Delaunay_graph&, const Vertex_handle&) const {
    return false;
  }
};

//=========================================================================

struct Null_functor
{
  Null_functor() {}
  template<typename T> Null_functor(T t) {}
};

//=========================================================================

template<class Functor>
struct Functor_exists
{
  typedef Tag_true  Value;
};

template<>
struct Functor_exists<Null_functor>
{
  typedef Tag_false Value;
};

//=========================================================================

template<class DG, class ET, class FT, class AS, class CDP,
	 class SI, class NS>
class Default_Voronoi_traits_2
{
 private:
  typedef Default_Voronoi_traits_2<DG,ET,FT,AS,CDP,SI,NS>   Self;

 public:
  typedef DG   Delaunay_graph;
  typedef ET   Edge_degeneracy_tester;
  typedef FT   Face_degeneracy_tester;
  typedef AS   Access_site_2;
  typedef CDP  Construct_dual_point_2;

  typedef NS   Nearest_site_2;
  typedef SI   Site_inserter;

  typedef typename Functor_exists<Nearest_site_2>::Value  Has_nearest_site_2;
  typedef typename Functor_exists<Site_inserter>::Value   Has_site_inserter;

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

  Access_site_2 access_site_2_object() const {
    return Access_site_2();
  }

  Construct_dual_point_2 construct_dual_point_2_object() const {
    return Construct_dual_point_2();
  }

  Nearest_site_2 nearest_site_2_object() const {
    return Nearest_site_2();
  }

  Site_inserter site_inserter_object() const {
    return Site_inserter();
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================

template<class DG, class ET, class FT, class NS>
class Default_cached_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_cached_Voronoi_traits_2<DG,ET,FT,NS>  Self;

 public:
  typedef DG           Delaunay_graph;
  typedef NS           Nearest_site_2;

  typedef Cached_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  //  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester_base,
  //					Edge_degeneracy_tester>
  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester_base>
  Face_degeneracy_tester;


 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

  Nearest_site_2 nearest_site_2_object() const {
    return Nearest_site_2();
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================


template<class DG, class ET, class FT, class NS>
class Default_ref_counted_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_ref_counted_Voronoi_traits_2<DG,ET,FT,NS>  Self;

 public:
  typedef DG           Delaunay_graph;
  typedef NS           Nearest_site_2;

  typedef Ref_counted_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  //  typedef Ref_counted_face_degeneracy_tester<Face_degeneracy_tester_base,
  //					     Edge_degeneracy_tester>
  typedef Ref_counted_face_degeneracy_tester<Face_degeneracy_tester_base>
  Face_degeneracy_tester;

 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

  Nearest_site_2 nearest_site_2_object() const {
    return Nearest_site_2();
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
