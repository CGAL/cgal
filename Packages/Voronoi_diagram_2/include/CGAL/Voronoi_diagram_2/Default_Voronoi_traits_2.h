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

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Cached_degeneracy_testers.h>
#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================
#if 0
template<class Vertex_handle, class Voronoi_vertex,
	 class Voronoi_edge, bool Use_site>
class Vertex_edge_maker;

template<class Vertex_handle_t, class Voronoi_vertex, class Voronoi_edge>
struct Vertex_edge_maker<Vertex_handle_t,Voronoi_vertex,Voronoi_edge,true>
{
  typedef Vertex_handle_t   Vertex_handle;
  typedef Voronoi_vertex    Voronoi_vertex_2;
  typedef Voronoi_edge      Voronoi_edge_2;

  static Voronoi_vertex_2 make_vertex(const Vertex_handle& v1,
				      const Vertex_handle& v2,
				      const Vertex_handle& v3) {
    Voronoi_vertex_2 vv;
    vv.set_sites(v1->site(), v2->site(), v3->site());
    return vv;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->site(), v2->site());
    return ve;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2,
				  const Vertex_handle& v3,
				  bool is_src) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->site(), v2->site(), v3->site(), is_src);
    return ve;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2,
				  const Vertex_handle& v3,
				  const Vertex_handle& v4) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->site(), v2->site(), v3->site(), v4->site());
    return ve;
  }
};
#endif

//=========================================================================

template<class DG, class ET, class FT, class PL>
class Default_Voronoi_traits_2
{
 private:
  typedef Default_Voronoi_traits_2<DG,ET,FT,PL>   Self;

 public:
  typedef DG  Dual_graph;
  typedef ET  Edge_degeneracy_tester;
  typedef FT  Face_degeneracy_tester;
  typedef PL  Point_locator;

  typedef Tag_true Has_point_locator;

  typedef typename Dual_graph::Vertex_handle  Vertex_handle;

  typedef CGAL::Object  Object;

  struct Construct_object_object {
    template<class T>
    Object operator()(const T& t) const {
      return CGAL::make_object(t);
    }
  };

  struct Assign {
    template<class T>
    bool operator()(T& t, const Object& o) const {
      return CGAL::assign(t, o);
    }
  };

  Assign assign_object() const {
    return Assign();
  }

  Construct_object_object construct_object_object() const {
    return Construct_object_object();
  }

  Point_locator point_locator_object() const {
    return Point_locator();
  }

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


template<class DG, class ET, class FT>
class Default_Voronoi_traits_2<DG,ET,FT,Tag_false>
{
 private:
  typedef Default_Voronoi_traits_2<DG,ET,FT,Tag_false>   Self;

 public:
  typedef DG  Dual_graph;
  typedef ET  Edge_degeneracy_tester;
  typedef FT  Face_degeneracy_tester;

  typedef CGAL::Object  Object;

  typedef Tag_false Has_point_locator;

  typedef typename Dual_graph::Vertex_handle  Vertex_handle;

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};




//=========================================================================
//=========================================================================

template<class DG, class ET, class FT>
class Default_cached_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_cached_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

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

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================


template<class DG, class ET, class FT>
class Default_ref_counted_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_ref_counted_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

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

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
