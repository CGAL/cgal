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

#ifndef CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Locate_result.h>
#include <cstdlib>
#include <algorithm>
#include <CGAL/Triangulation_utils_2.h>

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS
#include <CGAL/Voronoi_diagram_2/Identity_Voronoi_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class DG>
class AG_Nearest_site_2
{
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Locate_result_accessor<DG,false> Accessor;

 public:
  typedef DG                                          Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Edge               Edge;
  typedef typename Delaunay_graph::Point_2            Point_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Locate_result<DG,false> Query_result;

  typedef Arity_tag<2>    Arity;
  typedef Query_result   return_type;

 private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;
  typedef typename Delaunay_graph::Site_2             Site_2;

 public:
  Query_result operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    typename DG::Geom_traits::Oriented_side_of_bisector_2 side_of_bisector =
      dg.geom_traits().oriented_side_of_bisector_2_object();

    Vertex_handle v = dg.nearest_neighbor(p);
    if ( dg.dimension() == 0 ) {
      return Accessor::make_locate_result(v);
    }

    if ( dg.dimension() == 1 ) {
      Edge e = *dg.finite_edges_begin();
      Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
      Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

      Oriented_side os = side_of_bisector(v1->site(), v2->site(), p);
      
      if ( os == ON_ORIENTED_BOUNDARY ) {
	return Accessor::make_locate_result(e);
      } else {
	return Accessor::make_locate_result(v);
      }
    }

    CGAL_assertion( dg.dimension() == 2 );

    typename DG::Face_circulator fc_start = dg.incident_faces(v);
    typename DG::Face_circulator fc = fc_start;

    // first check if the point lies on a Voronoi vertex
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex(CW_CCW_2::ccw(index));
      Vertex_handle v2 = fc->vertex(CW_CCW_2::cw(index) );

      Oriented_side os1 = ON_POSITIVE_SIDE, os2 = ON_POSITIVE_SIDE;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), p);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(v->site(), v2->site(), p);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );
      CGAL_assertion( os2 != ON_NEGATIVE_SIDE );

      if ( os1 == ON_ORIENTED_BOUNDARY && os2 == ON_ORIENTED_BOUNDARY ) {
	Face_handle f(fc);
	return Accessor::make_locate_result(f);
      }

      ++fc;
    } while ( fc != fc_start );

    // now check if the point lies on a Voronoi edge
    fc_start = dg.incident_faces(v);
    fc = fc_start;
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex(CW_CCW_2::ccw(index));
      Vertex_handle v2 = fc->vertex(CW_CCW_2::cw(index) );

      Oriented_side os1 = ON_POSITIVE_SIDE, os2 = ON_POSITIVE_SIDE;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), p);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(v->site(), v2->site(), p);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );
      CGAL_assertion( os2 != ON_NEGATIVE_SIDE );
      CGAL_assertion( os1 != ON_ORIENTED_BOUNDARY ||
		      os2 != ON_ORIENTED_BOUNDARY );

      if ( os1 == ON_ORIENTED_BOUNDARY ) {
	Face_handle f(fc);
	Edge e(f, CW_CCW_2::cw(index));
	return Accessor::make_locate_result(e);
      } else if ( os2 == ON_ORIENTED_BOUNDARY ) {
	Face_handle f(fc);
	Edge e(f, CW_CCW_2::ccw(index));
	return Accessor::make_locate_result(e);
      }

      ++fc;
    } while ( fc != fc_start );

    return Accessor::make_locate_result(v);
  }
};


//=========================================================================
//=========================================================================

template<class DG>
class AG_Edge_degeneracy_tester
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                       Delaunay_graph;

  typedef typename DG::Edge                        Edge;
  typedef typename DG::Face_handle                 Face_handle;
  typedef typename DG::Edge_circulator             Edge_circulator;
  typedef typename DG::All_edges_iterator          All_edges_iterator;
  typedef typename DG::Finite_edges_iterator       Finite_edges_iterator;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

 private:
  typedef Triangulation_cw_ccw_2                     CW_CCW_2;

  typedef AG_Edge_degeneracy_tester<Delaunay_graph>  Self;

  typedef typename Delaunay_graph::Geom_traits       Geom_traits;

  typedef typename Delaunay_graph::Vertex_handle     Vertex_handle;

  typedef typename Geom_traits::Site_2               Site_2;

 public:
  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const
  {
    if ( dual.dimension() == 1 ) { return false; }

    if ( dual.is_infinite(f, i) ) { return false; }

    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = dual.data_structure().mirror_vertex(f, i);

    if ( dual.is_infinite(v3) || dual.is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( CW_CCW_2::ccw(i) );
    Vertex_handle v2 = f->vertex( CW_CCW_2::cw(i) );

    Site_2 s1 = v1->site();
    Site_2 s2 = v2->site();
    Site_2 s3 = v3->site();
    Site_2 s4 = v4->site();
    return dual.geom_traits().is_degenerate_edge_2_object()(s1,s2,s3,s4);
  }

  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    return operator()(dual, e.first, e.second);
  }

  bool operator()(const Delaunay_graph& dual,
		  const All_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Delaunay_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Delaunay_graph& dual,
		  const Edge_circulator& ec) const {
    return operator()(dual, *ec);
  }
};

//=========================================================================
//=========================================================================

template<class AG2>
class Apollonius_graph_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <AG2, AG_Edge_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   AG_Nearest_site_2<AG2> >
{
 private:
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;
  typedef AG_Nearest_site_2<AG2>                      AG_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>
  Face_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <AG2,Edge_tester,Face_tester,AG_Nearest_site_2>
  Base;

  typedef Apollonius_graph_Voronoi_traits_2<AG2>  Self;

 public:
  typedef typename AG2::Point_2                   Point_2;
  typedef typename AG2::Site_2                    Site_2;
  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Face_handle              Face_handle;

  typedef Tag_true                                Has_get_conflicts;
  typedef Tag_true                                Has_insert;
  typedef Tag_true                                Has_remove;

  struct Get_site_2 {
    typedef const Site_2&   result_type;
    typedef Vertex_handle   value_type;
    typedef Arity_tag<1>    Arity;

    result_type operator()(const Vertex_handle& v) const {
      return v->site();
    }
  };

  Get_site_2 get_site_2_object() const { return Get_site_2(); }

  struct Get_point_2 {
    typedef Point_2       result_type;
    typedef Face_handle   value_type;
    typedef Arity_tag<1>  Arity;

    Point_2 operator()(const Face_handle& f) const {
      typedef typename Base::Delaunay_graph::Geom_traits Geom_traits;

      return Geom_traits().construct_Apollonius_vertex_2_object()
	(f->vertex(0)->site(), f->vertex(1)->site(), f->vertex(2)->site());
    }
  };

  Get_point_2 get_point_2_object() const { return Get_point_2(); }
};


//=========================================================================
//=========================================================================

template<class AG2>
class Apollonius_graph_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <AG2, AG_Edge_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   AG_Nearest_site_2<AG2> >
{
 private:
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>
  Face_tester;

  typedef AG_Nearest_site_2<AG2>                      AG_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <AG2,Edge_tester,Face_tester,AG_Nearest_site_2>
  Base;

  typedef Apollonius_graph_cached_Voronoi_traits_2<AG2>  Self;
};

//=========================================================================
//=========================================================================

template<class AG2>
class Apollonius_graph_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <AG2, AG_Edge_degeneracy_tester<AG2>, 
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   AG_Nearest_site_2<AG2> >
{
 private:
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>
  Face_tester;

  typedef AG_Nearest_site_2<AG2>                      AG_Nearest_site_2;

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <AG2,Edge_tester,Face_tester,AG_Nearest_site_2>
  Base;

  typedef Apollonius_graph_ref_counted_Voronoi_traits_2<AG2>  Self;
};

//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt,class Agds,class LTag> class Apollonius_graph_2;
template<class Gt,class Agds,class LTag> class Apollonius_graph_hierarchy_2;

template<class Gt, class Agds, class LTag>
class Identity_Voronoi_traits_2< Apollonius_graph_2<Gt,Agds,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Apollonius_graph_2<Gt,Agds,LTag>,
    Apollonius_graph_Voronoi_traits_2< Apollonius_graph_2<Gt,Agds,LTag> >
  >
{};

template<class Gt, class Agds, class LTag>
class Identity_Voronoi_traits_2< Apollonius_graph_hierarchy_2<Gt,Agds,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Apollonius_graph_hierarchy_2<Gt,Agds,LTag>,
    Apollonius_graph_Voronoi_traits_2
    < Apollonius_graph_hierarchy_2<Gt,Agds,LTag> >
  >
{};


#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
