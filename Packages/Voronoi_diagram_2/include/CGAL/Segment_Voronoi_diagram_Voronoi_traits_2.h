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

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Locate_result.h>
#include <cstdlib>
#include <algorithm>

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS
#include <CGAL/Voronoi_diagram_2/Identity_Voronoi_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class DG>
class SVD_Nearest_site_2
{
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Locate_result_accessor<DG,false> Accessor;

 public:
  typedef DG                                          Delaunay_graph;

  typedef typename Delaunay_graph::Point_2            Point_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Locate_result<DG,false> Query_result;

  typedef Arity_tag<2>    Arity;
  typedef Query_result    return_type;

 private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;
  typedef typename Delaunay_graph::Site_2             Site_2;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Edge               Edge;

 public:
  Query_result operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    typename DG::Geom_traits::Oriented_side_of_bisector_2 side_of_bisector =
      dg.geom_traits().oriented_side_of_bisector_2_object();

    typename DG::Geom_traits::Equal_2 is_equal =
      dg.geom_traits().equal_2_object();

    Site_2 sp = Site_2::construct_site_2(p);

    Vertex_handle v = dg.nearest_neighbor(p);
    if ( dg.dimension() == 0 ) {
      return Accessor::make_locate_result(v);
    }

    if ( dg.dimension() == 1 ) {
      Edge e = *dg.finite_edges_begin();
      Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
      Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

      Oriented_side os = side_of_bisector(v1->site(), v2->site(), sp);
      
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

      // check if the query point is identical to an endpoint of a
      // segment that has a Voronoi face with zero area.
      if ( !dg.is_infinite(v1) && !dg.is_infinite(v2) ) {
	if ( v->is_point() && is_equal(v->site(), sp) &&
	     v1->is_segment() && v2->is_segment() ) {
	  bool b1 =
	    is_equal(v->site(), v1->site().source_site()) ||
	    is_equal(v->site(), v1->site().target_site());
	  bool b2 =
	    is_equal(v->site(), v2->site().source_site()) ||
	    is_equal(v->site(), v2->site().target_site());

	  if ( b1 && b2 ) { 
	    Face_handle f(fc);
	    return Accessor::make_locate_result(f);
	  }
	}
      }

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), sp);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(v->site(), v2->site(), sp);
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

      // check if the query point is identical to an endpoint of a
      // segment that has a Voronoi face with zero area.
      if ( !dg.is_infinite(v1) && !dg.is_infinite(v2) ) {
	if ( v->is_point() && is_equal(v->site(), sp) &&
	     v1->is_segment() && v2->is_segment() ) {
	  bool b1 =
	    is_equal(v->site(), v1->site().source_site()) ||
	    is_equal(v->site(), v1->site().target_site());
	  bool b2 =
	    is_equal(v->site(), v2->site().source_site()) ||
	    is_equal(v->site(), v2->site().target_site());

	  CGAL_assertion( !b1 || !b2 );

	  if ( b1 ) {
	    Face_handle f(fc);
	    Edge e(f, CW_CCW_2::cw(index));
	    return Accessor::make_locate_result(e);
	  } else if ( b2 ) {
	    Face_handle f(fc);
	    Edge e(f, CW_CCW_2::ccw(index));
	    return Accessor::make_locate_result(e);
	  }
	}
      }

      // check if the query point is lies on the bisector between a
      // segment and its endpoint
      if ( !dg.is_infinite(v1) ) {
	if ( v->is_point() && is_equal(v->site(), sp) && v1->is_segment() ) {
	  bool b =
	    is_equal(v->site(), v1->site().source_site()) ||
	    is_equal(v->site(), v1->site().target_site());

	  if ( b ) {
	    Face_handle f(fc);
	    Edge e(f, CW_CCW_2::cw(index));
	    return Accessor::make_locate_result(e);
	  }
	}
      }

      if ( !dg.is_infinite(v2) ) {
	if ( v->is_point() && is_equal(v->site(), sp) && v2->is_segment() ) {
	  bool b =
	    is_equal(v->site(), v2->site().source_site()) ||
	    is_equal(v->site(), v2->site().target_site());

	  if ( b ) {
	    Face_handle f(fc);
	    Edge e(f, CW_CCW_2::ccw(index));
	    return Accessor::make_locate_result(e);
	  }
	}
      }

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), sp);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(v->site(), v2->site(), sp);
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
class SVD_Edge_degeneracy_tester
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                           Delaunay_graph;

  typedef typename Delaunay_graph::Edge                Edge;
  typedef typename Delaunay_graph::Face_handle         Face_handle;
  typedef typename Delaunay_graph::Edge_circulator     Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator  All_edges_iterator;

  typedef typename Delaunay_graph::Finite_edges_iterator
  Finite_edges_iterator;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

 private:
  typedef SVD_Edge_degeneracy_tester<Delaunay_graph>   Self;

  typedef typename Delaunay_graph::Geom_traits         Geom_traits;

  typedef typename Delaunay_graph::Vertex_handle       Vertex_handle;

  typedef typename Delaunay_graph::Site_2              Site_2;

  typedef typename Geom_traits::Equal_2                Equal_2;
  typedef typename Geom_traits::Orientation_2          Orientation_2;

 private:
  bool is_degenerate_infinite_edge(const Delaunay_graph& dual,
				   const Face_handle& f, int i) const
  {
    CGAL_precondition( dual.is_infinite(f, i) );

    Vertex_handle v = f->vertex( dual.ccw(i) );
    Vertex_handle v_inf = f->vertex( dual.cw(i) );

    if ( dual.is_infinite(v) ) {
      std::swap(v, v_inf);
    }

    if ( v->storage_site().is_segment() ) { return false; }

    Vertex_handle vv[2];

    vv[0] = f->vertex(i);
    vv[1] = dual.data_structure().mirror_vertex(f, i);

    if ( vv[0] == vv[1] ) { return false; }

    Site_2 s_end[2];
    for (int i = 0; i < 2; i++) {
      if ( vv[i]->storage_site().is_point() ) { return false; }

      Equal_2 are_equal = dual.geom_traits().equal_2_object();
      if ( !are_equal(v->site(),vv[i]->site().source_site()) &&
	   !are_equal(v->site(),vv[i]->site().target_site()) ) {
	return false;
      }
      if ( are_equal(v->site(),vv[i]->site().source_site()) ) {
	s_end[i] = vv[i]->site().target_site();
      } else {
	s_end[i] = vv[i]->site().source_site();
      }
    }

    Orientation_2 orientation = dual.geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }

 public:
  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const
  {
    if ( dual.dimension() == 1 ) { return false; }

    if ( dual.is_infinite(f, i) ) {
      return is_degenerate_infinite_edge(dual, f, i);
    }

    Vertex_handle v3 = f->vertex(i);
    Vertex_handle v4 = dual.data_structure().mirror_vertex(f, i);

    if ( dual.is_infinite(v3) || dual.is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( dual.ccw(i) );
    Vertex_handle v2 = f->vertex( dual.cw(i) );

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

template<class DG>
class SVD_Face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef DG                                       Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle   Vertex_handle;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

 private:
  typedef SVD_Face_degeneracy_tester<Delaunay_graph>  Self;
  typedef SVD_Edge_degeneracy_tester<Delaunay_graph>  Edge_tester;

  typedef typename Delaunay_graph::Geom_traits      Geom_traits;
  typedef typename Delaunay_graph::Edge             Edge;
  typedef typename Delaunay_graph::Edge_circulator  Edge_circulator;
  typedef typename Delaunay_graph::Face_handle      Face_handle;
  typedef typename Delaunay_graph::size_type        size_type;
  typedef typename Delaunay_graph::Site_2           Site_2;

  typedef typename Geom_traits::Equal_2             Equal_2;
  typedef typename Geom_traits::Orientation_2       Orientation_2;

 public:
  bool operator()(const Delaunay_graph& dual, const Vertex_handle& v) const
  {
    if ( dual.dimension() < 2 ) { return false; }

    if ( dual.is_infinite(v) ) { return false; }

    // THIS TEST NEEDS TO USE GEOMETRY; I CANNOT DO IT IN AN ENTIRELY
    // COMBINATORIAL MANNER

    // SEGMENT SPECIFIC TEST
    if ( v->site().is_segment() ) { return false; }

    // THIS WORKS ONLY FOR SEGMENTS (OR MAYBE NOT...)
    Edge_circulator ec_start(v);
    Edge_circulator ec = ec_start;
    size_type deg = 0;       // vertex degree
    size_type n_degen = 0;   // number of degenerate/non-infinite edges
    size_type n_inf = 0;     // number of infinite edges
    // number of non-degenerate/non-infinite edges
    size_type n_non_degen = 0;
      
    Edge e[2];
    Edge_tester e_tester;
    do {
      if ( e_tester(dual, ec) ) { ++n_degen; }
      else if ( dual.is_infinite(ec) ) { ++n_inf; }
      else { 
	if ( !dual.is_infinite(ec) ) {
	  if ( n_non_degen < 2 ) {
	    e[n_non_degen] = *ec;
	  }
	  n_non_degen++;
	}
      }
      deg++;
      ++ec;
    } while ( ec != ec_start );

    if ( deg == n_degen ) { return true; }
    if ( n_non_degen != 2 ) { return false; }

    Vertex_handle vv[2];
    Site_2 s_end[2];
    for (int i = 0; i < 2; i++) {
      CGAL_assertion( !dual.is_infinite(e[i]) );
      CGAL_expensive_assertion( !e_tester(dual, e[i]) );

      Vertex_handle v1 = e[i].first->vertex( dual.ccw(e[i].second) );
      Vertex_handle v2 = e[i].first->vertex( dual.cw(e[i].second) );
      vv[i] = (v1 == v) ? v2 : v1;

      CGAL_assertion( v == v1 || v == v2 );

      if ( vv[i]->storage_site().is_point() ) { return false; }

      Equal_2 are_equal = dual.geom_traits().equal_2_object();
      if ( !are_equal(v->site(),vv[i]->site().source_site()) &&
	   !are_equal(v->site(),vv[i]->site().target_site()) ) {
	return false;
      }
      if ( are_equal(v->site(),vv[i]->site().source_site()) ) {
	s_end[i] = vv[i]->site().target_site();
      } else {
	s_end[i] = vv[i]->site().source_site();
      }
    }

    Orientation_2 orientation = dual.geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }

};


//=========================================================================
//=========================================================================

template<class SVD2>
class Segment_Voronoi_diagram_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <SVD2, SVD_Edge_degeneracy_tester<SVD2>, SVD_Face_degeneracy_tester<SVD2>,
   SVD_Nearest_site_2<SVD2> >
{
 private:
  typedef SVD_Edge_degeneracy_tester<SVD2>              Edge_tester;
  typedef SVD_Face_degeneracy_tester<SVD2>              Face_tester;
  typedef SVD_Nearest_site_2<SVD2>                      SVD_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <SVD2,Edge_tester,Face_tester,SVD_Nearest_site_2>
  Base;

  typedef Segment_Voronoi_diagram_Voronoi_traits_2<SVD2>  Self;

 public:
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;
  typedef typename SVD2::Vertex_handle             Vertex_handle;
  typedef typename SVD2::Face_handle               Face_handle;

  typedef Tag_false                                Has_get_conflicts;
  typedef Tag_false                                Has_insert;
  typedef Tag_false                                Has_remove;

  struct Get_site_2 {
    typedef Site_2          result_type;
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

      return Geom_traits().construct_svd_vertex_2_object()
	(f->vertex(0)->site(), f->vertex(1)->site(), f->vertex(2)->site());
    }
  };

  Get_point_2 get_point_2_object() const { return Get_point_2(); }
};



//=========================================================================
//=========================================================================

template<class SVD2>
class Segment_Voronoi_diagram_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <SVD2, SVD_Edge_degeneracy_tester<SVD2>, SVD_Face_degeneracy_tester<SVD2>,
   SVD_Nearest_site_2<SVD2> >
{
 private:
  typedef SVD_Edge_degeneracy_tester<SVD2>            Edge_tester;
  typedef SVD_Face_degeneracy_tester<SVD2>            Face_tester;
  typedef SVD_Nearest_site_2<SVD2>                    SVD_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <SVD2,Edge_tester,Face_tester,SVD_Nearest_site_2>
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;
};

//=========================================================================
//=========================================================================


template<class SVD2>
class Segment_Voronoi_diagram_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <SVD2, SVD_Edge_degeneracy_tester<SVD2>, SVD_Face_degeneracy_tester<SVD2>,
   SVD_Nearest_site_2<SVD2> >
{
 private:
  typedef SVD_Edge_degeneracy_tester<SVD2>             Edge_tester;
  typedef SVD_Face_degeneracy_tester<SVD2>             Face_tester;
  typedef SVD_Nearest_site_2<SVD2>                     SVD_Nearest_site_2;

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <SVD2,Edge_tester,Face_tester,SVD_Nearest_site_2>
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;
};

//=========================================================================
//=========================================================================


//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt, class STag, class DS, class LTag >
class Segment_Voronoi_diagram_hierarchy_2;

template<class Gt, class DS, class LTag>
class Segment_Voronoi_diagram_2;


template<class Gt, class DS, class LTag>
class Identity_Voronoi_traits_2< Segment_Voronoi_diagram_2<Gt,DS,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Segment_Voronoi_diagram_2<Gt,DS,LTag>,
    Segment_Voronoi_diagram_Voronoi_traits_2
    < Segment_Voronoi_diagram_2<Gt,DS,LTag> >
  >
{};

template<class Gt, class STag, class DS, class LTag>
class Identity_Voronoi_traits_2<
  Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>,
    Segment_Voronoi_diagram_Voronoi_traits_2
    < Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
  >
{};

#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
