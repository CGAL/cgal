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

#ifndef CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_vertex_base_2.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_edge_base_2.h>
#include <CGAL/Voronoi_diagram_2/Locate_result.h>
#include <cstdlib>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class DG>
class RT_Nearest_site_2
{
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Locate_result_accessor<DG,false> Accessor;

 public:
  typedef DG                                          Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Edge               Edge;

 private:
  typedef typename Delaunay_graph::Geom_traits        Geom_traits;

 public:
  typedef typename Geom_traits::Weighted_point_2      Weighted_point_2;
  typedef typename Geom_traits::Point_2               Point_2;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Locate_result<DG,false> Query_result;

  typedef Arity_tag<2>    Arity;
  typedef Query_result    return_type;

 private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;
  typedef typename Delaunay_graph::Vertex_circulator  Vertex_circulator;
  typedef typename Delaunay_graph::Face_circulator    Face_circulator;
  typedef typename Delaunay_graph::Edge_circulator    Edge_circulator;

 public:
  Query_result operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    typename DG::Geom_traits::Compare_power_distance_2 cmp_power_distance =
      dg.geom_traits().compare_power_distance_2_object();

    Vertex_handle v = dg.nearest_power_vertex(p);

    if ( dg.dimension() == 0 ) {
      return Accessor::make_locate_result(v);
    }

    if ( dg.dimension() == 1 ) {
      Edge_circulator ec = dg.incident_edges(v);
      Edge_circulator ec_start = ec;
      Comparison_result cr;

      do {
	Edge e = *ec;
	Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
	Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

	if ( v == v1 ) {
	  if ( !dg.is_infinite(v2) ) {
	    cr = cmp_power_distance(p, v2->point(), v->point());
	    CGAL_assertion( cr != SMALLER );
	    if ( cr == EQUAL ) {
	      return Accessor::make_locate_result( e );
	    }
	  }
	} else {
	  CGAL_assertion( v == v2 );
	  if ( !dg.is_infinite(v1) ) {
	    cr = cmp_power_distance(p, v1->point(), v->point());
	    CGAL_assertion( cr != SMALLER );
	    if ( cr == EQUAL ) {
	      return Accessor::make_locate_result( e );
	    }
	  }
	}
	++ec;
      } while ( ec != ec_start );

      return Accessor::make_locate_result(v);
    }

    CGAL_assertion( dg.dimension() == 2 );

    Face_circulator fc_start = dg.incident_faces(v);
    Face_circulator fc = fc_start;

    // first check if the point lies on a Voronoi vertex
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex(CW_CCW_2::ccw(index));
      Vertex_handle v2 = fc->vertex(CW_CCW_2::cw(index) );

      Comparison_result cr1 = LARGER, cr2 = LARGER;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	cr1 = cmp_power_distance(p, v1->point(), v->point());
      }
      if ( !dg.is_infinite(v2) ) {
	cr2 = cmp_power_distance(p, v2->point(), v->point());
      }

      CGAL_assertion( cr1 != SMALLER );
      CGAL_assertion( cr2 != SMALLER );

      if ( cr1 == EQUAL && cr2 == EQUAL ) {
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

      Comparison_result cr1 = LARGER, cr2 = LARGER;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	cr1 = cmp_power_distance(p, v1->point(), v->point());
      }
      if ( !dg.is_infinite(v2) ) {
	cr2 = cmp_power_distance(p, v2->point(), v->point());
      }

      CGAL_assertion( cr1 != SMALLER );
      CGAL_assertion( cr2 != SMALLER );
      CGAL_assertion( cr1 != EQUAL || cr2 != EQUAL );

      if ( cr1 == EQUAL ) {
	Face_handle f(fc);
	Edge e(f, CW_CCW_2::cw(index));
	return Accessor::make_locate_result(e);
      } else if ( cr2 == EQUAL ) {
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
class RT_Edge_degeneracy_tester
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
  typedef RT_Edge_degeneracy_tester<Delaunay_graph>  Self;

  typedef typename Delaunay_graph::Geom_traits       Geom_traits;

  typedef typename Delaunay_graph::Vertex_handle     Vertex_handle;

  typedef typename Geom_traits::Weighted_point_2     Site_2;

 public:
  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const
  {
    if ( dual.dimension() == 1 ) { return false; }

    if ( dual.is_infinite(f, i) ) { return false; }

    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = dual.tds().mirror_vertex(f, i);

    if ( dual.is_infinite(v3) || dual.is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( dual.ccw(i) );
    Vertex_handle v2 = f->vertex( dual.cw(i) );

    Site_2 s1 = v1->point();
    Site_2 s2 = v2->point();
    Site_2 s3 = v3->point();
    Site_2 s4 = v4->point();
    Oriented_side os =
      dual.geom_traits().power_test_2_object()(s1,s2,s3,s4);
    return os == ON_ORIENTED_BOUNDARY;
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

template<class DG> class Regular_triangulation_Voronoi_traits_2;
template<class DG> class RT_Voronoi_edge_2;

template<class DG>
class RT_Voronoi_vertex_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Voronoi_vertex_base_2
  <DG, typename DG::Geom_traits::Point_2,
   typename DG::Geom_traits::Weighted_point_2, RT_Voronoi_vertex_2<DG> >
{
  friend class Regular_triangulation_Voronoi_traits_2<DG>;
  friend class RT_Voronoi_edge_2<DG>;
#ifndef CGAL_CFG_NESTED_CLASS_FRIEND_DECLARATION_BUG
  friend class RT_Voronoi_edge_2<DG>::Base;
#else
  friend class
  CGAL_VORONOI_DIAGRAM_2_NS::Voronoi_edge_base_2
  <DG,
   typename DG::Geom_traits::Point_2,
   typename DG::Geom_traits::Weighted_point_2,
   RT_Voronoi_edge_2<DG>,
   RT_Voronoi_vertex_2<DG> >;
#endif

 private:
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Voronoi_vertex_base_2
  <DG, typename DG::Geom_traits::Point_2,
   typename DG::Geom_traits::Weighted_point_2,
   RT_Voronoi_vertex_2<DG> >
  Base;

 public:
  operator typename Base::Point_2() const {
    return
      typename Base::Geom_traits().construct_weighted_circumcenter_2_object()
      (this->s_[0], this->s_[1], this->s_[2]);
  }
};

//=========================================================================
  
template<class DG>
class RT_Voronoi_edge_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Voronoi_edge_base_2
  <DG, typename DG::Geom_traits::Point_2,
   typename DG::Geom_traits::Weighted_point_2,
   RT_Voronoi_edge_2<DG>, RT_Voronoi_vertex_2<DG> >
{
  friend class Regular_triangulation_Voronoi_traits_2<DG>;
  friend class RT_Voronoi_vertex_2<DG>;

 private:
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Voronoi_edge_base_2
  <DG, typename DG::Geom_traits::Point_2,
   typename DG::Geom_traits::Weighted_point_2,
   RT_Voronoi_edge_2<DG>, RT_Voronoi_vertex_2<DG> >
  Base;
};


//=========================================================================


template<class RT2>
class Regular_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <RT2, RT_Edge_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>,
   RT_Nearest_site_2<RT2> >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>              Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>
  Face_tester;

  typedef RT_Nearest_site_2<RT2>                      RT_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester,RT_Nearest_site_2>
  Base;

  typedef Regular_triangulation_Voronoi_traits_2<RT2>  Self;

 public:
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;
  typedef typename RT2::Vertex_handle                  Vertex_handle;

  typedef RT_Voronoi_vertex_2<RT2>                Voronoi_vertex_2;
  typedef RT_Voronoi_edge_2<RT2>                  Voronoi_edge_2;
  typedef Voronoi_edge_2                          Curve;

  typedef Tag_true                                Has_get_conflicts;
  typedef Tag_true                                Has_insert;

  static const Site_2& site(const Vertex_handle& v) {
    return v->point();
  }

  static Voronoi_vertex_2 make_vertex(const Vertex_handle& v1,
				      const Vertex_handle& v2,
				      const Vertex_handle& v3) {
    Voronoi_vertex_2 vv;
    vv.set_sites(v1->point(), v2->point(), v3->point());
    return vv;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->point(), v2->point());
    return ve;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2,
				  const Vertex_handle& v3,
				  bool is_src) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->point(), v2->point(), v3->point(), is_src);
    return ve;
  }

  static Voronoi_edge_2 make_edge(const Vertex_handle& v1,
				  const Vertex_handle& v2,
				  const Vertex_handle& v3,
				  const Vertex_handle& v4) {
    Voronoi_edge_2 ve;
    ve.set_sites(v1->point(), v2->point(), v3->point(), v4->point());
    return ve;
  }
};


//=========================================================================
//=========================================================================


template<class RT2>
class Regular_triangulation_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_cached_Voronoi_traits_2
  <RT2, RT_Edge_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>,
   RT_Nearest_site_2<RT2> >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>              Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>
  Face_tester;

  typedef RT_Nearest_site_2<RT2>                      RT_Nearest_site_2;

  typedef
  CGAL_VORONOI_DIAGRAM_2_NS::Default_cached_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester,RT_Nearest_site_2>
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

//=========================================================================
//=========================================================================

template<class RT2>
class Regular_triangulation_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <RT2, RT_Edge_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>,
   RT_Nearest_site_2<RT2> >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>               Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_face_degeneracy_tester<RT2>
  Face_tester;

  typedef RT_Nearest_site_2<RT2>                       RT_Nearest_site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester,RT_Nearest_site_2>
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

CGAL_END_NAMESPACE

//=========================================================================
//=========================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>
#include <CGAL/Hyperbola_2.h>

CGAL_BEGIN_NAMESPACE

template<class RT2>
Qt_widget& operator<<(Qt_widget& qt_w,
		      const RT_Voronoi_edge_2<RT2>& ve)
{
  typedef typename RT2::Geom_traits          Geom_traits;
  typedef typename Geom_traits::Assign_2     Assign_2;
  typedef typename Geom_traits::Point_2      Point_2;
  typedef typename Geom_traits::Segment_2    Segment_2;
  //  typedef typename Geom_traits::Ray_2        Ray_2;
  typedef typename Geom_traits::Line_2       Line_2;

  //  Assign_2 assign = Geom_traits().assign_2_object();
  //  Segment_2 s;
  //  Ray_2 r;
  //  Line_2 l;

  if ( ve.has_source() && ve.has_target() ) {
    typename Geom_traits::Construct_weighted_circumcenter_2 circumcenter;
    Point_2 c1 = circumcenter(ve.down(), ve.up(), ve.left());
    Point_2 c2 = circumcenter(ve.up(), ve.down(), ve.right());
    qt_w << Segment_2(c1, c2);
  } else if ( ve.has_source() && !ve.has_target() ) {
    typename Geom_traits::Construct_weighted_circumcenter_2 circumcenter;
    typename Geom_traits::Construct_radical_axis_2     c_bis;
    typename Geom_traits::Construct_ray_2          c_ray;
    Point_2 c = circumcenter(ve.down(), ve.up(), ve.left());
    Line_2 l = c_bis(ve.up(), ve.down());
    qt_w << c_ray(c, l);
  } else if ( !ve.has_source() && ve.has_target() ) {
    typename Geom_traits::Construct_weighted_circumcenter_2 circumcenter;
    typename Geom_traits::Construct_radical_axis_2     c_bis;
    typename Geom_traits::Construct_ray_2          c_ray;
    Point_2 c = circumcenter(ve.up(), ve.down(), ve.right());
    Line_2 l = c_bis(ve.down(), ve.up());
    qt_w << c_ray(c, l);
  } else {
    CGAL_assertion( !ve.has_source() && !ve.has_target() );
    typename Geom_traits::Construct_radical_axis_2 c_bis;
    qt_w << c_bis(ve.up(), ve.down());
  }

  return qt_w;
}



CGAL_END_NAMESPACE

#endif



#endif // CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
