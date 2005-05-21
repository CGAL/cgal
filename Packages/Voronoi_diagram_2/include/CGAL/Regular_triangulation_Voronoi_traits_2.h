#ifndef CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Default_Voronoi_traits_2.h>
#include <cstdlib>
#include <algorithm>

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class DG>
class RT_Edge_degeneracy_tester
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                       Dual_graph;

  typedef typename Dual_graph::Edge                Edge;
  typedef typename Dual_graph::Face_handle         Face_handle;
  typedef typename Dual_graph::Edge_circulator     Edge_circulator;
  typedef typename Dual_graph::All_edges_iterator  All_edges_iterator;
  typedef typename Dual_graph::Finite_edges_iterator Finite_edges_iterator;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

 private:
  typedef RT_Edge_degeneracy_tester<Dual_graph>    Self;

  typedef typename Dual_graph::Geom_traits         Geom_traits;

  typedef typename Dual_graph::Vertex_handle       Vertex_handle;

  typedef typename Geom_traits::Weighted_point_2   Site_2;

 public:
  bool operator()(const Dual_graph& dual, const Face_handle& f, int i) const {
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

  bool operator()(const Dual_graph& dual, const Edge& e) const {
    return operator()(dual, e.first, e.second);
  }

  bool operator()(const Dual_graph& dual,
		  const All_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Dual_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Dual_graph& dual, const Edge_circulator& ec) const {
    return operator()(dual, *ec);
  }
};

//=========================================================================
//=========================================================================

template<class DG, class Edge_tester>
class RT_Face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef DG                                      Dual_graph;
  typedef typename Dual_graph::Vertex_handle      Vertex_handle;
  typedef typename Dual_graph::Vertex_circulator  Vertex_circulator;
  typedef typename Dual_graph::Edge               Edge;
  typedef Edge_tester                             Edge_degeneracy_tester;

  typedef typename Dual_graph::All_vertices_iterator
  All_vertices_iterator;

  typedef typename Dual_graph::Finite_vertices_iterator
  Finite_vertices_iterator;

  typedef bool           result_type;
  typedef Arity_tag<3>   Arity;

 public:
  template<class A>
  bool operator()(const Dual_graph&, const Edge_tester&, const A& a) const {
    return false;
  }
};


//=========================================================================
//=========================================================================


template<class RT2>
class Regular_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <RT2,
   RT_Edge_degeneracy_tester<RT2>,
   RT_Face_degeneracy_tester<RT2,RT_Edge_degeneracy_tester<RT2> >
  >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>              Edge_tester;
  typedef RT_Face_degeneracy_tester<RT2,Edge_tester>  Face_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester>
  Base;

  typedef Regular_triangulation_Voronoi_traits_2<RT2>  Self;
};


//=========================================================================
//=========================================================================

template<class RT2>
class Regular_triangulation_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_cached_Voronoi_traits_2
  <RT2, RT_Edge_degeneracy_tester<RT2>,
   RT_Face_degeneracy_tester
   <RT2,CGAL_VORONOI_DIAGRAM_2_NS::Cached_edge_degeneracy_tester
    <RT_Edge_degeneracy_tester<RT2> >
   >
  >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>              Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Cached_edge_degeneracy_tester<Edge_tester>
  Cached_edge_tester;

  typedef RT_Face_degeneracy_tester<RT2,Cached_edge_tester>
  Face_tester;

  typedef
  CGAL_VORONOI_DIAGRAM_2_NS::Default_cached_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester>
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

//=========================================================================
//=========================================================================

template<class RT2>
class Regular_triangulation_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <RT2, RT_Edge_degeneracy_tester<RT2>,
   RT_Face_degeneracy_tester
   <RT2,CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
    <RT_Edge_degeneracy_tester<RT2> >
   >
  >
{
 private:
  typedef RT_Edge_degeneracy_tester<RT2>               Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
  <Edge_tester>
  Ref_counted_edge_tester;

  typedef RT_Face_degeneracy_tester<RT2,Ref_counted_edge_tester>
  Face_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <RT2,Edge_tester,Face_tester>
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

//=========================================================================
//=========================================================================


CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
