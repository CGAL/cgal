#ifndef CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Projector_classes.h>
#include <cstdlib>
#include <algorithm>

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class DG>
class DT_Edge_degeneracy_tester
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
  typedef Arity_tag<1>   Arity;

 private:
  typedef DT_Edge_degeneracy_tester<Dual_graph>    Self;

  typedef typename Dual_graph::Geom_traits         Geom_traits;

  typedef typename Dual_graph::Vertex_handle       Vertex_handle;

  typedef typename Geom_traits::Point_2            Point_2;

 public:
  DT_Edge_degeneracy_tester(const Dual_graph* dual = NULL) : dual_(dual) {}

  bool operator()(const Face_handle& f, int i) const {
    if ( dual_->is_infinite(f, i) ) { return false; }

    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = dual_->tds().mirror_vertex(f, i);

    if ( dual_->is_infinite(v3) || dual_->is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( dual_->ccw(i) );
    Vertex_handle v2 = f->vertex( dual_->cw(i) );

    Point_2 p1 = v1->point();
    Point_2 p2 = v2->point();
    Point_2 p3 = v3->point();
    Point_2 p4 = v4->point();
    Oriented_side os =
      dual_->geom_traits().side_of_oriented_circle_2_object()(p1,p2,p3,p4);
    return os == ON_ORIENTED_BOUNDARY;
  }
 
  bool operator()(const Edge& e) const {
    return operator()(e.first, e.second);
  }

  bool operator()(const All_edges_iterator& eit) const {
    return operator()(*eit);
  }

  bool operator()(const Finite_edges_iterator& eit) const {
    return operator()(*eit);
  }

  bool operator()(const Edge_circulator& ec) const {
    return operator()(*ec);
  }

 private:
  const Dual_graph* dual_;
};

//=========================================================================
//=========================================================================

template<class DG, class Edge_tester>
class DT_Face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef DG                                   Dual_graph;
  typedef typename Dual_graph::Vertex_handle   Vertex_handle;
  typedef typename Dual_graph::Edge            Edge;
  typedef Edge_tester                          Edge_degeneracy_tester;

  typedef typename Dual_graph::Finite_vertices_iterator
  Finite_vertices_iterator;

  typedef typename Dual_graph::All_vertices_iterator  All_vertices_iterator;
  typedef typename Dual_graph::Vertex_circulator      Vertex_circulator;

  typedef bool           result_type;
  typedef Arity_tag<1>   Arity;

 public:
  DT_Face_degeneracy_tester(const Dual_graph* dual = NULL) {}

  DT_Face_degeneracy_tester(const Dual_graph* dual,
			    const Edge_degeneracy_tester* e_tester) {}

  bool operator()(const Vertex_handle& v) const {
    return false;
  }
 
  bool operator()(const Vertex_circulator& vc) const {
    return false;
  }

#ifndef CGAL_T2_USE_ITERATOR_AS_HANDLE
  bool operator()(const All_vertices_iterator& vit) const {
    return false;
  }
#endif

  bool operator()(const Finite_vertices_iterator& vit) const {
    return false;
  }
};


//=========================================================================
//=========================================================================


template<class DT2>
class Delaunay_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <DT2,
   DT_Edge_degeneracy_tester<DT2>,
   DT_Face_degeneracy_tester<DT2,DT_Edge_degeneracy_tester<DT2> >
  >
{
 public:
  typedef DT2                                         Dual_graph;
  typedef DT_Edge_degeneracy_tester<DT2>              Edge_tester;
  typedef DT_Face_degeneracy_tester<DT2,Edge_tester>  Face_tester;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <DT2,Edge_tester,Face_tester>
  Base;

  typedef Delaunay_triangulation_Voronoi_traits_2<DT2>  Self;

  typedef typename Dual_graph::Triangulation_data_structure
  Dual_graph_data_structure;

  Delaunay_triangulation_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->tds();
  }
};


//=========================================================================
//=========================================================================
#if 0
template<class DT2>
class Delaunay_triangulation_cached_Voronoi_traits_2
  : public Default_cached_Voronoi_traits_2
  <DT2, DT_Edge_degeneracy_tester<DT2>,
   DT_Face_degeneracy_tester
   <DT2,Cached_edge_degeneracy_tester<DT_Edge_degeneracy_tester<DT2>,
				      Tds_project<DT2> >
   >,
   Tds_project<DT2>
  >
{
 public:
  typedef DT2  Dual_graph;
  typedef Tds_project<DT2>                            TDS_Project;
  typedef DT_Edge_degeneracy_tester<DT2>              Edge_tester;

  typedef Cached_edge_degeneracy_tester<Edge_tester,TDS_Project>
  Cached_edge_tester;

  typedef DT_Face_degeneracy_tester<DT2,Cached_edge_tester>
  Face_tester;

  typedef
  Default_cached_Voronoi_traits_2<DT2,Edge_tester,Face_tester,TDS_Project>
  Base;

  typedef Delaunay_triangulation_cached_Voronoi_traits_2<DT2>  Self;

  typedef typename Dual_graph::Triangulation_data_structure
  Dual_graph_data_structure;

  Delaunay_triangulation_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->tds();
  }
};
#endif

//=========================================================================
//=========================================================================

template<class DT2>
class Delaunay_triangulation_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <DT2, DT_Edge_degeneracy_tester<DT2>,
   DT_Face_degeneracy_tester
   <DT2,CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
    <DT_Edge_degeneracy_tester<DT2>,
     CGAL_VORONOI_DIAGRAM_2_NS::Tds_project<DT2> >
   >,
   CGAL_VORONOI_DIAGRAM_2_NS::Tds_project<DT2>
  >
{
 public:
  typedef DT2                                          Dual_graph;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Tds_project<DT2>  TDS_Project;
  typedef DT_Edge_degeneracy_tester<DT2>               Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
  <Edge_tester,TDS_Project>
  Ref_counted_edge_tester;

  typedef DT_Face_degeneracy_tester<DT2,Ref_counted_edge_tester>
  Face_tester;

  typedef
  CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <DT2,Edge_tester,Face_tester,TDS_Project>
  Base;

  typedef Delaunay_triangulation_cached_Voronoi_traits_2<DT2>  Self;

  typedef typename Dual_graph::Triangulation_data_structure
  Dual_graph_data_structure;

  Delaunay_triangulation_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->tds();
  }
};


//=========================================================================
//=========================================================================


CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H
