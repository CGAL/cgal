#ifndef CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Projector_classes.h>
#include <cstdlib>
#include <algorithm>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class DG>
class AG_Edge_degeneracy_tester
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                       Dual_graph;

  typedef typename Dual_graph::Edge                Edge;
  typedef typename Dual_graph::Face_handle         Face_handle;
  typedef typename Dual_graph::Edge_circulator     Edge_circulator;
  typedef typename Dual_graph::All_edges_iterator  All_edges_iterator;
  typedef typename Dual_graph::Finite_edges_iterator Finite_edges_iterator;

 private:
  typedef Triangulation_cw_ccw_2                   CW_CCW_2;

  typedef AG_Edge_degeneracy_tester<Dual_graph>    Self;

  typedef typename Dual_graph::Geom_traits         Geom_traits;

  typedef typename Dual_graph::Vertex_handle       Vertex_handle;

  typedef typename Geom_traits::Site_2             Site_2;

 public:
  AG_Edge_degeneracy_tester(const Dual_graph* dual = NULL) : dual_(dual) {}

  bool operator()(const Face_handle& f, int i) const {
    if ( dual_->is_infinite(f, i) ) { return false; }

    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = dual_->data_structure().mirror_vertex(f, i);

    if ( dual_->is_infinite(v3) || dual_->is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( CW_CCW_2::ccw(i) );
    Vertex_handle v2 = f->vertex( CW_CCW_2::cw(i) );

    Site_2 s1 = v1->site();
    Site_2 s2 = v2->site();
    Site_2 s3 = v3->site();
    Site_2 s4 = v4->site();
    return
      dual_->geom_traits().is_degenerate_edge_2_object()(s1,s2,s3,s4);
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
class AG_Face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef DG                                   Dual_graph;
  typedef typename Dual_graph::Vertex_handle   Vertex_handle;
  typedef typename Dual_graph::Edge            Edge;
  typedef Edge_tester                          Edge_degeneracy_tester;

  typedef typename Dual_graph::Finite_vertices_iterator
  Finite_vertices_iterator;

 public:
  AG_Face_degeneracy_tester(const Dual_graph* dual = NULL) {}

  AG_Face_degeneracy_tester(const Dual_graph* dual,
			    const Edge_degeneracy_tester* e_tester) {}

  bool operator()(const Vertex_handle& v) const {
    return false;
  }
 
  bool operator()(const Finite_vertices_iterator& vit) const {
    return false;
  }
};


//=========================================================================
//=========================================================================


template<class AG2>
class Apollonius_graph_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <AG2,
   AG_Edge_degeneracy_tester<AG2>,
   AG_Face_degeneracy_tester<AG2,AG_Edge_degeneracy_tester<AG2> >
  >
{
 public:
  typedef AG2                                         Dual_graph;
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;
  typedef AG_Face_degeneracy_tester<AG2,Edge_tester>  Face_tester;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <AG2,Edge_tester,Face_tester>
  Base;

  typedef Apollonius_graph_Voronoi_traits_2<AG2>  Self;

  typedef typename Dual_graph::Data_structure  Dual_graph_data_structure;

  Apollonius_graph_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};


//=========================================================================
//=========================================================================

#if 0
template<class AG2>
class Apollonius_graph_cached_Voronoi_traits_2
  : public Default_cached_Voronoi_traits_2
  <AG2, AG_Edge_degeneracy_tester<AG2>,
   AG_Face_degeneracy_tester
   <AG2,Cached_edge_degeneracy_tester<AG_Edge_degeneracy_tester<AG2>,
				      Data_structure_project<AG2> >
   >,
   Data_structure_project<AG2>
  >
{
 public:
  typedef AG2  Dual_graph;
  typedef Data_structure_project<AG2>                 DG_Project;
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;

  typedef Cached_edge_degeneracy_tester<Edge_tester,DG_Project>
  Cached_edge_tester;

  typedef AG_Face_degeneracy_tester<AG2,Cached_edge_tester>
  Face_tester;

  typedef
  Default_cached_Voronoi_traits_2<AG2,Edge_tester,Face_tester,DG_Project>
  Base;

  typedef Apollonius_graph_cached_Voronoi_traits_2<AG2>  Self;

  typedef typename Dual_graph::Data_structure  Dual_graph_data_structure;

  Apollonius_graph_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};
#endif

//=========================================================================
//=========================================================================


template<class AG2>
class Apollonius_graph_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <AG2, AG_Edge_degeneracy_tester<AG2>,
   AG_Face_degeneracy_tester
   <AG2,CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
    <AG_Edge_degeneracy_tester<AG2>,
     CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<AG2> >
   >,
   CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<AG2>
  >
{
 public:
  typedef AG2                                         Dual_graph;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<AG2>  DG_Project;
  typedef AG_Edge_degeneracy_tester<AG2>              Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
  <Edge_tester,DG_Project>
  Ref_counted_edge_tester;

  typedef AG_Face_degeneracy_tester<AG2,Ref_counted_edge_tester>
  Face_tester;

  typedef
  CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <AG2,Edge_tester,Face_tester,DG_Project>
  Base;

  typedef Apollonius_graph_cached_Voronoi_traits_2<AG2>  Self;

  typedef typename Dual_graph::Data_structure  Dual_graph_data_structure;

  Apollonius_graph_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};


//=========================================================================
//=========================================================================


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
