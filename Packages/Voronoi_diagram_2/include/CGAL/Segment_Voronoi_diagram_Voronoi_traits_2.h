#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Projector_classes.h>
#include <cstdlib>
#include <algorithm>

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class DG>
class SVD_Edge_degeneracy_tester
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                       Dual_graph;

  typedef typename Dual_graph::Edge                Edge;
  typedef typename Dual_graph::Face_handle         Face_handle;
  typedef typename Dual_graph::Edge_circulator     Edge_circulator;
  typedef typename Dual_graph::All_edges_iterator  All_edges_iterator;

  typedef typename Dual_graph::Finite_edges_iterator
  Finite_edges_iterator;

 private:
  typedef SVD_Edge_degeneracy_tester<Dual_graph>   Self;

  typedef typename Dual_graph::Geom_traits         Geom_traits;

  typedef typename Dual_graph::Vertex_handle       Vertex_handle;

  typedef typename Dual_graph::Site_2              Site_2;

  typedef typename Geom_traits::Equal_2            Equal_2;

 private:
  bool is_degenerate_infinite_edge(const Face_handle& f, int i) const
  {
    CGAL_precondition( dual_->is_infinite(f, i) );

    Vertex_handle v = f->vertex( dual_->ccw(i) );
    Vertex_handle v_inf = f->vertex( dual_->cw(i) );

    if ( dual_->is_infinite(v) ) {
      std::swap(v, v_inf);
    }

    if ( v->storage_site().is_segment() ) { return false; }

    Vertex_handle vv[2];

    vv[0] = f->vertex(i);
    vv[1] = dual_->data_structure().mirror_vertex(f, i);

    if ( vv[0] == vv[1] ) { return false; }

    Site_2 s_end[2];
    for (int i = 0; i < 2; i++) {
      if ( vv[i]->storage_site().is_point() ) { return false; }

      Equal_2 are_equal = dual_->geom_traits().equal_2_object();
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

    typename Geom_traits::Orientation_2
      orientation = dual_->geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }

 public:
  SVD_Edge_degeneracy_tester(const Dual_graph* dual = NULL) : dual_(dual) {}

  bool operator()(const Face_handle& f, int i) const {
    if ( dual_->is_infinite(f, i) ) {
      return is_degenerate_infinite_edge(f, i);
    }

    Vertex_handle v3 = f->vertex(i);
    Vertex_handle v4 = dual_->data_structure().mirror_vertex(f, i);

    if ( dual_->is_infinite(v3) || dual_->is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( dual_->ccw(i) );
    Vertex_handle v2 = f->vertex( dual_->cw(i) );

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
class SVD_Face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef DG                                    Dual_graph;
  typedef typename Dual_graph::Vertex_handle    Vertex_handle;

  typedef typename Dual_graph::Finite_vertices_iterator
  Finite_vertices_iterator;

  typedef Edge_tester                           Edge_degeneracy_tester;

 private:
  typedef SVD_Face_degeneracy_tester<Dual_graph,Edge_degeneracy_tester> Self;

  typedef typename Dual_graph::Geom_traits      Geom_traits;
  typedef typename Dual_graph::Edge             Edge;
  typedef typename Dual_graph::Face_handle      Face_handle;
  typedef typename Dual_graph::size_type        size_type;

  typedef typename Dual_graph::Edge_circulator  Edge_circulator;

  typedef typename Dual_graph::All_vertices_iterator
  All_vertices_iterator;

  typedef typename Dual_graph::Site_2           Site_2;

 public:
  SVD_Face_degeneracy_tester(const Dual_graph* dual = NULL)
    : dual_(dual), edge_tester(NULL) {}

  SVD_Face_degeneracy_tester(const Dual_graph* dual,
			     const Edge_degeneracy_tester* e_tester)
    : dual_(dual), edge_tester(e_tester) {}

  bool operator()(const Vertex_handle& v) const
  {
    if ( dual_->is_infinite(v) ) { return false; }

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
    do {
      if ( (*edge_tester)(ec) ) { ++n_degen; }
      else if ( dual_->is_infinite(ec) ) { ++n_inf; }
      else { 
	if ( !dual_->is_infinite(ec) ) {
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
      CGAL_assertion( !dual_->is_infinite(e[i]) );
      CGAL_assertion( !(*edge_tester)(e[i]) );

      Vertex_handle v1 = e[i].first->vertex( dual_->ccw(e[i].second) );
      Vertex_handle v2 = e[i].first->vertex( dual_->cw(e[i].second) );
      vv[i] = (v1 == v) ? v2 : v1;

      CGAL_assertion( v == v1 || v == v2 );

      if ( vv[i]->storage_site().is_point() ) { return false; }

      typename Geom_traits::Equal_2
	are_equal = dual_->geom_traits().equal_2_object();
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

    typename Geom_traits::Orientation_2
      orientation = dual_->geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }
 
  bool operator()(const Finite_vertices_iterator& vit) const {
    return operator()(Vertex_handle(vit));
  }

 private:
  const Dual_graph* dual_;
  const Edge_degeneracy_tester* edge_tester;
};


//=========================================================================
//=========================================================================


template<class SVD2>
struct Segment_Voronoi_diagram_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <SVD2,
   SVD_Edge_degeneracy_tester<SVD2>,
   SVD_Face_degeneracy_tester<SVD2,SVD_Edge_degeneracy_tester<SVD2> >
  >
{
  typedef SVD2                                          Dual_graph;
  typedef SVD_Edge_degeneracy_tester<SVD2>              Edge_tester;
  typedef SVD_Face_degeneracy_tester<SVD2,Edge_tester>  Face_tester;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Default_Voronoi_traits_2
  <SVD2,Edge_tester,Face_tester>
  Base;

  typedef Segment_Voronoi_diagram_Voronoi_traits_2<SVD2>  Self;

  typedef typename Dual_graph::Data_structure Dual_graph_data_structure;

  Segment_Voronoi_diagram_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};


//=========================================================================
//=========================================================================

#if 0
template<class SVD2>
struct Segment_Voronoi_diagram_cached_Voronoi_traits_2
  : public Default_cached_Voronoi_traits_2
  <SVD2,
   SVD_Edge_degeneracy_tester<SVD2>,
   SVD_Face_degeneracy_tester
   <SVD2,Cached_edge_degeneracy_tester<SVD_Edge_degeneracy_tester<SVD2>,
				       Data_structure_project<SVD2> >
   >,
   Data_structure_project<SVD2>
  >
{
  typedef SVD2    Dual_graph;
  typedef Data_structure_project<SVD2>                DG_Project;
  typedef SVD_Edge_degeneracy_tester<SVD2>            Edge_tester;

  typedef Cached_edge_degeneracy_tester<Edge_tester,DG_Project>
  Cached_edge_tester;
  typedef SVD_Face_degeneracy_tester<SVD2,Cached_edge_tester>
  Face_tester;

  typedef
  Default_cached_Voronoi_traits_2<SVD2,Edge_tester,Face_tester,DG_Project>
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;

  typedef typename Dual_graph::Data_structure Dual_graph_data_structure;

  Segment_Voronoi_diagram_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};
#endif

//=========================================================================
//=========================================================================


template<class SVD2>
struct Segment_Voronoi_diagram_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <SVD2,
   SVD_Edge_degeneracy_tester<SVD2>,
   SVD_Face_degeneracy_tester
   <SVD2,CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
    <SVD_Edge_degeneracy_tester<SVD2>,
     CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<SVD2> >
   >,
   CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<SVD2>
  >
{
  typedef SVD2                                         Dual_graph;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ds_project<SVD2>  DG_Project;
  typedef SVD_Edge_degeneracy_tester<SVD2>             Edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ref_counted_edge_degeneracy_tester
  <Edge_tester,DG_Project>
  Ref_counted_edge_tester;

  typedef SVD_Face_degeneracy_tester<SVD2,Ref_counted_edge_tester>
  Face_tester;

  typedef
  CGAL_VORONOI_DIAGRAM_2_NS::Default_ref_counted_Voronoi_traits_2
  <SVD2,Edge_tester,Face_tester,DG_Project>
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;

  typedef typename Dual_graph::Data_structure Dual_graph_data_structure;

  Segment_Voronoi_diagram_cached_Voronoi_traits_2(const Dual_graph* dg = NULL)
    : Base(dg) {}

  const Dual_graph_data_structure& dual_graph_data_structure() const {
    return this->dg_->data_structure();
  }
};


//=========================================================================
//=========================================================================


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
