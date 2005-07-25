#ifndef CGAL_VORONOI_DIAGRAM_2_IDENTITY_VORONOI_TRAITS_2_H
#define CGAL_VORONOI_DIAGRAM_2_IDENTITY_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=============================================================================
//=============================================================================

template<class DG>
struct Identity_edge_degeneracy_tester
{
  typedef DG                                             Delaunay_graph;

  typedef typename Delaunay_graph::Edge                  Edge;
  typedef typename Delaunay_graph::Face_handle           Face_handle;
  typedef typename Delaunay_graph::Edge_circulator       Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator    All_edges_iterator;
  typedef typename Delaunay_graph::Finite_edges_iterator Finite_edges_iterator;

  typedef bool           result_type;
  typedef Arity_tag<2>   Arity;

  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const All_edges_iterator& eit) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const Edge_circulator& ec) const {
    return false;
  }
};

//=============================================================================

template<class DG>
struct Identity_face_degeneracy_tester
  : public Default_face_degeneracy_tester<DG>
{};

//=============================================================================
//=============================================================================

template<class DG, class VT>
class Identity_Voronoi_traits_2_base
  : public VT
{
 public:
  typedef Identity_edge_degeneracy_tester<DG>  Edge_degeneracy_tester;
  typedef Identity_face_degeneracy_tester<DG>  Face_degeneracy_tester;

  const Edge_degeneracy_tester&  edge_degeneracy_tester_object() const {
    return edge_degeneracy_tester_;
  }

  const Face_degeneracy_tester&  face_degeneracy_tester_object() const {
    return face_degeneracy_tester_;
  }

 private:
  Edge_degeneracy_tester  edge_degeneracy_tester_;
  Face_degeneracy_tester  face_degeneracy_tester_;
};

//=============================================================================
//=============================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

//=============================================================================

template<class DG> class Identity_Voronoi_traits_2;

//=============================================================================

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_IDENTITY_VORONOI_TRAITS_2_H
