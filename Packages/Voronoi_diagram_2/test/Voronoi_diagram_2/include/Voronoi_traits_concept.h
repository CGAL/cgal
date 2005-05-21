#ifndef CGAL_VORONOI_TRAITS_CONCEPT_H
#define CGAL_VORONOI_TRAITS_CONCEPT_H 1

#include <CGAL/basic.h>


CGAL_BEGIN_NAMESPACE


template<class DG>
class Voronoi_traits_concept
{
 public:
  typedef DG                               Dual_graph;

  struct Edge_degeneracy_tester
  {
    typedef DG  Dual_graph;

    typedef typename DG::Edge                   Edge;
    typedef typename DG::Face_handle            Face_handle;
    typedef typename DG::Edge_circulator        Edge_circulator;
    typedef typename DG::All_edges_iterator     All_edges_iterator;
    typedef typename DG::Finite_edges_iterator  Finite_edges_iterator;

    typedef bool           result_type;
    typedef Arity_tag<1>   Arity;

    bool operator()(const Dual_graph&, const Edge&) const {
      return false;
    }

    bool operator()(const Dual_graph&, const Face_handle&, int) const {
      return false;
    }

    bool operator()(const Dual_graph&, const Edge_circulator&) const {
      return false;
    }

    bool operator()(const Dual_graph&, const All_edges_iterator&) const {
      return false;
    } 
  };

  struct Face_degeneracy_tester
  {
    typedef DG  Dual_graph;

    typedef typename DG::Vertex_handle             Vertex_handle;
    typedef typename DG::Vertex_circulator         Vertex_circulator;
    typedef typename DG::All_vertices_iterator     All_vertices_iterator;
    typedef typename DG::Finite_vertices_iterator  Finite_vertices_iterator;

    typedef bool           result_type;
    typedef Arity_tag<1>   Arity;

    bool operator()(const Dual_graph&, const Edge_degeneracy_tester&,
		    const Vertex_handle&) const {
      return false;
    }

    bool operator()(const Dual_graph&, const Edge_degeneracy_tester&,
		    const Vertex_circulator&) const {
      return false;
    }
  };

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 private:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_TRAITS_CONCEPT_H
