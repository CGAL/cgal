#ifndef CGAL_VORONOI_TRAITS_CONCEPT_H
#define CGAL_VORONOI_TRAITS_CONCEPT_H 1

#include <CGAL/basic.h>


CGAL_BEGIN_NAMESPACE


template<class DG>
class Voronoi_traits_concept
{
 public:
  typedef DG                               Dual_graph;
  typedef typename DG::Edge                Edge;
  typedef typename DG::Edge_circulator     Edge_circulator;
  typedef typename DG::All_edges_iterator  All_edges_iterator;
  typedef typename DG::Face_handle         Face_handle;
  typedef typename DG::Vertex_handle       Vertex_handle;

  typedef typename DG::Data_structure      Dual_graph_data_structure;

  Voronoi_traits_concept(const Dual_graph* = NULL) {}

  struct Edge_degeneracy_tester
  {
    typedef DG  Dual_graph;

    typedef typename DG::Edge                   Edge;
    typedef typename DG::Edge_circulator        Edge_circulator;
    typedef typename DG::All_edges_iterator     All_edges_iterator;
    typedef typename DG::Finite_edges_iterator  Finite_edges_iterator;

    typedef bool           result_type;
    typedef Arity_tag<1>   Arity;

    Edge_degeneracy_tester(const DG* = NULL) {}

    bool operator()(const Edge&) const {
      return false;
    }

    bool operator()(const Face_handle&, int) const {
      return false;
    }

    bool operator()(const Edge_circulator&) const {
      return false;
    }

    bool operator()(const All_edges_iterator&) const {
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

    Face_degeneracy_tester(const DG* = NULL) {}
    Face_degeneracy_tester(const DG*, const Edge_degeneracy_tester*) {}

    bool operator()(const Vertex_handle&) const {
      return false;
    }

    bool operator()(const Vertex_circulator&) const {
      return false;
    }
  };

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

  static const Dual_graph_data_structure& dual_graph_data_structure()
  {
    static Dual_graph_data_structure dual_graph_data_structure_static;
    return dual_graph_data_structure_static;
  }

 private:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_TRAITS_CONCEPT_H
