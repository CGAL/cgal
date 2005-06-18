#ifndef CGAL_VDA_TEST_CONCEPT_H
#define CGAL_VDA_TEST_CONCEPT_H 1


template<class DG>
void test_dual_graph_concept(const DG& dg)
{
  typedef typename DG::size_type                      size_type;
  typedef typename DG::Geom_traits                    Geom_traits;
  typedef typename DG::Triangulation_data_structure   Tds;

  //  typedef typename DG::Point_2                        Point_2;
  //  typedef typename DG::Site_2                         Site_2;

  typedef typename DG::Vertex                         Vertex;
  typedef typename DG::Face                           Face;
  typedef typename DG::Edge                           Edge;

  typedef typename DG::Vertex_handle                  Vertex_handle;
  typedef typename DG::Face_handle                    Face_handle;

  typedef typename DG::All_edges_iterator             All_edges_iterator;
  typedef typename DG::Finite_edges_iterator          Finite_edges_iterator;

  typedef typename DG::All_vertices_iterator          All_vertices_iterator;
  typedef typename DG::Finite_vertices_iterator       Finite_vertices_iterator;

  typedef typename DG::All_faces_iterator             All_faces_iterator;
  typedef typename DG::Finite_faces_iterator          Finite_faces_iterator;

  typedef typename DG::Face_circulator                Face_circulator;
  typedef typename DG::Edge_circulator                Edge_circulator;
  typedef typename DG::Vertex_circulator              Vertex_circulator;
}



template<class VT>
void test_voronoi_traits_concept(const VT& vt)
{
  typedef typename VT::Dual_graph              DG;
  typedef typename VT::Edge_degeneracy_tester  EDT;
  typedef typename VT::Face_degeneracy_tester  FDT;

  typedef typename VT::Vertex_handle           Vertex_handle;
  typedef typename VT::Point_2                 Point_2;
  typedef typename VT::Site_2                  Site_2;

  typedef typename VT::Voronoi_vertex_2        VV2;
  typedef typename VT::Voronoi_edge_2          VE2;
  typedef typename VT::Curve                   Curve;

  test_edt_concept( vt.edge_degeneracy_tester_object() );
  test_fdt_concept( vt.face_degeneracy_tester_object() );
}


template<class EDT>
void test_edt_concept(const EDT& edt)
{
  typedef typename EDT::result_type               result_type;
  typedef typename EDT::Arity                     Arity;

  typedef typename EDT::Dual_graph                DG;
  typedef typename EDT::Edge                      Edge;
  typedef typename EDT::Face_handle               Face_handle;
  typedef typename EDT::Edge_circulator           Edge_circulator;
  typedef typename EDT::All_edges_iterator        All_edges_iterator;
  typedef typename EDT::Finite_edges_iterator     Finite_edges_iterator;
}

template<class FDT>
void test_fdt_concept(const FDT& fdt)
{
  typedef typename FDT::result_type               result_type;
  typedef typename FDT::Arity                     Arity;

  typedef typename FDT::Dual_graph                DG;
  typedef typename FDT::Vertex_handle             Vertex_handle;
}




#endif // CGAL_VDA_TEST_CONCEPT_H
