#ifndef CGAL_KDS_KINETIC_REGULAR_CELL_BASE_3_H
#define CGAL_KDS_KINETIC_REGULAR_CELL_BASE_3_H

#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Delaunay_triangulation_cell_base_3.h>

CGAL_KDS_BEGIN_NAMESPACE
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits, class Cell_base= CGAL::Triangulation_cell_base_3<typename SimulationTraits::Instantaneous_kernel> >
class Regular_triangulation_cell_base_3: public Delaunay_triangulation_cell_base_3<SimulationTraits, Cell_base>
{
private:
  typedef typename Cell_base::Triangulation_data_structure   TDS;
public:
  typedef TDS	                         Triangulation_data_structure;
  typedef typename TDS::Cell_handle     Cell_handle;
  typedef typename TDS::Vertex_handle   Vertex_handle;
  typedef typename Cell_base::Geom_traits Traits;
  
  typedef typename SimulationTraits::Simulator::Event_key Edge_label;
  typedef Edge_label Facet_label;
  Regular_triangulation_cell_base_3(): Cell_base() {
  }
  
  Regular_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1, 
				    Vertex_handle v2, Vertex_handle v3): Cell_base(v0, v1, v2, v3){
  }
  
  Regular_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1, 
				    Vertex_handle v2, Vertex_handle v3,
				    Cell_handle f0, Cell_handle f1,
				    Cell_handle f2, Cell_handle f3): Cell_base(v0,v1,v2, v3, f0,f1,f2, f3){
  }
  
  template < typename TDS3 >
  struct Rebind_TDS {
    typedef typename Cell_base::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Delaunay_triangulation_cell_base_3<SimulationTraits, Cb3>  Other;
  };
};



CGAL_KDS_END_NAMESPACE


#endif
