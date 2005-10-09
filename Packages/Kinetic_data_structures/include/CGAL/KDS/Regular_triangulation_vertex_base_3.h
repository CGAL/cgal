#ifndef CGAL_KDS_KINETIC_REGULAR_VERTEX_BASE_3_H
#define CGAL_KDS_KINETIC_REGULAR_VERTEX_BASE_3_H

#include <CGAL/KDS/basic.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

CGAL_KDS_BEGIN_NAMESPACE
//! A class to track labels of edges of faces in a triangulation
template <class SimulationTraits, 
	  class Vertex_base= CGAL::Triangulation_vertex_base_3<typename SimulationTraits::Instantaneous_kernel> >
class Regular_triangulation_vertex_base_3: 
  public CGAL::Triangulation_vertex_base_with_info_3<typename SimulationTraits::Simulator::Event_key, 
				       typename SimulationTraits::Instantaneous_kernel,
				       Vertex_base>
{
private:
  typedef CGAL::Triangulation_vertex_base_with_info_3<typename SimulationTraits::Simulator::Event_key, 
						      typename SimulationTraits::Instantaneous_kernel,
						      Vertex_base> P;
  typedef typename Vertex_base::Triangulation_data_structure   TDS;
public:
  typedef TDS	                         Triangulation_data_structure;
  typedef typename TDS::Cell_handle     Cell_handle;
  typedef typename TDS::Vertex_handle   Vertex_handle;
  typedef typename Vertex_base::Geom_traits Traits;
  
  typedef typename SimulationTraits::Simulator::Event_key Label;
 
  Regular_triangulation_vertex_base_3(): P() {
  }
  
  Regular_triangulation_vertex_base_3(Cell_handle f): P(f){
  }
  
  template < typename TDS3 >
  struct Rebind_TDS {
    typedef typename Vertex_base::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Regular_triangulation_vertex_base_3<SimulationTraits, Cb3>  Other;
  };
};



CGAL_KDS_END_NAMESPACE


#endif
