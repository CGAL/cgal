#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_HV_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_HV_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/basic_hv.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Predicates_hv_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Constructions_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Traits_base_2.h>

namespace CGAL {

//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------

template<class R, class MTag, class ITag>
class Segment_Delaunay_graph_Linf_hv_traits_base_2
 : public Segment_Delaunay_graph_Linf_traits_base_2<R, MTag, ITag>
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  
  typedef Segment_Delaunay_graph_Linf_traits_base_2<R, MTag, ITag> Base;

  typedef typename Base::Kernel Kernel;
  typedef Kernel   K;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HV_2_NS::Vertex_conflict_C2<K,MTag>
  Vertex_conflict_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HV_2_NS::Finite_edge_interior_conflict_C2<K,MTag>
  Finite_edge_interior_conflict_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HV_2_NS::Infinite_edge_interior_conflict_C2<K,MTag>
  Infinite_edge_interior_conflict_2;

public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  Vertex_conflict_2
  vertex_conflict_2_object() const {
    return Vertex_conflict_2();
  }

  Finite_edge_interior_conflict_2
  finite_edge_interior_conflict_2_object() const {
    return Finite_edge_interior_conflict_2();
  }

  Infinite_edge_interior_conflict_2
  infinite_edge_interior_conflict_2_object() const {
    return Infinite_edge_interior_conflict_2();
  }

};

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_TRAITS_BASE_HV_2_H
