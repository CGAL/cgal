#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>

namespace CGAL {

template < class Gt,
	   class ST = Segment_Delaunay_graph_storage_traits_2<Gt>,
	   class STag = Tag_false,
	   class D_S = Triangulation_data_structure_2<
              Segment_Delaunay_graph_hierarchy_vertex_base_2<
		Segment_Delaunay_graph_vertex_base_2<ST> >,
              Segment_Delaunay_graph_face_base_2<Gt> >,
	   class LTag = Tag_false>
class Segment_Delaunay_graph_Linf_hierarchy_2
  : public Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag,D_S,LTag,
             Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag> >
{
};


} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H
