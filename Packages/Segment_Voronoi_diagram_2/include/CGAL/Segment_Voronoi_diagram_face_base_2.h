#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H

#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>

CGAL_BEGIN_NAMESPACE


template <class Gt,
	  class Fb = Triangulation_ds_face_base_2<> >
class Segment_Voronoi_diagram_face_base_2
  : public Apollonius_graph_face_base_2<Gt,Fb>
{};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
