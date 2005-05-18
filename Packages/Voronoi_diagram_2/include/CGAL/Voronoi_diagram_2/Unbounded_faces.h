#ifndef CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
#define CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================

template<class VDA, class Base_it>
class Bounded_face_tester
{
 private:
  // this class returns true if the face is bounded

  // this is essentially VDA::Non_degenerate_faces_iterator
  typedef Base_it                                      Base_iterator;
  typedef typename VDA::Dual_graph::Vertex_circulator  Dual_vertex_circulator;
  typedef typename VDA::Dual_graph::Vertex_handle      Dual_vertex_handle;

 public:
  Bounded_face_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    Dual_vertex_handle v = it.base();

    Dual_vertex_circulator vc = vda_->dual().incident_vertices(v);
    Dual_vertex_circulator vc_start = vc;
    do {
      if ( vda_->dual().is_infinite(vc) ) { return false; }
      ++vc;
    } while ( vc != vc_start );
    return true;
  }

 private:
  const VDA* vda_;
};

//=========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
