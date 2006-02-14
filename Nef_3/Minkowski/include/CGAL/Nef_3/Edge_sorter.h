#ifndef CGAL_NEF3_EDGE_SORTER_H
#define CGAL_NEF3_EDGE_SORTER_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_io_parser.h>

CGAL_BEGIN_NAMESPACE



template<typename Nef_>
class Edge_sorter : public Modifier_base<typename Nef_::SNC_structure> {
  
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;

 public:
  Edge_sorter() {}
      
  void operator()(SNC_structure& snc) {
    std::sort(snc.halfedges_begin(), snc.halfedges_end(), sort_edges<SNC_structure>(snc));
  }
};
  
CGAL_END_NAMESPACE
#endif //CGAL_NEF3_EDGE_SORTER_H
