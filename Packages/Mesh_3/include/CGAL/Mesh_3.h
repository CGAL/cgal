#ifndef CGAL_MESH_3_H
#define CGAL_MESH_3_H

#include <CGAL/basic.h>

namespace CGAL {

// Tr: a conform Delaunay triangulation
template <class Tr>
class Mesh_3 : public Tr {
public:
  typedef Tr Triangulation;
  typedef Mesh_3<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  void fill_facette_map();
};

#endif
