#ifndef MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

CGAL_BEGIN_NAMESPACE 

template <class Triangulation_3, class HalfedgeDS, class Converter >
class Marching_tetrahedra_traits_skin_surface_3 {
public:
  typedef typename Triangulation_3::Vertex_handle            Vertex_handle;
  typedef typename Triangulation_3::Edge                     Edge;
  typedef typename Triangulation_3::Cell_handle              Cell_handle;
  typedef typename Triangulation_3::Geom_traits::Point_3     Triang_Point;

  typedef typename HalfedgeDS::Traits                 HDS_K;
  typedef typename HDS_K::RT                   HDS_rt;
  typedef typename HDS_K::Point_3              HDS_point;

  Marching_tetrahedra_traits_skin_surface_3(HDS_rt iso_value=0)
    : iso_value(iso_value) {
  }

  // These two functions are required by the marching tetrahedra algorithm
  Sign sign(Cell_handle ch, int i) const {
    return CGAL_NTS sign(
       ch->surf->value(converter(ch->vertex(i)->point())) - iso_value);
  }
  HDS_point intersection(Cell_handle ch, int i, int j) const {
    // Precondition: ch is not an infinite cell: their surface is not set
    return ch->surf->to_surface(
      converter(ch->vertex(i)->point()),
      converter(ch->vertex(j)->point()));
  }

private:
  Converter converter;
  HDS_rt iso_value;
};

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
