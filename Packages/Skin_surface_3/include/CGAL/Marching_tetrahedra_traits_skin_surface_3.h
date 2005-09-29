#ifndef MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

// #include <CGAL/Polyhedron_3_incremental_builder_3.h>
// #include <CGAL/Modifier_base.h>
// #include <CGAL/Cartesian_converter.h>

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
  
  Sign sign(const Vertex_handle vh) const {
    return CGAL::sign(
      vh->cell()->surf->value(converter(vh->point())) - iso_value);
  }
  HDS_point intersection(const Edge &e) const {
    // Precondition: e.first is not an infinite cell: they have not surface set
    return e.first->surf->to_surface(
      converter(e.first->vertex(e.second)->point()),
      converter(e.first->vertex(e.third)->point()));
  }

  // Additional functions, not belonging to the traits concept:
  HDS_rt value(const Cell_handle &ch, const HDS_point &p) const {
    return ch->surf->value(p);
  }
  HDS_rt value(const Cell_handle &ch, const Triang_Point &p) const {
    return ch->surf->value(converter(p));
  }

  
private:
  Converter converter;
  HDS_rt iso_value;
};

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
