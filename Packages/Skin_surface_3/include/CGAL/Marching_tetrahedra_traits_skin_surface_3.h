#ifndef MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

// #include <CGAL/Polyhedron_incremental_builder_3.h>
// #include <CGAL/Modifier_base.h>
// #include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE 

template <class T_, class HDS_, class Converter_ >
class Marching_tetrahedra_traits_skin_surface_3 {
public:
  typedef T_         T;
  typedef HDS_       HDS;
  typedef Converter_ Converter;

  typedef typename T::Vertex_handle            Sc_vertex_handle;
  typedef typename T::Edge                     Sc_edge;
  typedef typename T::Facet                    Sc_facet;
  typedef typename T::Cell_handle              Sc_cell_handle;
  typedef typename T::Geom_traits::Point_3     Sc_point;

  typedef typename HDS::Traits                 Mesh_K;
  typedef typename Mesh_K::RT                  Mesh_rt;
  typedef typename Mesh_K::Point_3             Mesh_point;

  Marching_tetrahedra_traits_skin_surface_3(Mesh_rt iso_value=0)
    : iso_value(iso_value) {
  }
  
  Sign sign(Sc_vertex_handle const vh) {
    return CGAL::sign(
      vh->cell()->surf->value(converter(vh->point())) - iso_value);
  }
  Mesh_rt value(Sc_cell_handle const &ch, Mesh_point const &p) {
    return ch->surf->value(p);
  }
  Mesh_rt value(Sc_cell_handle const &ch, Sc_point const &p) {
    return ch->surf->value(converter(p));
  }
  Mesh_point intersection(Sc_edge const& e) {
    // Precondition: e.first is not an infinite cell: they have not surface set
    return e.first->surf->to_surface(
      converter(e.first->vertex(e.second)->point()),
      converter(e.first->vertex(e.third)->point()));
  }
  
  Converter converter;
  Mesh_rt iso_value;
};

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
