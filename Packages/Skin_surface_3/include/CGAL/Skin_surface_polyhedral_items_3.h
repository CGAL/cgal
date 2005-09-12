#ifndef SKIN_SURFACE_POLYHEDRAL_ITEMS_3_H
#define SKIN_SURFACE_POLYHEDRAL_ITEMS_3_H

#include <CGAL/HalfedgeDS_face_base.h>
#include <set>

template <class Refs,  class  T, class PolyT>
struct Skin_surface_polyhedral_vertex
  : public CGAL::HalfedgeDS_vertex_base<Refs, T,
					typename PolyT::Point_3> {
	
  typedef typename PolyT::Point_3               Point;
  typedef typename PolyT::RT                    RT;
  typedef CGAL::HalfedgeDS_vertex_base<Refs, T, Point> Base;
	
  Skin_surface_polyhedral_vertex() : Base() {}
  Skin_surface_polyhedral_vertex(const Skin_surface_polyhedral_vertex &v)
    : Base(v) {}
  Skin_surface_polyhedral_vertex(const Point &p) : Base(p) {}
};

// Faces have references to the tetrahedra that contain them:
template <class Refs, class TetraComplex>
struct Skin_surface_polyhedral_face : public CGAL::HalfedgeDS_face_base<Refs> {
  Skin_surface_polyhedral_face() : CGAL::HalfedgeDS_face_base<Refs>(), tetras() {}
  Skin_surface_polyhedral_face(const Skin_surface_polyhedral_face &f)
    : CGAL::HalfedgeDS_face_base<Refs>(f), tetras(f.tetras) {}

  std::set<typename TetraComplex::Cell_handle> tetras;
};

// An items type using my face.
template <class TetraComplex>
struct Skin_surface_polyhedral_items_3 : public CGAL::Polyhedron_items_3 {
  template <class Refs, class Traits>
  struct Vertex_wrapper {
    typedef typename Traits::Point_3 Point;
    typedef Skin_surface_polyhedral_vertex<Refs, CGAL::Tag_true, Traits> Vertex;
  };
  template <class Refs, class Traits>
  struct Face_wrapper {
    typedef Skin_surface_polyhedral_face<Refs, TetraComplex> Face;
  };
};

#endif // SKIN_SURFACE_POLYHEDRAL_ITEMS_3_H
