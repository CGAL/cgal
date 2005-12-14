#ifndef SKIN_SURFACE_H
#define SKIN_SURFACE_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulated_mixed_complex_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Marching_tetrahedra_observer_skin_surface_3.h>
#include <CGAL/marching_tetrahedra_3.h>

CGAL_BEGIN_NAMESPACE

template <class Regular_3, class SkinSurfaceTraits_3>
void skin_surface_construct_bounding_box_3(
  Regular_3 &reg,
  SkinSurfaceTraits_3 const& traits) {
  typedef typename Regular_3::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular_3::Geom_traits     GT;
  typedef typename GT::Bare_point             Point;
  typedef typename GT::Point_3                Weighted_point;
  typedef typename GT::RT                     RT;
  
  Finite_vertices_iterator vit = reg.finite_vertices_begin();
  if (vit != reg.finite_vertices_end()) {
    Bbox_3 bbox = vit->point().bbox();
    RT max_weight=vit->point().weight();
    while (++vit != reg.finite_vertices_end()) {
      bbox = bbox + vit->point().bbox();
      if (max_weight < vit->point().weight())
	max_weight = vit->point().weight();
    }

    // add a bounding octahedron:
    RT dx = bbox.xmax() - bbox.xmin();
    RT dy = bbox.ymax() - bbox.ymin();
    RT dz = bbox.zmax() - bbox.zmin();
  
    Point mid(bbox.xmin() + dx/2, bbox.ymin() + dy/2, bbox.zmin() + dz/2);
    RT dr = sqrt(CGAL::to_double(max_weight)) + .001;
  
    reg.insert(Weighted_point(
      Point(bbox.xmax()+(dy+dz+dr)/traits.shrink_factor(),mid.y(),mid.z()),-1));
    reg.insert(Weighted_point(
      Point(bbox.xmin()-(dy+dz+dr)/traits.shrink_factor(),mid.y(),mid.z()),-1));
    reg.insert(Weighted_point(
      Point(mid.x(),bbox.ymax()+(dx+dz+dr)/traits.shrink_factor(),mid.z()),-1));
    reg.insert(Weighted_point(
      Point(mid.x(),bbox.ymin()-(dx+dz+dr)/traits.shrink_factor(),mid.z()),-1));
    reg.insert(Weighted_point(
      Point(mid.x(),mid.y(),bbox.zmax()+(dx+dy+dr)/traits.shrink_factor()),-1));
    reg.insert(Weighted_point(
      Point(mid.x(),mid.y(),bbox.zmin()-(dx+dy+dr)/traits.shrink_factor()),-1));
  }
}

template <class InputIterator, class Polyhedron_3, class SkinSurfaceTraits_3>
void skin_surface_3(InputIterator first, InputIterator last,
  Polyhedron_3 &polyhedron, const SkinSurfaceTraits_3 &skin_surface_traits,
  bool verbose = false) {
  if (first == last) {
    std::cout << " No input balls" << std::endl;
    return;
  }

  // Types
  typedef SkinSurfaceTraits_3                              Skin_surface_traits;
  typedef typename Skin_surface_traits::Regular_traits     Regular_traits;
  typedef typename Regular_traits::Bare_point              Reg_point;
  typedef typename Regular_traits::Weighted_point          Reg_weighted_point;

  typedef Regular_triangulation_3<Regular_traits> Regular;
  typedef Triangulated_mixed_complex_3<SkinSurfaceTraits_3>
                                                  Triangulated_mixed_complex;
  typedef Marching_tetrahedra_traits_skin_surface_3<
    Triangulated_mixed_complex,
    Polyhedron_3,
    typename SkinSurfaceTraits_3::T2P_converter>  Marching_tetrahedra_traits;
  typedef Marching_tetrahedra_observer_skin_surface_3<
    Triangulated_mixed_complex, Polyhedron_3>     Marching_tetrahedra_observer;
    
  // Code
  Regular regular;
  Triangulated_mixed_complex triangulated_mixed_complex;

  while (first != last) {
    regular.insert((*first));
    first++;
  }

  skin_surface_construct_bounding_box_3(regular,skin_surface_traits);
  
  if (verbose) {
    std::cerr << "Triangulation ready" << std::endl;
  }

  // Construct the triangulated mixed complex:
  triangulate_mixed_complex_3(
    regular, triangulated_mixed_complex, skin_surface_traits);

  CGAL_assertion(triangulated_mixed_complex.is_valid());
  if (verbose) {
    std::cerr << "Triangulated mixed complex ready" << std::endl;
  }

  // Extract the coarse mesh using marching_tetrahedra
  Marching_tetrahedra_traits marching_traits;
  Marching_tetrahedra_observer marching_observer;
  
  marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits, marching_observer);

  if (verbose) {
    std::cerr << "Mesh ready" << std::endl;
  }
  
}


CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_H
