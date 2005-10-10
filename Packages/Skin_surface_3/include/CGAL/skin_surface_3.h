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

template <class InputIterator, class Polyhedron_3, class SkinSurfaceTraits_3>
void skin_surface_3(InputIterator first, InputIterator last,
  Polyhedron_3 &polyhedron, const SkinSurfaceTraits_3 &skin_surface_traits) {

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
  
  if (first != last) {
    // Construct regular triangulation ...
    Bbox_3 bbox = (*first).bbox();
    double max_weight=1;
    while (first != last) {
      max_weight = std::max(max_weight, (*first).weight());
      bbox = bbox + (*first).bbox();
      regular.insert((*first));
      first++;
    }

    // add a bounding octahedron:
    Reg_point mid((bbox.xmin() + bbox.xmax())/2,
                  (bbox.ymin() + bbox.ymax())/2,
                  (bbox.zmin() + bbox.zmax())/2);
    double size = 1.5*((bbox.xmax() - bbox.xmin() +
		        bbox.ymax() - bbox.ymin() +
		        bbox.zmax() - bbox.zmin())/2 + max_weight);
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x()+size,mid.y(),mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x()-size,mid.y(),mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y()+size,mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y()-size,mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()+size),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()-size),-1));
  }

  // Construct the triangulated mixed complex:
  triangulate_mixed_complex_3(
    regular, triangulated_mixed_complex, skin_surface_traits);

  // Extract the coarse mesh using marching_tetrahedra
  Marching_tetrahedra_traits marching_traits;
  Marching_tetrahedra_observer marching_observer;
  
  marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits, marching_observer);
}


CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_H
