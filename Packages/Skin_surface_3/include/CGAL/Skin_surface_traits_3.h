#ifndef CGAL_SKIN_SURFACE_TRAITS_3_H
#define CGAL_SKIN_SURFACE_TRAITS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
//#include <CGAL/Regular_triangulation_3.h>

//#include <CGAL/Skin_surface_simplicial_complex_3.h>

//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Skin_surface_polyhedral_items_3.h>

CGAL_BEGIN_NAMESPACE 

template <class Kernel>
class Construct_anchor_point_3 {
public:
  typedef typename Kernel::RT            RT;
  typedef typename Kernel::Point_3       Point;
  
  Construct_anchor_point_3(RT &shrink) : shrink(shrink) {}
  
  Point operator()(const Point &anchor_del, const Point &anchor_vor) {
    return anchor_del + shrink*(anchor_vor - anchor_del);
  }
  RT shrink;
};

// The following four converters convert geometric objects between the
//   datastructures mentioned above:
// - RegToSimplConverter
// - RegToMeshConverter 
// - SimplToMeshConverter 
// - MeshToSimplConverter 

template <
  class RegularKernel = Exact_predicates_inexact_constructions_kernel,
  class TriangulatedMixedComplexKernel =
    Exact_predicates_exact_constructions_kernel,
  class PolyhedronKernel = Simple_cartesian<double> >  
class Skin_surface_traits_3 {
public:
  typedef Skin_surface_traits_3<RegularKernel, TriangulatedMixedComplexKernel,
				PolyhedronKernel>                         Self;
  // Kernel definitions and converters
  typedef RegularKernel                      Regular_kernel;
  typedef TriangulatedMixedComplexKernel     Triangulated_mixed_complex_kernel;
  typedef PolyhedronKernel                   Polyhedron_kernel;

  typedef Regular_triangulation_euclidean_traits_3<Regular_kernel>
                                             Regular_traits;

  typedef Regular_triangulation_euclidean_traits_3<
           Triangulated_mixed_complex_kernel> Triangulated_mixed_complex_traits;

  typedef typename Triangulated_mixed_complex_traits::Construct_weighted_circumcenter_3
                                              Construct_weighted_circumcenter_3;

  typedef Weighted_converter_3 <
    Cartesian_converter<Regular_kernel,Triangulated_mixed_complex_kernel> >
                                              R2T_converter;
  typedef Weighted_converter_3 <
    Cartesian_converter<Regular_kernel, Polyhedron_kernel> >
                                              R2P_converter;
  typedef Weighted_converter_3 <
    Cartesian_converter<Triangulated_mixed_complex_kernel,Polyhedron_kernel> >
                                              T2P_converter; 
  typedef Weighted_converter_3 <
    Cartesian_converter<Polyhedron_kernel,Triangulated_mixed_complex_kernel> >
                                              P2T_converter; 

  Skin_surface_traits_3(typename Regular_kernel::RT shrink = .5) :
    shrink(shrink) {}
  
  // Function objects:
  R2T_converter r2t_converter_object() const {
    return R2T_converter();
  }
  R2P_converter r2p_converter_object() const {
    return R2P_converter();
  }
  T2P_converter t2p_converter_object() const {
    return T2P_converter();
  }
  P2T_converter p2t_converter_object() const {
    return P2T_converter();
  }

//   template <class Kernel>
//   Construct_weighted_circumcenter_3<
//     Regular_triangulation_euclidean_traits_3<Kernel> >
//   construct_weighted_circumcenter_3_object() const {
//     return Regular_triangulation_euclidean_traits_3<Kernel>().
//       construct_weighted_circumcenter_3_object();
//   }

  template <class Kernel>
  In_smallest_orthogonal_sphere_3<
    Regular_triangulation_euclidean_traits_3<Kernel> >
  in_smallest_orthogonal_sphere_3_object() const {
    return Regular_triangulation_euclidean_traits_3<Kernel>().
      in_smallest_orthogonal_sphere_3_object();
  }

  Construct_anchor_point_3<Triangulated_mixed_complex_kernel>
  construct_weighted_circumcenter_3_object() const {
    return Construct_anchor_point_3<Triangulated_mixed_complex_kernel>(shrink);
  }
  typename Regular_kernel::RT shrink_factor() const {
    return shrink;
  }

  typename Regular_kernel::RT shrink;
};

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_TRAITS_3_H
