#ifndef CGAL_SKIN_SURFACE_TRAITS_3_H
#define CGAL_SKIN_SURFACE_TRAITS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

CGAL_BEGIN_NAMESPACE 

// NGHK: CHANGE COMMENT
// Data structures:
// RegularTriangulation3 is the triangulation of the balls for which the
//   skin surface is computed
// SimplicialComplex3 is the triangulation that will contain the tetrahedral
//   decomposition of the space into tetrahedra containing patches of degree 2
// Mesh3 will be used to store the final mesh
//
// The following four converters convert geometric objects between the
//   datastructures mentioned above:
// - RegToSimplConverter
// - RegToMeshConverter 
// - SimplToMeshConverter 
// - MeshToSimplConverter 
template <
  class RegularTriangulationTraits = Exact_predicates_inexact_constructions_kernel,
  class SimplicialComplexTraits = Exact_predicates_exact_constructions_kernel,
  class MeshTraits = Simple_cartesian<double>,
  class RegToSimplConverter = Weighted_converter_3 < Cartesian_converter<
    RegularTriangulationTraits, SimplicialComplexTraits > >,
  class SimplToMeshConverter = Weighted_converter_3 < Cartesian_converter<
    SimplicialComplexTraits, MeshTraits,
    To_double<typename SimplicialComplexTraits::FT> > >,
  class MeshToSimplConverter = Weighted_converter_3 < Cartesian_converter<
    MeshTraits, SimplicialComplexTraits > > >
class Skin_surface_traits_3 {
public:
  typedef RegularTriangulationTraits Regular_K;
  typedef SimplicialComplexTraits    Simplicial_K;
  typedef MeshTraits                 Mesh_K;
	
  typedef RegToSimplConverter   R2S_converter;
  typedef SimplToMeshConverter  S2M_converter; 
  typedef MeshToSimplConverter  M2S_converter; 

  // Do you want to store the centers into the simplicial complex, such that
  // you can update the simplicial complex for another the shrink factor.
  typedef Tag_true              Store_centers_in_simplicial_complex;
};

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_TRAITS_3_H
