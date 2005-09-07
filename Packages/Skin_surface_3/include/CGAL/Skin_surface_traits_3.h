#ifndef CGAL_SKIN_SURFACE_TRAITS_3_H
#define CGAL_SKIN_SURFACE_TRAITS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/Skin_surface_simplicial_complex_3.h>

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
  class Regular_K_ = Exact_predicates_inexact_constructions_kernel,
  class Simplicial_K_ = Exact_predicates_exact_constructions_kernel,
  class Mesh_K_ = Simple_cartesian<double> >  
class Skin_surface_traits_3 {
public:
  typedef Skin_surface_traits_3<Regular_K_, Simplicial_K_, Mesh_K_> Self;
  // Kernel definitions and converters
  typedef Regular_K_                                       Regular_K;
  typedef Simplicial_K_                                    Simplicial_K;
  typedef Mesh_K_                                          Mesh_K;
	
  typedef Weighted_converter_3 <Cartesian_converter<Regular_K, Simplicial_K> >
                                                           R2S_converter;
  typedef Weighted_converter_3 <Cartesian_converter<Simplicial_K,Mesh_K,
    To_double<typename Simplicial_K::FT> > >               S2M_converter; 
  typedef Weighted_converter_3 <Cartesian_converter<Mesh_K,Simplicial_K > >
                                                           M2S_converter; 

  // Regular triangulation
  typedef CGAL::Regular_triangulation_euclidean_traits_3<Regular_K>
                                                           Regular_traits;
  typedef CGAL::Regular_triangulation_3<Regular_traits>    Regular;

  // Simplicial complex
  typedef Skin_surface_simplicial_vertex_base_3<Self>      Simplicial_vertex;
  typedef Skin_surface_simplicial_cell_base_3<Self>        Simplicial_cell;
  typedef Triangulation_data_structure_3 <Simplicial_vertex,Simplicial_cell>
                                                           Simplicial_TDS;
  typedef Skin_surface_simplicial_complex_3<Simplicial_K, Simplicial_TDS>
                                                           Simplicial;

  // Polyhedral mesh

  
  // Tags
  typedef Tag_false                                        Cache_anchor;

  // Function objects:
  R2S_converter r2s_converter_object() const
  { return R2S_converter(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_TRAITS_3_H
