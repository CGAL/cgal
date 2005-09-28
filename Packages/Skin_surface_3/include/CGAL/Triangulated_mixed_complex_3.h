#ifndef TRIANGULATED_MIXED_COMPLEX_3
#define TRIANGULATED_MIXED_COMPLEX_3

#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulated_mixed_complex_cell_3.h>

CGAL_BEGIN_NAMESPACE

template <
  class SkinSurfaceTraits_3,
  class GT=typename SkinSurfaceTraits_3::Triangulated_mixed_complex_kernel,
  class PolyhedronKernel_3=typename SkinSurfaceTraits_3::Polyhedron_kernel,
  class Tds = Triangulation_data_structure_3 <
    Triangulation_vertex_base_3<GT>,
    Triangulated_mixed_complex_cell_3<GT, PolyhedronKernel_3> > >
class Triangulated_mixed_complex_3 : public Triangulation_3<GT, Tds> {
};

CGAL_END_NAMESPACE

#endif // TRIANGULATED_MIXED_COMPLEX_3
