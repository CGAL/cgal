#include "Polyhedron_type_fwd.h"
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#  include "Image_type_fwd.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#  include "implicit_functions/Implicit_function_interface.h"
#endif

#include "Kernel_type.h"
#include "Meshing_thread.h"

struct Mesh_parameters;

namespace CGAL { namespace Three {
  class Scene_item;
  class Scene_interface;
}}

typedef std::list<std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> > Polylines_container;

Meshing_thread* cgal_code_mesh_3(const Polyhedron* pMesh,
                                 const Polylines_container&,
                                 QString filename,
                                 const double angle,
                                 const double facet_sizing,
                                 const double approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 bool protect_features,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene);

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
Meshing_thread* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene);
#endif

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
Meshing_thread* cgal_code_mesh_3(const CGAL::Image_3* pImage,
                                 const Polylines_container& polylines,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 bool protect_features,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene);
#endif
