#include "config.h"

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS

#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include "implicit_functions/Implicit_function_interface.h"
#include "Meshing_thread.h"
#include "Mesh_function.h"

#include <CGAL/Bbox_3.h>


Meshing_thread* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape)
{
  typedef Mesh_function<Function_mesh_domain> Mesh_function;
  
  if( pfunction == NULL ) { return NULL; }
  
  CGAL::Bbox_3 domain_bbox (pfunction->bbox().xmin,
                            pfunction->bbox().ymin,
                            pfunction->bbox().zmin,
                            pfunction->bbox().xmax,
                            pfunction->bbox().ymax,
                            pfunction->bbox().zmax);
  
  Function_mesh_domain* p_domain =
    new Function_mesh_domain(Function_wrapper(*pfunction), domain_bbox, 1e-7);
  
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item();
  
  Mesh_parameters param;
  param.protect_features = false;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(), p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}

#endif // CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
