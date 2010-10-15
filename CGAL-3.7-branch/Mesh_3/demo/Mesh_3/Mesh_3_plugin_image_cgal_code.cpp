#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include "Image_type.h"
#include "Meshing_thread.h"
#include "Mesh_function.h"

#include <CGAL/Bbox_3.h>


Meshing_thread* cgal_code_mesh_3(const Image* pImage,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape)
{
  typedef Mesh_function<Image_mesh_domain> Mesh_function;
  
  if( NULL == pImage ) { return NULL; }
  
  Image_mesh_domain* p_domain = new Image_mesh_domain(*pImage, 1e-6);
  
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item();
  
  Mesh_parameters param;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(), p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}
