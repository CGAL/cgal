#include "config_mesh_3.h"
#include "Mesh_3_plugin_cgal_code.h"

#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Bbox_3.h>

#include <Polyhedron_type.h>
#include <C3t3_type.h>

#include <Scene_c3t3_item.h>

#include "Mesh_function.h"
#include "Facet_extra_criterion.h"

#include <CGAL/Timer.h>
using namespace CGAL::Three;

typedef Tr::Point Point_3;

Meshing_thread* cgal_code_mesh_3(const Polyhedron* pMesh,
                                 const Polylines_container& polylines,
                                 QString filename,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 bool protect_features,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene)
{
  if(!pMesh) return 0;

  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << facet_angle << std::endl
            << "  facets size bound: " << facet_sizing << std::endl
            << "  approximation bound: " << facet_approx << std::endl
            << "  tetrahedra size bound: " << tet_sizing << std::endl;
  std::cerr << "Build AABB tree...";
  CGAL::Timer timer;
  timer.start();
  // Create domain
  Polyhedral_mesh_domain* p_domain = new Polyhedral_mesh_domain(*pMesh);
  if(polylines.empty() && protect_features) {
      p_domain->detect_features();
  }
  if(! polylines.empty()){
    p_domain->add_features(polylines.begin(), polylines.end());
    protect_features = true; // so that it will be passed in make_mesh_3
  }

  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->set_scene(scene);

  Mesh_parameters param;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.manifold = manifold;
  param.protect_features = protect_features;

  typedef ::Mesh_function<Polyhedral_mesh_domain> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS

Meshing_thread* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene)
{
  if (pfunction == NULL) { return NULL; }

  CGAL::Bbox_3 domain_bbox(pfunction->bbox().xmin,
                           pfunction->bbox().ymin,
                           pfunction->bbox().zmin,
                           pfunction->bbox().xmax,
                           pfunction->bbox().ymax,
                           pfunction->bbox().zmax);

  Function_mesh_domain* p_domain =
    new Function_mesh_domain(Function_wrapper(*pfunction), domain_bbox, 1e-7);

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->set_scene(scene);

  Mesh_parameters param;
  param.protect_features = false;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.manifold = manifold;

  typedef ::Mesh_function<Function_mesh_domain> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}
#endif // CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS


#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES

Meshing_thread* cgal_code_mesh_3(const Image* pImage,
                                 const Polylines_container& polylines,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
                                 bool protect_features,
                                 const int manifold,
                                 CGAL::Three::Scene_interface* scene)
{
  if (NULL == pImage) { return NULL; }

  Image_mesh_domain* p_domain = new Image_mesh_domain(*pImage, 1e-6);


  if(protect_features && polylines.empty()){
    std::vector<std::vector<Point_3> > polylines_on_bbox;
    CGAL::polylines_to_protect<Point_3>(*pImage, polylines_on_bbox);
    p_domain->add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());
  }
  if(! polylines.empty()){
    // Insert edge in domain
    p_domain->add_features(polylines.begin(), polylines.end());
    protect_features = true; // so that it will be passed in make_mesh_3
  }
  CGAL::Timer timer;
  timer.start();
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item;
  p_new_item->set_scene(scene);

  Mesh_parameters param;
  param.protect_features = false;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tet_sizing;
  param.tet_shape = tet_shape;
  param.manifold = manifold;
  
  typedef ::Mesh_function<Image_mesh_domain> Mesh_function;
  Mesh_function* p_mesh_function = new Mesh_function(p_new_item->c3t3(),
                                                     p_domain, param);
  return new Meshing_thread(p_mesh_function, p_new_item);
}

#endif //CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES


//#include "Polyhedron_demo_mesh_3_plugin_cgal_code.moc"
//#include "Scene_c3t3_item.moc" //Check this one, it's strange moc include.

