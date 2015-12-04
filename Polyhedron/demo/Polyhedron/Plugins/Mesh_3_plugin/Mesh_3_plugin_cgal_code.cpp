#include "config_mesh_3.h"

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Bbox_3.h>

#include <Polyhedron_type.h>
#include <C3t3_type.h>

#include <Scene_c3t3_item.h>
#include <CGAL/Three/Scene_item.h>

#include "Mesh_function.h"

#include <CGAL/Timer.h>
using namespace CGAL::Three;
// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Edge_criteria Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             QString filename,
                             const double angle,
                             const double facet_sizing,
                             const double approx,
                             const double tet_sizing,
                             const double tet_shape,
                             const bool protect_features,
                             CGAL::Three::Scene_interface* scene)
{
  if(!pMesh) return 0;

  // remesh

  // Set mesh criteria
  Edge_criteria edge_criteria(facet_sizing);
  Facet_criteria facet_criteria(angle, facet_sizing, approx); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape, tet_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  CGAL::Timer timer;
  timer.start();
  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << angle << std::endl
            << "  facets size bound: " << facet_sizing << std::endl
            << "  approximation bound: " << approx << std::endl
            << "  tetrahedra size bound: " << tet_sizing << std::endl;
  std::cerr << "Build AABB tree...";
  // Create domain
  Polyhedral_mesh_domain domain(*pMesh);
  if(protect_features) {
      domain.detect_features();
  }
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // Meshing
  std::cerr << "Mesh...";
  CGAL::parameters::internal::Features_options features =
          protect_features ?
              CGAL::parameters::features(domain) :
              CGAL::parameters::no_features();

  Scene_c3t3_item* new_item = new Scene_c3t3_item(
    CGAL::make_mesh_3<C3t3>(domain, criteria,
      features,
      CGAL::parameters::no_perturb(),
      CGAL::parameters::no_exude()));

  new_item->set_scene(scene);
  std::cerr << "done (" << timer.time() << " ms, " << new_item->c3t3().triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(new_item->c3t3().triangulation().number_of_vertices() > 0)
  {

    std::ofstream medit_out("out.mesh");
    new_item->c3t3().output_to_medit(medit_out);

    const CGAL::Three::Scene_item::Bbox& bbox = new_item->bbox();
    new_item->setPosition((float)(bbox.xmin + bbox.xmax)/2.f,
                          (float)(bbox.ymin + bbox.ymax)/2.f,
                          (float)(bbox.zmin + bbox.zmax)/2.f);
    return new_item;
  }
  else {
    delete new_item;
    return 0;
  }
}

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS

Scene_item* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tet_sizing,
                                 const double tet_shape,
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

  // Set mesh criteria
  Edge_criteria edge_criteria(facet_sizing);
  Facet_criteria facet_criteria(facet_angle, facet_sizing, facet_approx); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape, tet_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  Scene_c3t3_item* p_new_item = new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(*p_domain,
                                                                           criteria,
                                                                           CGAL::parameters::no_perturb(),
                                                                           CGAL::parameters::no_exude()));

  CGAL::Timer timer;
  p_new_item->set_scene(scene);
  std::cerr << "done (" << timer.time() << " ms, "
    << p_new_item->c3t3().triangulation().number_of_vertices() << " vertices)"
    << std::endl;

  if (p_new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    const Scene_item::Bbox& bbox = p_new_item->bbox();
    p_new_item->setPosition((float)(bbox.xmin + bbox.xmax) / 2.f,
                            (float)(bbox.ymin + bbox.ymax) / 2.f,
                            (float)(bbox.zmin + bbox.zmax) / 2.f);
    return p_new_item;
  }
  else {
    delete p_new_item;
    return 0;
  }
}
#endif // CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS


#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES

Scene_item* cgal_code_mesh_3(const Image* pImage,
                             const double facet_angle,
                             const double facet_sizing,
                             const double facet_approx,
                             const double tet_sizing,
                             const double tet_shape,
                             CGAL::Three::Scene_interface* scene)
{

  if (NULL == pImage) { return NULL; }

  Image_mesh_domain* p_domain
    = new Image_mesh_domain(*pImage, 1e-6);

  // Set mesh criteria
  Edge_criteria edge_criteria(facet_sizing);
  Facet_criteria facet_criteria(facet_angle, facet_sizing, facet_approx); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape, tet_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);
  CGAL::Timer timer;
  timer.start();
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(*p_domain,
                                                                            criteria,
                                                                            CGAL::parameters::no_perturb(),
                                                                            CGAL::parameters::no_exude())
);

  p_new_item->set_scene(scene);

  std::cerr << "done (" << timer.time() << " ms, "
    << p_new_item->c3t3().triangulation().number_of_vertices() << " vertices)"
    << std::endl;

  if (p_new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    const Scene_item::Bbox& bbox = p_new_item->bbox();
    p_new_item->setPosition((float)(bbox.xmin + bbox.xmax) / 2.f,
                            (float)(bbox.ymin + bbox.ymax) / 2.f,
                            (float)(bbox.zmin + bbox.zmax) / 2.f);
    return p_new_item;
  }
  else {
    delete p_new_item;
    return 0;
  }

}

#endif //CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES


//#include "Polyhedron_demo_mesh_3_plugin_cgal_code.moc"
//#include "Scene_c3t3_item.moc" //Check this one, it's strange moc include.

