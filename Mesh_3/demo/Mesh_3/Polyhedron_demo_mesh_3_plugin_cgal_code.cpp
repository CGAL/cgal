//#define CGAL_MESH_3_VERBOSE
#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <fstream>

#include <CGAL/Timer.h>


// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

Scene_c3t3_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                                  const QString filename,
                                  const double angle,
                                  const double sizing,
                                  const double approx,
                                  const double tets_sizing,
                                  const double tet_shape,
                                  const bool lloyd,
                                  const bool odt,
                                  const bool perturb,
                                  const bool exude)
{
  if(!pMesh) return 0;

  // remesh

  typedef Tr::Geom_traits GT;

  // Set mesh criteria
  Facet_criteria facet_criteria(angle, sizing, approx); // angle, size, approximation
  Cell_criteria cell_criteria(4, tets_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  CGAL::Timer timer;
  timer.start();
  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << angle << std::endl
            << "  facets size bound: " << sizing << std::endl
            << "  approximation bound: " << approx << std::endl
            << "  tetrahedra size bound: " << tets_sizing << std::endl
            << "  tetrahedron radius-edge: " << tet_shape << std::endl;
  std::cerr << "Build AABB tree...";
  // Create domain
  Mesh_domain domain(*pMesh);
  std::cerr << "done (" << timer.time() << " s)" << std::endl;

  // Meshing
  timer.stop();timer.reset();timer.start();
  std::cerr << "Mesh...";
  
  namespace cgp = CGAL::parameters;
  namespace cgpi = cgp::internal;
  
  cgpi::Lloyd_options lloyd_obj = lloyd ? cgp::lloyd() : cgp::no_lloyd();
  cgpi::Odt_options odt_obj = odt ? cgp::odt() : cgp::no_odt();
  cgpi::Perturb_options perturb_obj = perturb ? cgp::perturb() : cgp::no_perturb();
  cgpi::Exude_options exude_obj = exude ? cgp::exude() : cgp::no_exude();
  
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      lloyd_obj, odt_obj, perturb_obj, exude_obj);
  
  std::cerr << "done (" << timer.time() << " s, " << c3t3.triangulation().number_of_vertices() << " vertices)" << std::endl;  
  Scene_c3t3_item* new_item = new Scene_c3t3_item(c3t3);


  if(new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    std::ofstream medit_out("out.mesh");
    new_item->c3t3().output_to_medit(medit_out);
    const Scene_item::Bbox& bbox = new_item->bbox();
    new_item->setPosition((bbox.xmin + bbox.xmax)/2.f,
                          (bbox.ymin + bbox.ymax)/2.f,
                          (bbox.zmin + bbox.zmax)/2.f);
    return new_item;
  }
  else {
    delete new_item;
    return 0;
  }
}
