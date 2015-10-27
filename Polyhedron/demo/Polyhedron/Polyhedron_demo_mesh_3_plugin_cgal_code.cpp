#include <CGAL/AABB_intersections.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <Polyhedron_type.h>
#include <C3t3_type.h>

#include <Scene_polyhedron_item.h>
#include <Scene_polygon_soup_item.h>
#include <Scene_c3t3_item.h>
#include <Scene_item.h>

#include <fstream>
#include <sstream>

#include <CGAL/Timer.h>

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

  Scene_c3t3_item* new_item =
          new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(domain, criteria, features));
  new_item->set_scene(scene);
  std::cerr << "done (" << timer.time() << " ms, " << new_item->c3t3().triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    std::ofstream medit_out("out.mesh");
    new_item->c3t3().output_to_medit(medit_out);

    const Scene_item::Bbox& bbox = new_item->bbox();
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

#include "Polyhedron_demo_mesh_3_plugin_cgal_code.moc"
//#include "Scene_c3t3_item.moc" //Check this one, it's strange moc include.

