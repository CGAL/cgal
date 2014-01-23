#include <CGAL/AABB_intersections.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"
#include "C3t3_type.h"
#include "Scene_c3t3_item.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <fstream>

#include <CGAL/Timer.h>

// Mesh Criteria
typedef C3t3::Triangulation Tr;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

#include "Scene_item.h"
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             const QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing)
{
  if(!pMesh) return 0;

  // remesh

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
            << "  tetrahedra size bound: " << tets_sizing << std::endl;
  std::cerr << "Build AABB tree...";
  // Create domain
  Mesh_domain domain(*pMesh);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // Meshing
  std::cerr << "Mesh...";
  Scene_c3t3_item* new_item = 
    new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(domain, criteria));

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

#include "Scene_c3t3_item.moc"
