#define CGAL_SURFACE_MESHER_VERBOSE 0
#undef CGAL_SURFACE_MESHER_VERBOSE

#include <QApplication>

#include <CGAL/AABB_tree/AABB_polyhedral_oracle.h>
#include <CGAL/AABB_tree/AABB_tree.h>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"
#include <fstream>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>

#include <QTime>
#include <QInputDialog>



void MainWindow::on_actionRemeshing_triggered()
{
  if(onePolygonIsSelected())
  {
    int index = getSelectedPolygonIndex();

    // get active polyhedron
    Polyhedron* pMesh = scene->polyhedron(index);

    if(!pMesh) return;

    // remesh

    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;

    Tr tr; // 3D-Delaunay triangulation
    C2t3 c2t3(tr); // 2D-complex in 3D-Delaunay triangulation

    // TODO: get parameters using ONE dialog box
    // sizing and approximation parameters should be expressed as ratio of 
    // scene bbox diagonal.

    double diag = scene->len_diagonal();

    bool ok;
    const double angle = 
      QInputDialog::getDouble(this, "Min triangle angle",
      "Angle:",
      25, // default value
      1, // min
      30, // max
      2, // decimals
      &ok);
    if(!ok) return;

    const double sizing = 
      QInputDialog::getDouble(this, "Sizing",
      "Size:",
      diag * 0.05, // default value
      diag * 10e-6, // min
      diag, // max
      4, // decimals
      &ok);
    if(!ok) return;

    const double approx = 
      QInputDialog::getDouble(this, "Approximation error",
      "Error:",
      diag * 0.005, // default value
      diag * 10e-7, // min
      diag, // max
      6, // decimals
      &ok);
    if(!ok) return;

    // meshing parameters
    CGAL::Surface_mesh_default_criteria_3<Tr> facets_criteria(angle,sizing,approx);

    QApplication::setOverrideCursor(Qt::WaitCursor);

    // AABB tree
    QTime time;
    time.start();
    std::cout << "Build AABB tree...";
    typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel; 
    typedef CGAL::AABB_tree<Simple_cartesian_kernel,Polyhedron::Facet_handle,Polyhedron> Tree;
    Tree tree;
    tree.build_faces(*pMesh);
    std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

    // input surface
    typedef CGAL::AABB_polyhedral_oracle<Polyhedron,GT,Simple_cartesian_kernel> Input_surface;
    Input_surface input(&tree);

    // initial point set
    time.start();
    std::cout << "Insert initial point set...";
    unsigned int nb_initial_points = 10;
    Polyhedron::Point_iterator it;
    typedef CGAL::Cartesian_converter<Kernel,GT> Converter;
    Converter convert;
    unsigned int i = 0;
    for(it = pMesh->points_begin();
        it != pMesh->points_end(), i < nb_initial_points;
	it++, i++)
      tr.insert(convert(*it));
    std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

    // remesh
    time.start();
    std::cout << "Remesh...";
    CGAL::make_surface_mesh(c2t3, input, facets_criteria, CGAL::Manifold_with_boundary_tag());
    std::cout << "done (" << time.elapsed() << " ms, " << tr.number_of_vertices() << " vertices)" << std::endl;

    if(tr.number_of_vertices() > 0)
    {
      // add remesh as new polyhedron
      Polyhedron *pRemesh = new Polyhedron;
      CGAL::Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron> builder(c2t3);
      pRemesh->delegate(builder);
      scene->addPolyhedron(pRemesh,
			   QObject::tr("%1 remeshed (%2 %3 %4)")
			   .arg(scene->polyhedronName(index))
			   .arg(angle)
			   .arg(sizing)
			   .arg(approx),
			   Qt::magenta,
			   true,
			   scene->polyhedronRenderingMode(index));
    }

    QApplication::restoreOverrideCursor();
  }
}
