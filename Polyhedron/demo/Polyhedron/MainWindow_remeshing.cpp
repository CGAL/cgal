#include "MainWindow.h"
#include "Scene.h"
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/Collisions/AABB_polyhedral_oracle.h>
#include <CGAL/Collisions/AABB_tree.h>

#define CGAL_C2T3_USE_POLYHEDRON
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <QTime>
#include <QInputDialog>

void MainWindow::on_actionRemeshing_triggered()
{
  if(onePolygonIsSelected())
  {
    int index = getSelectedPolygonIndex();

    // get active polyhedron
    Polyhedron* pMesh = scene->polyhedron(index);

    // remesh
    double init = clock();

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
    typedef CGAL::Cartesian<double> Cartesian_kernel; 
    typedef CGAL::AABB_tree<GT,Polyhedron::Facet_handle,Polyhedron> Tree;
    Tree tree;
    tree.build_faces(*pMesh);
    std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

    // input surface
    typedef CGAL::AABB_polyhedral_oracle<Polyhedron,GT,GT> Input_surface;
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
    CGAL::make_surface_mesh(c2t3, input, facets_criteria, CGAL::Manifold_tag());
    std::cout << "done (" << time.elapsed() << " ms, " << tr.number_of_vertices() << " vertices)" << std::endl;

    if(tr.number_of_vertices() > 0)
    {
      // add remesh as new polyhedron
      Polyhedron *pRemesh = new Polyhedron;

      // TODO: dump output to polyhedron

      // for now...
      std::ofstream out("d:\\remesh.off");
      CGAL::output_surface_facets_to_off(out, c2t3);

      // NO IDEA WHY THIS LINE DOES NOT COMPILE
 //     scene->addPolyhedron(pRemesh,
	//tr("%1 (remesh)").arg(scene->polyhedronName(index)),
	//Qt::magenta,
	//scene->isPolyhedronActivated(index),
	//scene->polyhedronRenderingMode(index));
    }

    QApplication::restoreOverrideCursor();
  }
}

/*


int main() {

// defining the surface
std::ifstream file_input("data/triceratops.off");
Polyhedral_surface surface(file_input);

// defining meshing criteria
CGAL::Surface_mesh_default_criteria_3<Tr> 
facets_criteria(30.,  // angular bound
0.5,  // radius bound
0.5); // distance bound
CGAL::Surface_mesh_default_edges_criteria_3<Tr>
edges_criteria(0.5,  // radius bound
0.5); // distance bound

// meshing surface
CGAL::make_piecewise_smooth_surface_mesh(c2t3, surface,
facets_criteria, edges_criteria,
CGAL::Manifold_tag());

std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}

*/