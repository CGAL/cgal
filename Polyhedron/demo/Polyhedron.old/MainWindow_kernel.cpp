#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"
#include <CGAL/Polyhedron_kernel.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Dualizer.h>
#include <CGAL/translate.h>

void MainWindow::on_actionKernel_triggered()
{
  if(onePolygonIsSelected())
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    int index = getSelectedPolygonIndex();

    // get active polyhedron
    Polyhedron* pMesh = scene->polyhedron(index);

    typedef CGAL::Gmpzf ET; // choose exact integral type
    typedef Polyhedron_kernel<Kernel,ET> Polyhedron_kernel;

    // get triangles from polyhedron
    std::list<Triangle> triangles;
    get_triangles(*pMesh,std::back_inserter(triangles));

    // solve LP 
    std::cout << "Solve linear program...";
    Polyhedron_kernel kernel;
    if(!kernel.solve(triangles.begin(),triangles.end()))
    {
      std::cout << "done (empty kernel)" << std::endl;
      QApplication::restoreOverrideCursor();
      return;
    }
    std::cout << "done" << std::endl;

    // add kernel as new polyhedron
    Polyhedron *pKernel = new Polyhedron;

    // get inside point
    Point inside_point = kernel.inside_point();
    Vector translate = inside_point - CGAL::ORIGIN;

    // compute dual of translated polyhedron w.r.t. inside point.
    std::cout << "Compute dual of translated polyhedron...";
    std::list<Point> dual_points;
    std::list<Triangle>::iterator it;
    for(it = triangles.begin();
      it != triangles.end();
      it++)
    {
      const Triangle& triangle = *it;
      const Point p0 = triangle[0] - translate;
      const Point p1 = triangle[1] - translate;
      const Point p2 = triangle[2] - translate;
      Plane plane(p0,p1,p2); 
      Vector normal = plane.orthogonal_vector();
      normal = normal / std::sqrt(normal*normal);
      // compute distance to origin (do not use plane.d())
      FT distance_to_origin = std::sqrt(CGAL::squared_distance(Point(CGAL::ORIGIN),plane));
      Point dual_point = CGAL::ORIGIN + normal / distance_to_origin;
      dual_points.push_back(dual_point);
    }
    std::cout << "ok" << std::endl;

    // compute convex hull in dual space
    std::cout << "convex hull in dual space...";
    Polyhedron convex_hull;
    CGAL::convex_hull_3(dual_points.begin(),dual_points.end(),convex_hull);
    std::cout << "ok" << std::endl;

    // dualize and translate back to get final kernel
    Dualizer<Polyhedron,Kernel> dualizer;
    dualizer.run(convex_hull,*pKernel);
    ::translate<Polyhedron,Kernel>(*pKernel,translate);
    pKernel->inside_out();

    scene->addPolyhedron(pKernel,
      tr("%1 (kernel)").arg(scene->polyhedronName(index)),
      Qt::magenta,
      scene->isPolyhedronActivated(index),
      scene->polyhedronRenderingMode(index));

    QApplication::restoreOverrideCursor();
  }
}
