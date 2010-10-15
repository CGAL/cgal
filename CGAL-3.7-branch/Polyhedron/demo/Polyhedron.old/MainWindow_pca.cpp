#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/centroid.h>
#include <CGAL/bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>


#include <CGAL/Make_quad_soup.h> // output for plane fitting
#include <CGAL/Make_bar.h> // output for line fitting

// specify this to compile on windows
#undef min
#undef max

void MainWindow::on_actionFitPlane_triggered()
{
  if(onePolygonIsSelected())
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // get active polyhedron
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);

    // get triangles from the mesh
    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator f;
    for(f = pMesh->facets_begin();
      f != pMesh->facets_end();
      ++f)
    {
      const Point& a = f->halfedge()->vertex()->point();
      const Point& b = f->halfedge()->next()->vertex()->point();
      const Point& c = f->halfedge()->prev()->vertex()->point();
      triangles.push_back(Triangle(a,b,c));
    }

    // fit plane to triangles
    Plane plane;
    std::cout << "Fit plane...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());

    // compute bounding box diagonal
    Iso_cuboid bbox = CGAL::bounding_box(pMesh->points_begin(),pMesh->points_end());

    // compute scale (for rendering) using diagonal of bbox
    Point cmin = (bbox.min)();
    Point cmax = (bbox.max)();
    FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));
    Vector u1 = plane.base1();
    u1 = u1 / std::sqrt(u1*u1);
    u1 = u1 * 0.7 * diag;
    Vector u2 = plane.base2();
    u2 = u2 / std::sqrt(u2*u2);
    u2 = u2 * 0.7 * diag;
    std::list<Point> points;
    points.push_back(center_of_mass + u1);
    points.push_back(center_of_mass + u2);
    points.push_back(center_of_mass - u1);
    points.push_back(center_of_mass - u2);

    // add best fit plane as new polyhedron
    Polyhedron *pFit = new Polyhedron;
    typedef std::list<Point>::iterator Iterator;
    Make_quad_soup<Polyhedron,Kernel,Iterator> quad;
    quad.run(points.begin(),points.end(),*pFit);

    scene->addPolyhedron(pFit,
      tr("%1 (plane fit)").arg(scene->polyhedronName(index)),
      Qt::magenta,
      scene->isPolyhedronActivated(index),
      scene->polyhedronRenderingMode(index));

    QApplication::restoreOverrideCursor();
  }
}

void MainWindow::on_actionFitLine_triggered()
{
  if(onePolygonIsSelected())
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // get active polyhedron
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);

    // get triangles from the mesh
    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator f;
    for(f = pMesh->facets_begin();
      f != pMesh->facets_end();
      ++f)
    {
      const Point& a = f->halfedge()->vertex()->point();
      const Point& b = f->halfedge()->next()->vertex()->point();
      const Point& c = f->halfedge()->prev()->vertex()->point();
      triangles.push_back(Triangle(a,b,c));
    }

    // fit line to triangles
    Line line;
    std::cout << "Fit line...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());

    // compute bounding box diagonal
    Iso_cuboid bbox = CGAL::bounding_box(pMesh->points_begin(),pMesh->points_end());

    // compute scale for rendering using diagonal of bbox
    Point cmin = bbox.min();
    Point cmax = bbox.max();
    FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));

    // construct a 3D bar
    Vector u = line.to_vector();
    u = u / std::sqrt(u*u);

    Point a = center_of_mass + u * diag;
    Point b = center_of_mass - u * diag;

    Plane plane_a = line.perpendicular_plane(a);

    Vector u1 = plane_a.base1();
    u1 = u1 / std::sqrt(u1*u1);
    u1 = u1 * 0.01 * diag;
    Vector u2 = plane_a.base2();
    u2 = u2 / std::sqrt(u2*u2);
    u2 = u2 * 0.01 * diag;

    Point points[8];

    points[0] = a + u1;
    points[1] = a + u2;
    points[2] = a - u1;
    points[3] = a - u2;

    points[4] = b + u1;
    points[5] = b + u2;
    points[6] = b - u1;
    points[7] = b - u2;

    // add best fit line as new polyhedron bar
    Polyhedron *pFit = new Polyhedron;
    Make_bar<Polyhedron,Kernel> bar;
    bar.run(points,*pFit);

    scene->addPolyhedron(pFit,
      tr("%1 (line fit)").arg(scene->polyhedronName(index)),
      Qt::magenta,
      scene->isPolyhedronActivated(index),
      scene->polyhedronRenderingMode(index));

    QApplication::restoreOverrideCursor();
  }
}
