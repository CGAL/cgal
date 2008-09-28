#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/self_intersect.h>
#include <CGAL/Make_triangle_soup.h>

void MainWindow::on_actionSelfIntersection_triggered()
{
  if(onePolygonIsSelected())
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // get active polyhedron
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);

    // compute self-intersections
    typedef std::list<Triangle>::iterator Iterator;
    typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;
    std::list<Triangle> triangles; // intersecting triangles
    typedef std::back_insert_iterator<std::list<Triangle> > OutputIterator;
    std::cout << "Self-intersect...";
    ::self_intersect<Polyhedron,Kernel,OutputIterator>(*pMesh,std::back_inserter(triangles));
    std::cout << "ok (" << triangles.size() << " triangle(s))" << std::endl;

    // add intersecting triangles as a new polyhedron, i.e., a triangle soup.
    if(triangles.size() != 0)
    {
      Polyhedron *pSoup = new Polyhedron;
      Make_triangle_soup<Polyhedron,Kernel,Iterator> soup_builder;
      soup_builder.run(triangles.begin(),triangles.end(),*pSoup);
      scene->addPolyhedron(pSoup,
	tr("%1 (intersecting triangles)").arg(scene->polyhedronName(index)),
	Qt::magenta,
	scene->isPolyhedronActivated(index),
	scene->polyhedronRenderingMode(index));
    }

    QApplication::restoreOverrideCursor();
  }
}

