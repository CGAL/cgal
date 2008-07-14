#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Make_soup.h>

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/self_intersect.h>



void MainWindow::on_actionSelf_intersection_triggered()
{
	if(onePolygonIsSelected())
	{
		QApplication::setOverrideCursor(Qt::WaitCursor);

		// get active polyhedron
		int index = getSelectedPolygonIndex();
		Polyhedron* pMesh = scene->polyhedron(index);

		// types
		typedef std::list<Triangle>::iterator Iterator;
		typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

		// compute self-intersections
		std::list<Triangle> intersecting_triangles;
		typedef std::back_insert_iterator<std::list<Triangle> > OutputIterator;
		self_intersect<Polyhedron,Kernel,OutputIterator>(*pMesh,std::back_inserter(intersecting_triangles));

		// add intersecting triangles as a new polyhedron (a triangle soup)
		Polyhedron *pSoup = new Polyhedron;
		Make_soup<Polyhedron,Kernel,Iterator> soup;
		soup.run(intersecting_triangles.begin(),intersecting_triangles.end(),*pSoup);

		scene->addPolyhedron(pSoup,
			tr("%1 (intersecting triangles)").arg(scene->polyhedronName(index)),
			scene->polyhedronColor(index), // PA: to be changed to red 
			scene->isPolyhedronActivated(index),
			scene->polyhedronRenderingMode(index));

		QApplication::setOverrideCursor(Qt::ArrowCursor);
	}
}

