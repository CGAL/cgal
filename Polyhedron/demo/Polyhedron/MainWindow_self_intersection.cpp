#include "MainWindow.h"
#include "Scene.h"

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/self_intersect.h>
#include <CGAL/Make_triangle_soup.h>

void MainWindow::on_actionSelf_intersection_triggered()
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
		std::list<Triangle> intersecting_triangles;
		typedef std::back_insert_iterator<std::list<Triangle> > OutputIterator;
		std::cout << "Self-intersect...";
		self_intersect<Polyhedron,Kernel,OutputIterator>(*pMesh,std::back_inserter(intersecting_triangles));
		std::cout << "ok (" << intersecting_triangles.size() << " triangle(s))" << std::endl;

		// add intersecting triangles as a new polyhedron (a triangle soup)
		if(intersecting_triangles.size() != 0)
		{
			Polyhedron *pSoup = new Polyhedron;
			Make_triangle_soup<Polyhedron,Kernel,Iterator> soup_builder;
			soup_builder.run(intersecting_triangles.begin(),intersecting_triangles.end(),*pSoup);
			scene->addPolyhedron(pSoup,
				tr("%1 (intersecting triangles)").arg(scene->polyhedronName(index)),
				Qt::red,
				scene->isPolyhedronActivated(index),
				scene->polyhedronRenderingMode(index));
		}

		QApplication::setOverrideCursor(Qt::ArrowCursor);
	}
}

