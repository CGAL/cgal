#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Monge_via_jet_fitting.h>

void MainWindow::on_actionEstimateCurvature_triggered()
{
	if(onePolygonIsSelected())
	{
		// wait cursor
		QApplication::setOverrideCursor(Qt::WaitCursor);

		// get active polyhedron
		int index = getSelectedPolygonIndex();
		Polyhedron* pMesh = scene->polyhedron(index);

		// types
		typedef CGAL::Monge_via_jet_fitting<Kernel> Fitting;
		typedef Fitting::Monge_form Monge_form;

		// store curvature directions as point quadruplets (little ribbons)
		std::vector<Point> min_directions;
		std::vector<Point> max_directions;

		Polyhedron::Vertex_iterator v;
		for(v = pMesh->vertices_begin();
		    v != pMesh->vertices_end();
		    v++)
		{
			std::vector<Point> points;

			// pick central point
			points.push_back(v->point());

			// and its neighbors
			Polyhedron::Halfedge_around_vertex_circulator he = v->vertex_begin();
			Polyhedron::Halfedge_around_vertex_circulator end = he;
			CGAL_For_all(he,end)
			{
				const Point& p = he->opposite()->vertex()->point();
				points.push_back(p);
			}

			// estimate curvature by fitting
			Fitting monge_fit;
			const int dim_monge = 2;
			const int dim_fitting = 2;
			Monge_form monge_form = monge_fit(points.begin(),points.end(),dim_fitting,dim_monge);

			Vector umin = monge_form.minimal_principal_direction();
			Vector umax = monge_form.maximal_principal_direction();
		}

		// add principal curvature directions as new polyhedron
		Polyhedron *pMin_curvature_directions = new Polyhedron;
		Polyhedron *pMax_curvature_directions = new Polyhedron;

		scene->addPolyhedron(pMin_curvature_directions,
			tr("%1 (min curvatures)").arg(scene->polyhedronName(index)),
			scene->polyhedronColor(index),
			scene->isPolyhedronActivated(index),
			scene->polyhedronRenderingMode(index));
		scene->addPolyhedron(pMax_curvature_directions,
			tr("%1 (max curvatures)").arg(scene->polyhedronName(index)),
			scene->polyhedronColor(index),
			scene->isPolyhedronActivated(index),
			scene->polyhedronRenderingMode(index));

		// default cursor
		QApplication::setOverrideCursor(Qt::ArrowCursor);
	}
}
