#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Make_quad_soup.h>

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
		std::list<Point> min_directions;
		std::list<Point> max_directions;

		Polyhedron::Vertex_iterator v;
		for(v = pMesh->vertices_begin();
		    v != pMesh->vertices_end();
		    v++)
		{
			std::vector<Point> points;

			// pick central point
			const Point& central_point = v->point();
			points.push_back(central_point);

			// compute min edge len around central vertex
			// to scale the ribbons used to display the directions
			FT min_edge_len = 1e38;

			// and its neighbors
			Polyhedron::Halfedge_around_vertex_circulator he = v->vertex_begin();
			Polyhedron::Halfedge_around_vertex_circulator end = he;
			CGAL_For_all(he,end)
			{
				const Point& p = he->opposite()->vertex()->point();
				points.push_back(p);
				FT edge_len = std::sqrt(CGAL::squared_distance(central_point,p));
				min_edge_len = edge_len < min_edge_len ? edge_len : min_edge_len;
			}

			if(points.size() > 5)
			{
				// estimate curvature by fitting
				Fitting monge_fit;
				const int dim_monge = 2;
				const int dim_fitting = 2;
				Monge_form monge_form = monge_fit(points.begin(),points.end(),dim_fitting,dim_monge);

				// make monge form comply with vertex normal (to get correct orientation)
				Vector n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
				monge_form.comply_wrt_given_normal(n);

				Vector normal = min_edge_len * monge_form.normal_direction();
				Vector umin = min_edge_len * monge_form.minimal_principal_direction();
				Vector umax = min_edge_len * monge_form.maximal_principal_direction();

				Point lifted_point = central_point + 0.1 * normal;

				min_directions.push_back(lifted_point +  0.2 * umin + 0.02 * umax);
				min_directions.push_back(lifted_point -  0.2 * umin + 0.02 * umax);
				min_directions.push_back(lifted_point -  0.2 * umin - 0.02 * umax);
				min_directions.push_back(lifted_point +  0.2 * umin - 0.02 * umax);

				max_directions.push_back(lifted_point + 0.02 * umin +  0.2 * umax);
				max_directions.push_back(lifted_point - 0.02 * umin +  0.2 * umax);
				max_directions.push_back(lifted_point - 0.02 * umin -  0.2 * umax);
				max_directions.push_back(lifted_point + 0.02 * umin -  0.2 * umax);
			}
		}

		// add principal curvature directions as new polyhedron
		Polyhedron *pMin_curvature_directions = new Polyhedron;
		Polyhedron *pMax_curvature_directions = new Polyhedron;

		typedef std::list<Point>::iterator Iterator;
		Make_quad_soup<Polyhedron,Kernel,Iterator> min_soup;
		min_soup.run(min_directions.begin(),min_directions.end(),*pMin_curvature_directions);
		scene->addPolyhedron(pMin_curvature_directions,
			tr("%1 (min curvatures)").arg(scene->polyhedronName(index)),
			Qt::red,
			scene->isPolyhedronActivated(index),
			scene->polyhedronRenderingMode(index));

		Make_quad_soup<Polyhedron,Kernel,Iterator> max_soup;
		max_soup.run(max_directions.begin(),max_directions.end(),*pMax_curvature_directions);
		scene->addPolyhedron(pMax_curvature_directions,
			tr("%1 (max curvatures)").arg(scene->polyhedronName(index)),
			Qt::blue,
			scene->isPolyhedronActivated(index),
			scene->polyhedronRenderingMode(index));

		// default cursor
		QApplication::setOverrideCursor(Qt::ArrowCursor);
	}
}
