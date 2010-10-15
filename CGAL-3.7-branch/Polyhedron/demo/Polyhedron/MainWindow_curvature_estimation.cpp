#ifdef CGAL_LAPACK_ENABLED

#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Make_quad_soup.h>
#include <CGAL/compute_normal.h>

void MainWindow::on_actionEstimateCurvature_triggered()
{
  if(onePolygonIsSelected())
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // get active polyhedron
    int index = getSelectedSceneItemIndex();
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
      Polyhedron::Halfedge_around_vertex_circulator he = v->vertex_begin();
      Polyhedron::Halfedge_around_vertex_circulator end = he;
      CGAL_For_all(he,end)
      {
	const Point& p = he->opposite()->vertex()->point();
	points.push_back(p);
	FT edge_len = std::sqrt(CGAL::squared_distance(central_point,p));
	min_edge_len = edge_len < min_edge_len ? edge_len : min_edge_len; // avoids #undef min
      }

      const double du = 0.2;
      const double dv = 0.02;
      const double lift = 0.1;
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

	Point lifted_point = central_point + lift * normal;

	min_directions.push_back(lifted_point +  du * umin + dv * umax);
	min_directions.push_back(lifted_point -  du * umin + dv * umax);
	min_directions.push_back(lifted_point -  du * umin - dv * umax);
	min_directions.push_back(lifted_point +  du * umin - dv * umax);

	max_directions.push_back(lifted_point + dv * umin +  du * umax);
	max_directions.push_back(lifted_point - dv * umin +  du * umax);
	max_directions.push_back(lifted_point - dv * umin -  du * umax);
	max_directions.push_back(lifted_point + dv * umin -  du * umax);
      }
    }

    // add principal curvature directions as new polyhedron
    Polyhedron *pMin = new Polyhedron;
    Polyhedron *pMax = new Polyhedron;

    typedef std::list<Point>::iterator Iterator;
    Make_quad_soup<Polyhedron,Kernel,Iterator> min_soup;
    min_soup.run(min_directions.begin(),min_directions.end(),*pMin);
    scene->addPolyhedron(pMin,
      tr("%1 (min curvatures)").arg(scene->polyhedronName(index)),
      Qt::red,
      scene->isPolyhedronVisible(index),
      scene->polyhedronRenderingMode(index));

    Make_quad_soup<Polyhedron,Kernel,Iterator> max_soup;
    max_soup.run(max_directions.begin(),max_directions.end(),*pMax);
    scene->addPolyhedron(pMax,
      tr("%1 (max curvatures)").arg(scene->polyhedronName(index)),
      Qt::blue,
      scene->isPolyhedronVisible(index),
      scene->polyhedronRenderingMode(index));

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#else //  CGAL_LAPACK_ENABLED

#include <QMessageBox>
void MainWindow::on_actionEstimateCurvature_triggered()
{
  QMessageBox::warning(this, "Function not available", 
		       "This function is not available."
		       " You need to configure LAPACK support.");
}

#endif // CGAL_LAPACK_ENABLED
