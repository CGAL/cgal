#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QGLWidget>
#include <QMessageBox>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QTime>
#include <QApplication>
#include <QInputDialog>

#include "render_edges.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>

Scene::Scene()
{
	m_pPolyhedron = NULL;

	// view options
	m_view_points = true;
	m_view_segments = true;
	m_view_polyhedron = true;
}

Scene::~Scene()
{
	delete m_pPolyhedron;
}

int
Scene::open(QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"...").arg(filename);

  QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

  QFileInfo fileinfo(filename);
  std::ifstream in(filename.toUtf8());

  if(!in || !fileinfo.isFile() || ! fileinfo.isReadable())
  {
	std::cerr << "cannot open file" << std::endl;
    QApplication::restoreOverrideCursor();
    return -1;
  }

  if(m_pPolyhedron != NULL)
	  delete m_pPolyhedron;

  // allocate new polyhedron
  m_pPolyhedron = new Polyhedron;
  in >> *m_pPolyhedron;
  if(!in)
  {
	std::cerr << "file is not a valid OFF file" << std::endl;
    QApplication::restoreOverrideCursor();

    delete m_pPolyhedron;
	m_pPolyhedron = NULL;

    return -1;
  }

  QApplication::restoreOverrideCursor();

  cerr << " Ok.\n";
  return 0;
}

void Scene::draw()
{
	// draw black edges
	if(m_pPolyhedron != NULL && m_view_polyhedron)
	{
		::glDisable(GL_LIGHTING);
		::glColor3ub(0,0,0);
		::glLineWidth(1.0f);
		gl_render_edges(*m_pPolyhedron);
	}

	// draw red points
	if(m_points.size() != 0 && m_view_points)
	{
		::glDisable(GL_LIGHTING);
		::glColor3ub(180,0,0);
		::glPointSize(2.0f);
		::glBegin(GL_POINTS);
		std::list<Point>::iterator it;
		for(it = m_points.begin(); it != m_points.end(); it++)
		{
			const Point& p = *it;
			::glVertex3d(p.x(),p.y(),p.z());
		}
		::glEnd();
	}

	// draw green segments
	if(m_segments.size() != 0 && m_view_segments)
	{
		::glDisable(GL_LIGHTING);
		::glColor3ub(0,100,0);
		::glLineWidth(2.0f);
		::glBegin(GL_LINES);
		std::list<Segment>::iterator it;
		for(it = m_segments.begin(); it != m_segments.end(); it++)
		{
			const Segment& s = *it;
			const Point& p = s.source();
			const Point& q = s.target();
			::glVertex3d(p.x(),p.y(),p.z());
			::glVertex3d(q.x(),q.y(),q.z());
		}
		::glEnd();
	}
}

Point Scene::random_point()
{
		FT x = (double)rand() / RAND_MAX - 0.5;
		FT y = (double)rand() / RAND_MAX - 0.5;
		FT z = (double)rand() / RAND_MAX - 0.5;
		return Point(x,y,z);
}

Vector Scene::random_vector()
{
		FT x = (double)rand() / RAND_MAX;
		FT y = (double)rand() / RAND_MAX;
		FT z = (double)rand() / RAND_MAX;
		return Vector(x,y,z);
}

void Scene::generate_inside_points(const unsigned int nb_trials)
{
	typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;

	std::cout << "Construct AABB tree...";
	Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done." << std::endl;

    QTime time;
    time.start();
	std::cout << "Generate inside points from " << nb_trials << " trials: ";

	unsigned int i;
	Vector vec = random_vector();
	for(i=0;i<nb_trials;i++)
	{
		Point p = random_point();
		Ray ray(p,vec);
		int nb = (int)tree.number_of_intersected_primitives(ray);
		if(nb % 2 != 0)
			m_points.push_back(p);
	}
	std::cout << m_points.size() << " points, " << time.elapsed() << " ms." << std::endl;
}

void Scene::generate_boundary_segments(const unsigned int nb_slices)
{
	typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;
	typedef Tree::Object_and_primitive_id Object_and_primitive_id;

	std::cout << "Construct AABB tree...";
	Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done." << std::endl;

    QTime time;
    time.start();
	std::cout << "Generate boundary segments from " << nb_slices << " slices: ";

	Vector normal((FT)0.0,(FT)0.0,(FT)1.0);
	unsigned int i;
	for(i=0;i<nb_slices;i++)
	{
		FT z = -0.5 + (FT)i / (FT)nb_slices;
		Point p((FT)0.0, (FT)0.0, z);
		Plane plane(p,normal);

		std::list<Object_and_primitive_id> intersections;
		tree.all_intersections(plane,std::back_inserter(intersections));

		std::list<Object_and_primitive_id>::iterator it;
		for(it = intersections.begin();
			it != intersections.end();
			it++)
		{
			Object_and_primitive_id op = *it;
			CGAL::Object object = op.first;
			Segment segment;
			if(CGAL::assign(segment,object))
				m_segments.push_back(segment);
		}
	}
	std::cout << m_segments.size() << " segments, " << time.elapsed() << " ms." << std::endl;
}

void Scene::generate_boundary_points(const unsigned int nb_points)
{
	typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;
	typedef Tree::Object_and_primitive_id Object_and_primitive_id;

	std::cout << "Construct AABB tree...";
	Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done." << std::endl;

    QTime time;
    time.start();
	std::cout << "Generate boundary points: ";

	unsigned int nb = 0;
	unsigned int nb_lines = 0;
	while(nb < nb_points)
	{
		Point p = random_point();
		Point q = random_point();
		Line line(p,q);

		std::list<Object_and_primitive_id> intersections;
		tree.all_intersections(line,std::back_inserter(intersections));
		nb_lines++;

		std::list<Object_and_primitive_id>::iterator it;
		for(it = intersections.begin();
			it != intersections.end();
			it++)
		{
			Object_and_primitive_id op = *it;
			CGAL::Object object = op.first;
			Point point;
			if(CGAL::assign(point,object))
			{
				m_points.push_back(point);
				nb++;
			}
		}
	}
	std::cout << nb_lines << " lines launched, " << time.elapsed() << " ms." << std::endl;

}

void Scene::generate_edge_points(const unsigned int nb_points)
{
	typedef CGAL::AABB_polyhedron_segment_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;
	typedef Tree::Object_and_primitive_id Object_and_primitive_id;

	std::cout << "Construct AABB tree...";
	Tree tree(m_pPolyhedron->edges_begin(),m_pPolyhedron->edges_end());
	std::cout << "done." << std::endl;

    QTime time;
    time.start();
	std::cout << "Generate edge points: ";

	unsigned int nb = 0;
	unsigned int nb_planes = 0;
	while(nb < nb_points)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);

		std::list<Object_and_primitive_id> intersections;
		tree.all_intersections(plane,std::back_inserter(intersections));
		nb_planes++;

		std::list<Object_and_primitive_id>::iterator it;
		for(it = intersections.begin();
			it != intersections.end();
			it++)
		{
			Object_and_primitive_id op = *it;
			CGAL::Object object = op.first;
			Point point;
			if(CGAL::assign(point,object))
			{
				m_points.push_back(point);
				nb++;
			}
		}
	}
	std::cout << nb_planes << " planes launched, " << time.elapsed() << " ms." << std::endl;
}

void Scene::benchmark_distances()
{
	std::cout << "to be implemented" << std::endl;
}

void Scene::benchmark_intersections()
{
	typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;

	std::cout << "Construct AABB tree...";
	Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done." << std::endl;

    QTime time;
    time.start();
	std::cout << "Benchmark do_intersect" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Ray ray(p,q);
		tree.do_intersect(ray);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;

	// with line
	nb = 0;
    time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Line line(p,q);
		tree.do_intersect(line);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;

	// with segment
	nb = 0;
    time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Segment segment(p,q);
		tree.do_intersect(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;
}

void Scene::toggle_view_poyhedron()
{
	m_view_polyhedron = !m_view_polyhedron;
}

void Scene::toggle_view_segments()
{
	m_view_segments = !m_view_segments;
}

void Scene::toggle_view_points()
{
	m_view_points = !m_view_points;
}


