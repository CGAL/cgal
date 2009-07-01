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

Scene::Scene()
{
	m_pPolyhedron = NULL;

	// view options
	m_view_points = true;
	m_view_segments = true;
	m_view_polyhedron = true;

	// distance functions
	m_max_distance_function = (FT)0.0;
	m_signed_distance_function = false;

	m_red_ramp.build_red();
	m_blue_ramp.build_blue();
}

Scene::~Scene()
{
	delete m_pPolyhedron;
}

int
Scene::open(QString filename)
{
	QTextStream cerr(stderr);
	cerr << QString("Opening file \"%1\"...\n").arg(filename);
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

	std::cerr << " Ok.\n";
	return 0;
}

void Scene::draw()
{
	if(m_view_polyhedron)
		draw_polyhedron();

	if(m_view_points)
		draw_points();

	if(m_view_segments)
		draw_segments();

	if(m_signed_distance_function)
		draw_signed_distance_function();
	else
		draw_unsigned_distance_function();
}

void Scene::draw_polyhedron()
{
	// draw black edges
	if(m_pPolyhedron != NULL)
	{
		::glDisable(GL_LIGHTING);
		::glColor3ub(0,0,0);
		::glLineWidth(1.0f);
		gl_render_edges(*m_pPolyhedron);
	}
}

void Scene::draw_segments()
{
	if(m_segments.size() != 0)
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

void Scene::draw_points()
{
	// draw red points
	if(m_points.size() != 0)
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
}

void Scene::draw_unsigned_distance_function()
{
	if(m_max_distance_function == (FT)0.0)
		return;

	::glDisable(GL_LIGHTING);
	::glShadeModel(GL_SMOOTH);
	::glBegin(GL_QUADS);
	int i,j;
	const int nb_quads = 99;
	for(i=0;i<nb_quads;i++)
		for(j=0;j<nb_quads;j++)
		{
			Point_distance& pd00 = m_distance_function[i][j];
			Point_distance& pd01 = m_distance_function[i][j+1];
			Point_distance& pd11 = m_distance_function[i+1][j+1];
			Point_distance& pd10 = m_distance_function[i+1][j];
			Point& p00 = pd00.first;
			Point& p01 = pd01.first;
			Point& p11 = pd11.first;
			Point& p10 = pd10.first;
			FT& d00 = pd00.second;
			FT& d01 = pd01.second;
			FT& d11 = pd11.second;
			FT& d10 = pd10.second;
			unsigned int i00 = 255-(unsigned int)(255.0 * d00 / m_max_distance_function);
			unsigned int i01 = 255-(unsigned int)(255.0 * d01 / m_max_distance_function);
			unsigned int i11 = 255-(unsigned int)(255.0 * d11 / m_max_distance_function);
			unsigned int i10 = 255-(unsigned int)(255.0 * d10 / m_max_distance_function);
			::glColor3ub(m_thermal_ramp.r(i00),m_thermal_ramp.g(i00),m_thermal_ramp.b(i00));
			::glVertex3d(p00.x(),p00.y(),p00.z());
			::glColor3ub(m_thermal_ramp.r(i01),m_thermal_ramp.g(i01),m_thermal_ramp.b(i01));
			::glVertex3d(p01.x(),p01.y(),p01.z());
			::glColor3ub(m_thermal_ramp.r(i11),m_thermal_ramp.g(i11),m_thermal_ramp.b(i11));
			::glVertex3d(p11.x(),p11.y(),p11.z());
			::glColor3ub(m_thermal_ramp.r(i10),m_thermal_ramp.g(i10),m_thermal_ramp.b(i10));
			::glVertex3d(p10.x(),p10.y(),p10.z());
		}
		::glEnd();
}

void Scene::draw_signed_distance_function()
{
	if(m_max_distance_function == (FT)0.0)
		return;

	::glDisable(GL_LIGHTING);
	::glShadeModel(GL_SMOOTH);
	::glBegin(GL_QUADS);
	int i,j;
	const int nb_quads = 99;
	for(i=0;i<nb_quads;i++)
		for(j=0;j<nb_quads;j++)
		{
			Point_distance& pd00 = m_distance_function[i][j];
			Point_distance& pd01 = m_distance_function[i][j+1];
			Point_distance& pd11 = m_distance_function[i+1][j+1];
			Point_distance& pd10 = m_distance_function[i+1][j];
			Point& p00 = pd00.first;
			Point& p01 = pd01.first;
			Point& p11 = pd11.first;
			Point& p10 = pd10.first;
			FT& d00 = pd00.second;
			FT& d01 = pd01.second;
			FT& d11 = pd11.second;
			FT& d10 = pd10.second;

			// determines grey level
			unsigned int i00 = 255-(unsigned)(255.0 * (double)std::fabs(d00) / m_max_distance_function);
			unsigned int i01 = 255-(unsigned)(255.0 * (double)std::fabs(d01) / m_max_distance_function);
			unsigned int i11 = 255-(unsigned)(255.0 * (double)std::fabs(d11) / m_max_distance_function);
			unsigned int i10 = 255-(unsigned)(255.0 * (double)std::fabs(d10) / m_max_distance_function);

			// assembles one quad
			if(d00 > 0.0)
				::glColor3ub(m_red_ramp.r(i00),m_red_ramp.g(i00),m_red_ramp.b(i00));
			else
				::glColor3ub(m_blue_ramp.r(i00),m_blue_ramp.g(i00),m_blue_ramp.b(i00));
			::glVertex3d(p00.x(),p00.y(),p00.z());

			if(d01 > 0.0)
				::glColor3ub(m_red_ramp.r(i01),m_red_ramp.g(i01),m_red_ramp.b(i01));
			else
				::glColor3ub(m_blue_ramp.r(i01),m_blue_ramp.g(i01),m_blue_ramp.b(i01));
			::glVertex3d(p01.x(),p01.y(),p01.z());

			if(d11 > 0)
				::glColor3ub(m_red_ramp.r(i11),m_red_ramp.g(i11),m_red_ramp.b(i11));
			else
				::glColor3ub(m_blue_ramp.r(i11),m_blue_ramp.g(i11),m_blue_ramp.b(i11));
			::glVertex3d(p11.x(),p11.y(),p11.z());

			if(d10 > 0)
				::glColor3ub(m_red_ramp.r(i10),m_red_ramp.g(i10),m_red_ramp.b(i10));
			else
				::glColor3ub(m_blue_ramp.r(i10),m_blue_ramp.g(i10),m_blue_ramp.b(i10));
			::glVertex3d(p10.x(),p10.y(),p10.z());
		}
		::glEnd();
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
	std::cout << nb_lines << " line queries, " << time.elapsed() << " ms." << std::endl;
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
	std::cout << nb_planes << " plane queries, " << time.elapsed() << " ms." << std::endl;
}

void Scene::benchmark_distances()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	// TODO: check what the KD-tree is doing in there for large models.
	std::cout << "One single call to initialize KD-tree" << std::endl;
	tree.closest_point(CGAL::ORIGIN);

	std::cout << "------- Now, the real benchmark -------" << std::endl;
	bench_closest_point(tree);
	bench_squared_distance(tree);
	bench_closest_point_and_primitive(tree);
}

void Scene::benchmark_intersections()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	bench_do_intersect(tree);
	bench_nb_intersections(tree);
	bench_any_intersection(tree);
	bench_all_intersections(tree);
	bench_all_intersected_primitives(tree);
}

void Scene::unsigned_distance_function()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	m_max_distance_function = (FT)0.0;
	int i,j;
	for(i=0;i<100;i++)
	{
		FT x = -0.5 + (FT)i/100.0;
		for(j=0;j<100;j++)
		{
			FT y = -0.5 + (FT)j/100.0;
			Point query(x,y,0.0);
			FT sq_distance = tree.squared_distance(query);
			FT distance = std::sqrt(sq_distance);
			m_distance_function[i][j] = Point_distance(query,distance);
			m_max_distance_function = distance > m_max_distance_function ?
				distance : m_max_distance_function;
		}
	}
	m_signed_distance_function = false;
}

void Scene::unsigned_distance_function_to_edges()
{
	typedef CGAL::AABB_polyhedron_segment_primitive<Kernel,Polyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Edge_tree;

	QTime time;
	time.start();
	std::cout << "Construct AABB tree from edges...";
	Edge_tree tree(m_pPolyhedron->edges_begin(),m_pPolyhedron->edges_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	m_max_distance_function = (FT)0.0;
	int i,j;
	for(i=0;i<100;i++)
	{
		FT x = -0.5 + (FT)i/100.0;
		for(j=0;j<100;j++)
		{
			FT y = -0.5 + (FT)j/100.0;
			Point query(x,y,0.0);
			FT sq_distance = tree.squared_distance(query);
			FT distance = std::sqrt(sq_distance);
			m_distance_function[i][j] = Point_distance(query,distance);
			m_max_distance_function = distance > m_max_distance_function ?
distance : m_max_distance_function;
		}
	}
	m_signed_distance_function = false;
}

void Scene::signed_distance_function()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	m_max_distance_function = (FT)0.0;
	Vector vec = random_vector();
	int i,j;
	for(i=0;i<100;i++)
	{
		FT x = -0.5 + (FT)i/100.0;
		for(j=0;j<100;j++)
		{
			// compute distance
			FT y = -0.5 + (FT)j/100.0;
			Point query(x,y,0.0);
			FT sq_distance = tree.squared_distance(query);
			FT unsigned_distance = std::sqrt(sq_distance);

			// get sign through ray casting (random vector)
			Ray ray(query,vec);
			unsigned int nbi = tree.number_of_intersected_primitives(ray);
			FT sign = nbi%2 == 0 ? (FT)1.0 : (FT)-1.0;
			FT signed_distance = sign * unsigned_distance;

			m_distance_function[i][j] = Point_distance(query,signed_distance);
			m_max_distance_function = unsigned_distance > m_max_distance_function ?
unsigned_distance : m_max_distance_function;
		}
	}
	m_signed_distance_function = true;
}

void Scene::bench_squared_distance(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark squared distance" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.squared_distance(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}


void Scene::bench_closest_point(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark closest point" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.closest_point(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}

void Scene::bench_closest_point_and_primitive(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark closest point and primitive" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.closest_point_and_primitive(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}


void Scene::bench_all_intersected_primitives(Facet_tree& tree)
{
	std::list<Primitive_id> primitive_ids;

	QTime time;
	time.start();
	std::cout << "Benchmark all_intersected_primitives" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Ray ray(p,q);
		tree.all_intersected_primitives(ray,std::back_inserter(primitive_ids));
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;
	primitive_ids.clear();

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Line line(p,q);
		tree.all_intersected_primitives(line,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;
	primitive_ids.clear();

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Segment segment(p,q);
		tree.all_intersected_primitives(segment,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);
		tree.all_intersected_primitives(plane,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}


void Scene::bench_do_intersect(Facet_tree& tree)
{
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

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);
		tree.do_intersect(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_nb_intersections(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark number_of_intersected_primitives" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Ray ray(p,q);
		tree.number_of_intersected_primitives(ray);
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
		tree.number_of_intersected_primitives(line);
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
		tree.number_of_intersected_primitives(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);
		tree.number_of_intersected_primitives(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_any_intersection(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark any_intersection" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Ray ray(p,q);
		tree.any_intersection(ray);
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
		tree.any_intersection(line);
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
		tree.any_intersection(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);
		tree.any_intersection(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_all_intersections(Facet_tree& tree)
{
	std::list<Object_and_primitive_id> intersections;

	QTime time;
	time.start();
	std::cout << "Benchmark all_intersections" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Ray ray(p,q);
		tree.all_intersections(ray,std::back_inserter(intersections));
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;
	intersections.clear();

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Point q = random_point();
		Line line(p,q);
		tree.all_intersections(line,std::back_inserter(intersections));
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
		tree.all_intersections(segment,std::back_inserter(intersections));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Point p = random_point();
		Vector vec = random_vector();
		Plane plane(p,vec);
		tree.all_intersections(plane,std::back_inserter(intersections));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
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


