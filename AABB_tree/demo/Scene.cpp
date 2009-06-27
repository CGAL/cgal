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

#include "render.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

Scene::Scene()
{
	m_pPolyhedron = NULL;
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
	if(m_pPolyhedron != NULL)
	{
		::glDisable(GL_LIGHTING);
		::glColor3ub(0,0,0);
		gl_render_edges(*m_pPolyhedron);
	}

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

	m_points.clear();

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


