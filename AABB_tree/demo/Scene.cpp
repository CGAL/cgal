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
#include <QApplication>

#include "render.h"
#include <CGAL/IO/Polyhedron_iostream.h>


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
	if(m_pPolyhedron != NULL)
		gl_render_facets(*m_pPolyhedron);
}


