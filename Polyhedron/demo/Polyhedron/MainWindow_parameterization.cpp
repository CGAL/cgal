#include <QApplication>
#include <QTime>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>

#include <iostream>

void MainWindow::on_actionMVC_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  // get active polyhedron
  int index = getSelectedPolygonIndex();
  Polyhedron* pMesh = scene->polyhedron(index);

  // parameterize
  QTime time;
  time.start();
  typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron> Adaptor;
  Adaptor adaptor(*pMesh);  
  std::cerr << "Parameterize...";
  typedef CGAL::Parameterizer_traits_3<Adaptor> Parameterizer;
  Parameterizer::Error_code err = CGAL::parameterize(adaptor);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  QApplication::setOverrideCursor(Qt::ArrowCursor);
}
