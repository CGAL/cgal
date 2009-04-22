#include <QTime>
#include <QApplication>

#include "MainWindow.h"
#include "Scene.h"

#include "Polyhedron_type.h"
#include <CGAL/Subdivision_method_3.h>

void MainWindow::on_actionLoop_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  QTime time;
  time.start();
  std::cout << "Loop subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Loop_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  scene->polyhedronChanged(poly);
}

void MainWindow::on_actionCatmullClark_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  QTime time;
  time.start();
  std::cout << "Catmull-Clark subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::CatmullClark_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  scene->polyhedronChanged(poly);
}

void MainWindow::on_actionSqrt3_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  QTime time;
  time.start();
  std::cout << "Sqrt3 subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Sqrt3_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  scene->polyhedronChanged(poly);
}
