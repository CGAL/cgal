#include "MainWindow.h"
#include "Scene.h"

#include <CGAL/Subdivision_method_3.h>

void MainWindow::on_actionCatmullClark_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  CGAL::Subdivision_method_3::CatmullClark_subdivision(*poly, 1);
  poly->compute_normals();
  viewer->updateGL();
}

void MainWindow::on_actionSqrt3_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  CGAL::Subdivision_method_3::Sqrt3_subdivision(*poly, 1);
  poly->compute_normals();
  viewer->updateGL();
}
