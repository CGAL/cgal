#include "MainWindow.h"
#include "Scene.h"

#include "Polyhedron_type.h"
#include <CGAL/Subdivision_method_3.h>

void MainWindow::on_actionLoop_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
	// PA: does not compile if uncomment this
  // CGAL::Subdivision_method_3::Loop_subdivision(*poly, 1);
  scene->polyhedronChanged(poly);
}

void MainWindow::on_actionCatmullClark_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  //CGAL::Subdivision_method_3::CatmullClark_subdivision(*poly, 1);
  scene->polyhedronChanged(poly);
}

void MainWindow::on_actionSqrt3_triggered()
{
  Polyhedron* poly = getSelectedPolygon();
  if(!poly) return;
  //CGAL::Subdivision_method_3::Sqrt3_subdivision(*poly, 1);
  scene->polyhedronChanged(poly);
}



