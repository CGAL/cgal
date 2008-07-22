#include "MainWindow.h"
#include "Scene.h"

void MainWindow::on_actionInsideOut_triggered()
{
  if(onePolygonIsSelected())
  {
    // get selected polyhedron
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);

		// inside out
		pMesh->inside_out();

    // update scene
    scene->polyhedronChanged(index);
    QApplication::restoreOverrideCursor();
  }
}
