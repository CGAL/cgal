#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

void MainWindow::on_actionInsideOut_triggered()
{
  if(onePolygonIsSelected())
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

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
