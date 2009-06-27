#include <QApplication>
#include "MainWindow.h"
#include "Scene.h"
#include "types.h"

void MainWindow::on_actionInsideOut_triggered()
{
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // get selected polyhedron
    Polyhedron* pMesh = m_pScene->polyhedron();

    // inside out
    pMesh->inside_out();

    // default cursor
    QApplication::restoreOverrideCursor();
}
