#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <QApplication>
#include <QInputDialog>
#include <QTime>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

void MainWindow::on_actionSimplify_triggered()
{
  if(onePolygonIsSelected())
  {
    // get selected polyhedron
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);

    // get option (#edges)
    bool ok;
    const unsigned int nb_edges = 
      QInputDialog::getInteger(this, tr("Stop condition"),
      tr("Number of edges:"),
      pMesh->size_of_halfedges () / 4, // default value: current #edges / 2 
      3, // min = one triangle
      pMesh->size_of_halfedges(), // max #edges
      1, // step for the spinbox
      &ok);

    // check user cancellation
    if(!ok)
      return;

    // simplify
    QTime time;
    time.start();
    std::cout << "Simplify...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    namespace SMS = CGAL::Surface_mesh_simplification;
    SMS::Count_stop_predicate< Polyhedron > stop(nb_edges); // target #edges
    SMS::edge_collapse( *pMesh, stop,
      CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,*pMesh))
      .edge_index_map(boost::get(CGAL::edge_external_index,*pMesh)));
    std::cout << "ok (" << time.elapsed() << " ms, " 
      << pMesh->size_of_halfedges() / 2 << " edges)" << std::endl;

    // update scene
    scene->polyhedronChanged(index);
    QApplication::restoreOverrideCursor();
  }
}
