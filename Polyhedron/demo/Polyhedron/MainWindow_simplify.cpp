#include "MainWindow.h"
#include "Scene.h"

#include <QInputDialog>

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

		// simplify
                bool ok;
		const unsigned int nb_edges = 
                  QInputDialog::getInteger(this, tr("Stop condition"),
                                           tr("Number of edges:"),
                                           pMesh->size_of_halfedges () / 4, 
                                           // current value: nb of edges /2 
                                           10, // min
                                           pMesh->size_of_halfedges(), // max
                                           1, // step for the spinbox
                                           &ok);
                if(!ok) return;

1000; // TODO: should be an option 
		namespace SMS = CGAL::Surface_mesh_simplification;

		// wait cursor
		QApplication::setOverrideCursor(Qt::WaitCursor);

		SMS::Count_stop_predicate< Polyhedron > stop(nb_edges); // target # edges
		SMS::edge_collapse( *pMesh, stop,
		                    CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,*pMesh))
				                       .edge_index_map(boost::get(CGAL::edge_external_index,*pMesh)));

		// update scene
		scene->polyhedronChanged(index);

		// default cursor
		QApplication::restoreOverrideCursor();
	}
}
