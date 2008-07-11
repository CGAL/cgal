#include "MainWindow.h"
#include "Scene.h"

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

void MainWindow::on_actionSimplify_triggered()
{
  if(onePolygonIsSelected())
  {
    int index = getSelectedPolygonIndex();
    Polyhedron* pMesh = scene->polyhedron(index);
    // simplify
    //unsigned int nb_edges = 1000; // TODO: should be an option 
    //namespace SMS = CGAL::Surface_mesh_simplification ;
    //SMS::Count_stop_predicate< Polyhedron > stop(nb_edges); // target # edges
    //SMS::edge_collapse( *pMesh, stop,
    //                     CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,*pMesh))
    //		                       .edge_index_map(boost::get(CGAL::edge_external_index,*pMesh)));

    // recompute normals
    pMesh->compute_normals();

    // Tell the scene that polyhedron #index has been changed
    scene->polyhedronChanged(index);
  }
}
