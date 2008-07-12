#include "MainWindow.h"
#include "Scene.h"

#include <CGAL/convex_hull_3.h>

void MainWindow::on_actionConvexHull_triggered()
{
  if(onePolygonIsSelected())
  {
    int index = getSelectedPolygonIndex();

    // get active polyhedron
    Polyhedron* pMesh = scene->polyhedron(index);

    // add convex hull as new polyhedron
    Polyhedron *pConvex_hull = new Polyhedron;
    CGAL::convex_hull_3(pMesh->points_begin(),pMesh->points_end(),*pConvex_hull);
    //pConvex_hull->compute_normals();

    scene->addPolyhedron(pConvex_hull,
                         tr("%1 (convex hull)").arg(scene->polyhedronName(index)),
                         scene->polyhedronColor(index),
                         scene->isPolyhedronActivated(index),
                         scene->polyhedronRenderingMode(index));
  }
}
