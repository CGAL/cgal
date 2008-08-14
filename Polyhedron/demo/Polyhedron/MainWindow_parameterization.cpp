#include <QApplication>
#include <QTime>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>

#include <iostream>

void MainWindow::parameterize(const Parameterization_method method)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  // get active polyhedron
  int index = getSelectedPolygonIndex();
  Polyhedron* pMesh = scene->polyhedron(index);

  // parameterize
  QTime time;
  time.start();
  typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron> Adaptor;
  Adaptor adaptor(*pMesh);  

  bool success;
  switch(method)
  {
    case PARAM_MVC:
    {
      std::cerr << "Parameterize (MVC)...";
      typedef CGAL::Mean_value_coordinates_parameterizer_3<Adaptor> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(adaptor,Parameterizer());
      success = err == Parameterizer::OK;
      break;
    }
    case PARAM_DCP:
    {
      std::cerr << "Parameterize (DCP)...";
      typedef CGAL::Discrete_conformal_map_parameterizer_3<Adaptor> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(adaptor,Parameterizer());
      success = err == Parameterizer::OK;
    }
  }

  if(success)
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  else
  {
    std::cout << "failure" << std::endl;
    QApplication::setOverrideCursor(Qt::ArrowCursor);
    return;
  }

  // add parameterized mesh 
  Polyhedron *pParameterization = new Polyhedron(*pMesh);
  Polyhedron::Vertex_iterator it1, it2;
  for(it1 = pMesh->vertices_begin(), 
      it2 = pParameterization->vertices_begin();
      it1 != pMesh->vertices_end(),
      it2 != pParameterization->vertices_end();
      it1++, it2++)
  {
    // (u,v) pair is stored in any halfedge
    FT u = adaptor.info(it1->halfedge())->uv().x();
    FT v = adaptor.info(it1->halfedge())->uv().y();
    it2->point() = Point(u-0.5,v-0.5,0.0);
  }

  scene->addPolyhedron(pParameterization,
    tr("%1 (parameterization)").arg(scene->polyhedronName(index)),
    Qt::magenta,
    scene->isPolyhedronActivated(index),
    scene->polyhedronRenderingMode(index));

  QApplication::setOverrideCursor(Qt::ArrowCursor);
}

void MainWindow::on_actionMVC_triggered()
{
  parameterize(PARAM_MVC);
}

void MainWindow::on_actionDCP_triggered()
{
  std::cerr << "DCP...";
  parameterize(PARAM_DCP);
}
