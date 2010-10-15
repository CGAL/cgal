#ifdef CGAL_POLYHEDRON_DEMO_USE_PARAMETRIZATION
#ifdef CGAL_TAUCS_ENABLED

#include <QApplication>
#include <QTime>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"
#include "Textured_polyhedron_type.h"

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>

#include <CGAL/Textured_polyhedron_builder.h>

#include <iostream>

void MainWindow::parameterize(const Parameterization_method method)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  // get active polyhedron
  int index = getSelectedSceneItemIndex();
  Polyhedron* pMesh = scene->polyhedron(index);
  if(pMesh == NULL)
  {
    QApplication::restoreOverrideCursor();
    return;
  }

  // parameterize
  QTime time;
  time.start();
  typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron> Adaptor;
  Adaptor adaptor(*pMesh);  

  bool success = false;
  switch(method)
  {
  case PARAM_MVC:
    {
      std::cout << "Parameterize (MVC)...";
      typedef CGAL::Mean_value_coordinates_parameterizer_3<Adaptor> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(adaptor,Parameterizer());
      success = err == Parameterizer::OK;
      break;
    }
  case PARAM_DCP:
    {
      std::cout << "Parameterize (DCP)...";
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
    QApplication::restoreOverrideCursor();
    return;
  }

  // add textured polyhedon to the scene
  Textured_polyhedron *pTex_polyhedron = new Textured_polyhedron();
  Textured_polyhedron_builder<Polyhedron,Textured_polyhedron,Kernel> builder;
  builder.run(*pMesh,*pTex_polyhedron);
  pTex_polyhedron->compute_normals();

  Polyhedron::Vertex_iterator it1;
  Textured_polyhedron::Vertex_iterator it2;
  for(it1 = pMesh->vertices_begin(), 
    it2 = pTex_polyhedron->vertices_begin();
    it1 != pMesh->vertices_end(),
    it2 != pTex_polyhedron->vertices_end();
  it1++, it2++)
  {
    // (u,v) pair is stored per halfedge
    FT u = adaptor.info(it1->halfedge())->uv().x();
    FT v = adaptor.info(it1->halfedge())->uv().y();
    it2->u() = u;
    it2->v() = v;
  }

  scene->addTexPolyhedron(pTex_polyhedron,
    tr("%1 (parameterized)").arg(scene->polyhedronName(index)),
    Qt::white,
    scene->isPolyhedronVisible(index),
    scene->polyhedronRenderingMode(index));

  QApplication::restoreOverrideCursor();
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

#else // #ifdef CGAL_TAUCS_ENABLED 

#include "MainWindow.h"
#include <QMessageBox>
void MainWindow::on_actionMVC_triggered()
{
  QMessageBox::warning(this, "Function not available", "This function is not available. "
		       "You need to configure TAUCS support.");
}

void MainWindow::on_actionDCP_triggered()
{
  QMessageBox::warning(this, "Function not available", "This function is not available. "
		       "You need to configure TAUCS support.");
}

#endif // #ifdef CGAL_TAUCS_ENABLED 
#endif // CGAL_POLYHEDRON_DEMO_USE_PARAMETRIZATION
