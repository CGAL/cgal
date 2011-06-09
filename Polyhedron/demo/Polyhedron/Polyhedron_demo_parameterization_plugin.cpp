#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_textured_polyhedron_item.h"
#include "Textured_polyhedron_type.h"
#include "Polyhedron_type.h"

#include <QTime>

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>

#include <CGAL/Textured_polyhedron_builder.h>

#include <iostream>

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

typedef Kernel::FT FT;

class Polyhedron_demo_parameterization_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMVC"
                         << "actionDCP";
  }

public slots:
  void on_actionMVC_triggered();
  void on_actionDCP_triggered();

protected:
  enum Parameterization_method { PARAM_MVC, PARAM_DCP };
  void parameterize(Parameterization_method method);
}; // end Polyhedron_demo_parameterization_plugin



void Polyhedron_demo_parameterization_plugin::parameterize(const Parameterization_method method)
{
  // get active polyhedron
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if(!poly_item)
    return;

  Polyhedron* pMesh = poly_item->polyhedron();
  if(!pMesh)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

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

  Scene_item* new_item = new Scene_textured_polyhedron_item(pTex_polyhedron);

  new_item->setName(tr("%1 (parameterized)").arg(poly_item->name()));
  new_item->setColor(Qt::white);
  new_item->setRenderingMode(poly_item->renderingMode());

  scene->addItem(new_item);

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_parameterization_plugin::on_actionMVC_triggered()
{
  parameterize(PARAM_MVC);
}

void Polyhedron_demo_parameterization_plugin::on_actionDCP_triggered()
{
  std::cerr << "DCP...";
  parameterize(PARAM_DCP);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_parameterization_plugin, Polyhedron_demo_parameterization_plugin)

#include "Polyhedron_demo_parameterization_plugin.moc"
