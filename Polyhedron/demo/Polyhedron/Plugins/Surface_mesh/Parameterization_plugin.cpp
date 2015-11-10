#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_textured_polyhedron_item.h"
#include "Textured_polyhedron_type.h"
#include "Polyhedron_type.h"

#include <QTime>

#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>

#include <CGAL/Textured_polyhedron_builder.h>

#include <iostream>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

typedef Kernel::FT FT;
using namespace CGAL::Three;
class Polyhedron_demo_parameterization_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMVC"
                         << "actionDCP"
                         << "actionLSC";
  }

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface)
  {
      mw = mainWindow;
      scene = scene_interface;
      actions_map["actionMVC"] = new QAction("MVC", mw);
      actions_map["actionMVC"]->setProperty("subMenuName",
        "Planar Parameterization of Triangulated Surface Meshes");

      actions_map["actionDCP"] = new QAction ("DCP", mw);
      actions_map["actionDCP"]->setProperty("subMenuName",
        "Planar Parameterization of Triangulated Surface Meshes");

      actions_map["actionLSC"] = new QAction("LSC", mw);
      actions_map["actionLSC"]->setProperty("subMenuName",
        "Planar Parameterization of Triangulated Surface Meshes");

      connect(actions_map["actionMVC"], SIGNAL(triggered()),
              this, SLOT(on_actionMVC_triggered()));
      connect(actions_map["actionDCP"], SIGNAL(triggered()),
              this, SLOT(on_actionDCP_triggered()));
      connect(actions_map["actionLSC"], SIGNAL(triggered()),
              this, SLOT(on_actionLSC_triggered()));
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionMVC_triggered();
  void on_actionDCP_triggered();
  void on_actionLSC_triggered();

protected:
  enum Parameterization_method { PARAM_MVC, PARAM_DCP, PARAM_LSC };
  void parameterize(Parameterization_method method);
}; // end Polyhedron_demo_parameterization_plugin



void Polyhedron_demo_parameterization_plugin::parameterize(const Parameterization_method method)
{
  // get active polyhedron
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
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
  case PARAM_LSC:
    {
      std::cout << "Parameterize (LSC)...";
      typedef CGAL::LSCM_parameterizer_3<Adaptor> Parameterizer;
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
    it1 != pMesh->vertices_end() &&
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

  poly_item->setVisible(false);
  scene->itemChanged(index);
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

void Polyhedron_demo_parameterization_plugin::on_actionLSC_triggered()
{
  std::cerr << "LSC...";
  parameterize(PARAM_LSC);
}

#include "Parameterization_plugin.moc"
