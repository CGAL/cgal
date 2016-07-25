#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_textured_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Textured_polyhedron_type.h"
#include "Polyhedron_type.h"
#include "Messages_interface.h"

#include <QTime>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>

#include <CGAL/Textured_polyhedron_builder.h>

#include <iostream>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/boost/graph/properties.h>

typedef Kernel::FT FT;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
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
  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m_i)
  {
      mw = mainWindow;
      scene = scene_interface;
      message = m_i;
      QAction* actionMVC = new QAction("Mean Value Coordinates", mw);
      QAction* actionDCP = new QAction ("Discrete Conformal Map", mw);
      QAction* actionLSC = new QAction("Least Square Conformal Map", mw);
      actionMVC->setObjectName("actionMVC");
      actionDCP->setObjectName("actionDCP");
      actionLSC->setObjectName("actionLSC");

      _actions << actionMVC
               << actionDCP
               << actionLSC;
      autoConnectActions();
      Q_FOREACH(QAction *action, _actions)
        action->setProperty("subMenuName",
                            "Triangulated Surface Mesh Parameterization");


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
private:
  QList<QAction*> _actions;
  Messages_interface* message;
}; // end Polyhedron_demo_parameterization_plugin



void Polyhedron_demo_parameterization_plugin::parameterize(const Parameterization_method method)
{
  typedef Kernel::Point_2                                               Point_2;
  typedef boost::graph_traits<Polyhedron>::edge_descriptor              P_edge_descriptor;
  typedef boost::graph_traits<Polyhedron>::halfedge_descriptor          P_halfedge_descriptor;
  typedef boost::graph_traits<Polyhedron>::vertex_descriptor            P_vertex_descriptor;


  typedef CGAL::Unique_hash_map<P_halfedge_descriptor,Point_2>          UV_uhm;
  typedef CGAL::Unique_hash_map<P_edge_descriptor,bool>                 Seam_edge_uhm;
  typedef CGAL::Unique_hash_map<vertex_descriptor,bool>                 Seam_vertex_uhm;

  typedef boost::associative_property_map<UV_uhm>                       UV_pmap;
  typedef boost::associative_property_map<Seam_edge_uhm>                Seam_edge_pmap;
  typedef boost::associative_property_map<Seam_vertex_uhm>              Seam_vertex_pmap;


  typedef CGAL::Seam_mesh<Polyhedron, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;
  typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor           halfedge_descriptor;

  // get active polyhedron
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if(!poly_item)
    return;
  poly_item->polyhedron()->normalize_border();
  if(poly_item->polyhedron()->size_of_border_halfedges()==0)
    message->warning("The polyhedron has no border, therefore the Parameterization cannot apply.");
  Polyhedron* pMesh = poly_item->polyhedron();
  if(!pMesh)
    return;
  Scene_polyhedron_selection_item* sel_item = NULL;
  bool is_seamed = false;
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices())
  {
    sel_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(id));
    if(!sel_item)
      continue;
    if(sel_item->selected_edges.empty()
       || !sel_item->selected_facets.empty()
       || !sel_item->selected_vertices.empty())
      continue;
    is_seamed = true;
  }
  // Two property maps to store the seam edges and vertices
  Seam_edge_uhm seam_edge_uhm(false);
  Seam_edge_pmap seam_edge_pm(seam_edge_uhm);

  Seam_vertex_uhm seam_vertex_uhm(false);
  Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  ///////////////////////////////////

  // parameterize
  QTime time;
  time.start();
  //fill pmaps

  halfedge_descriptor smhd;
  if(!is_seamed)
  {
    for(Polyhedron::Halfedge_iterator hit = pMesh->halfedges_begin(); hit != pMesh->halfedges_end(); ++hit)
    {
      P_vertex_descriptor svd(source(hit, *pMesh)), tvd(target(hit, *pMesh));
      P_edge_descriptor ed = edge(svd, tvd,*pMesh).first;
      if(is_border(ed,*pMesh)){
        if(smhd == boost::graph_traits<Polyhedron>::null_halfedge()){
          smhd = halfedge(edge(svd,tvd,*pMesh).first,*pMesh);
          if (smhd == boost::graph_traits<Polyhedron>::null_halfedge())
            smhd = opposite(smhd, *pMesh);
        }
        break;
      }
    }
  }
  if(is_seamed)
  {
    BOOST_FOREACH(P_edge_descriptor ed, sel_item->selected_edges)
    {
      P_halfedge_descriptor hd = halfedge(ed, *pMesh);
      P_vertex_descriptor svd(source(hd, *pMesh)), tvd(target(hd, *pMesh));
      if(! is_border(ed,*pMesh)){
        put(seam_edge_pm, ed, true);
        put(seam_vertex_pm, svd, true);
        put(seam_vertex_pm, tvd, true);
        if(smhd == boost::graph_traits<Polyhedron>::null_halfedge()){
          smhd = hd;
        }
      }
    }
  }
  Seam_mesh *sMesh = new Seam_mesh(*pMesh, seam_edge_pm, seam_vertex_pm);
  UV_uhm uv_uhm;
  UV_pmap uv_pm(uv_uhm);

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,*sMesh); // a halfedge on the virtual border
  bool success = false;

  switch(method)
  {
  case PARAM_MVC:
    {
      std::cout << "Parameterize (MVC)...";
      typedef CGAL::Mean_value_coordinates_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == Parameterizer::OK;
      break;
    }
  case PARAM_DCP:
    {
      std::cout << "Parameterize (DCP)...";
      typedef CGAL::Discrete_conformal_map_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == Parameterizer::OK;
      break;
    }
  case PARAM_LSC:
    {
      std::cout << "Parameterize (LSC)...";
      typedef CGAL::LSCM_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == Parameterizer::OK;
      break;
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
    FT u = uv_pm[it1->halfedge()].x();
    FT v = uv_pm[it1->halfedge()].y();
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
  std::cerr << "MVC...";
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
