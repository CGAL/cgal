#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_textured_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Textured_polyhedron_type.h"
#include "Polyhedron_type.h"

#include "Scene.h"
#include <QTime>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QPen>
#include <QDockWidget>

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
#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "ui_Parameterization_widget.h"


class Navigation : public CGAL::Qt::GraphicsViewNavigation
{
public:
  Navigation()
    :CGAL::Qt::GraphicsViewNavigation(),
      prev_pos(QPoint(0,0))
  {
  }
protected:
  bool eventFilter(QObject *obj, QEvent *ev)
  {
    QGraphicsView* v = qobject_cast<QGraphicsView*>(obj);
    if(v == NULL) {
      QWidget* viewport = qobject_cast<QWidget*>(obj);
      if(viewport == NULL) {
        return false;
      }
      v = qobject_cast<QGraphicsView*>(viewport->parent());
      if(v == NULL) {
        return false;
      }
    }
    switch(ev->type())
    {
    case QEvent::MouseMove: {
      QMouseEvent* me = static_cast<QMouseEvent*>(ev);
      if(is_dragging)
      {
        qreal dir[2] = {v->mapToScene(me->pos()).x() - prev_pos.x(),
                        v->mapToScene(me->pos()).y() - prev_pos.y()};

        v->translate(dir[0],dir[1]);
        v->update();
      }
      prev_pos = v->mapToScene(me->pos());
      break;
    }

    case QEvent::MouseButtonPress: {
      is_dragging = true;
      break;
    }
    case QEvent::MouseButtonRelease: {
      is_dragging = false;
      break;
    }
    case QEvent::Wheel: {
      QWheelEvent* event = static_cast<QWheelEvent*>(ev);
      QPointF old_pos = v->mapToScene(event->pos());
      if(event->delta() <0)
        v->scale(1.2, 1.2);
      else
        v->scale(0.8, 0.8);
      QPointF new_pos = v->mapToScene(event->pos());
      QPointF delta = new_pos - old_pos;
      v->translate(delta.x(), delta.y());
      v->update();
      break;
    }

    case QEvent::MouseButtonDblClick: {
      v->fitInView(v->scene()->itemsBoundingRect(), Qt::KeepAspectRatio);
      break;
    }
    default:
      CGAL::Qt::GraphicsViewNavigation::eventFilter(obj, ev);
    }
    return false;
  }
private:
  bool is_dragging;
  QPointF prev_pos;
};
typedef Kernel::FT                                                    FT;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor            vertex_descriptor;
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

class UVItem : public QGraphicsItem
{
public :
  UVItem(Textured_polyhedron* t_m, QRectF brect)
    :QGraphicsItem(),
       texMesh(t_m), bounding_rect(brect)
  {
  }

  QRectF boundingRect() const
  {
    return bounding_rect;
  }

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
  {
    QPen pen;
    pen.setColor(Qt::black);
    pen.setWidth(0);
    painter->setPen(pen);
    for( Textured_polyhedron::Edge_iterator
         ei = texMesh->edges_begin();
         ei != texMesh->edges_end();
         ++ei)
    {
      QPointF source(ei->vertex()->u(), ei->vertex()->v()), target(ei->opposite()->vertex()->u(), ei->opposite()->vertex()->v());
      painter->drawLine(source, target);
    }
  }

private:
  Textured_polyhedron* texMesh;
  QRectF bounding_rect;

};

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
      Scene* true_scene = static_cast<Scene*>(scene);
      connect(true_scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
              this, SLOT(destroyPolyline(CGAL::Three::Scene_item*)));
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
      dock_widget = new QDockWidget("UVMapping", mw);
      ui_widget.setupUi(dock_widget);
      graphics_scene = new QGraphicsScene(dock_widget);
      ui_widget.graphicsView->setScene(graphics_scene);
      ui_widget.graphicsView->setRenderHints(QPainter::Antialiasing);
      navigation = new Navigation();
      ui_widget.graphicsView->installEventFilter(navigation);
      ui_widget.graphicsView->viewport()->installEventFilter(navigation);
      addDockWidget(dock_widget);
      dock_widget->setVisible(false);

  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void closure()
  {
    dock_widget->hide();
  }
public Q_SLOTS:
  void on_actionMVC_triggered();
  void on_actionDCP_triggered();
  void on_actionLSC_triggered();
  void replacePolyline()
  {
    int id = scene->mainSelectionIndex();
    Q_FOREACH(UVItem* pl, projections)
    {
      if(pl==NULL
         || pl != projections[scene->item(id)])
        continue;
      if(!graphics_scene->items().empty())
        graphics_scene->removeItem(graphics_scene->items().first());
      graphics_scene->addItem(pl);
      ui_widget.graphicsView->fitInView(pl->boundingRect(), Qt::KeepAspectRatio);
      dock_widget->setWindowTitle(tr("UVMapping for %1").arg(scene->item(id)->name()));
      break;
    }
  }
  void destroyPolyline(CGAL::Three::Scene_item* item)
  {
    Q_FOREACH(UVItem* pli, projections)
    {
      if(projections.key(pli) != item)
        continue;
      graphics_scene->removeItem(pli);
      delete pli;
      projections.remove(item);
      break;
    }
    dock_widget->setWindowTitle(tr("UVMapping"));
  }

protected:
  enum Parameterization_method { PARAM_MVC, PARAM_DCP, PARAM_LSC };
  void parameterize(Parameterization_method method);
private:
  QList<QAction*> _actions;
  QDockWidget* dock_widget;
  Ui::Parameterization ui_widget;
  QGraphicsScene *graphics_scene;
  Navigation* navigation;
  QMap<Scene_item*, UVItem*> projections;
}; // end Polyhedron_demo_parameterization_plugin



void Polyhedron_demo_parameterization_plugin::parameterize(const Parameterization_method method)
{

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
QString new_item_name;
  switch(method)
  {
  case PARAM_MVC:
    {
      std::cout << "Parameterize (MVC)...";
      new_item_name = tr("%1 (parameterized (MVC))").arg(poly_item->name());
      typedef CGAL::Mean_value_coordinates_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == Parameterizer::OK;
      break;
    }
  case PARAM_DCP:
    {
      new_item_name = tr("%1 (parameterized (DCP))").arg(poly_item->name());
      std::cout << "Parameterize (DCP)...";
      typedef CGAL::Discrete_conformal_map_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer::Error_code err = CGAL::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == Parameterizer::OK;
      break;
    }
  case PARAM_LSC:
    {
      new_item_name = tr("%1 (parameterized (LSC))").arg(poly_item->name());
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

  Polyhedron::Edge_iterator it1;
  Textured_polyhedron::Edge_iterator it2;
  QPointF min(FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX);
  for(it1 = pMesh->edges_begin(),
      it2 = pTex_polyhedron->edges_begin();
      it1 != pMesh->edges_end()&&
      it2 != pTex_polyhedron->edges_end();
      it1++, it2++)
  {
    FT u1 = uv_pm[halfedge(source(it1, *sMesh), *sMesh)].x();
    FT v1 = uv_pm[halfedge(source(it1, *sMesh), *sMesh)].y();
    source(it2, *pTex_polyhedron)->u() = u1;
    source(it2, *pTex_polyhedron)->v() = v1;
    FT u2 = uv_pm[halfedge(target(it1, *sMesh), *sMesh)].x();
    FT v2 = uv_pm[halfedge(target(it1, *sMesh), *sMesh)].y();
    target(it2, *pTex_polyhedron)->u() = u2;
    target(it2, *pTex_polyhedron)->v() = v2;

    if(u1<min.x())
      min.setX(u1);
    if(u1>max.x())
      max.setX(u1);
    if(v1<min.y())
      min.setY(v1);
    if(v1>max.y())
      max.setY(v1);

    if(u2<min.x())
      min.setX(u2);
    if(u2>max.x())
      max.setX(u2);
    if(v2<min.y())
      min.setY(v2);
    if(v2>max.y())
      max.setY(v2);
  }
  UVItem *projection
      = new UVItem(pTex_polyhedron,QRectF(min, max));
  Scene_textured_polyhedron_item* new_item = new Scene_textured_polyhedron_item(pTex_polyhedron);

  new_item->setName(new_item_name);
  new_item->setColor(Qt::white);
  new_item->setRenderingMode(poly_item->renderingMode());
  connect( new_item, SIGNAL(selectionChanged()),
           this, SLOT(replacePolyline()) );
  poly_item->setVisible(false);
  scene->itemChanged(index);
  scene->addItem(new_item);

  if(!graphics_scene->items().empty())
    graphics_scene->removeItem(graphics_scene->items().first());
  graphics_scene->addItem(projection);
  projections[new_item] = projection;
  if(dock_widget->isHidden())
    dock_widget->setVisible(true);
  dock_widget->setWindowTitle(tr("UVMapping for %1").arg(new_item->name()));
  ui_widget.graphicsView->fitInView(projection->boundingRect(), Qt::KeepAspectRatio);
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
