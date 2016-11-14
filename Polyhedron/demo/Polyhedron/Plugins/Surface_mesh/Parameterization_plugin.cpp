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
#include <Messages_interface.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Textured_polyhedron_builder.h>

#include <iostream>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include "ui_Parameterization_widget.h"
#include "ui_ARAP_dialog.h"

namespace SMP = CGAL::Surface_mesh_parameterization;

typedef boost::unordered_set<Textured_polyhedron::Base::Facet_handle> Component;
typedef std::vector<Component> Components;
struct Is_selected_property_map{
  typedef boost::graph_traits<Textured_polyhedron::Base>::edge_descriptor edge_descriptor;
  std::vector<bool>* is_selected_ptr;
  Is_selected_property_map()
    : is_selected_ptr(NULL) {}
  Is_selected_property_map(std::vector<bool>& is_selected)
    : is_selected_ptr( &is_selected) {}

  std::size_t id(edge_descriptor ed) { return ed.halfedge()->id()/2; }

  friend bool get(Is_selected_property_map map, edge_descriptor ed)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    return (*map.is_selected_ptr)[map.id(ed)];
  }

  friend void put(Is_selected_property_map map, edge_descriptor ed, bool b)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    (*map.is_selected_ptr)[map.id(ed)]=b;
  }
};

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


typedef boost::graph_traits<Textured_polyhedron::Base>::
                                         edge_descriptor              T_edge_descriptor;
typedef boost::graph_traits<Textured_polyhedron::Base>::
                                         halfedge_descriptor          T_halfedge_descriptor;
typedef boost::graph_traits<Textured_polyhedron::Base>::
                                         vertex_descriptor            T_vertex_descriptor;


typedef CGAL::Unique_hash_map<T_halfedge_descriptor,Point_2>          UV_uhm;
typedef CGAL::Unique_hash_map<T_edge_descriptor,bool>                 Seam_edge_uhm;
typedef CGAL::Unique_hash_map<T_vertex_descriptor,bool>               Seam_vertex_uhm;

typedef boost::associative_property_map<UV_uhm>                       UV_pmap;
typedef boost::associative_property_map<Seam_edge_uhm>                Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm>              Seam_vertex_pmap;


typedef CGAL::Seam_mesh<Textured_polyhedron::Base, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;
typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor           halfedge_descriptor;

class UVItem : public QGraphicsItem
{
public :
  UVItem(Textured_polyhedron* t_m,
         Components* components,
         std::vector<std::vector<float> >uv_borders,
         QRectF brect)
    :QGraphicsItem(),
      texMesh(t_m),
      bounding_rect(brect),
      components(components),
      m_borders(uv_borders),
      m_current_component(0)
  {
  }

  ~UVItem()
  {
    delete components;
  }

  std::vector<float> borders() { return m_borders[m_current_component]; }

  QRectF boundingRect() const
  {
    return bounding_rect;
  }
  QString item_name()const{ return texMesh_name; }
  void set_item_name(QString s){ texMesh_name = s;}

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
  {
    QPen pen;
    pen.setColor(Qt::black);
    pen.setWidth(0);
    painter->setPen(pen);
    for( Component::iterator
         fi = components->at(m_current_component).begin();
         fi != components->at(m_current_component).end();
         ++fi)
    {
      Textured_polyhedron::Facet_handle f(*fi);
      QPointF pt_A(f->halfedge()->u(), f->halfedge()->v());
      QPointF pt_B(f->halfedge()->next()->u(), f->halfedge()->next()->v());
      QPointF pt_C(f->halfedge()->next()->next()->u(), f->halfedge()->next()->next()->v());

      painter->drawLine(pt_A, pt_B);
      painter->drawLine(pt_B, pt_C);
      painter->drawLine(pt_C, pt_A);
    }
  }
  int number_of_components()const{return components->size();}
  int current_component()const{return m_current_component;}
  void set_current_component(int n){m_current_component = n;}
private:
  Textured_polyhedron* texMesh;
  QString texMesh_name;
  QRectF bounding_rect;
  Components* components;
  std::vector<std::vector<float> > m_borders;
  int m_current_component;
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
            Messages_interface* msg)
  {
    mw = mainWindow;
    scene = scene_interface;
    messages = msg;
    Scene* true_scene = static_cast<Scene*>(scene);
    connect(true_scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
            this, SLOT(destroyPolyline(CGAL::Three::Scene_item*)));
    QAction* actionMVC = new QAction("Mean Value Coordinates", mw);
    QAction* actionDCP = new QAction ("Discrete Conformal Map", mw);
    QAction* actionLSC = new QAction("Least Square Conformal Map", mw);
    QAction* actionDAP = new QAction("Discrete Authalic", mw);
    QAction* actionARAP = new QAction("As Rigid As Possible", mw);
    actionMVC->setObjectName("actionMVC");
    actionDCP->setObjectName("actionDCP");
    actionLSC->setObjectName("actionLSC");
    actionDAP->setObjectName("actionDAP");
    actionARAP->setObjectName("actionARAP");

    _actions << actionMVC
             << actionDCP
             << actionLSC
             << actionDAP
             << actionARAP;
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
    ui_widget.component_numberLabel->setText("Component : 1");
    connect(ui_widget.prevButton, &QPushButton::clicked, this, &Polyhedron_demo_parameterization_plugin::on_prevButton_pressed);
    connect(ui_widget.nextButton, &QPushButton::clicked, this, &Polyhedron_demo_parameterization_plugin::on_nextButton_pressed);
    addDockWidget(dock_widget);
    dock_widget->setVisible(false);
    current_uv_item = NULL;

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
  void on_actionDAP_triggered();
  void on_actionARAP_triggered();
  void on_prevButton_pressed();
  void on_nextButton_pressed();
  void replacePolyline()
  {
    if(current_uv_item)
      qobject_cast<Scene_textured_polyhedron_item*>(projections.key(current_uv_item))->add_border_edges(std::vector<float>(0));
    int id = scene->mainSelectionIndex();
    Q_FOREACH(UVItem* pl, projections)
    {
      if(pl==NULL
         || pl != projections[scene->item(id)])
        continue;
      current_uv_item = pl;
      break;
    }
    if(!current_uv_item)
    {
      dock_widget->setWindowTitle(tr("UVMapping"));
      ui_widget.component_numberLabel->setText(QString("Component :"));
    }
    else
    {
      if(!graphics_scene->items().empty())
        graphics_scene->removeItem(graphics_scene->items().first());
      graphics_scene->addItem(current_uv_item);
      ui_widget.graphicsView->fitInView(current_uv_item->boundingRect(), Qt::KeepAspectRatio);
      ui_widget.component_numberLabel->setText(QString("Component : %1/%2").arg(current_uv_item->current_component()+1).arg(current_uv_item->number_of_components()));
      dock_widget->setWindowTitle(tr("UVMapping for %1").arg(current_uv_item->item_name()));
      qobject_cast<Scene_textured_polyhedron_item*>(projections.key(current_uv_item))->add_border_edges(current_uv_item->borders());

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
    if(projections.empty() || projections.first() == NULL)
    {
      current_uv_item = NULL;
      dock_widget->setWindowTitle(tr("UVMapping"));
      ui_widget.component_numberLabel->setText(QString("Component :"));
    }
    else
      current_uv_item = projections.first();
  }

protected:
  enum Parameterization_method { PARAM_MVC, PARAM_DCP, PARAM_LSC, PARAM_DAP, PARAM_ARAP};
  void parameterize(Parameterization_method method);
private:
  Messages_interface *messages;
  QList<QAction*> _actions;
  QDockWidget* dock_widget;
  Ui::Parameterization ui_widget;
  QGraphicsScene *graphics_scene;
  Navigation* navigation;
  QMap<Scene_item*, UVItem*> projections;
  UVItem* current_uv_item;
}; // end Polyhedron_demo_parameterization_plugin

void Polyhedron_demo_parameterization_plugin::on_prevButton_pressed()
{
  int id = scene->mainSelectionIndex();
  Q_FOREACH(UVItem* pl, projections)
  {
    if(pl==NULL
       || pl != projections[scene->item(id)])
      continue;

    current_uv_item = pl;
    break;
  }
  if(current_uv_item == NULL)
    return;
  current_uv_item->set_current_component((std::max)(0,current_uv_item->current_component()-1));
  replacePolyline();
}

void Polyhedron_demo_parameterization_plugin::on_nextButton_pressed()
{
  int id = scene->mainSelectionIndex();
  Q_FOREACH(UVItem* pl, projections)
  {
    if(pl==NULL
       || pl != projections[scene->item(id)])
      continue;

    current_uv_item = pl;
    break;
  }
  if(current_uv_item == NULL)
    return;
  current_uv_item->set_current_component((std::min)(current_uv_item->number_of_components()-1,current_uv_item->current_component()+1));
  ui_widget.component_numberLabel->setText(QString("Component : %1/%2").arg(current_uv_item->current_component()+1).arg(current_uv_item->number_of_components()));
  replacePolyline();
}
void Polyhedron_demo_parameterization_plugin::parameterize(const Parameterization_method method)
{

  // get active polyhedron
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if(!poly_item)
  {
    messages->error("Selected item is not of the right type.");
    return;
  }

  Polyhedron* pMesh = poly_item->polyhedron();
  if(!pMesh)
  {
    messages->error("Selected item has no valid polyhedron.");
    return;
  }
  pMesh->normalize_border();
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
  if(!is_seamed && pMesh->is_closed())
  {
    messages->error("The selected mesh has no border and is not seamed.");
    return;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);

  ///////////////////////////////////
  ////////// PARAMETERIZE ///////////
  ///////////////////////////////////

  QTime time;
  time.start();
  // add textured polyhedon to the scene
  Textured_polyhedron *tpMesh = new Textured_polyhedron();
  Textured_polyhedron_builder<Polyhedron,Textured_polyhedron,Kernel> builder;
  builder.run(*pMesh,*tpMesh);
  tpMesh->compute_normals();
  tpMesh->normalize_border();
  Textured_polyhedron::Base tMesh = static_cast<Textured_polyhedron::Base&>(*tpMesh);

  CGAL::set_halfedgeds_items_id(tMesh);
  std::vector<bool> mark(tpMesh->size_of_halfedges()/2,false);
  std::vector<T_edge_descriptor> seam_edges;
  if(is_seamed)
  {
    //create a textured_polyhedron edges selection from the ids of the corresponding vertices
    BOOST_FOREACH(P_edge_descriptor ed, sel_item->selected_edges)
    {
      Polyhedron::Vertex_handle a(source(ed, *pMesh)), b(target(ed, *pMesh));

      for(Textured_polyhedron::Edge_iterator it =
          tMesh.edges_begin(); it != tMesh.edges_end();
          ++it)
      {
        Textured_polyhedron::Vertex_handle ta(source(it, tMesh)), tb(target(it, tMesh));
        if((ta->id() == a->id() && tb->id() == b->id())
           ||
           (ta->id() == b->id() && tb->id() == a->id()))
        {
          sel_item->selected_vertices.insert(a);
          sel_item->selected_vertices.insert(b);
          T_edge_descriptor ted(it);
          seam_edges.push_back(ted);
          break;
        }
      }

    }
    qDebug()<<sel_item->selected_edges.size()<<", "<<seam_edges.size();
    //fill seam mesh pmaps
    BOOST_FOREACH(T_edge_descriptor ed, seam_edges)
    {
      T_halfedge_descriptor hd = halfedge(ed, tMesh);
      T_vertex_descriptor svd(source(hd, tMesh)), tvd(target(hd, tMesh));
      if(!is_border(ed, tMesh))
      {
        put(seam_edge_pm, ed, true);
        put(seam_vertex_pm, svd, true);
        put(seam_vertex_pm, tvd, true);
        mark[hd->id()/2] = true;
      }
    }
  }

  Seam_mesh *sMesh = new Seam_mesh(tMesh, seam_edge_pm, seam_vertex_pm);
  sMesh->set_seam_edges_number(seam_edges.size());

  UV_uhm uv_uhm;
  UV_pmap uv_pm(uv_uhm);

  QString new_item_name;
  //determine the different connected_components
  boost::property_map<Textured_polyhedron::Base, boost::face_external_index_t>::type fim
      = get(boost::face_external_index, tMesh);
  boost::vector_property_map<int,
      boost::property_map<Textured_polyhedron::Base, boost::face_external_index_t>::type>
      fccmap(fim);

  Is_selected_property_map edge_pmap(mark);

  int number_of_components =
      CGAL::Polygon_mesh_processing::connected_components(
        tMesh,
        fccmap,
        CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(
          edge_pmap));
  Components t_components(number_of_components);
  Textured_polyhedron::Base::Facet_iterator fit;
  for(fit = tMesh.facets_begin();
      fit != tMesh.facets_end();
      ++fit)
  {
    t_components.at(fccmap[fit]).insert(fit);
  }

  //once per component
  std::vector<std::vector<float> >uv_borders;
  uv_borders.resize(number_of_components);
  for(int current_component=0; current_component<number_of_components; ++current_component)
  {

    T_halfedge_descriptor phd;
    std::vector<T_halfedge_descriptor> border;
    PMP::border_halfedges(t_components.at(current_component),
                          tMesh,
                          std::back_inserter(border));

    BOOST_FOREACH(T_halfedge_descriptor hd, border)
    {
        uv_borders[current_component].push_back(source(hd, tMesh)->point().x());
        uv_borders[current_component].push_back(source(hd, tMesh)->point().y());
        uv_borders[current_component].push_back(source(hd, tMesh)->point().z());

        uv_borders[current_component].push_back(target(hd, tMesh)->point().x());
        uv_borders[current_component].push_back(target(hd, tMesh)->point().y());
        uv_borders[current_component].push_back(target(hd, tMesh)->point().z());
    }

    if(!border.empty())
      phd = opposite(border.front(), tMesh);
    //in case there are no border, take the first halfedge of the seam.
    if(phd == boost::graph_traits<Textured_polyhedron::Base>::null_halfedge())
    {
      phd = halfedge(*seam_edges.begin(), tMesh);
    }

    halfedge_descriptor bhd(phd);
    bhd = opposite(bhd,*sMesh); // a halfedge on the virtual border
    bool success = false;
    switch(method)
    {
    case PARAM_MVC:
    {
      std::cout << "Parameterize (MVC)...";
      new_item_name = tr("%1 (parameterized (MVC))").arg(poly_item->name());
      typedef SMP::Mean_value_coordinates_parameterizer_3<Seam_mesh> Parameterizer;
      SMP::Error_code err = SMP::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == SMP::OK;
      break;
    }
    case PARAM_DCP:
    {
      new_item_name = tr("%1 (parameterized (DCP))").arg(poly_item->name());
      std::cout << "Parameterize (DCP)...";
      typedef SMP::Discrete_conformal_map_parameterizer_3<Seam_mesh> Parameterizer;
      SMP::Error_code err = SMP::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = err == SMP::OK;
      break;
    }
    case PARAM_LSC:
    {
      new_item_name = tr("%1 (parameterized (LSC))").arg(poly_item->name());
      std::cout << "Parameterize (LSC)...";
      typedef SMP::Two_vertices_parameterizer_3<Seam_mesh> Border_parameterizer;
      typedef SMP::LSCM_parameterizer_3<Seam_mesh, Border_parameterizer> Parameterizer_with_border;
      typedef SMP::LSCM_parameterizer_3<Seam_mesh> Parameterizer;
      if(!border.empty())
      {
        double max_dist = 0;
        int indice_max = 0;
        Kernel::Point_3 a = target(border[0], tMesh)->point();
        for(std::size_t i=1; i<border.size(); ++i)
        {
          Kernel::Point_3 b = target(border[i], tMesh)->point();
          double dist = std::sqrt((b.x()-a.x())*(b.x()-a.x())+
                                  (b.y()-a.y())*(b.y()-a.y())+
                                  (b.z()-a.z())*(b.z()-a.z()));

          if(dist>max_dist)
          {
            max_dist = dist;
            indice_max = i;
          }
        }
        T_halfedge_descriptor phd2 = opposite(border[indice_max], tMesh);

        boost::graph_traits<Seam_mesh>::vertex_descriptor vp1 = target(halfedge_descriptor(phd),*sMesh);
        boost::graph_traits<Seam_mesh>::vertex_descriptor vp2 = target(halfedge_descriptor(phd2),*sMesh);
        SMP::Error_code err = SMP::parameterize(*sMesh,
                                                Parameterizer_with_border(Border_parameterizer(vp1, vp2)),
                                                bhd, uv_pm);
        success = err == SMP::OK;
      }
      else
      {
        SMP::Error_code err = SMP::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
        success = err == SMP::OK;
      }
      break;
    }
    case PARAM_DAP:
    {
      new_item_name = tr("%1 (parameterized (DAP))").arg(poly_item->name());
      std::cout << "Parameterize (DAP)...";
      typedef SMP::Discrete_authalic_parameterizer_3<Seam_mesh> Parameterizer;
      SMP::Error_code err = SMP::parameterize(*sMesh, Parameterizer(), bhd, uv_pm);
      success = (err == SMP::OK);
      break;
    }
    case PARAM_ARAP:
    {
      new_item_name = tr("%1 (parameterized (ARAP))").arg(poly_item->name());
      std::cout << "Parameterize (ARAP)...";

      QDialog dialog(mw);
      Ui::ARAP_dialog ui;
      ui.setupUi(&dialog);
      connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
      connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

      // Get values
      QApplication::restoreOverrideCursor();

      int i = dialog.exec();
      if (i == QDialog::Rejected)
        return;

      FT lambda = ui.lambdaSpinBox->value();
      QApplication::setOverrideCursor(Qt::WaitCursor);

      typedef SMP::ARAP_parameterizer_3<Seam_mesh> Parameterizer;
      SMP::Error_code err = SMP::parameterize(*sMesh, Parameterizer(lambda), bhd, uv_pm);
      success = (err == SMP::OK);
      break;
    }
    }//end switch

    QApplication::restoreOverrideCursor();

    if(success)
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    else
    {
      std::cout << "failure" << std::endl;
      return;
    }
  } //end for each component


  Textured_polyhedron::Base::Halfedge_iterator it1;
  Textured_polyhedron::Halfedge_iterator it2;
  QPointF min(FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX);
  for(it1 = tMesh.halfedges_begin(),
      it2 = tpMesh->halfedges_begin();
      it1 != tMesh.halfedges_end()&&
      it2 != tpMesh->halfedges_end();
      ++it1, ++it2)
  {
    Seam_mesh::halfedge_descriptor hd(it1);
    FT u = uv_pm[target(hd,*sMesh)].x();
    FT v = uv_pm[target(hd,*sMesh)].y();
    it2->u() = u;
    it2->v() = v;
    if(u<min.x())
      min.setX(u);
    if(u>max.x())
      max.setX(u);
    if(v<min.y())
      min.setY(v);
    if(v>max.y())
      max.setY(v);
  }
  Components* components = new Components(0);
  components->resize(number_of_components);
  Textured_polyhedron::Base::Facet_iterator bfit;
  Textured_polyhedron::Facet_iterator tfit;
  for(bfit = tMesh.facets_begin(),
      tfit = tpMesh->facets_begin();
      bfit != tMesh.facets_end()&&
      tfit != tpMesh->facets_end();
      ++bfit, ++tfit)
  {
    components->at(fccmap[bfit]).insert(tfit);
  }
  UVItem *projection
      = new UVItem(tpMesh, components, uv_borders, QRectF(min, max));
  projection->set_item_name(new_item_name);
  Scene_textured_polyhedron_item* new_item = new Scene_textured_polyhedron_item(tpMesh);
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
  if(current_uv_item)
    qobject_cast<Scene_textured_polyhedron_item*>(projections.key(current_uv_item))->add_border_edges(std::vector<float>(0));
  current_uv_item = projection;
  qobject_cast<Scene_textured_polyhedron_item*>(projections.key(current_uv_item))->add_border_edges(current_uv_item->borders());
  if(dock_widget->isHidden())
    dock_widget->setVisible(true);
  dock_widget->setWindowTitle(tr("UVMapping for %1").arg(new_item->name()));
  ui_widget.component_numberLabel->setText(QString("Component : %1/%2").arg(current_uv_item->current_component()+1).arg(current_uv_item->number_of_components()));
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

void Polyhedron_demo_parameterization_plugin::on_actionDAP_triggered()
{
  std::cerr << "DAP...";
  parameterize(PARAM_DAP);
}

void Polyhedron_demo_parameterization_plugin::on_actionARAP_triggered()
{
  std::cerr << "ARAP...";
  parameterize(PARAM_ARAP);
}

#include "Parameterization_plugin.moc"
