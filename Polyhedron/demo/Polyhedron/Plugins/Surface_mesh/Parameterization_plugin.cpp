#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>

#include "Scene_surface_mesh_item.h"
#include "Scene_textured_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "SMesh_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include "Scene.h"
#include <QElapsedTimer>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QPen>
#include <QDockWidget>
#include <Messages_interface.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/property_map.h>


#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/internal/orbifold_cone_helper.h>
#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>


#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>


#include <CGAL/boost/graph/properties.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "ui_Parameterization_widget.h"
#include "ui_OTE_dialog.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_textured_surface_mesh_item Scene_textured_facegraph_item;
typedef SMesh Textured_face_graph;
typedef SMesh Base_face_graph;
typedef SMesh Face_graph;
typedef EPICK Traits;

namespace SMP = CGAL::Surface_mesh_parameterization;

typedef boost::unordered_set<boost::graph_traits<Base_face_graph>::face_descriptor> Component;
typedef std::vector<Component> Components;

struct Is_selected_property_map{
  typedef boost::graph_traits<Base_face_graph>::edge_descriptor edge_descriptor;
  typedef boost::property_map<Base_face_graph, boost::halfedge_index_t>::type HIndexMap;
  std::vector<bool>* is_selected_ptr;
  Base_face_graph* graph;
  HIndexMap idmap;
  Is_selected_property_map()
    : is_selected_ptr(NULL), graph(NULL) {}
  Is_selected_property_map(std::vector<bool>& is_selected,
                           Base_face_graph* graph)
    : is_selected_ptr( &is_selected), graph(graph)
  {
    idmap = get(boost::halfedge_index, *graph);
  }

  std::size_t id(edge_descriptor ed) { return get(idmap, halfedge(ed, *graph))/2; }

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
  { }

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
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
      QPoint pos = event->pos();
#else
      QPointF pos = event->position();
#endif
      QPointF old_pos = v->mapToScene(pos.x(), pos.y());
      if(event->angleDelta().y() <0)
        v->scale(1.2, 1.2);
      else
        v->scale(0.8, 0.8);
      QPointF new_pos = v->mapToScene(pos.x(), pos.y());
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

namespace SMP = CGAL::Surface_mesh_parameterization;

typedef Traits::FT                                                  FT;
typedef boost::graph_traits<Face_graph>::vertex_descriptor          P_vertex_descriptor;
typedef Traits::Point_2                                             Point_2;
typedef boost::graph_traits<Face_graph>::edge_descriptor            P_edge_descriptor;
typedef boost::graph_traits<Face_graph>::halfedge_descriptor        P_halfedge_descriptor;

// Textured polyhedron
typedef boost::graph_traits<Base_face_graph>::
                                         edge_descriptor            T_edge_descriptor;
typedef boost::graph_traits<Base_face_graph>::
                                         halfedge_descriptor        T_halfedge_descriptor;
typedef boost::graph_traits<Base_face_graph>::
                                         vertex_descriptor          T_vertex_descriptor;

// Seam
typedef CGAL::Unique_hash_map<T_halfedge_descriptor,Point_2>        UV_uhm;
typedef CGAL::Unique_hash_map<T_edge_descriptor,bool>               Seam_edge_uhm;
typedef CGAL::Unique_hash_map<T_vertex_descriptor,bool>             Seam_vertex_uhm;

typedef boost::associative_property_map<UV_uhm>                     UV_pmap;
typedef boost::associative_property_map<Seam_edge_uhm>              Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm>            Seam_vertex_pmap;

typedef CGAL::Seam_mesh<Base_face_graph,
                        Seam_edge_pmap, Seam_vertex_pmap>           Seam_mesh;

typedef boost::graph_traits<Seam_mesh>::vertex_descriptor           s_vertex_descriptor;
typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor         s_halfedge_descriptor;
typedef boost::graph_traits<Seam_mesh>::face_descriptor             s_face_descriptor;

typedef boost::graph_traits<Seam_mesh>::edges_size_type             s_edges_size_type;

typedef boost::unordered_set<boost::graph_traits<Base_face_graph>::
face_descriptor>                                                    Component;
typedef std::vector<Component>                                      Components;

typedef boost::unordered_set<s_face_descriptor>                       SComponent;
typedef std::vector<SComponent>                                     SComponents;

class UVItem : public QGraphicsItem
{
public :
  UVItem(Components* components,
         Base_face_graph* graph,
         std::vector<std::vector<float> >uv_borders,
         QRectF brect)
    :
      QGraphicsItem(),
      bounding_rect(brect),
      components(components),
      graph(graph),
      m_borders(uv_borders),
      m_concatenated_borders(),
      m_current_component(0)
  {
    std::size_t total_border_size = 0;
    for(std::size_t i=0; i<m_borders.size(); ++i)
      total_border_size += m_borders[i].size();

    m_concatenated_borders.resize(total_border_size);
    for(std::size_t i=0; i<m_borders.size(); ++i)
    {
      const std::vector<float>& ith_border = m_borders[i];
      m_concatenated_borders.insert(m_concatenated_borders.end(),
                                    ith_border.begin(), ith_border.end());
    }
  }

  ~UVItem()
  {
    delete components;
  }

  const std::vector<std::vector<float> >& borders() const { return m_borders; }
  const std::vector<float>& concatenated_borders() const { return m_concatenated_borders; }

  QRectF boundingRect() const
  {
    return bounding_rect;
  }

  QString item_name()const{ return texMesh_name; }
  void set_item_name(QString s){ texMesh_name = s;}

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
  {
    QPen pen;
    QBrush brush;
    brush.setColor(QColor(100, 100, 255));
    brush.setStyle(Qt::SolidPattern);
    pen.setColor(Qt::black);
    pen.setWidth(0);
    painter->setPen(pen);
    painter->setBrush(brush);
    SMesh::Property_map<halfedge_descriptor,float> u,v;

    u = graph->add_property_map<halfedge_descriptor,float>("h:u", 0.0f).first;
    v = graph->add_property_map<halfedge_descriptor,float>("h:v", 0.0f).first;

    for( Component::iterator
         fi = components->at(m_current_component).begin();
         fi != components->at(m_current_component).end();
         ++fi)
    {
      boost::graph_traits<Base_face_graph>::face_descriptor f(*fi);

      QPointF points[3];
      boost::graph_traits<Base_face_graph>::halfedge_descriptor h = halfedge(f, *graph);;
      points[0] = QPointF(get(u, h), get(v, h));
      h = next(halfedge(f, *graph), *graph);
      points[1] = QPointF(get(u, h), get(v, h));
      h = next(next(halfedge(f, *graph), *graph), *graph);
      points[2] = QPointF(get(u, h), get(v, h));
      painter->drawPolygon(points,3);
    }
  }

  int number_of_components()const{return static_cast<int>(components->size());}
  int current_component()const{return m_current_component;}
  void set_current_component(int n){m_current_component = n;}

private:
  QString texMesh_name;
  QRectF bounding_rect;
  Components* components;
  Base_face_graph* graph;
  std::vector<std::vector<float> > m_borders;
  std::vector<float> m_concatenated_borders;
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
  QList<QAction*> actions() const
  {
    return _actions;
  }

  ~Polyhedron_demo_parameterization_plugin()
  {
    delete navigation;
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
    QAction* actionIAP = new QAction("Iterative Authalic", mw);
    QAction* actionARAP = new QAction("As Rigid As Possible", mw);
    QAction* actionOTE = new QAction("Orbifold Tutte Embedding", mw);
    QAction* actionBTP = new QAction("Tutte Barycentric", mw);
    actionMVC->setObjectName("actionMVC");
    actionDCP->setObjectName("actionDCP");
    actionLSC->setObjectName("actionLSC");
    actionDAP->setObjectName("actionDAP");
    actionIAP->setObjectName("actionIAP");
    actionARAP->setObjectName("actionARAP");
    actionOTE->setObjectName("actionOTE");
    actionBTP->setObjectName("actionBTP");

    _actions << actionARAP
             << actionBTP
             << actionDAP
             << actionIAP
             << actionDCP
             << actionLSC
             << actionMVC
             << actionOTE;
    autoConnectActions();
    Q_FOREACH(QAction *action, _actions)
      action->setProperty("subMenuName",
                          "Triangulated Surface Mesh Parameterization"
                          );
    dock_widget = new QDockWidget(
          "UVMapping "
          , mw);
    ui_widget.setupUi(dock_widget);
    dock_widget->setWindowTitle(tr(
                                  "UVMapping "
                                  ));
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

  bool applicable(QAction*) const
  {
    if (scene->selectionIndices().size() == 1)
    {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }

    Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices())
    {
      //if one facegraph is found in the selection, it's fine
      if (qobject_cast<Scene_facegraph_item*>(scene->item(id)))
        return true;
    }
    return false;
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
  void on_actionIAP_triggered();
  void on_actionARAP_triggered();
  void on_actionOTE_triggered();
  void on_actionBTP_triggered();
  void on_prevButton_pressed();
  void on_nextButton_pressed();

  void replacePolyline()
  {
    if(current_uv_item){
      Scene_textured_facegraph_item* t_item =
          qobject_cast<Scene_textured_facegraph_item*>(projections.key(current_uv_item));
     t_item->add_border_edges(std::vector<float>(0));
    }

    int id = scene->mainSelectionIndex();

    Q_FOREACH(UVItem* pl, projections)
    {
      if(pl==NULL || pl != projections[scene->item(id)])
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
      ui_widget.graphicsView->fitInView(current_uv_item->boundingRect(),
                                        Qt::KeepAspectRatio);
      ui_widget.component_numberLabel->setText(
            QString("Component : %1/%2").arg(current_uv_item->current_component()+1)
            .arg(current_uv_item->number_of_components()));
      dock_widget->setWindowTitle(tr("UVMapping for %1")
                                  .arg(current_uv_item->item_name()));
      Scene_textured_facegraph_item* t_item =
          qobject_cast<Scene_textured_facegraph_item*>(projections.key(current_uv_item));
     t_item->add_border_edges(
            current_uv_item->concatenated_borders());
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
  enum Parameterization_method { PARAM_MVC, PARAM_DCP, PARAM_LSC, PARAM_DAP,
                                 PARAM_IAP, PARAM_ARAP, PARAM_OTE, PARAM_BTP};
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
  Scene_facegraph_item* poly_item = NULL;
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices())
  {
    poly_item = qobject_cast<Scene_facegraph_item*>(scene->item(id));
    if(!poly_item)
    {
      continue;
    }
    else
    {
      index = id;
      break;
    }
  }

  if(!poly_item)
  {
    CGAL::Three::Three::error("Selected item is not of the right type.");
    return;
  }

  Face_graph* pMesh = poly_item->face_graph();
  if(!pMesh)
  {
    CGAL::Three::Three::error("Selected item has no valid polyhedron.");
    return;
  }
  Scene_polyhedron_selection_item* sel_item = NULL;
  bool is_seamed = false;
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices())
  {
    sel_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(id));
    if(!sel_item)
      continue;
    if(sel_item->selected_edges.empty())
      continue;
    if(method == PARAM_OTE && sel_item->selected_vertices.empty())
       continue;
    is_seamed = true;
  }

  if(method == PARAM_OTE &&
     (sel_item == NULL || sel_item->selected_vertices.empty())) {
    std::cerr << "\nError: no cones/seam selected; Aborting parameterization." << std::endl;
    return;
  }

  // Two property maps to store the seam edges and vertices
  Seam_edge_uhm seam_edge_uhm(false);
  Seam_edge_pmap seam_edge_pm(seam_edge_uhm);

  Seam_vertex_uhm seam_vertex_uhm(false);
  Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);

  if(!is_seamed && is_closed(*pMesh))
  {
    CGAL::Three::Three::error("The selected mesh has no (real or virtual) border.");
    return;
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);

  ///////////////////////////////////
  ////////// PARAMETERIZE ///////////
  ///////////////////////////////////

  QElapsedTimer time;
  time.start();
  // add textured polyhedon to the scene

  // \todo for surface_mesh
  Base_face_graph tMesh = *pMesh;
  std::vector<bool> mark(num_halfedges(tMesh)/2,false);
  std::vector<T_edge_descriptor> seam_edges;
  typedef boost::property_map<Base_face_graph, boost::vertex_index_t>::type VIDMap;
  VIDMap vidmap = get(boost::vertex_index, tMesh);
  if(is_seamed)
  {
    //create a textured_polyhedron edges selection from the ids of the corresponding vertices
    typedef boost::property_map<Base_face_graph, boost::halfedge_index_t>::type HIDMap;
    HIDMap hidmap = get(boost::halfedge_index, tMesh);
    for(P_edge_descriptor ed : sel_item->selected_edges)
    {
      boost::graph_traits<Face_graph>::vertex_descriptor a(source(ed, *pMesh)), b(target(ed, *pMesh));

      for(boost::graph_traits<Textured_face_graph>::edge_iterator it =
          edges(tMesh).begin(); it != edges(tMesh).end();
          ++it)
      {
        boost::graph_traits<Textured_face_graph>::vertex_descriptor ta(source(*it, tMesh)), tb(target(*it, tMesh));

        if((get(vidmap, ta) == get(vidmap, a) && get(vidmap,tb) == get(vidmap,b))
           ||
           (get(vidmap,ta) == get(vidmap,b) && get(vidmap,tb) == get(vidmap,a)))
        {
          T_edge_descriptor ted(*it);
          seam_edges.push_back(ted);
          break;
        }
      }

    }
    qDebug() << sel_item->selected_edges.size() << ", " << seam_edges.size();
    //fill seam mesh pmaps
    for(T_edge_descriptor ed : seam_edges)
    {
      T_halfedge_descriptor hd = halfedge(ed, tMesh);
      T_vertex_descriptor svd(source(hd, tMesh)), tvd(target(hd, tMesh));
      if(!is_border(ed, tMesh))
      {
        put(seam_edge_pm, ed, true);
        put(seam_vertex_pm, svd, true);
        put(seam_vertex_pm, tvd, true);
        mark[get(hidmap, hd)/2] = true;
      }
    }
  }

  // map the cones from the selection plugin to the textured polyhedron
  boost::unordered_set<T_vertex_descriptor> unordered_cones;
  if(method == PARAM_OTE) {
    for(P_vertex_descriptor vd : sel_item->selected_vertices) {
      boost::graph_traits<Face_graph>::vertex_descriptor pvd(vd);
      boost::graph_traits<Textured_face_graph>::vertex_iterator it = vertices(tMesh).begin(),
          end = vertices(tMesh).end();
      for(; it!=end; ++it) {
        boost::graph_traits<Textured_face_graph>::vertex_descriptor tvd(*it);
        if(get(vidmap, *it) == get(vidmap, pvd)) {
          unordered_cones.insert(tvd);
        }
      }
    }
  }
  Seam_mesh sMesh(tMesh, seam_edge_pm, seam_vertex_pm);
  sMesh.set_seam_edges_number(static_cast<s_edges_size_type>(seam_edges.size()));

  // The parameterized values
  UV_uhm uv_uhm;
  UV_pmap uv_pm(uv_uhm);

  QString new_item_name;
  //determine the different connected_components
  boost::container::flat_map<boost::graph_traits<Base_face_graph>::face_descriptor, int> face_component_map;
  boost::associative_property_map< boost::container::flat_map<s_face_descriptor, int> >
      fccmap(face_component_map);

  Is_selected_property_map edge_pmap(mark, &tMesh);

  int number_of_components =
      CGAL::Polygon_mesh_processing::connected_components(
        tMesh,
        fccmap,
        CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(
          edge_pmap));

  // Next is the gathering of the border halfedges of the connected component.
  // It is wrong to pass the underlying mesh tMesh: a sphere split in half does
  // not have any border if border_halfedges() is run with tMesh.
  //
  // The proper way would be to completely redesign the plugin to use Seam meshes
  // everywhere. But that's not worth it. Instead, we abuse the fact that faces
  // are the same in tMesh and sMesh.

  // the SEAM MESH faces of each connected component
  SComponents s_components(number_of_components);

  for(boost::graph_traits<Base_face_graph>::face_iterator fit = faces(tMesh).begin();
      fit != faces(tMesh).end(); ++fit) {
    s_components.at(fccmap[*fit]).insert(s_face_descriptor(*fit));
  }

  // once per component
  std::vector<std::vector<float> >uv_borders;
  uv_borders.resize(number_of_components);

  // to track whether the components are successfully parameterized
  SMP::Error_code status = SMP::OK;

  for(int current_component=0; current_component<number_of_components; ++current_component)
  {
    std::vector<s_halfedge_descriptor> border;
    PMP::border_halfedges(s_components.at(current_component),
                          sMesh, std::back_inserter(border));

    std::cout << sMesh.number_of_seam_edges() << " seams" << std::endl;
    std::cout << (s_components.at(current_component)).size() << " faces" << std::endl;
    std::cout << border.size() << " border halfedges" << std::endl;

    // find longest border in the connected component
    s_halfedge_descriptor bhd; // a halfedge on the (possibly virtual) border
    boost::unordered_set<s_halfedge_descriptor> visited;
    FT result_len = 0;
    for(s_halfedge_descriptor hd : border)
    {
      assert(is_border(hd, sMesh));

      if(visited.find(hd) == visited.end())
      {
        FT len = 0;
        for(s_halfedge_descriptor haf : halfedges_around_face(hd, sMesh))
        {
          len += PMP::edge_length(haf, sMesh);
          visited.insert(haf);
        }

        if(result_len < len)
        {
          result_len = len;
          bhd = hd;
        }
      }
    }
    CGAL_postcondition(bhd != s_halfedge_descriptor());
    CGAL_postcondition(is_border(bhd, sMesh));
    typedef boost::property_map<Base_face_graph, boost::vertex_point_t>::type VPMap;
    VPMap vpmap =get(boost::vertex_point, tMesh);

    // collect the border edges for that connected component
    for(s_halfedge_descriptor haf : halfedges_around_face(bhd, sMesh))
    {
        uv_borders[current_component].push_back(get(vpmap, source(haf, tMesh)).x());
        uv_borders[current_component].push_back(get(vpmap, source(haf, tMesh)).y());
        uv_borders[current_component].push_back(get(vpmap, source(haf, tMesh)).z());

        uv_borders[current_component].push_back(get(vpmap, target(haf, tMesh)).x());
        uv_borders[current_component].push_back(get(vpmap, target(haf, tMesh)).y());
        uv_borders[current_component].push_back(get(vpmap, target(haf, tMesh)).z());
    }

    switch(method)
    {
    case PARAM_MVC:
    {
      std::cout << "Parameterize (MVC)..." << std::endl;
      new_item_name = tr("%1 (parameterized (MVC))").arg(poly_item->name());
      typedef SMP::Mean_value_coordinates_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(), bhd, uv_pm);
      break;
    }
    case PARAM_DCP:
    {
      new_item_name = tr("%1 (parameterized (DCP))").arg(poly_item->name());
      std::cout << "Parameterize (DCP)..." << std::endl;
      typedef SMP::Discrete_conformal_map_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(), bhd, uv_pm);
      break;
    }
    case PARAM_LSC:
    {
      new_item_name = tr("%1 (parameterized (LSC))").arg(poly_item->name());
      std::cout << "Parameterize (LSC)..." << std::endl;
      typedef SMP::LSCM_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(), bhd, uv_pm);
      break;
    }
    case PARAM_DAP:
    {
      new_item_name = tr("%1 (parameterized (DAP))").arg(poly_item->name());
      std::cout << "Parameterize (DAP)..." << std::endl;
      typedef SMP::Discrete_authalic_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(), bhd, uv_pm);
      break;
    }
    case PARAM_IAP:
    {
      new_item_name = tr("%1 (parameterized (IAP))").arg(poly_item->name());
      std::cout << "Parameterize (IAP)..." << std::endl;
      typedef SMP::Iterative_authalic_parameterizer_3<Seam_mesh> Parameterizer;
      Parameterizer parameterizer;
      status = parameterizer.parameterize(sMesh, bhd, uv_pm, 15 /*iterations*/);
      break;
    }
    case PARAM_ARAP:
    {
      new_item_name = tr("%1 (parameterized (ARAP))").arg(poly_item->name());
      std::cout << "Parameterize (ARAP)..." << std::endl;
      FT lambda = 10000; // a big value to ensure the parameterization is ARAP (and not ASAP)
      typedef SMP::ARAP_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(lambda), bhd, uv_pm);
      break;
    }
    case PARAM_OTE:
    {
      new_item_name = tr("%1 (parameterized (OTE))").arg(poly_item->name());
      std::cout << "Parameterize (OTE)..." << std::endl;

      // OTE cannot handle multiple connected components right now
      // @todo (need to remove the assertions such as cones.size() == 4
      //        and check where and when cones are used (passed by ID, for ex.?))
      if(number_of_components != 1) {
        std::cerr << "Orbifold Tutte Embedding can only handle one connected component" << std::endl;
        status = SMP::ERROR_NO_TOPOLOGICAL_BALL;
        break;
      }

      typedef SMP::Orbifold_Tutte_parameterizer_3<Seam_mesh> Parameterizer;

      // Get orbifold type
      QDialog dialog(mw);
      Ui::OTE_dialog ui;
      ui.setupUi(&dialog);
      connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
      connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

      QApplication::restoreOverrideCursor();

      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Aborting parameterization" << std::endl;
        QApplication::restoreOverrideCursor();
        return;
      }

      SMP::Orbifold_type orb = static_cast<SMP::Orbifold_type>(ui.OrbComboBox->currentIndex());
      std::cout << "Selected orbifold type: " << ui.OrbComboBox->currentText().toStdString() << std::endl;

      if((unordered_cones.size() != 3 && unordered_cones.size() != 4) ||
         (unordered_cones.size() == 3 && orb == SMP::Parallelogram ) ||
         (unordered_cones.size() == 4 && orb != SMP::Parallelogram)) {
        std::cerr << "Error: incompatible orbifold type and number of cones" << std::endl;
        std::cerr << "Types I, II & III require 3 selected vertices" << std::endl;
        std::cerr << "Type IV requires 4 selected vertices" << std::endl;
        QApplication::restoreOverrideCursor();
        return;
      }

      // Now, parameterize
      Parameterizer parameterizer(orb);

      // Mark cones in the seam mesh
      boost::unordered_map<s_vertex_descriptor, SMP::Cone_type> cmap;
      if(!SMP::locate_unordered_cones(sMesh, unordered_cones.begin(), unordered_cones.end(), cmap))
      {
        std::cerr << "Error: invalid cone or seam selection" << std::endl;
        QApplication::restoreOverrideCursor();
        return;
      }

      QApplication::setOverrideCursor(Qt::WaitCursor);

      // Fill the index property map
      typedef boost::unordered_map<s_vertex_descriptor, int> Indices;
      Indices indices;
      CGAL::Polygon_mesh_processing::connected_component(
             face(opposite(bhd, sMesh), sMesh),
             sMesh,
             boost::make_function_output_iterator(
             SMP::internal::Index_map_filler<Seam_mesh, Indices>(sMesh, indices)));
      boost::associative_property_map<Indices> vimap(indices);

      // Call to parameterizer
      status = parameterizer.parameterize(sMesh, bhd, cmap, uv_pm, vimap);
      break;
    }
    case PARAM_BTP:
    {
      std::cout << "Parameterize (BTP)..." << std::endl;
      new_item_name = tr("%1 (parameterized (BTP))").arg(poly_item->name());
      typedef SMP::Barycentric_mapping_parameterizer_3<Seam_mesh> Parameterizer;
      status = SMP::parameterize(sMesh, Parameterizer(), bhd, uv_pm);
      break;
    }
    } //end switch

    std::cout << "Connected component " << current_component << ": ";
    if(status == SMP::OK) {
      std::cout << "success (in " << time.elapsed() << " ms)" << std::endl;
    } else {
      std::cerr << "failure: " << SMP::get_error_message(status) << std::endl;
      QApplication::restoreOverrideCursor();
      return;
    }

    if(status != SMP::OK)
      break;

  } //end for each component

  QApplication::restoreOverrideCursor();
  QPointF pmin(FLT_MAX, FLT_MAX), pmax(-FLT_MAX, -FLT_MAX);

  SMesh::Property_map<halfedge_descriptor, float> umap;
  SMesh::Property_map<halfedge_descriptor, float> vmap;

  umap = tMesh.add_property_map<halfedge_descriptor, float>("h:u", 0.0f).first;
  vmap = tMesh.add_property_map<halfedge_descriptor, float>("h:v", 0.0f).first;

  tMesh.property_stats(std::cerr);
  Base_face_graph::Halfedge_iterator it;
  for(it = tMesh.halfedges_begin();
      it != tMesh.halfedges_end();
      ++it)
  {
    Seam_mesh::halfedge_descriptor hd(*it);
    FT u = uv_pm[target(hd, sMesh)].x();
    FT v = uv_pm[target(hd, sMesh)].y();
    put(umap, *it, static_cast<float>(u));
    put(vmap, *it, static_cast<float>(v));
    if(u<pmin.x())
      pmin.setX(u);
    if(u>pmax.x())
      pmax.setX(u);
    if(v<pmin.y())
      pmin.setY(v);
    if(v>pmax.y())
      pmax.setY(v);
  }

  Components* components = new Components(0);
  components->resize(number_of_components);
  boost::graph_traits<Base_face_graph>::face_iterator bfit;

  for(bfit = faces(tMesh).begin();
      bfit != faces(tMesh).end();
      ++bfit)
  {
    components->at(fccmap[*bfit]).insert(*bfit);
  }

  Scene_textured_facegraph_item* new_item = new Scene_textured_facegraph_item(tMesh);
  UVItem *projection = new UVItem(components,new_item->textured_face_graph(), uv_borders, QRectF(pmin, pmax));
  projection->set_item_name(new_item_name);

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
  if(current_uv_item){
    Scene_textured_facegraph_item* t_item =
        qobject_cast<Scene_textured_facegraph_item*>(projections.key(current_uv_item));
   t_item->add_border_edges(std::vector<float>(0));
  }
  current_uv_item = projection;
  Scene_textured_facegraph_item* t_item =
      qobject_cast<Scene_textured_facegraph_item*>(projections.key(current_uv_item));
  t_item->add_border_edges(
        current_uv_item->concatenated_borders());
  if(dock_widget->isHidden()){
    dock_widget->setVisible(true);
    dock_widget->raise();
  }
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

void Polyhedron_demo_parameterization_plugin::on_actionIAP_triggered()
{
  std::cerr << "IAP...";
  parameterize(PARAM_IAP);
}

void Polyhedron_demo_parameterization_plugin::on_actionARAP_triggered()
{
  std::cerr << "ARAP...";
  parameterize(PARAM_ARAP);
}

void Polyhedron_demo_parameterization_plugin::on_actionOTE_triggered()
{
  std::cerr << "OTE...";
  parameterize(PARAM_OTE);
}

void Polyhedron_demo_parameterization_plugin::on_actionBTP_triggered()
{
  std::cerr << "BTP...";
  parameterize(PARAM_BTP);
}

#include "Parameterization_plugin.moc"
