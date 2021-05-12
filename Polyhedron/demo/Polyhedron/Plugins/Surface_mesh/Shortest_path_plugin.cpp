#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Three.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_shortest_path_item.h"
#include "Messages_interface.h"
#include "Scene.h"
#include "ui_Shortest_path_widget.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QElapsedTimer>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>

#include <map>
#include <algorithm>
#include <vector>

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;

using namespace CGAL::Three;
class Polyhedron_demo_shortest_path_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
private:

  typedef boost::property_map<FaceGraph, boost::vertex_index_t>::type VertexIndexMap;
  typedef boost::property_map<FaceGraph, CGAL::halfedge_index_t>::type HalfedgeIndexMap;
  typedef boost::property_map<FaceGraph, CGAL::face_index_t>::type FaceIndexMap;
  typedef boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VertexPointMap;

  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, FaceGraph> Surface_mesh_shortest_path_traits;
  typedef CGAL::Surface_mesh_shortest_path<Surface_mesh_shortest_path_traits, VertexIndexMap, HalfedgeIndexMap, FaceIndexMap, VertexPointMap> Surface_mesh_shortest_path;

  struct ShortestPathsPointsVisitor
  {
    typedef std::vector<Surface_mesh_shortest_path::Point_3> Container;
    Container& m_container;

    ShortestPathsPointsVisitor(Container& container)
      : m_container(container)
    {
    }

    void point(const Surface_mesh_shortest_path::Point_3& point)
    {
      std::cout << point << std::endl;
      m_container.push_back(point);
    }
  };

  typedef std::map<Scene_facegraph_item*, Scene_polyhedron_shortest_path_item* > Shortest_paths_map;

public:

  QList<QAction*> actions() const
  {
      return QList<QAction*>() << actionMakeShortestPaths;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* messages)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->m_messages = messages;

    dock_widget = new QDockWidget(
          "Shortest Path"
          , mw);

    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    dock_widget->setWindowTitle(tr(
                                  "Shortest Path "
                                  ));
    addDockWidget(dock_widget);

    connect(ui_widget.Selection_type_combo_box, SIGNAL(currentIndexChanged(int)),  this, SLOT(on_Selection_type_combo_box_changed(int)));
    connect(ui_widget.Primitives_type_combo_box, SIGNAL(currentIndexChanged(int)), this, SLOT(on_Primitives_type_combo_box_changed(int)));

    actionMakeShortestPaths = new QAction(
          "Make Shortest Path"
          , this->mw);
    actionMakeShortestPaths->setProperty("subMenuName", "Triangulated Surface Mesh Shortest Paths");
    actionMakeShortestPaths->setObjectName("actionMakeShortestPaths");

    Scene* trueScene = dynamic_cast<Scene*>(scene_interface);
    // This is for later
    if(trueScene) {
        connect(trueScene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(item_about_to_be_destroyed(CGAL::Three::Scene_item*)));
        connect(trueScene, SIGNAL(newItem(int)), this, SLOT(new_item(int)));
    }
    autoConnectActions();
  }

  virtual void closure()
  {
    dock_widget->hide();
  }
private:

  Scene_polyhedron_shortest_path_item::Selection_mode get_selection_mode() const;
  Scene_polyhedron_shortest_path_item::Primitives_mode get_primitives_mode() const;

  void check_and_set_ids(FaceGraph* polyhedron);

public Q_SLOTS:
  void on_actionMakeShortestPaths_triggered();
  void on_Selection_type_combo_box_changed(int index);
  void on_Primitives_type_combo_box_changed(int index);
  void new_item(int index);
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item);

private:
  Shortest_paths_map m_shortestPathsMap;

  Messages_interface* m_messages;
  QAction* actionMakeShortestPaths;
  QDockWidget* dock_widget;
  Ui::Shortest_path ui_widget;
};

Scene_polyhedron_shortest_path_item::Selection_mode Polyhedron_demo_shortest_path_plugin::get_selection_mode() const
{
  return (Scene_polyhedron_shortest_path_item::Selection_mode) ui_widget.Selection_type_combo_box->currentIndex();
}

Scene_polyhedron_shortest_path_item::Primitives_mode Polyhedron_demo_shortest_path_plugin::get_primitives_mode() const
{
  return (Scene_polyhedron_shortest_path_item::Primitives_mode) ui_widget.Primitives_type_combo_box->currentIndex();
}

void Polyhedron_demo_shortest_path_plugin::item_about_to_be_destroyed(CGAL::Three::Scene_item* sceneItem)
{
    // if polyhedron item
    Scene_facegraph_item* polyhedronItem = qobject_cast<Scene_facegraph_item*>(sceneItem);
    if(polyhedronItem)
    {
      Shortest_paths_map::iterator found = m_shortestPathsMap.find(polyhedronItem);

      if (found != m_shortestPathsMap.end())
      {
        Scene_polyhedron_shortest_path_item* shortestPathItem = found->second;
        m_shortestPathsMap.erase(found);
        scene->erase(scene->item_id(shortestPathItem));
      }
    }

    // if polyhedron selection item
    Scene_polyhedron_shortest_path_item* shortestPathItem = qobject_cast<Scene_polyhedron_shortest_path_item*>(sceneItem);
    if(shortestPathItem)
    {
      Scene_facegraph_item* polyhedronItem = shortestPathItem->polyhedron_item();
      Shortest_paths_map::iterator found = m_shortestPathsMap.find(polyhedronItem);

      if (found != m_shortestPathsMap.end())
      {
        m_shortestPathsMap.erase(found);
      }
    }
  }

void Polyhedron_demo_shortest_path_plugin::new_item(int itemIndex)
{
  Scene_polyhedron_shortest_path_item* item = qobject_cast<Scene_polyhedron_shortest_path_item*>(scene->item(itemIndex));

  if (!item)
  {
    return;
  }

  if(item->polyhedron_item() == NULL)
  {
    Scene_facegraph_item* polyhedronItem = getSelectedItem<Scene_facegraph_item>();

    if(!polyhedronItem)
    {
      CGAL_assertion(item->polyhedron_item() == NULL); // which means it is coming from selection_io loader
      CGAL::Three::Three::information(tr("Error: please select corresponding polyhedron item from Geometric Objects list."));
      scene->erase(itemIndex);
      return;
    }

    if(!item->deferred_load(polyhedronItem, this->scene, this->m_messages, this->mw))
    {
      CGAL::Three::Three::information("Error: loading selection item is not successful!");
      scene->erase(itemIndex);
      return;
    }
  }

  check_and_set_ids(item->polyhedron_item()->polyhedron());

  Scene_polyhedron_shortest_path_item::Selection_mode selectionMode = get_selection_mode(); // Scene_polyhedron_shortest_path_item::INSERT_POINTS_MODE;

  std::cout << "Selection mode: " << selectionMode << std::endl;

  item->set_selection_mode(selectionMode);

  Scene_polyhedron_shortest_path_item::Primitives_mode primitivesMode = get_primitives_mode(); // Scene_polyhedron_shortest_path_item::FACE_MODE;

  std::cout << "Primitives mode: " << primitivesMode << std::endl;

  item->set_primitives_mode(primitivesMode);

  item->setRenderingMode(Flat);

  if(item->name() == "unamed")
  {
    item->setName(tr("%1 (shortest path computation item)").arg(item->polyhedron_item()->name()));
  }

  m_shortestPathsMap.insert(std::make_pair(item->polyhedron_item(), item));
}

void Polyhedron_demo_shortest_path_plugin::on_actionMakeShortestPaths_triggered()
{
  Scene_facegraph_item* polyhedronItem = getSelectedItem<Scene_facegraph_item>();
  if (polyhedronItem)
  {
    if (m_shortestPathsMap.find(polyhedronItem) == m_shortestPathsMap.end())
    {
      dock_widget->show();
      dock_widget->raise();
      // The other parts of initialization will be handled by the 'new_item' callback
      Scene_group_item* group = new Scene_group_item(QString("%1 Shortest Path").arg(polyhedronItem->name()));
      Scene_polyhedron_shortest_path_item* sp_item =new Scene_polyhedron_shortest_path_item(polyhedronItem, this->scene, this->m_messages, this->mw);
      sp_item->setName(tr("Source Points for %1").arg(polyhedronItem->name()));
      scene->addItem(sp_item);
      scene->addItem(group);
      scene->changeGroup(sp_item, group);
    }
    else
    {
      CGAL::Three::Three::warning(tr("A shortest path item for this polyhedron already exists (only one allowed per for now)"));
    }
  }
  else
  {
    CGAL::Three::Three::warning("No polyhedron selected.");
  }
}

void Polyhedron_demo_shortest_path_plugin::on_Selection_type_combo_box_changed(int index)
{
  std::cout << "Selection mode changed: " << index << std::endl;

  for (Shortest_paths_map::iterator it = m_shortestPathsMap.begin(); it != m_shortestPathsMap.end(); ++it)
  {
    it->second->set_selection_mode(get_selection_mode());
  }
}

void Polyhedron_demo_shortest_path_plugin::on_Primitives_type_combo_box_changed(int index)
{
  std::cout << "Primitives mode changed: " << index << std::endl;

  for (Shortest_paths_map::iterator it = m_shortestPathsMap.begin(); it != m_shortestPathsMap.end(); ++it)
  {
    it->second->set_primitives_mode(get_primitives_mode());
  }
}

void Polyhedron_demo_shortest_path_plugin::check_and_set_ids(FaceGraph* polyhedron)
{
  typedef boost::graph_traits<FaceGraph>::vertex_iterator vertex_iterator;
  typedef boost::graph_traits<FaceGraph>::halfedge_iterator halfedge_iterator;
  typedef boost::graph_traits<FaceGraph>::face_iterator face_iterator;

  boost::property_map<Face_graph, boost::vertex_index_t>::type vimap
      = get(boost::vertex_index, *polyhedron);

  vertex_iterator testVertex1 = vertices(*polyhedron).begin();
  vertex_iterator testVertex2 = ++vertices(*polyhedron).begin();
  if(get(vimap, *testVertex1) == get(vimap, *testVertex2))
  {
    polyhedron->collect_garbage();
  }
  boost::property_map<Face_graph, boost::halfedge_index_t>::type himap
      = get(boost::halfedge_index, *polyhedron);

  halfedge_iterator testHalfedge1 = halfedges(*polyhedron).begin();
  halfedge_iterator testHalfedge2 = ++halfedges(*polyhedron).begin();

  if (get(himap, *testHalfedge1) == get(himap, *testHalfedge2))
  {
    polyhedron->collect_garbage();
  }

  face_iterator testFacet1 = faces(*polyhedron).begin();
  face_iterator testFacet2 = ++faces(*polyhedron).begin();
  boost::property_map<Face_graph, boost::face_index_t>::type fimap
      = get(boost::face_index, *polyhedron);

  if (get(fimap, *testFacet1) == get(fimap, *testFacet2))
  {
    polyhedron->collect_garbage();
  }
}



#include "Shortest_path_plugin.moc"
