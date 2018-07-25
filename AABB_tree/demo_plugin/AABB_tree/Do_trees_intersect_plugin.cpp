#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include "Messages_interface.h"
#include <CGAL/Three/Three.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_movable_sm_item.h"
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Rigid_mesh_collision_detection.h>
#include "Scene.h"

class DoTreesIntersectplugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    Q_FOREACH(Scene::Item_id i, scene->selectionIndices())
    {
      if(! qobject_cast<Scene_surface_mesh_item*>(scene->item(i)))
        return false;
    }
    return items.empty();
  }
  
  QList<QAction*> actions() const Q_DECL_OVERRIDE
  {
    return _actions;
  }
  
  void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface* mi) Q_DECL_OVERRIDE
  {
    this->messageInterface = mi;
    this->scene = sc;
    this->mw = mw;
    QAction *actionCreateTrees= new QAction(QString("Start Intersection Tests"), mw);
    actionCreateTrees->setProperty("submenuName", "AABB_tree");
    if(actionCreateTrees) {
      connect(actionCreateTrees, SIGNAL(triggered()),
              this, SLOT(start()));
      _actions << actionCreateTrees;
    }
  }
private Q_SLOTS:
  void start()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Q_FOREACH(Scene::Item_id i, scene->selectionIndices())
    {
      Scene_surface_mesh_item* item=qobject_cast<Scene_surface_mesh_item*>(scene->item(i));
      connect(item, &Scene_surface_mesh_item::aboutToBeDestroyed,
              this, &DoTreesIntersectplugin::cleanup);     
      
      CGAL::qglviewer::Vec pos((item->bbox().min(0) + item->bbox().max(0))/2.0, 
                               (item->bbox().min(1) + item->bbox().max(1))/2.0, 
                               (item->bbox().min(2) + item->bbox().max(2))/2.0);
      
      Scene_movable_sm_item* mov_item = new Scene_movable_sm_item(pos,item->face_graph(),"");
      connect(mov_item, &Scene_movable_sm_item::aboutToBeDestroyed,
              this, &DoTreesIntersectplugin::cleanup);
      connect(mov_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
              this, &DoTreesIntersectplugin::update_trees);
      mov_item->setName(item->name());
      mov_item->setRenderingMode(Wireframe);
      item->setVisible(false);
      items.push_back(mov_item);
      scene->setSelectedItem(scene->addItem(mov_item));
      mov_item->redraw();
    }
    connect(static_cast<Scene*>(scene), &Scene::itemIndexSelected,
            this, &DoTreesIntersectplugin::update_trees);
    Q_FOREACH(Scene_movable_sm_item* item, items)
    {
      meshes.push_back(*item->getFaceGraph());
    }
    col_det = new CGAL::Rigid_mesh_collision_detection<SMesh, EPICK>(meshes);
    update_trees();
    items.back()->setRenderingMode(Wireframe);
    
    QApplication::restoreOverrideCursor();
  }
  
  
public Q_SLOTS:
  void update_trees()
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first());
    
    Scene_movable_sm_item* sel_item = qobject_cast<Scene_movable_sm_item*>(scene->item(scene->mainSelectionIndex()));
    if(!sel_item)
      return;    
   
    std::size_t mesh_id = 0;
    std::size_t sel_id = 0;
    Q_FOREACH(Scene_movable_sm_item* item, items)
    {
      if(item == sel_item)
      {
        sel_id = mesh_id;
        break;
      }
      ++mesh_id;
    }
    mesh_id = 0;
    Q_FOREACH(Scene_movable_sm_item* item, items)
    {
      if(mesh_id == sel_id)
      {
        ++mesh_id;
        item->setColor(QColor(255,184,61));
        continue;
      }
      item->setColor(QColor(Qt::green));
      const double* matrix = item->manipulatedFrame()->matrix();
      item->setFMatrix(matrix);
      EPICK::Aff_transformation_3 translation(CGAL::TRANSLATION, -EPICK::Vector_3(item->center().x,
                                                                                  item->center().y,
                                                                                  item->center().z));
      EPICK::Aff_transformation_3 rota(
            matrix[0], matrix[4], matrix[8],matrix[12],
          matrix[1], matrix[5], matrix[9],matrix[13],
          matrix[2], matrix[6], matrix[10],matrix[14]);
      EPICK::Aff_transformation_3 transfo = 
          rota*translation;
      
      col_det->set_transformation(mesh_id++, transfo);
      
      item->setRenderingMode(Wireframe);
      item->itemChanged();
    }
    const double* matrix = sel_item->manipulatedFrame()->matrix();
    sel_item->setFMatrix(matrix);
    EPICK::Aff_transformation_3 translation(CGAL::TRANSLATION, -EPICK::Vector_3(sel_item->center().x,
                                                                                sel_item->center().y,
                                                                                sel_item->center().z));
    EPICK::Aff_transformation_3 rota(
          matrix[0], matrix[4], matrix[8],matrix[12],
        matrix[1], matrix[5], matrix[9],matrix[13],
        matrix[2], matrix[6], matrix[10],matrix[14]);
    EPICK::Aff_transformation_3 transfo = 
        rota*translation;
    std::vector<std::pair<std::size_t, bool> > inter_and_incl
        = col_det->
        set_transformation_and_all_intersections_and_all_inclusions(sel_id, transfo);
    for(std::size_t i=0; i<inter_and_incl.size(); ++i)
    {
      std::size_t id = inter_and_incl[i].first;
      bool including = inter_and_incl[i].second;
      if(including)
        items[id]->setColor(QColor(Qt::blue));
      else
        items[id]->setColor(QColor(Qt::red));
    }
    sel_item->setRenderingMode(Wireframe);
    sel_item->itemChanged();
    viewer->update();
  }
  
  void cleanup()
  {
    Q_FOREACH(Scene_movable_sm_item* item, items)
      if(item)
      {
        scene->erase(scene->item_id(item));
        item = nullptr;
      }
    items.clear();
    delete col_det;
    col_det = nullptr;
  }
  
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  CGAL::Rigid_mesh_collision_detection<SMesh, EPICK> *col_det;
  std::vector<Scene_movable_sm_item*> items;
  std::vector<SMesh> meshes;
};
#include "Do_trees_intersect_plugin.moc"
