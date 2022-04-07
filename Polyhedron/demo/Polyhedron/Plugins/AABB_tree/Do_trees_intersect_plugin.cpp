#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include "Messages_interface.h"
#include <CGAL/Three/Three.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_movable_sm_item.h"
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Rigid_triangle_mesh_collision_detection.h>
#include "Scene.h"

class DoTreesIntersectplugin:
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:

  bool eventFilter(QObject *, QEvent *event) Q_DECL_OVERRIDE
  {
    if(event->type() != QEvent::KeyPress)
      return false;
    QKeyEvent * e = static_cast<QKeyEvent*>(event);
    if (e->key()==Qt::Key_W){
      do_transparency = !do_transparency;
      change_display();
      return true;
    }
    return false;
  }

  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    if(scene->selectionIndices().size() <2)
      return false;
    Q_FOREACH(Scene::Item_id i, scene->selectionIndices())
    {
      if(! qobject_cast<Scene_surface_mesh_item*>(scene->item(i)))
        return false;
    }
    return (! group_item);
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
    QAction *actionCreateTrees= new QAction(QString("Collision Detection"), mw);
    actionCreateTrees->setProperty("subMenuName", "Polygon Mesh Processing");
    actionCreateTrees->setProperty("submenuName", "AABB_tree");
    if(actionCreateTrees) {
      connect(actionCreateTrees, SIGNAL(triggered()),
              this, SLOT(start()));
      _actions << actionCreateTrees;
    }
    do_transparency = false;
    group_item = nullptr;
  }
private Q_SLOTS:
  void start()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Q_FOREACH(Scene::Item_id i, scene->selectionIndices())
    {
      Scene_surface_mesh_item* item=qobject_cast<Scene_surface_mesh_item*>(scene->item(i));
      if (!CGAL::is_triangle_mesh(*item->face_graph()))
      {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(mw, "Error", QString("%1 is not pure triangle. Aborting.").arg(item->name()));
        return;
      }
    }
    group_item = new Scene_group_item("Test Items");
    connect(group_item, &Scene_group_item::aboutToBeDestroyed,
            this, [this](){
      items.clear();
      if(col_det)
        delete col_det;
      col_det = nullptr;
      group_item = nullptr;});

    scene->addItem(group_item);
    Q_FOREACH(Scene::Item_id i, scene->selectionIndices())
    {
      Scene_surface_mesh_item* item=qobject_cast<Scene_surface_mesh_item*>(scene->item(i));
      connect(item, &Scene_surface_mesh_item::aboutToBeDestroyed,
              this, &DoTreesIntersectplugin::cleanup);

      CGAL::qglviewer::Vec pos(((item->bbox().min)(0) + (item->bbox().max)(0))/2.0,
                               ((item->bbox().min)(1) + (item->bbox().max)(1))/2.0,
                               ((item->bbox().min)(2) + (item->bbox().max)(2))/2.0);

      Scene_movable_sm_item* mov_item = new Scene_movable_sm_item(pos,item->face_graph(),"");
      connect(mov_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
              this, &DoTreesIntersectplugin::update_trees);
      mov_item->setName(item->name());
      if(do_transparency)
      {
        mov_item->setRenderingMode(Flat);
        mov_item->setAlpha(120);
      }
      else
      {
        mov_item->setRenderingMode(Wireframe);
      }
      item->setVisible(false);
      items.push_back(mov_item);
      scene->addItem(mov_item);
      scene->changeGroup(mov_item, group_item);
      group_item->lockChild(mov_item);
      mov_item->redraw();
    }
    scene->setSelectedItem(group_item->getChildren().last());
    connect(static_cast<Scene*>(scene), &Scene::itemIndexSelected,
            this, &DoTreesIntersectplugin::update_trees);
    col_det = new CGAL::Rigid_triangle_mesh_collision_detection<SMesh>();
    col_det->reserve(items.size());
    Q_FOREACH(Scene_movable_sm_item* item, items)
    {
      col_det->add_mesh(*item->getFaceGraph());
    }
    init_trees();
    static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first())->installEventFilter(this);
    QApplication::restoreOverrideCursor();
    CGAL::Three::Three::information("Press `W` to switch between Wireframe and Transparency mode.");
  }

public Q_SLOTS:
  void init_trees()
  {
    if(items.empty())
      return;
    Q_FOREACH(Scene_movable_sm_item* item, items)
      item->setColor(QColor(Qt::green));

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
        prev_ids.push_back(sel_id);
        continue;
      }
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

      if(do_transparency)
      {
        item->setRenderingMode(Flat);
        item->setAlpha(120);
      }
      else
      {
        item->setRenderingMode(Wireframe);
      }
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
    col_det->set_transformation(sel_id, transfo);
    std::vector<std::pair<std::size_t, bool> > inter_and_incl
        = col_det->get_all_intersections_and_inclusions(sel_id);
    for(std::size_t i=0; i<inter_and_incl.size(); ++i)
    {
      std::size_t id = inter_and_incl[i].first;
      bool including = inter_and_incl[i].second;
      if(including)
        items[id]->setColor(QColor(Qt::blue));
      else
        items[id]->setColor(QColor(Qt::red));
      prev_ids.push_back(id);
    }
    if(do_transparency)
    {
      sel_item->setRenderingMode(Flat);
      sel_item->setAlpha(120);
    }
    else
    {
      sel_item->setRenderingMode(Wireframe);
    }
    sel_item->itemChanged();
    viewer->update();
  }

  void update_trees()
  {
    if(items.empty())
      return;
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
        ++mesh_id;
        item->setColor(QColor(255,184,61));
        break;
      }
      ++mesh_id;
    }
    for(std::size_t i = 0; i< prev_ids.size(); ++i)
    {
      std::size_t id = prev_ids[i];
      if(id == sel_id)
      {
        continue;
      }
      Scene_movable_sm_item* item = items[id];
      item->setColor(QColor(Qt::green));
      item->itemChanged();
    }
    prev_ids.clear();
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
    col_det->set_transformation(sel_id, transfo);
    std::vector<std::pair<std::size_t, bool> > inter_and_incl
        = col_det->get_all_intersections_and_inclusions(sel_id);
    for(std::size_t i=0; i<inter_and_incl.size(); ++i)
    {
      std::size_t id = inter_and_incl[i].first;
      bool including = inter_and_incl[i].second;
      if(including)
        items[id]->setColor(QColor(Qt::blue));
      else
        items[id]->setColor(QColor(Qt::red));
      prev_ids.push_back(id);
    }
    prev_ids.push_back(sel_id);
    sel_item->itemChanged();
    viewer->update();
  }

  void cleanup()
  {
    if(!group_item)
      return;
    scene->erase(scene->item_id(group_item));
    group_item = nullptr;
    items.clear();
    prev_ids.clear();
    delete col_det;
    col_det = nullptr;
  }

  //switch transparent/wireframe.
  void change_display()
  {
    for(Scene_movable_sm_item* item : items)
    {
      if(do_transparency)
      {
        item->setRenderingMode(Flat);
        item->setAlpha(120);
      }
      else
      {
        item->setRenderingMode(Wireframe);
      }
    }
  }

private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  CGAL::Rigid_triangle_mesh_collision_detection<SMesh> *col_det;
  std::vector<Scene_movable_sm_item*> items;
  std::vector<std::size_t> prev_ids;
  Scene_group_item* group_item;
  bool do_transparency;
};
#include "Do_trees_intersect_plugin.moc"
