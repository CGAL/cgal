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
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_do_intersect_transform_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include "Scene.h"

typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_do_intersect_transform_traits<AABB_triangle_traits, EPICK> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

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
    const int indexA = scene->selectionAindex();
    const int indexB = scene->selectionBindex();
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(indexA))
        && qobject_cast<Scene_surface_mesh_item*>(scene->item(indexB))
        && base_item == NULL
        && query_item == NULL;
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
    QAction *actionCreateTrees= new QAction(QString("Start Intersection Tests(A/B)"), mw);
    actionCreateTrees->setProperty("submenuName", "AABB_tree");
    if(actionCreateTrees) {
      connect(actionCreateTrees, SIGNAL(triggered()),
              this, SLOT(start()));
      _actions << actionCreateTrees;
    }
    query_item = 0;
    base_item = 0;
  }
private Q_SLOTS:
  void start()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const int indexA = scene->selectionAindex();
    const int indexB = scene->selectionBindex();
    Scene_surface_mesh_item* itemA = qobject_cast<Scene_surface_mesh_item*>(scene->item(indexA));
    Scene_surface_mesh_item* itemB = qobject_cast<Scene_surface_mesh_item*>(scene->item(indexB));
    connect(itemA, &Scene_surface_mesh_item::aboutToBeDestroyed,
            this, &DoTreesIntersectplugin::cleanup);
    connect(itemB, &Scene_surface_mesh_item::aboutToBeDestroyed,
            this, &DoTreesIntersectplugin::cleanup);
    
    
    SMesh* tm = itemA->face_graph();
    SMesh* tm2 = itemB->face_graph();
    trees[tm] = new Tree (tm->faces_begin(), tm->faces_end(), *tm);
    trees[tm2] = new Tree (tm2->faces_begin(), tm2->faces_end(), *tm2);
    CGAL::qglviewer::Vec pos((itemB->bbox().min(0) + itemB->bbox().max(0))/2.0, 
                          (itemB->bbox().min(1) + itemB->bbox().max(1))/2.0, 
                          (itemB->bbox().min(2) + itemB->bbox().max(2))/2.0);
    
    query_item = new Scene_movable_sm_item(pos,tm2,"");
    connect(query_item, &Scene_movable_sm_item::aboutToBeDestroyed,
               this, &DoTreesIntersectplugin::cleanup);
    connect(query_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, &DoTreesIntersectplugin::update_trees);
    query_item->setRenderingMode(Flat);
    query_item->setName(itemB->name());
    itemB->setVisible(false);
    itemA->setVisible(false);
    scene->setSelectedItem(scene->addItem(query_item));
    pos = CGAL::qglviewer::Vec((itemA->bbox().min(0) + itemA->bbox().max(0))/2.0, 
                              (itemA->bbox().min(1) + itemA->bbox().max(1))/2.0, 
                              (itemA->bbox().min(2) + itemA->bbox().max(2))/2.0);
    
    base_item = new Scene_movable_sm_item(pos,tm,"");
    connect(base_item, &Scene_movable_sm_item::aboutToBeDestroyed,
            this, &DoTreesIntersectplugin::cleanup);
    connect(base_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, &DoTreesIntersectplugin::update_trees);
    base_item->setRenderingMode(Wireframe);
    base_item->setName(itemA->name());
    scene->addItem(base_item);
    update_trees();
    query_item->redraw();
    base_item->redraw();
    
    connect(static_cast<Scene*>(scene), &Scene::itemIndexSelected,
            this, &DoTreesIntersectplugin::update_trees);
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
    Scene_movable_sm_item* other_item = (base_item == sel_item) 
        ? query_item
        : base_item;
    
    const double* matrixB = base_item->manipulatedFrame()->matrix();
    base_item->setFMatrix(matrixB);
    EPICK::Aff_transformation_3 translationB(CGAL::TRANSLATION, -EPICK::Vector_3(base_item->center().x,
                                                                                 base_item->center().y,
                                                                                 base_item->center().z));
    EPICK::Aff_transformation_3 rotaB(
          matrixB[0], matrixB[4], matrixB[8],matrixB[12],
        matrixB[1], matrixB[5], matrixB[9],matrixB[13],
        matrixB[2], matrixB[6], matrixB[10],matrixB[14]);
    EPICK::Aff_transformation_3 transfoB = 
        rotaB*translationB;
    trees[base_item->getFaceGraph()]->traits().set_transformation(transfoB);
    
    const double* matrixA = query_item->manipulatedFrame()->matrix();
    query_item->setFMatrix(matrixA);
    EPICK::Aff_transformation_3 translationA(CGAL::TRANSLATION, -EPICK::Vector_3(query_item->center().x,
                                                                                 query_item->center().y,
                                                                                 query_item->center().z));
    EPICK::Aff_transformation_3 rotaA(
          matrixA[0], matrixA[4], matrixA[8],matrixA[12],
        matrixA[1], matrixA[5], matrixA[9],matrixA[13],
        matrixA[2], matrixA[6], matrixA[10],matrixA[14]);
    EPICK::Aff_transformation_3 transfoA = 
        rotaA*translationA;
    trees[query_item->getFaceGraph()]->traits().set_transformation(transfoA);
    
    if(trees[sel_item->getFaceGraph()]->do_intersect(*trees[other_item->getFaceGraph()]))
      sel_item->setColor(QColor(Qt::red));
    else
    {
#if 0
      typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type VPM;
      VPM vpm2 = get(CGAL::vertex_point, *query_item->getFaceGraph());
      CGAL::Side_of_triangle_mesh<SMesh, EPICK,
          VPM, Tree> sotm1(*t1);
      if(sotm1((transfoA).transform(vpm2[*query_item->getFaceGraph()->vertices().begin()])) != CGAL::ON_UNBOUNDED_SIDE)
      {        
        query_item->setColor(QColor(Qt::blue));
      }
      else
      {
        CGAL::Side_of_triangle_mesh<SMesh, EPICK,
            VPM, Tree> sotm2(*t2);
        if(sotm2(base_item->face_graph()->point(*base_item->face_graph()->vertices().begin())) != CGAL::ON_UNBOUNDED_SIDE)
          query_item->setColor(QColor(Qt::blue));
        else
#endif
          sel_item->setColor(QColor(Qt::green));
      #if 0
      }
#endif
        
    }
    sel_item->setRenderingMode(Flat);
    other_item->setRenderingMode(Wireframe);
    sel_item->itemChanged();
    other_item->itemChanged();
    viewer->update();
  }
  
  void cleanup()
  {
    if(query_item)
    {
      delete trees[query_item->getFaceGraph()];
      trees[query_item->getFaceGraph()] = NULL;
      scene->erase(scene->item_id(query_item));
      query_item = NULL;
    }
    if(base_item)
    {
      delete trees[base_item->getFaceGraph()];
      trees[base_item->getFaceGraph()] = NULL;
      scene->erase(scene->item_id(base_item));
      base_item = NULL;
    }
  }
  
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  std::map<SMesh*, Tree*> trees;
  Scene_movable_sm_item* query_item;
  Scene_movable_sm_item* base_item;
};
#include "Do_trees_intersect_plugin.moc"
