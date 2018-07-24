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
    t1 = 0;
    t2 = 0;
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
            [this](){
      if(query_item)
        scene->erase(scene->item_id(query_item));
      cleanup();
    });
    connect(itemB, &Scene_surface_mesh_item::aboutToBeDestroyed,
            [this](){
      if(query_item)
        scene->erase(scene->item_id(query_item));
      cleanup();
    });
    
    
    SMesh* tm = itemA->face_graph();
    SMesh* tm2 = itemB->face_graph();
    t1 = new Tree (tm->faces_begin(), tm->faces_end(), *tm);
    t2 = new Tree (tm2->faces_begin(), tm2->faces_end(), *tm2);
    CGAL::qglviewer::Vec pos((itemB->bbox().min(0) + itemB->bbox().max(0))/2.0, 
                          (itemB->bbox().min(1) + itemB->bbox().max(1))/2.0, 
                          (itemB->bbox().min(2) + itemB->bbox().max(2))/2.0);
    
    query_item = new Scene_movable_sm_item(pos,tm2,"");
    connect(query_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
            [this](){
      cleanup();
    });
    connect(query_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, &DoTreesIntersectplugin::update_trees);
    query_item->setRenderingMode(Flat);
    query_item->setName(itemB->name());
    itemB->setVisible(false);
    itemA->setVisible(false);
    scene->setSelectedItem(scene->addItem(query_item));
    base_item = new Scene_movable_sm_item(pos,tm,"");
    connect(base_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
            [this](){
      cleanup();
    });
    connect(base_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, &DoTreesIntersectplugin::update_trees);
    base_item->setRenderingMode(Wireframe);
    base_item->setName(itemA->name());
    scene->addItem(base_item);
    update_trees();
    query_item->redraw();
    base_item->redraw();
    QApplication::restoreOverrideCursor();
  }
public Q_SLOTS:
  void update_trees()
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first());
    
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
    t1->traits().set_transformation(transfoB);
    
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
    t2->traits().set_transformation(transfoA);
    if(t2->do_intersect(*t1))
      query_item->setColor(QColor(Qt::red));
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
          query_item->setColor(QColor(Qt::green));
      #if 0
      }
#endif
        
    }
    
    viewer->update();
  }
  
  void cleanup()
  {
    if(t1)
      delete t1;
    t1 = NULL;
    if(t2)
      delete t2;
    t2 = NULL;
    query_item = NULL;
    base_item = NULL;
  }
  
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  Tree *t1, *t2;
  Scene_movable_sm_item* query_item;
  Scene_movable_sm_item* base_item;
};
#include "Do_trees_intersect_plugin.moc"
