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
        && fixed_item == NULL
        && moving_item == NULL;
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
    moving_item = 0;
    fixed_item = 0;
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
      if(moving_item)
        scene->erase(scene->item_id(moving_item));
      cleanup();
    });
    connect(itemB, &Scene_surface_mesh_item::aboutToBeDestroyed,
            [this](){
      if(moving_item)
        scene->erase(scene->item_id(moving_item));
      cleanup();
    });
    
    
    SMesh* tm = itemA->face_graph();
    SMesh* tm2 = itemB->face_graph();
    t1 = new Tree (tm->faces_begin(), tm->faces_end(), *tm);
    t2 = new Tree (tm2->faces_begin(), tm2->faces_end(), *tm2);
    CGAL::qglviewer::Vec pos((itemB->bbox().min(0) + itemB->bbox().max(0))/2.0, 
                          (itemB->bbox().min(1) + itemB->bbox().max(1))/2.0, 
                          (itemB->bbox().min(2) + itemB->bbox().max(2))/2.0);
    
    moving_item = new Scene_movable_sm_item(pos,tm2,"");
    connect(moving_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
            [this](){
      cleanup();
    });
    fixed_item = itemA;
    connect(moving_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, &DoTreesIntersectplugin::update_trees);
    moving_item->setRenderingMode(Flat);
    moving_item->setName(itemB->name());
    itemB->setVisible(false);
    itemA->setRenderingMode(Wireframe);
    scene->addItem(moving_item);
    update_trees();
    moving_item->redraw();
    QApplication::restoreOverrideCursor();
  }
public Q_SLOTS:
  void update_trees()
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first());
    const double* matrix = moving_item->manipulatedFrame()->matrix();
    moving_item->setFMatrix(matrix);
    EPICK::Aff_transformation_3 translation(CGAL::TRANSLATION, -EPICK::Vector_3(moving_item->center().x,
                                                                               moving_item->center().y,
                                                                               moving_item->center().z));
    EPICK::Aff_transformation_3 rota(
          matrix[0], matrix[4], matrix[8],matrix[12],
        matrix[1], matrix[5], matrix[9],matrix[13],
        matrix[2], matrix[6], matrix[10],matrix[14]);
    EPICK::Aff_transformation_3 transfo = 
        rota*translation;
    t2->traits().set_transformation(transfo);
    if(t2->do_intersect(*t1))
      moving_item->setColor(QColor(Qt::red));
    else
    {
      typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type VPM;
      VPM vpm2 = get(CGAL::vertex_point, *moving_item->getFaceGraph());
      CGAL::Side_of_triangle_mesh<SMesh, EPICK,
          VPM, Tree> sotm1(*t1);
      if(sotm1((transfo).transform(vpm2[*moving_item->getFaceGraph()->vertices().begin()])) != CGAL::ON_UNBOUNDED_SIDE)
      {        
        moving_item->setColor(QColor(Qt::blue));
      }
      else
      {
        CGAL::Side_of_triangle_mesh<SMesh, EPICK,
            VPM, Tree> sotm2(*t2);
        if(sotm2(fixed_item->face_graph()->point(*fixed_item->face_graph()->vertices().begin())) != CGAL::ON_UNBOUNDED_SIDE)
          moving_item->setColor(QColor(Qt::blue));
        else
          moving_item->setColor(QColor(Qt::green));
      }
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
    moving_item = NULL;
    fixed_item = NULL;
  }
  
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  Tree *t1, *t2;
  Scene_movable_sm_item* moving_item;
  Scene_surface_mesh_item* fixed_item;
};
#include "Do_trees_intersect_plugin.moc"
