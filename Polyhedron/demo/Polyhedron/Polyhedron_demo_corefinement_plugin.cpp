#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_item.h"
#include "Scene_combinatorial_map_item.h"
#include "Polyhedron_type.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polylines_item.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

//#define PRINT_EACH_VOLUME

class Polyhedron_demo_corefinement_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:

  bool applicable() const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionPolyhedronCorefinement_3;
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionPolyhedronCorefinement_3 = new QAction("Polyhedra corefinement", mw);
    if(actionPolyhedronCorefinement_3) {
      connect(actionPolyhedronCorefinement_3, SIGNAL(triggered()),
              this, SLOT(corefinement()));
    }
  }

private:

  QAction*  actionPolyhedronCorefinement_3;

public slots:
  void corefinement();

}; // end class Polyhedron_demo_corefinement_plugin


struct Is_on_polyline{
  bool operator()(Polyhedron::Halfedge_handle he) const { 
    return he->is_feature_edge();
  }
};

struct Set_vertex_corner{
  template <typename Info>
  void add_info_to_node(int, Polyhedron*,const Info&) {
  }

  void operator()(Polyhedron::Vertex_handle v, int, Polyhedron*) {
    ++v->nb_of_feature_edges;
  }
};

void Polyhedron_demo_corefinement_plugin::corefinement()
{
  int indexA = scene->selectionAindex();
  int indexB = scene->selectionBindex();

  Scene_polyhedron_item* itemA = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(indexA));
  Scene_polyhedron_item* itemB = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(indexB));
  if(!itemA || !itemB || itemA == itemB)
  {
    Q_FOREACH(int index, scene->selectionIndices()) {
      Scene_polyhedron_item* item = 
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
      if(!item)
        return;
    }
    if(scene->selectionIndices().size() == 2) {
      indexA = scene->selectionIndices()[0];
      indexB = scene->selectionIndices()[1];
      itemA = 
        qobject_cast<Scene_polyhedron_item*>(scene->item(indexA));
      itemB = 
        qobject_cast<Scene_polyhedron_item*>(scene->item(indexB));
    }
  }
  std::vector<Polyhedron*> poly_ptrs;
  Q_FOREACH(int index, scene->selectionIndices()) {
    Scene_polyhedron_item* item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(!item)
      return;
    else if(item != itemA) {
      poly_ptrs.push_back(item->polyhedron());
      if(poly_ptrs.back() == 0) return;
    }
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if(itemA && itemB && itemA != itemB) {
  
   // perform Boolean operation
    QTime time;
    time.start();

    //since the visitor modify input polyhedra, we make copies.
    Polyhedron A( *itemA->polyhedron() );
    Polyhedron B( *itemB->polyhedron() );
    
    itemA->setVisible(false);
    itemB->setVisible(false);
    
    Scene_polylines_item* new_item = new Scene_polylines_item();
    
  #ifdef _COREFINEMENT_OUTPUT_IS_POLYHEDRON
    typedef CGAL::Polyhedron_corefinement<Polyhedron> Corefinement;
    Corefinement corefinement;

    typedef std::list<std::pair<Polyhedron*,int> > Decomposition;
    Decomposition decomposition;
    
    #ifndef PRINT_EACH_VOLUME
    int features=Corefinement::Join_tag+Corefinement::Intersection_tag+Corefinement::P_minus_Q_tag+Corefinement::Q_minus_P_tag;
    #else
    int features=Corefinement::Decomposition_tag;
    #endif
    corefinement(A, B, std::back_inserter(new_item->polylines),std::back_inserter(decomposition),features);
    
    
    for (Decomposition::iterator it=decomposition.begin();it!=decomposition.end();++it)
    {
      Polyhedron* new_poly=it->first;
      Scene_polyhedron_item* new_item_bool= new Scene_polyhedron_item(new_poly);
      new_item_bool->setName(
        QString::fromStdString(
          corefinement.get_type_str(itemA->name().toStdString(),itemB->name().toStdString(),it->second)
        )
      );
      scene->addItem(new_item_bool);
    }
  #else  
    Scene_combinatorial_map_item* cmap_item = new Scene_combinatorial_map_item(scene,static_cast<void*>(&A)); 
    typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron> Split_visitor; 
    cmap_item->m_combinatorial_map=new Combinatorial_map_3();
    Split_visitor visitor(cmap_item->m_combinatorial_map);
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections(visitor);
    polyline_intersections(A, B, std::back_inserter(new_item->polylines));
    cmap_item->setName(QString("%1_and_%2_corefined").arg(itemA->name()).arg(itemB->name()));
    scene->addItem(cmap_item);
  #endif
    new_item->setName(tr("boundary intersection"));
    new_item->setColor(Qt::green);
    new_item->setRenderingMode(Wireframe);
    scene->addItem(new_item);
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
      
  }

  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_corefinement_plugin, Polyhedron_demo_corefinement_plugin)

#include "Polyhedron_demo_corefinement_plugin.moc"
