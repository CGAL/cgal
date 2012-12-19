#define CGAL_USE_SEGMENT_APPROACH
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef CGAL_USE_SEGMENT_APPROACH
#include <CGAL/intersection_of_Polyhedra_3.h>
#else
#include <CGAL/intersection_Polyhedron_3_Polyhedron_3.h>
#endif
#include <CGAL/iterator.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polylines_item.h"

#include <boost/foreach.hpp>

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

class Polyhedron_demo_intersection_plugin :
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
    return QList<QAction*>() << actionPolyhedronIntersection_3;
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionPolyhedronIntersection_3 = new QAction("Intersect polyhedra", mw);
    if(actionPolyhedronIntersection_3) {
      connect(actionPolyhedronIntersection_3, SIGNAL(triggered()),
              this, SLOT(intersection()));
    }
  }

private:

  QAction*  actionPolyhedronIntersection_3;

public slots:
  void intersection();

}; // end class Polyhedron_demo_intersection_plugin


#ifdef CGAL_USE_SEGMENT_APPROACH
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
#endif

void Polyhedron_demo_intersection_plugin::intersection()
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
#ifndef CGAL_USE_SEGMENT_APPROACH
  if(!itemA || !itemB || itemA == itemB)
    return;
#else
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
#endif

  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_polylines_item* new_item = new Scene_polylines_item();
 // perform Boolean operation
  QTime time;
  time.start();

#ifdef CGAL_USE_SEGMENT_APPROACH
  typedef CGAL::Node_visitor_for_polyline_split<Polyhedron,
    Is_on_polyline, Set_vertex_corner> Split_visitor;

  CGAL::Intersection_of_Polyhedra_3<Polyhedron,
    Kernel,
    Split_visitor> polyline_intersections;

  typedef std::pair<Polyhedron::Facet_handle,
                    Polyhedron::Facet_handle> Pair_of_facet_handles;
  typedef std::vector<Pair_of_facet_handles> Polyline_info;
  typedef std::vector<Polyline_info> Polylines_infos;
  Polylines_infos polylines_infos;

  typedef Scene_polylines_item::Polylines_container Polylines_container;
  typedef std::back_insert_iterator<Polylines_container> To_container;
  typedef std::back_insert_iterator<Polylines_infos> To_infos;
  To_container to_container(new_item->polylines);
  To_infos to_infos(polylines_infos);
  CGAL::Dispatch_output_iterator<
    CGAL::cpp11::tuple<Polylines_container::value_type,
                       Polyline_info>,
    CGAL::cpp11::tuple<To_container, 
                       To_infos> > out_iterator(to_container,
                                                to_infos);

  if(itemA && itemB && itemA != itemB) {
    Polyhedron* A = itemA->polyhedron();
    Polyhedron* B = itemB->polyhedron();
    polyline_intersections(*A, *B, out_iterator);
  } else {
    if(itemA) {
      Polyhedron* A = itemA->polyhedron();
      polyline_intersections(*A, 
                             poly_ptrs.begin(),
                             poly_ptrs.end(),
                             out_iterator,
                             0);
    } else {
      polyline_intersections(poly_ptrs.begin(),
                             poly_ptrs.end(),
                             out_iterator,
                             0);
    }
  }
  QStringList polylines_metadata;
  BOOST_FOREACH(Polyline_info& info, polylines_infos) {
    std::set<int> indices;
    BOOST_FOREACH(Pair_of_facet_handles p, info)
    {
      indices.insert(p.first->patch_id());
      indices.insert(p.second->patch_id());
    }
    QString metadata;
    BOOST_FOREACH(int index, indices) {
      metadata = metadata + QString(" %1").arg(index);
    }
    std::cerr << "new polyline metadata: " << qPrintable(metadata) << "\n";
    polylines_metadata << metadata;
  }
  new_item->setProperty("polylines metadata", polylines_metadata);
  new_item->setName(tr("intersection"));
#else
  Polyhedron* A = itemA->polyhedron();
  Polyhedron* B = itemB->polyhedron();
  CGAL::intersection_Polyhedron_3_Polyhedron_3(*A, *B, std::back_inserter(new_item->polylines));

  QString name = tr("%1 intersection %2");
   
  new_item->setName(name.arg(itemA->name(), itemB->name()));
  itemA->setRenderingMode(Wireframe);
  itemB->setRenderingMode(Wireframe);
  scene->itemChanged(indexA);
  scene->itemChanged(indexB);
#endif
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  new_item->setColor(Qt::green);
  new_item->setRenderingMode(Wireframe);
  scene->addItem(new_item);

  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_intersection_plugin, Polyhedron_demo_intersection_plugin)

#include "Polyhedron_demo_intersection_plugin.moc"
