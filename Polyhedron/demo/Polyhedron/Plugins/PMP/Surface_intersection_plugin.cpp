#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"
#include "Messages_interface.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QElapsedTimer>
#include <QMessageBox>

typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_polylines_item::Polyline Polyline_3;
typedef boost::graph_traits<Scene_face_graph_item::Face_graph>::face_descriptor face_descriptor;

using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;

class Polyhedron_demo_intersection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "surface_intersection_plugin.json")

public:

  bool applicable(QAction* action) const {

    if(scene->selectionIndices().size() != 2)
      return false;
    if(action == actionPolyhedronIntersection_3)
    {
      return
          qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().first())) &&
          qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().last()));
    }
    else if(action == actionSurfacePolylineIntersection)
    {
      if((qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().first())) &&
              qobject_cast<Scene_polylines_item*>(scene->item(scene->selectionIndices().last()))
              )
          ||
          (qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().last())) &&
           qobject_cast<Scene_polylines_item*>(scene->item(scene->selectionIndices().first()))
           ))
        return true;
      else
        return false;
    }
    else //if (action == actionPolylinesIntersection)
    {
      return  qobject_cast<Scene_polylines_item*>(scene->item(scene->selectionIndices().first())) &&
          qobject_cast<Scene_polylines_item*>(scene->item(scene->selectionIndices().last())) ;
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionPolyhedronIntersection_3
                             << actionSurfacePolylineIntersection
                             << actionPolylinesIntersection
                                ;
  }

  void init(QMainWindow* mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface* mi) {
    this->scene = scene_interface;
    this->mi = mi;
    actionPolyhedronIntersection_3 = new QAction("Surface Intersection", mw);
    actionPolyhedronIntersection_3->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionPolyhedronIntersection_3, SIGNAL(triggered()),
            this, SLOT(intersectionSurfaces()));
    actionSurfacePolylineIntersection = new QAction("Surface-Polyline Intersection", mw);
    actionSurfacePolylineIntersection->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionSurfacePolylineIntersection, SIGNAL(triggered()),
            this, SLOT(intersectionSurfacePolyline()));

    actionPolylinesIntersection = new QAction("Polylines Intersection", mw);
    actionPolylinesIntersection->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionPolylinesIntersection, SIGNAL(triggered()),
            this, SLOT(intersectionPolylines()));
  }

private:

  QAction*  actionPolyhedronIntersection_3;
  QAction*  actionSurfacePolylineIntersection;
  QAction*  actionPolylinesIntersection;
  Scene_interface *scene;
  Messages_interface* mi;

public Q_SLOTS:
  void intersectionSurfaces();
  void intersectionSurfacePolyline();
  void intersectionPolylines();

}; // end class Polyhedron_demo_intersection_plugin

void Polyhedron_demo_intersection_plugin::intersectionSurfaces()
{
  Scene_face_graph_item* itemA = NULL;
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_face_graph_item* itemB =
      qobject_cast<Scene_face_graph_item*>(scene->item(index));

    if(itemB)
    {
      if (itemA==NULL)
      {
        itemA = itemB;
        continue;
      }
      if(!is_triangle_mesh(*itemA->face_graph())
         || !is_triangle_mesh(*itemB->face_graph()))
      {
        CGAL::Three::Three::error("The two meshes must be triangle meshes.");
      }
      QApplication::setOverrideCursor(Qt::WaitCursor);

      Scene_polylines_item* new_item = new Scene_polylines_item();
     // perform Boolean operation
      QElapsedTimer time;
      time.start();

      try{
        PMP::surface_intersection(*itemA->polyhedron(),
                                  *itemB->polyhedron(),
                                  std::back_inserter(new_item->polylines),
                                  PMP::parameters::throw_on_self_intersection(true));
      }
      catch(CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception)
      {
        QMessageBox::warning((QWidget*)NULL,
          tr("Self-intersections Found"),
          tr("Some self-intersections were found amongst intersecting facets"));
        delete new_item;
        QApplication::restoreOverrideCursor();
        return;
      }

      QString name = tr("%1 intersection %2");

      new_item->setName(name.arg(itemA->name(), itemB->name()));
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      if (new_item->polylines.empty())
        delete new_item;
      else{
        new_item->setColor(Qt::green);
        new_item->setRenderingMode(Wireframe);
        scene->addItem(new_item);
        new_item->invalidateOpenGLBuffers();
      }

      QApplication::restoreOverrideCursor();
    }
  }
}
#include <CGAL/intersections.h>
void Polyhedron_demo_intersection_plugin::intersectionPolylines()
{
  typedef std::pair<
      std::pair<std::size_t, std::size_t>,
      std::pair<std::size_t, std::size_t> >  Poly_intersection;

  Scene_polylines_item* itemA = NULL;
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polylines_item* itemB =
      qobject_cast<Scene_polylines_item*>(scene->item(index));

    if(itemB)
    {
      if (itemA==NULL)
      {
        itemA = itemB;
        continue;
      }

      QApplication::setOverrideCursor(Qt::WaitCursor);

      Scene_points_with_normal_item* new_point_item = new Scene_points_with_normal_item();
      Scene_polylines_item* new_pol_item = new Scene_polylines_item();
     // perform Boolean operation
      QElapsedTimer time;
      time.start();
      std::vector<Polyline_3> polyA, polyB;
      Q_FOREACH(const Polyline_3& poly, itemA->polylines)
      {
        polyA.push_back(poly);
      }
      Q_FOREACH(const Polyline_3& poly, itemB->polylines)
      {
        polyB.push_back(poly);
      }

      std::vector<Poly_intersection> poly_intersections;
      PMP::internal::compute_polylines_polylines_intersection(polyA, polyB,std::back_inserter(poly_intersections) , Kernel());

      Q_FOREACH(const Poly_intersection& inter, poly_intersections)
      {
        Kernel::Segment_3 segA(polyA[inter.first.first][inter.first.second], polyA[inter.first.first][inter.first.second +1]);
        Kernel::Segment_3 segB(polyB[inter.second.first][inter.second.second], polyB[inter.second.first][inter.second.second+1]);

        CGAL::cpp11::result_of<Kernel::Intersect_3(Kernel::Segment_3, Kernel::Segment_3)>::type
            result = CGAL::intersection(segA, segB);
        if (result) {
          if (const Kernel::Segment_3* s = boost::get<Kernel::Segment_3>(&*result)) {
            Polyline_3 p;
            p.push_back(s->point(0));
            p.push_back(s->point(1));
            new_pol_item->polylines.push_back(p);
          } else {
            const Kernel::Point_3* p = boost::get<Kernel::Point_3 >(&*result);
            new_point_item->point_set()->insert(*p);
          }
        }
      }
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;



      if (new_point_item->point_set()->empty())
        delete new_point_item;
      else{
        QString point_name = tr("%1 intersection points %2");

        new_point_item->setName(point_name.arg(itemA->name(), itemB->name()));
        new_point_item->setColor(Qt::green);
        new_point_item->setRenderingMode(Points);
        new_point_item->setPointSize(5);
        scene->addItem(new_point_item);
        new_point_item->invalidateOpenGLBuffers();
      }

      if (new_pol_item->polylines.empty())
        delete new_pol_item;
      else{
        QString polyline_name = tr("%1 intersection segments %2");
        new_pol_item->setName(polyline_name.arg(itemA->name(), itemB->name()));
        new_pol_item->setColor(Qt::green);
        new_pol_item->setRenderingMode(Wireframe);
        scene->addItem(new_pol_item);
        new_pol_item->invalidateOpenGLBuffers();
      }

      QApplication::restoreOverrideCursor();
    }
  }
}

void Polyhedron_demo_intersection_plugin::intersectionSurfacePolyline()
{
  typedef std::pair<std::size_t,
      std::pair<std::size_t, std::size_t> >Poly_intersection;

  Scene_item* it1 = scene->item(scene->selectionIndices().first());
  Scene_item* it2 = scene->item(scene->selectionIndices().last());;

  Scene_face_graph_item* itemA = qobject_cast<Scene_face_graph_item*>(it1);
  if(!itemA)
    itemA = qobject_cast<Scene_face_graph_item*>(it2);
  if(!itemA)
    return;

  Scene_polylines_item* itemB = qobject_cast<Scene_polylines_item*>(it2);
  if(!itemB)
    itemB = qobject_cast<Scene_polylines_item*>(it1);
  if(!itemB)
    return;

  if(!is_triangle_mesh(*itemA->face_graph()))
  {
    CGAL::Three::Three::error("The mesh must be a triangle mesh.");
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_points_with_normal_item* new_point_item = new Scene_points_with_normal_item();
  Scene_polylines_item* new_pol_item = new Scene_polylines_item();
  // perform Boolean operation
  QElapsedTimer time;
  time.start();
  Scene_face_graph_item::Face_graph tm = *itemA->face_graph();
  std::vector<face_descriptor> Afaces;
  Q_FOREACH(face_descriptor f, faces(tm))
  {
    Afaces.push_back(f);
  }
  boost::property_map<Scene_face_graph_item::Face_graph ,CGAL::vertex_point_t>::type
  vpm = get(CGAL::vertex_point,tm);

  std::vector<Polyline_3> polylines;
  //Polyline_3 polyline;
  Q_FOREACH(const Polyline_3 &pol, itemB->polylines)
  {
    polylines.push_back(pol);
  }
  std::vector<Poly_intersection> poly_intersections;
  PMP::internal::compute_face_polylines_intersection(faces(*itemA->face_graph()),
                                                     polylines,
                                                     *itemA->face_graph(),
                                                     std::back_inserter(poly_intersections),
                                                     CGAL::Polygon_mesh_processing::parameters::all_default());

  Q_FOREACH(const Poly_intersection& inter, poly_intersections)
  {
    Kernel::Segment_3 segment(polylines[inter.second.first][inter.second.second], polylines[inter.second.first][inter.second.second+1]);

    Kernel::Triangle_3 triangle(
          get(vpm, target(halfedge(Afaces[inter.first],tm), tm)),
        get(vpm, target(next(halfedge(Afaces[inter.first], tm), tm), tm)),
        get(vpm, target(next(next(halfedge(Afaces[inter.first], tm), tm), tm), tm))
        );

    CGAL::cpp11::result_of<Kernel::Intersect_3(Kernel::Segment_3, Kernel::Segment_3)>::type
        result = CGAL::intersection(triangle, segment);
    if (result) {
      if (const Kernel::Segment_3* s = boost::get<Kernel::Segment_3>(&*result)) {
        Polyline_3 p;
        p.push_back(s->point(0));
        p.push_back(s->point(1));
        new_pol_item->polylines.push_back(p);
      } else {
        const Kernel::Point_3* p = boost::get<Kernel::Point_3 >(&*result);
        new_point_item->point_set()->insert(*p);
      }
    }
  }
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;



  if (new_point_item->point_set()->empty())
    delete new_point_item;
  else{
    QString point_name = tr("%1 intersection points %2");

    new_point_item->setName(point_name.arg(itemA->name(), itemB->name()));
    new_point_item->setColor(Qt::green);
    new_point_item->setRenderingMode(Points);
    new_point_item->setPointSize(5);
    scene->addItem(new_point_item);
    new_point_item->invalidateOpenGLBuffers();
  }

  if (new_pol_item->polylines.empty())
    delete new_pol_item;
  else{
    QString polyline_name = tr("%1 intersection segments %2");
    new_pol_item->setName(polyline_name.arg(itemA->name(), itemB->name()));
    new_pol_item->setColor(Qt::green);
    new_pol_item->setRenderingMode(Wireframe);
    scene->addItem(new_pol_item);
    new_pol_item->invalidateOpenGLBuffers();
  }

  QApplication::restoreOverrideCursor();
}

#include "Surface_intersection_plugin.moc"
