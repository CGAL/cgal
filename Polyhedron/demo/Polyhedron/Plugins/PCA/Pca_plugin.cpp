#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_plane_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/centroid.h>
#include <CGAL/bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/boost/graph/helpers.h>

#include "Kernel_type.h"
typedef Kernel::Plane_3 Plane;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Line_3 Line;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;

template <typename TriangleMesh, typename OutputIterator>
CGAL::Bbox_3 triangles(const TriangleMesh& mesh,
                         OutputIterator out)
{
  CGAL::Bbox_3 bb;
  typename boost::property_map<TriangleMesh,CGAL::vertex_point_t>::const_type vpm =
      get(CGAL::vertex_point, mesh);
  BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor fd, faces(mesh)){
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd = halfedge(fd,mesh);
    Triangle t(get(vpm,source(hd,mesh)),
               get(vpm,target(hd,mesh)),
               get(vpm,target(next(hd,mesh),mesh)));
    *out++ = t;
    bb = bb + t.bbox();
  }
  return bb;
}


using namespace CGAL::Three;
class Polyhedron_demo_pca_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  QList<QAction*> actions() const {
    return _actions;
  }


  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      scene = scene_interface;
      QAction *actionFitPlane = new QAction("Fit Plane", mw);
      QAction *actionFitLine = new QAction("Fit Line", mw);

      connect(actionFitPlane, SIGNAL(triggered()),
              this, SLOT(on_actionFitPlane_triggered()));
      connect(actionFitLine, SIGNAL(triggered()),
              this, SLOT(on_actionFitLine_triggered()));
      _actions << actionFitPlane
               << actionFitLine;
      Q_FOREACH(QAction* action, _actions)
        action->setProperty("subMenuName", "Principal Component Analysis");


  }


  bool applicable(QAction*) const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }


public Q_SLOTS:
  void on_actionFitPlane_triggered();
  void on_actionFitLine_triggered();

private:
  Scene_interface* scene;
  QList<QAction*> _actions;
}; // end Polyhedron_demo_pca_plugin

void Polyhedron_demo_pca_plugin::on_actionFitPlane_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  QApplication::setOverrideCursor(Qt::WaitCursor);

    std::list<Triangle> triangles;
  if(item){
    Polyhedron* pMesh = item->polyhedron();
    ::triangles(*pMesh,std::back_inserter(triangles));
  } else if(sm_item){
    SMesh* pMesh = sm_item->polyhedron();
    ::triangles(*pMesh,std::back_inserter(triangles));
    }
  if(! triangles.empty()){
    QString item_name = (item)?item->name():sm_item->name();
    // fit plane to triangles
    Plane plane;
    std::cout << "Fit plane...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());
    Scene_plane_item* new_item = new Scene_plane_item(this->scene);
    new_item->setPosition(center_of_mass.x(),
                          center_of_mass.y(),
                          center_of_mass.z());
    const Vector& normal = plane.orthogonal_vector();
    new_item->setNormal(normal.x(), normal.y(), normal.z());
    new_item->setName(tr("%1 (plane fit)").arg(item_name));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode((item)?item->renderingMode():sm_item->renderingMode());
    scene->addItem(new_item);

  }
  else
    {
      Scene_points_with_normal_item* item =
        qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

      if (item)
        {
          Point_set* points = item->point_set();

          // fit plane to triangles
          Plane plane;
          Point center_of_mass;
          std::cout << "Fit plane...";
         CGAL::linear_least_squares_fitting_3
           (points->points().begin(),points->points().end(),plane, center_of_mass,
            CGAL::Dimension_tag<0>());
          std::cout << "ok" << std::endl;

          // compute centroid
          Scene_plane_item* new_item = new Scene_plane_item(this->scene);
          new_item->setPosition(center_of_mass.x(),
                                center_of_mass.y(),
                                center_of_mass.z());
          const Vector& normal = plane.orthogonal_vector();
          new_item->setNormal(normal.x(), normal.y(), normal.z());
          new_item->setName(tr("%1 (plane fit)").arg(item->name()));
          new_item->setColor(Qt::magenta);
          new_item->setRenderingMode(item->renderingMode());
          scene->addItem(new_item);

        }
    }
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_pca_plugin::on_actionFitLine_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

    QApplication::setOverrideCursor(Qt::WaitCursor);

  CGAL::Bbox_3 bb;

    std::list<Triangle> triangles;
   if(item){
    Polyhedron* pMesh = item->polyhedron();
    bb = ::triangles(*pMesh,std::back_inserter(triangles));
  } else if(sm_item){
    SMesh* pMesh = sm_item->polyhedron();
    bb = ::triangles(*pMesh,std::back_inserter(triangles));
  }

  if(! triangles.empty()){
    QString item_name = (item)?item->name():sm_item->name();
    // fit line to triangles
    Line line;
    std::cout << "Fit line...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());

    // compute bounding box diagonal
    Iso_cuboid bbox(bb);

    // compute scale for rendering using diagonal of bbox
    Point cmin = (bbox.min)();
    Point cmax = (bbox.max)();
    FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));

    // construct a 3D bar
    Vector u = line.to_vector();
    u = u / std::sqrt(u*u);

    Point a = center_of_mass + u * diag;
    Point b = center_of_mass - u * diag;

    Plane plane_a = line.perpendicular_plane(a);

    Vector u1 = plane_a.base1();
    u1 = u1 / std::sqrt(u1*u1);
    u1 = u1 * 0.01 * diag;
    Vector u2 = plane_a.base2();
    u2 = u2 / std::sqrt(u2*u2);
    u2 = u2 * 0.01 * diag;

    Point points[8];

    points[0] = a + u1;
    points[1] = a + u2;
    points[2] = a - u1;
    points[3] = a - u2;

    points[4] = b + u1;
    points[5] = b + u2;
    points[6] = b - u1;
    points[7] = b - u2;

    // add best fit line as new polyhedron bar
    SMesh* pFit = new SMesh;
    CGAL::make_hexahedron(points[0],
                          points[1],
                          points[2],
                          points[3],
                          points[4],
                          points[5],
                          points[6],
                          points[7],
                          *pFit);
    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(pFit);
    new_item->setName(tr("%1 (line fit)").arg(item_name));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode((item)?item->renderingMode(): sm_item->renderingMode());
    scene->addItem(new_item);
  }
  else
  {
    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

    if (item)
      {
        Point_set* point_set = item->point_set();
        Line line;
        Point center_of_mass;
        std::cout << "Fit line...";
        CGAL::linear_least_squares_fitting_3(point_set->points().begin(),
                                             point_set->points().end(),
                                             line, center_of_mass,
                                             CGAL::Dimension_tag<0>());
        std::cout << "ok" << std::endl;

        // compute bounding box diagonal
        Iso_cuboid bbox = CGAL::bounding_box(point_set->points().begin(),
                                             point_set->points().end());

        // compute scale for rendering using diagonal of bbox
        Point cmin = (bbox.min)();
        Point cmax = (bbox.max)();
        FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));

        // construct a 3D bar
        Vector u = line.to_vector();
        u = u / std::sqrt(u*u);

        Point a = center_of_mass + u * diag;
        Point b = center_of_mass - u * diag;

        Plane plane_a = line.perpendicular_plane(a);

        Vector u1 = plane_a.base1();
        u1 = u1 / std::sqrt(u1*u1);
        u1 = u1 * 0.01 * diag;
        Vector u2 = plane_a.base2();
        u2 = u2 / std::sqrt(u2*u2);
        u2 = u2 * 0.01 * diag;

        Point points[8];

        points[0] = a + u1;
        points[1] = a + u2;
        points[2] = a - u1;
        points[3] = a - u2;

        points[4] = b + u1;
        points[5] = b + u2;
        points[6] = b - u1;
        points[7] = b - u2;

        SMesh* pFit = new SMesh;
        CGAL::make_hexahedron(points[0],
                              points[1],
                              points[2],
                              points[3],
                              points[4],
                              points[5],
                              points[6],
                              points[7],
                              *pFit);
        Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(pFit);
        new_item->setName(tr("%1 (line fit)").arg(item->name()));
        new_item->setColor(Qt::magenta);
        new_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(new_item);
      }
  }

  QApplication::restoreOverrideCursor();
}

#include "Pca_plugin.moc"
