#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QApplication>

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/optimal_bounding_box.h>

using namespace CGAL::Three;

typedef Scene_surface_mesh_item Scene_facegraph_item;

class Create_obb_mesh_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const;

  bool applicable(QAction*) const
  {
    if(scene->mainSelectionIndex() != -1
       && scene->item(scene->mainSelectionIndex())->isFinite())
      return true;
    return false;
  }

protected:
  void gather_mesh_points(std::vector<Point_3>& points);
  void obb();

public Q_SLOTS:
  void createObb()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    obb();
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionObb;
}; // end Create_obb_mesh_plugin class

void Create_obb_mesh_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionObb = new QAction(tr("Create &Optimal Bbox Mesh"), mainWindow);
  actionObb->setObjectName("createObbMeshAction");
  connect(actionObb, SIGNAL(triggered()), this, SLOT(createObb()));
}

QList<QAction*> Create_obb_mesh_plugin::actions() const
{
  return QList<QAction*>() << actionObb;
}

void Create_obb_mesh_plugin::gather_mesh_points(std::vector<Point_3>& points)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_facegraph_item* item = qobject_cast<Scene_facegraph_item*>(scene->item(index));

  Scene_polyhedron_selection_item* selection_item =
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item || selection_item)
  {
    typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type PointPMap;
    typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

    std::vector<vertex_descriptor> selected_vertices;

    if(item != NULL)
    {
      FaceGraph& pmesh = *item->polyhedron();
      selected_vertices.assign(vertices(pmesh).begin(), vertices(pmesh).end());
      PointPMap pmap = get(CGAL::vertex_point, pmesh);
      for(vertex_descriptor v : selected_vertices)
        points.push_back(get(pmap, v));

    }
    else if(selection_item != NULL) // using selection of faces
    {
      FaceGraph& pmesh = *selection_item->polyhedron();
      for(face_descriptor f : selection_item->selected_facets)
      {
        for(vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
          selected_vertices.push_back(v);
      }

      PointPMap pmap = get(CGAL::vertex_point, pmesh);
      for(vertex_descriptor v : selected_vertices)
        points.push_back(get(pmap, v));
    }

    CGAL_assertion(points.size() >= 3);
  }

  if(point_set_item)
  {
    Point_set* points_set = point_set_item->point_set();
    if(points_set == NULL)
        return;

    std::cout << "points_set->size()= " << points_set->size() << std::endl;
    for(const Point_3& p : points_set->points())
      points.push_back(p);
  }
}

void Create_obb_mesh_plugin::obb()
{
  // gather point coordinates
  std::vector<Point_3> points;
  gather_mesh_points(points);

  // find obb
  std::array<Point_3, 8> obb_points;
  CGAL::oriented_bounding_box(points, obb_points);

  Scene_facegraph_item* item;
  SMesh* p = new SMesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);
  item = new Scene_facegraph_item(p);

  item->setName("Optimal bbox mesh");
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
}

#include "Create_obb_mesh_plugin.moc"
