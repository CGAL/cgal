#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QApplication>

#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/optimal_bounding_box.h>

using namespace CGAL::Three;

typedef Scene_surface_mesh_item Scene_facegraph_item;

class Create_obb_mesh_plugin
  : public QObject,
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
    bool at_least_one_non_empty = false;
    for(int index : scene->selectionIndices())
    {
      Scene_item* item = scene->item(index);
      if(!item->isFinite())
        return false;

      Scene_facegraph_item* sm_item = qobject_cast<Scene_facegraph_item*>(item);
      Scene_polygon_soup_item* ps_item = qobject_cast<Scene_polygon_soup_item*>(item);
      Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(item);
      Scene_points_with_normal_item* pts_item = qobject_cast<Scene_points_with_normal_item*>(item);
      if(!sm_item && !ps_item && !selection_item && !pts_item)
        return false;

      if(!item->isEmpty())
        at_least_one_non_empty = true;
    }

    return at_least_one_non_empty;
  }

protected:
  void gather_mesh_points(std::vector<Point_3>& points);
  int dimensionality(const std::vector<Point_3>& points);
  bool obb();

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
};

void
Create_obb_mesh_plugin::
init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionObb = new QAction(tr("Create &Optimal Bounding Box"), mainWindow);
  actionObb->setObjectName("createObbMeshAction");
  connect(actionObb, SIGNAL(triggered()), this, SLOT(createObb()));
}

QList<QAction*>
Create_obb_mesh_plugin::
actions() const
{
  return QList<QAction*>() << actionObb;
}

int
Create_obb_mesh_plugin::
dimensionality(const std::vector<Point_3>& points)
{
  if(points.empty())
    return -1;

  int d = 0;
  Point_3 p0 = points[0], p1, p2;

  auto it = std::cbegin(points), end = std::cend(points);
  while(it != end)
  {
    if(p0 != *it)
    {
      p1 = *it;
      d = 1;
      break;
    }
    ++it;
  }

  while(it != end)
  {
    if(!collinear(p0, p1, *it))
    {
      p2 = *it;
      d = 2;
      break;
    }
    ++it;
  }

  while(it != end)
  {
    if(!coplanar(p0, p1, p2, *it))
    {
      d = 3;
      break;
    }
    ++it;
  }

  return d;
}

void
Create_obb_mesh_plugin::
gather_mesh_points(std::vector<Point_3>& points)
{
  for(int index : scene->selectionIndices())
  {
    Scene_item* item = scene->item(index);

    // Surface Mesh
    Scene_facegraph_item* sm_item = qobject_cast<Scene_facegraph_item*>(item);
    if(sm_item)
    {
      FaceGraph* sm_ptr = sm_item->polyhedron();
      if(sm_ptr == nullptr)
        continue;

      for(auto v : vertices(*sm_ptr))
        points.push_back(get(CGAL::vertex_point, *sm_ptr, v));

      continue;
    }

    // Polygon soup
    Scene_polygon_soup_item* ps_item = qobject_cast<Scene_polygon_soup_item*>(item);
    if(ps_item)
    {
      for(const Point_3& p : ps_item->points())
        points.push_back(p);

      continue;
    }

    // Selection
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(item);
    if(selection_item != nullptr)
    {
      FaceGraph* sm_ptr = selection_item->polyhedron();
      if(sm_ptr == nullptr)
        continue;

      auto vpm = get(CGAL::vertex_point, *sm_ptr);

      for(auto f : selection_item->selected_facets)
      {
        // @todo avoid duplication
        for(auto v : vertices_around_face(halfedge(f, *sm_ptr), *sm_ptr))
          points.push_back(get(vpm, v));
      }

      for(auto e : selection_item->selected_edges)
      {
        points.push_back(get(vpm, source(e, *sm_ptr)));
        points.push_back(get(vpm, target(e, *sm_ptr)));
      }

      for(auto v : selection_item->selected_vertices)
        points.push_back(get(vpm, v));

      continue;
    }

    // Point set
    Scene_points_with_normal_item* pts_item = qobject_cast<Scene_points_with_normal_item*>(item);
    if(pts_item)
    {
      Point_set* pts_ptr = pts_item->point_set();
      if(pts_ptr == nullptr)
          return;

      for(const Point_3& p : pts_ptr->points())
        points.push_back(p);

      continue;
    }
  }
}

bool
Create_obb_mesh_plugin::
obb()
{
  // gather point coordinates
  std::vector<Point_3> points;
  gather_mesh_points(points);

  const int d = dimensionality(points);
  if(d != 3)
  {
    std::cerr << "Dimensionality of the point set is: " << d << std::endl;
    QMessageBox::warning(mw, "Error", "At least 4 non-coplanar points are required to compute an OBB.");
    return false;
  }

  // compute the OBB
  std::array<Point_3, 8> obb_points;
  CGAL::oriented_bounding_box(points, obb_points);

  Scene_facegraph_item* item;
  SMesh* p = new SMesh;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);

  std::cout << "Optimal bounding box: " << obb_points[0]
            << "\n                      " << obb_points[7]
            << std::endl;

  QString name;
  for(int index : scene->selectionIndices())
  {
    Scene_item* item = scene->item(index);
    if(!item->isEmpty())
    {
      if(name.size() > 0)
      {
        name = name + " and others";
        break;
      }
      else
      {
        name = item->name();
      }
    }
  }

  item = new Scene_facegraph_item(p);
  item->setName(name + " (OBB)");
  item->setRenderingMode(Wireframe);
  item->setColor(Qt::black);
  scene->addItem(item);

  return true;
}

#include "Create_obb_mesh_plugin.moc"
