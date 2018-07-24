#include <QtCore/qglobal.h>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Scene_interface.h>
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Messages_interface.h"

#include <QAction>
#include <QApplication>
#include <QMainWindow>
#include <QObject>
#include <QInputDialog>
#include <limits>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/internal/repair_extra.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

using namespace CGAL::Three;
class Polyhedron_demo_repair_polyhedron_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;

    actionRemoveIsolatedVertices = new QAction(tr("Remove Isolated Vertices"), mw);
    actionRemoveDegenerateFaces = new QAction(tr("Remove Degenerate Faces"), mw);
    actionRemoveSelfIntersections = new QAction(tr("Remove Self-Intersections"), mw);
    actionStitchCloseBorderHalfedges = new QAction(tr("Stitch Close Border Halfedges"), mw);
    actionDuplicateNMVertices = new QAction(tr("Duplicate Non-Manifold Vertices"), mw);
    actionAutorefine = new QAction(tr("Autorefine Mesh"), mw);
    actionAutorefineAndRMSelfIntersections = new QAction(tr("Autorefine and Remove Self-Intersections"), mw);

    actionRemoveIsolatedVertices->setObjectName("actionRemoveIsolatedVertices");
    actionRemoveDegenerateFaces->setObjectName("actionRemoveDegenerateFaces");
    actionRemoveSelfIntersections->setObjectName("actionRemoveSelfIntersections");
    actionStitchCloseBorderHalfedges->setObjectName("actionStitchCloseBorderHalfedges");
    actionDuplicateNMVertices->setObjectName("actionDuplicateNMVertices");
    actionAutorefine->setObjectName("actionAutorefine");
    actionAutorefineAndRMSelfIntersections->setObjectName("actionAutorefineAndRMSelfIntersections");

    actionRemoveDegenerateFaces->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionStitchCloseBorderHalfedges->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionRemoveSelfIntersections->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionRemoveIsolatedVertices->setProperty("subMenuName", "Polygon Mesh Processing/Repair");
    actionDuplicateNMVertices->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionAutorefine->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionAutorefineAndRMSelfIntersections->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");

    autoConnectActions();
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionRemoveDegenerateFaces
                             << actionRemoveIsolatedVertices
                             << actionRemoveSelfIntersections
                             << actionStitchCloseBorderHalfedges
                             << actionDuplicateNMVertices
                             << actionAutorefine
                             << actionAutorefineAndRMSelfIntersections;
  }

  bool applicable(QAction*) const
  {
    int item_id = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(scene->item(item_id)) ||
           qobject_cast<Scene_surface_mesh_item*>(scene->item(item_id));
  }
  template <typename Item>
  void on_actionRemoveIsolatedVertices_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionRemoveDegenerateFaces_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionRemoveSelfIntersections_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionStitchCloseBorderHalfedges_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionDuplicateNMVertices_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionAutorefine_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionAutorefineAndRMSelfIntersections_triggered(Scene_interface::Item_id index);

public Q_SLOTS:
  void on_actionRemoveIsolatedVertices_triggered();
  void on_actionRemoveDegenerateFaces_triggered();
  void on_actionRemoveSelfIntersections_triggered();
  void on_actionStitchCloseBorderHalfedges_triggered();
  void on_actionDuplicateNMVertices_triggered();
  void on_actionAutorefine_triggered();
  void on_actionAutorefineAndRMSelfIntersections_triggered();

private:
  QAction* actionRemoveIsolatedVertices;
  QAction* actionRemoveDegenerateFaces;
  QAction* actionRemoveSelfIntersections;
  QAction* actionStitchCloseBorderHalfedges;
  QAction* actionDuplicateNMVertices;
  QAction* actionAutorefine;
  QAction* actionAutorefineAndRMSelfIntersections;

  Messages_interface* messages;
}; // end Polyhedron_demo_repair_polyhedron_plugin

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveIsolatedVertices_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    std::size_t nbv =
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(
        *poly_item->polyhedron());
    messages->information(tr(" %1 isolated vertices have been removed.")
      .arg(nbv));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveIsolatedVertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionRemoveIsolatedVertices_triggered<Scene_polyhedron_item>(index);
  on_actionRemoveIsolatedVertices_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveDegenerateFaces_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    std::size_t nbv =
      CGAL::Polygon_mesh_processing::remove_degenerate_faces(
      *poly_item->polyhedron());
    messages->information(tr(" %1 degenerate faces have been removed.")
      .arg(nbv));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveDegenerateFaces_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionRemoveDegenerateFaces_triggered<Scene_polyhedron_item>(index);
  on_actionRemoveDegenerateFaces_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveSelfIntersections_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    bool solved =
      CGAL::Polygon_mesh_processing::remove_self_intersections(
      *poly_item->polyhedron());
    if (!solved)
      messages->information(tr("Some self-intersection could not be fixed"));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveSelfIntersections_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionRemoveSelfIntersections_triggered<Scene_polyhedron_item>(index);
  on_actionRemoveSelfIntersections_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefine_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    try{
      CGAL::Polygon_mesh_processing::experimental::autorefine(*poly_item->polyhedron());
    }
    catch(CGAL::Polygon_mesh_processing::Corefinement::Triple_intersection_exception)
    {
      messages->warning(tr("The result of the requested operation is not handled (triple intersection)."));
    }
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefine_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionAutorefine_triggered<Scene_polyhedron_item>(index);
  on_actionAutorefine_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefineAndRMSelfIntersections_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    try{
      bool solved = 
        CGAL::Polygon_mesh_processing::experimental::
          autorefine_and_remove_self_intersections(*poly_item->polyhedron());
      if (!solved)
        messages->information(tr("Self-intersection could not be removed due to non-manifold edges in the output"));
    }
    catch(CGAL::Polygon_mesh_processing::Corefinement::Triple_intersection_exception)
    {
      messages->warning(tr("The result of the requested operation is not handled (triple intersection)."));
    }
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefineAndRMSelfIntersections_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionAutorefineAndRMSelfIntersections_triggered<Scene_polyhedron_item>(index);
  on_actionAutorefineAndRMSelfIntersections_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionStitchCloseBorderHalfedges_triggered(Scene_interface::Item_id index)
{
  typedef typename boost::graph_traits<typename Item::Face_graph>::halfedge_descriptor halfedge_descriptor;
  namespace PMP =   CGAL::Polygon_mesh_processing;

  if (Item* poly_item = qobject_cast<Item*>(scene->item(index)))
  {
    double epsilon = QInputDialog::getDouble(mw,
                                             QString("Choose Epsilon"),
                                             QString("Snapping distance for endpoints"),
                                             0,
                                             -(std::numeric_limits<double>::max)(),
                                             (std::numeric_limits<double>::max)(), 10);

    std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > halfedges_to_stitch;
    PMP::collect_close_stitchable_boundary_edges(*poly_item->polyhedron(), epsilon,
                                                 get(boost::vertex_point, *poly_item->polyhedron()),
                                                 halfedges_to_stitch);
    PMP::stitch_borders(*poly_item->polyhedron(), halfedges_to_stitch);
    messages->information(tr(" %1 pairs of halfedges stitched.").arg(halfedges_to_stitch.size()));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionStitchCloseBorderHalfedges_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionStitchCloseBorderHalfedges_triggered<Scene_polyhedron_item>(index);
  on_actionStitchCloseBorderHalfedges_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionDuplicateNMVertices_triggered(Scene_interface::Item_id index)
{
  namespace PMP =   CGAL::Polygon_mesh_processing;

  if (Item* poly_item = qobject_cast<Item*>(scene->item(index)))
  {
    std::size_t nb_vd = PMP::duplicate_non_manifold_vertices(*poly_item->polyhedron());
    messages->information(tr(" %1 vertices created").arg(nb_vd));
    if (nb_vd)
    {
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionDuplicateNMVertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionDuplicateNMVertices_triggered<Scene_polyhedron_item>(index);
  on_actionDuplicateNMVertices_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

#include "Repair_polyhedron_plugin.moc"
