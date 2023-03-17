#include <QtCore/qglobal.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include <QAction>
#include <QApplication>
#include <QMainWindow>
#include <QObject>
#include <QInputDialog>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/internal/repair_extra.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include "ui_RemoveNeedlesDialog.h"
#include "ui_SelfSnapDialog.h"
#include <cmath>
#include <limits>
#include <vector>

using namespace CGAL::Three;
class Polyhedron_demo_repair_polyhedron_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "repair_polyhedron_plugin.json")

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
    actionExtractNMVertices = new QAction(tr("Extract Non-Manifold Vertices"), mw);
    actionMergeDuplicatedVerticesOnBoundaryCycles = new QAction(tr("Merge Duplicated Vertices on Boundary Cycles"), mw);
    actionAutorefine = new QAction(tr("Autorefine Mesh"), mw);
    actionAutorefineAndRMSelfIntersections = new QAction(tr("Autorefine and Remove Self-Intersections"), mw);
    actionRemoveNeedlesAndCaps = new QAction(tr("Remove Needles And Caps"));
    actionSnapBorders = new QAction(tr("Snap Boundaries"));

    actionRemoveIsolatedVertices->setObjectName("actionRemoveIsolatedVertices");
    actionRemoveDegenerateFaces->setObjectName("actionRemoveDegenerateFaces");
    actionRemoveSelfIntersections->setObjectName("actionRemoveSelfIntersections");
    actionStitchCloseBorderHalfedges->setObjectName("actionStitchCloseBorderHalfedges");
    actionDuplicateNMVertices->setObjectName("actionDuplicateNMVertices");
    actionExtractNMVertices->setObjectName("actionExtractNMVertices");
    actionMergeDuplicatedVerticesOnBoundaryCycles->setObjectName("actionMergeDuplicatedVerticesOnBoundaryCycles");
    actionAutorefine->setObjectName("actionAutorefine");
    actionAutorefineAndRMSelfIntersections->setObjectName("actionAutorefineAndRMSelfIntersections");
    actionRemoveNeedlesAndCaps->setObjectName("actionRemoveNeedlesAndCaps");
    actionSnapBorders->setObjectName("actionSnapBorders");

    actionRemoveDegenerateFaces->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionStitchCloseBorderHalfedges->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionRemoveSelfIntersections->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionRemoveIsolatedVertices->setProperty("subMenuName", "Polygon Mesh Processing/Repair");
    actionDuplicateNMVertices->setProperty("subMenuName", "Polygon Mesh Processing/Repair");
    actionExtractNMVertices->setProperty("subMenuName", "Polygon Mesh Processing/Repair");
    actionMergeDuplicatedVerticesOnBoundaryCycles->setProperty("subMenuName", "Polygon Mesh Processing/Repair");
    actionAutorefine->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionAutorefineAndRMSelfIntersections->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");
    actionSnapBorders->setProperty("subMenuName", "Polygon Mesh Processing/Repair/Experimental");

    autoConnectActions();
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionRemoveDegenerateFaces
                             << actionRemoveIsolatedVertices
                             << actionRemoveSelfIntersections
                             << actionStitchCloseBorderHalfedges
                             << actionDuplicateNMVertices
                             << actionExtractNMVertices
                             << actionMergeDuplicatedVerticesOnBoundaryCycles
                             << actionAutorefine
                             << actionAutorefineAndRMSelfIntersections
                             << actionRemoveNeedlesAndCaps
                             << actionSnapBorders;
  }

  bool applicable(QAction*) const
  {
    int item_id = scene->mainSelectionIndex();
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(item_id));
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
  void on_actionExtractNMVertices_triggered(Scene_interface::Item_id index);
  template <typename Item>
  void on_actionMergeDuplicatedVerticesOnBoundaryCycles_triggered(Scene_interface::Item_id index);
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
  void on_actionExtractNMVertices_triggered();
  void on_actionMergeDuplicatedVerticesOnBoundaryCycles_triggered();
  void on_actionAutorefine_triggered();
  void on_actionAutorefineAndRMSelfIntersections_triggered();
  void on_actionRemoveNeedlesAndCaps_triggered();
  void on_actionSnapBorders_triggered();

private:
  QAction* actionRemoveIsolatedVertices;
  QAction* actionRemoveDegenerateFaces;
  QAction* actionRemoveSelfIntersections;
  QAction* actionStitchCloseBorderHalfedges;
  QAction* actionDuplicateNMVertices;
  QAction* actionExtractNMVertices;
  QAction* actionMergeDuplicatedVerticesOnBoundaryCycles;
  QAction* actionAutorefine;
  QAction* actionAutorefineAndRMSelfIntersections;
  QAction* actionRemoveNeedlesAndCaps;
  QAction* actionSnapBorders;

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
    CGAL::Three::Three::information(tr(" %1 isolated vertices have been removed.")
      .arg(nbv));
    poly_item->setNbIsolatedvertices(0);
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveIsolatedVertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionRemoveIsolatedVertices_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveNeedlesAndCaps_triggered()
{
  QCursor tmp_cursor(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(tmp_cursor);

  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
  {
    return;
  }

  QDialog dialog;
  Ui::NeedleDialog ui;
  ui.setupUi(&dialog);
  ui.collapseBox->setValue(sm_item->diagonalBbox()*0.01);
  if(dialog.exec() != QDialog::Accepted)
    return;
  CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(*sm_item->face_graph(),
                                                                CGAL::parameters::cap_threshold(std::cos((ui.capBox->value()/180.0) * CGAL_PI))
                                                                                 .needle_threshold(ui.needleBox->value())
                                                                                 .collapse_length_threshold(ui.collapseBox->value())
                                                                                 .flip_triangle_height_threshold(ui.flipBox->value()));
  sm_item->invalidateOpenGLBuffers();
  sm_item->itemChanged();
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionSnapBorders_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(!sm_item)
  {
    return;
  }


  QDialog dialog;
  Ui::SelfSnapDialog ui;
  ui.setupUi(&dialog);
  connect(ui.use_local_tolerance, SIGNAL(toggled(bool)),
          ui.tolerances, SLOT(setDisabled(bool)));

  if(dialog.exec() != QDialog::Accepted)
    return;

  QCursor tmp_cursor(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(tmp_cursor);

  typedef Scene_surface_mesh_item::Face_graph Face_graph;
  Face_graph& tm = *sm_item->face_graph();
  typedef boost::graph_traits<Face_graph>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<Face_graph>::vertex_descriptor vertex_descriptor;

  CGAL::Polygon_mesh_processing::stitch_borders(tm);
#if 1
/// detection of non-manifold parts
  std::map< std::pair<Kernel::Point_3, Kernel::Point_3>, std::vector<halfedge_descriptor> > edges;
  for(halfedge_descriptor h : halfedges(tm))
  {
    if (is_border(h,tm))
      edges[CGAL::make_sorted_pair(tm.point(target(h,tm)), tm.point(source(h,tm)))].push_back(h);
  }

  std::vector<int> fccs(num_faces(tm),-1);
  int nbcc = CGAL::Polygon_mesh_processing::connected_components(tm, CGAL::make_property_map(fccs));
  //this has to be done per cycle so as to keep 2 patches
  // remove the smallest CCs
  std::vector<int> cc_sizes(nbcc, 0);
  for(int i : fccs)
    cc_sizes[i]+=1;
  std::set<int> ccs_to_rm;
  for (auto p : edges)
    if (p.second.size() >= 2)
      for(halfedge_descriptor h : p.second)
      {
        int ccid = fccs[face(opposite(h, tm),tm)];
        if ( cc_sizes[ccid]<=4 )
          ccs_to_rm.insert(ccid);
      }
  std::cout << "removing " << ccs_to_rm.size() << " ccs\n";
  CGAL::Polygon_mesh_processing::remove_connected_components(tm, ccs_to_rm, CGAL::make_property_map(fccs));
  std::cout << "input is valid after cc removal:"<< CGAL::is_valid_polygon_mesh(tm) << "\n";
///
#endif

  if (ui.use_local_tolerance->isChecked())
  {
    CGAL::Polygon_mesh_processing::experimental::snap_borders(tm, CGAL::parameters::do_simplify_border(ui.do_simplify_border->isChecked()));
    CGAL::Polygon_mesh_processing::stitch_borders(tm);
    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(tm);
  }
  else
  {
    std::vector<double> tolerances/* = 0.005, 0.0125, 0.025, 0.05, 0.07 */;
    bool ok;
    Q_FOREACH(QString tol_text, ui.tolerances->text().split(","))
    {
      double d = tol_text.toDouble(&ok);
      if (ok)
        tolerances.push_back(d);
      else
        QMessageBox(QMessageBox::Warning,
          QString("Invalid value"),
          QString("\""+tol_text+"\" is not a valid double, ignored."),
          QMessageBox::Ok,
          this->mw).exec();
    }

    for (double tol : tolerances )
    {
      std::cout << "using tol = " << tol << "\n";
      CGAL::Constant_property_map<vertex_descriptor, double> tolerance_map(tol);
      CGAL::Polygon_mesh_processing::experimental::snap_borders(tm, tolerance_map, CGAL::parameters::do_simplify_border(ui.do_simplify_border->isChecked()));

      CGAL::Polygon_mesh_processing::stitch_borders(tm);
      CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(tm);

      // post processing
      std::vector<halfedge_descriptor> remaining_cycles;
      CGAL::Polygon_mesh_processing::extract_boundary_cycles(tm, std::back_inserter(remaining_cycles));

      for (halfedge_descriptor hc : remaining_cycles)
      {
        if (next(next(hc,tm),tm)==prev(hc,tm))
        {
          //get smallest halfedge
          halfedge_descriptor hm = hc;
          double min_l = CGAL::Polygon_mesh_processing::edge_length(hc, tm);

          double el = CGAL::Polygon_mesh_processing::edge_length(next(hc, tm), tm);
          if (el<min_l)
          {
            min_l=el;
            hm=next(hc, tm);
          }

          el = CGAL::Polygon_mesh_processing::edge_length(prev(hc, tm), tm);
          if (el<min_l)
          {
            min_l=el;
            hm=prev(hc, tm);
          }
          if (el>tol)
            continue;
          if (!CGAL::Euler::does_satisfy_link_condition(edge(hm, tm), tm))
          {
            // simply fill the face
            std::array<vertex_descriptor,3> vr = { source(hm, tm), target(hm, tm), target(next(hm, tm), tm) };
            CGAL::Euler::add_face(vr, tm);
            continue;
          }

          std::array<vertex_descriptor,3> vr = { source(hm, tm), target(hm, tm), target(next(hm, tm), tm) };
          CGAL::Euler::add_face(vr, tm);
          CGAL::Euler::collapse_edge(edge(hm, tm), tm);
        }
      }
    }
  }
  CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(tm);

  sm_item->invalidateOpenGLBuffers();
  sm_item->itemChanged();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveDegenerateFaces_triggered(Scene_interface::Item_id index)
{
  Item* poly_item =
    qobject_cast<Item*>(scene->item(index));
  if (poly_item)
  {
    std::size_t nbv = faces(*poly_item->polyhedron()).size();
      CGAL::Polygon_mesh_processing::remove_degenerate_faces(
      *poly_item->polyhedron());
    nbv -= faces(*poly_item->polyhedron()).size();
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
    CGAL::Three::Three::information(tr(" %1 degenerate faces have been removed.")
                          .arg(nbv));
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveDegenerateFaces_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
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
      CGAL::Polygon_mesh_processing::experimental::remove_self_intersections(
      *poly_item->polyhedron(), CGAL::parameters::preserve_genus(false));
    if (!solved)
      CGAL::Three::Three::information(tr("Some self-intersection could not be fixed"));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveSelfIntersections_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
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
    catch(const CGAL::Polygon_mesh_processing::Corefinement::Triple_intersection_exception&)
    {
      CGAL::Three::Three::warning(tr("The result of the requested operation is not handled (triple intersection)."));
    }
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefine_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
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
        CGAL::Three::Three::information(tr("Self-intersection could not be removed due to non-manifold edges in the output"));
    }
    catch(const CGAL::Polygon_mesh_processing::Corefinement::Triple_intersection_exception&)
    {
      CGAL::Three::Three::warning(tr("The result of the requested operation is not handled (triple intersection)."));
    }
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionAutorefineAndRMSelfIntersections_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
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
    CGAL::Three::Three::information(tr(" %1 pairs of halfedges stitched.").arg(halfedges_to_stitch.size()));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionStitchCloseBorderHalfedges_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
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
    CGAL::Three::Three::information(tr(" %1 vertices created").arg(nb_vd));
    if (nb_vd)
    {
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
  }
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionExtractNMVertices_triggered(Scene_interface::Item_id index)
{
  namespace PMP =   CGAL::Polygon_mesh_processing;
  typedef typename boost::graph_traits<typename Item::Face_graph>::halfedge_descriptor halfedge_descriptor;

  if (Item* poly_item = qobject_cast<Item*>(scene->item(index)))
  {
    std::vector<halfedge_descriptor> hds;
    PMP::non_manifold_vertices(*poly_item->face_graph(), std::back_inserter(hds));
    if(hds.empty())
    {
      CGAL::Three::Three::information(tr(" No NM vertex found."));
      return;
    }
    Scene_points_with_normal_item* pitem = new Scene_points_with_normal_item();
    pitem->setColor(Qt::red);
    pitem->setName(QString("%1 nm vertices").arg(poly_item->name()));
    for(const auto& h : hds)
    {
      pitem->point_set()->insert(poly_item->face_graph()->point(target(h, *poly_item->face_graph())));
    }
    scene->addItem (pitem);
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionDuplicateNMVertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionDuplicateNMVertices_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionExtractNMVertices_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionExtractNMVertices_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

template <typename Item>
void Polyhedron_demo_repair_polyhedron_plugin::on_actionMergeDuplicatedVerticesOnBoundaryCycles_triggered(Scene_interface::Item_id index)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(Item* poly_item = qobject_cast<Item*>(scene->item(index)))
  {
    const std::size_t old_nv = vertices(*poly_item->polyhedron()).size();
    PMP::merge_duplicated_vertices_in_boundary_cycles(*poly_item->polyhedron());
    const std::size_t new_nv = vertices(*poly_item->polyhedron()).size();
    CGAL::Three::Three::information(tr(" %1 vertices merged").arg(old_nv - new_nv));
    if(old_nv != new_nv)
    {
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionMergeDuplicatedVerticesOnBoundaryCycles_triggered()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  on_actionMergeDuplicatedVerticesOnBoundaryCycles_triggered<Scene_surface_mesh_item>(index);
  QApplication::restoreOverrideCursor();
}

#include "Repair_polyhedron_plugin.moc"
