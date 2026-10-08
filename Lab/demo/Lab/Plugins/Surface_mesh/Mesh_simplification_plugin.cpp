#include "ui_Mesh_simplification_dialog.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QElapsedTimer>
#include <QAction>

#include <QMessageBox>
#include <QButtonGroup>

#include <optional>
#include <type_traits>

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef EPICK::FT FT;

namespace SMS = ::CGAL::Surface_mesh_simplification;

using namespace CGAL::Three;

class CGAL_Lab_mesh_simplification_plugin
  : public QObject,
    public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

public:
  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    mw = mainWindow;
    scene = scene_interface;
    QAction *actionSimplify = new QAction("Simplification", mw);
    actionSimplify->setProperty("subMenuName",
                                "Triangulated Surface Mesh Simplification");
    connect(actionSimplify, SIGNAL(triggered()),
            this, SLOT(on_actionSimplify_triggered()));
    _actions <<actionSimplify;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSimplify_triggered();

private :
  Scene_interface *scene;
  QMainWindow *mw;
  QList<QAction*> _actions;
}; // end CGAL_Lab_mesh_simplification_plugin

void CGAL_Lab_mesh_simplification_plugin::on_actionSimplify_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_selection_item* selection_item =
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
  if(selection_item && selection_item->selected_edges.empty())
  {
    QMessageBox::warning(mw, "Empty Edges", "There are no selected edges. Aborting.");
    return;
  }

  Scene_facegraph_item* poly_item =
    qobject_cast<Scene_facegraph_item*>(scene->item(index));

  if (poly_item || selection_item)
  {
    FaceGraph& pmesh = (poly_item != nullptr) ? *poly_item->polyhedron()
                                              : *selection_item->polyhedron();

    // get option
    QDialog dialog(mw);
    Ui::Mesh_simplification_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));

    // Set up button groups for exclusive selection
    QButtonGroup *strategyGroup = new QButtonGroup(&dialog);
    strategyGroup->addButton(ui.radio_GH);
    strategyGroup->addButton(ui.radio_LT);
    strategyGroup->addButton(ui.radio_EL);
    ui.radio_GH->setChecked(true);

    QButtonGroup *stopGroup = new QButtonGroup(&dialog);
    stopGroup->addButton(ui.radio_EdgeLength);
    stopGroup->addButton(ui.radio_EdgeCount);
    stopGroup->addButton(ui.radio_FacetCount);
    ui.radio_EdgeCount->setChecked(true);

    ui.useBoundedPlacement->setChecked(true);

    Scene_interface::Bbox bbox = poly_item != nullptr ? poly_item->bbox()
                                                      : (selection_item != nullptr ? selection_item->bbox()
                                                                                   : scene->bbox());

    double diag_length = CGAL::sqrt((bbox.xmax() - bbox.xmin())*(bbox.xmax() - bbox.xmin()) +
                                    (bbox.ymax() - bbox.ymin())*(bbox.ymax() - bbox.ymin()) +
                                    (bbox.zmax() - bbox.zmin())*(bbox.zmax() - bbox.zmin()));

    ui.spin_EdgeCount->setMinimum(0);
    ui.spin_EdgeCount->setMaximum(static_cast<int>(num_edges(pmesh)));
    ui.spin_EdgeCount->setValue(static_cast<int>(num_edges(pmesh)/2));
    ui.spin_FacetCount->setMinimum(0);
    ui.spin_FacetCount->setMaximum(static_cast<int>(num_faces(pmesh)));
    ui.spin_FacetCount->setValue(static_cast<int>(num_faces(pmesh)/2));
    ui.spin_EdgeLength->setValue(diag_length * 0.05);

    ui.spin_EdgeLength->setEnabled(ui.radio_EdgeLength->isChecked());
    ui.spin_EdgeCount->setEnabled(ui.radio_EdgeCount->isChecked());
    ui.spin_FacetCount->setEnabled(ui.radio_FacetCount->isChecked());
    connect(ui.radio_EdgeLength, &QRadioButton::toggled,
            ui.spin_EdgeLength, &QWidget::setEnabled);
    connect(ui.radio_EdgeCount, &QRadioButton::toggled,
            ui.spin_EdgeCount, &QWidget::setEnabled);
    connect(ui.radio_FacetCount, &QRadioButton::toggled,
            ui.spin_FacetCount, &QWidget::setEnabled);

    // Disable Edge Length stop for GH and LT strategies
    ui.radio_EdgeLength->setEnabled(false); // GH is default
    connect(ui.radio_GH, &QRadioButton::toggled,
            [&ui](bool checked) {
              ui.radio_EdgeLength->setEnabled(!checked);
              if (checked &&  ui.radio_EdgeLength->isChecked()) {
                  ui.radio_EdgeCount->setChecked(true);
              }
            });
    connect(ui.radio_LT, &QRadioButton::toggled,
            [&ui](bool checked) {
              ui.radio_EdgeLength->setEnabled(!checked);
              if (checked && ui.radio_EdgeLength->isChecked()) {
                  ui.radio_EdgeCount->setChecked(true);
              }
            });
    connect(ui.radio_EL, &QRadioButton::toggled,
            [&ui](bool checked) {
              ui.radio_EdgeLength->setEnabled(checked);
              if (checked && !ui.radio_EdgeLength->isChecked() &&
                             !ui.radio_EdgeCount->isChecked() &&
                             !ui.radio_FacetCount->isChecked()) {
                  ui.radio_EdgeLength->setChecked(true);
              }
            });

    // check user cancellation
    if(dialog.exec() == QDialog::Rejected)
      return;

    // simplify
    QElapsedTimer time;
    time.start();
    std::cout << "Simplify...";

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();

    auto call_edge_collapse = [&](const auto& stop, const auto& cost, const auto& base_placement)
    {
      using Base_placement = std::decay_t<decltype(base_placement)>;

      if(selection_item)
      {
        SMS::Constrained_placement<
          Base_placement,
          Scene_polyhedron_selection_item::Is_constrained_map<
            Scene_polyhedron_selection_item::Selection_set_edge> > placement(selection_item->constrained_edges_pmap(), base_placement);

        if(ui.useBoundedPlacement->isChecked())
        {
          SMS::Bounded_normal_change_filter<> filter;
          SMS::edge_collapse(pmesh, stop,
                            CGAL::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
                                             .get_cost(cost)
                                             .filter(filter)
                                             .get_placement(placement));
        }
        else
        {
          SMS::edge_collapse(pmesh, stop,
                            CGAL::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
                                             .get_cost(cost)
                                             .get_placement(placement));
        }
      }
      else
      {
        if(ui.useBoundedPlacement->isChecked())
        {
          SMS::Bounded_normal_change_filter<> filter;
          SMS::edge_collapse(pmesh, stop,
                            CGAL::parameters::vertex_index_map(get(boost::vertex_index, pmesh))
                                             .get_cost(cost)
                                             .filter(filter)
                                             .get_placement(base_placement));
        }
        else
        {
          SMS::edge_collapse(pmesh, stop,
                            CGAL::parameters::vertex_index_map(get(boost::vertex_index, pmesh))
                                             .get_cost(cost)
                                             .get_placement(base_placement));
        }
      }
    };

    if (ui.radio_LT->isChecked())
    {
      if (ui.radio_EdgeCount->isChecked())
        call_edge_collapse(SMS::Edge_count_stop_predicate<FaceGraph>(ui.spin_EdgeCount->value()),
                           SMS::LindstromTurk_cost<FaceGraph>(),
                           SMS::LindstromTurk_placement<FaceGraph>());
      else
        call_edge_collapse(SMS::Face_count_stop_predicate<FaceGraph>(ui.spin_FacetCount->value()),
                           SMS::LindstromTurk_cost<FaceGraph>(),
                           SMS::LindstromTurk_placement<FaceGraph>());
    }
    else if (ui.radio_GH->isChecked())
    {
      SMS::GarlandHeckbert_policies<FaceGraph, EPICK> policies(pmesh);
      if (ui.radio_EdgeCount->isChecked())
        call_edge_collapse(SMS::Edge_count_stop_predicate<FaceGraph>(ui.spin_EdgeCount->value()),
                           policies, policies);
      else
        call_edge_collapse(SMS::Face_count_stop_predicate<FaceGraph>(ui.spin_FacetCount->value()),
                           policies, policies);
    }
    else if (ui.radio_EL->isChecked())
    {
      if (ui.radio_EdgeLength->isChecked())
        call_edge_collapse(SMS::Edge_length_stop_predicate<FT>(ui.spin_EdgeLength->value()),
                           SMS::Edge_length_cost<FaceGraph>(),
                           SMS::Midpoint_placement<FaceGraph>());
      else if (ui.radio_EdgeCount->isChecked())
        call_edge_collapse(SMS::Edge_count_stop_predicate<FaceGraph>(ui.spin_EdgeCount->value()),
                           SMS::Edge_length_cost<FaceGraph>(),
                           SMS::Midpoint_placement<FaceGraph>());
      else
        call_edge_collapse(SMS::Face_count_stop_predicate<FaceGraph>(ui.spin_FacetCount->value()),
                           SMS::Edge_length_cost<FaceGraph>(),
                           SMS::Midpoint_placement<FaceGraph>());
    }
    else
    {
      std::cerr << "Unknown simplification strategy selected." << std::endl;
      return;
    }

    std::cout << "ok (" << time.elapsed() << " ms, " << num_halfedges(pmesh) / 2 << " edges)" << std::endl;

    // update scene
    if (poly_item != nullptr)
    {
      poly_item->invalidateOpenGLBuffers();
      poly_item->polyhedron()->collect_garbage();
    }
    else
    {
      selection_item->polyhedron_item()->polyhedron()->collect_garbage();
      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
      selection_item->invalidateOpenGLBuffers();
    }

    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
  }
}

#include "Mesh_simplification_plugin.moc"
