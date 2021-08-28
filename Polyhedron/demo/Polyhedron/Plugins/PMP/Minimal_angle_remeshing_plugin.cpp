//#define CGAL_PMP_REMESHING_VERBOSE
//#define CGAL_PMP_REMESHING_DEBUG
//#define CGAL_PMP_REMESHING_VERY_VERBOSE
//#define CGAL_PMP_REMESHING_VERBOSE_PROGRESS

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"

#include "Scene_polyhedron_selection_item.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/minimal_angle_remeshing.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/utility.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/property_map.h>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>
#include <cmath>

#include "ui_Minimal_angle_remeshing_dialog.h"


typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

using namespace CGAL::Three;
class Polyhedron_demo_minimal_angle_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "minimal_angle_remeshing_plugin.json")

  typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

  struct Visitor
  {
    typedef typename Scene_polyhedron_selection_item::Selection_set_facet Container;
    Container& faces;

    Visitor(Container& container)
      : faces(container)
    {}

    void before_subface_creations(face_descriptor fd)
    {
      Container::iterator it = faces.find(fd);
      faces.erase(it);
    }
    void after_subface_created(face_descriptor fd)
    {
      faces.insert(fd);
    }
    void after_subface_creations(){}
  };

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionMinimalAngleRemeshing_ = new QAction("Minimal Angle Remeshing", mw);
    actionMinimalAngleRemeshing_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionMinimalAngleRemeshing_) {
      connect(actionMinimalAngleRemeshing_, SIGNAL(triggered()),
        this, SLOT(minimal_angle_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMinimalAngleRemeshing_;
  }

  bool applicable(QAction*) const
  {
    if (scene->selectionIndices().size() == 1)
    {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }

    bool ok(true), found_poly(false);

    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (!qobject_cast<Scene_facegraph_item*>(scene->item(index)))
        ok = false;
      else
        found_poly=true;
    }
    return ok && found_poly;
  }

public Q_SLOTS:
  void minimal_angle_remeshing()
  {
    if (scene->selectionIndices().size() > 1)
    {
      minimal_angle_remeshing_of_several_polyhedra();
      return;
    }

    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_facegraph_item* poly_item = qobject_cast<Scene_facegraph_item*>(scene->item(index));

    if (poly_item)
    {
      // get the named parameters
      NamedParameters np;
      bool suc = get_named_parameters(np);
      if (!suc) {
        std::cout << "Remeshing aborted" << std::endl;
      }
      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);
      // start timing
      QElapsedTimer time;
      time.start();
      // triangulate if necessary
      if (!CGAL::is_triangle_mesh(*poly_item->polyhedron()))
      {
        QApplication::restoreOverrideCursor();
        if (QMessageBox::Ok ==
          QMessageBox::question(mw, tr("Error - Triangulate Faces?"),
            tr("The input mesh is not a triangulated surface mesh.\n"
              "Do you wish to triangulate faces first, or cancel remeshing ?"),
            (QMessageBox::Ok | QMessageBox::Cancel), QMessageBox::Ok))
        {
          QApplication::setOverrideCursor(Qt::WaitCursor);
          CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
        }
        else
        {
          return;
        }
      }
      // perform minimal angle remeshing
      perform_minimal_angle_remeshing(*poly_item->polyhedron(), np);
      // change status
      poly_item->setItemIsMulticolor(false);
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
      // finish timing
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
      // restore to default cursor
      QApplication::restoreOverrideCursor();
    }
    else
    {
      std::cout << "Can't remesh that type of thing" << std::endl;
    }
  }

  void minimal_angle_remeshing_of_several_polyhedra()
  {
    // get the named parameters
    NamedParameters np;
    bool suc = get_named_parameters(np);
    if (!suc) {
      std::cout << "Remeshing aborted" << std::endl;
      return;
    }
    // collect the selections
    std::vector<Scene_facegraph_item*> selection;
    for(int index : scene->selectionIndices())
    {
      Scene_facegraph_item* poly_item = qobject_cast<Scene_facegraph_item*>(scene->item(index));
      if (poly_item == NULL)
      {
        std::cout << scene->item(index)->name().data()
                  << " is not a FaceGraph, remeshing skipped\n";
        continue;
      }
      else
      {
        selection.push_back(poly_item);
      }
    }
    //check non-triangulated surfaces
    for (Scene_facegraph_item* poly_item : selection)
    {
      if (!CGAL::is_triangle_mesh(*poly_item->polyhedron()))
      {
        if (QMessageBox::Ok == QMessageBox::question(mw,
              tr("Error - Triangulate Faces?"),
              tr("The input mesh ").append(poly_item->name())
               .append(tr(" is not a triangulated surface mesh.\n"
                "Do you wish to triangulate faces first, or cancel remeshing ?")),
              (QMessageBox::Ok | QMessageBox::Cancel),
              QMessageBox::Ok))
        {
          QApplication::setOverrideCursor(Qt::WaitCursor);
          CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
          QApplication::restoreOverrideCursor();
        }
        else
        {
          QApplication::restoreOverrideCursor();
          return;
        }
      }
    }
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    // perform minimal angle remeshing
    int total_time = 0;
    for(Scene_facegraph_item* poly_item : selection)
    {
      QElapsedTimer time;
      time.start();
      std::cout << "Minimal angle remeshing of "
                << poly_item->name().toStdString() << " started..." << std::endl;
      perform_minimal_angle_remeshing(*poly_item->polyhedron(), np);
      std::cout << "Minimal angle remeshing of "
                << poly_item->name().toStdString() << " done." << std::endl;
      total_time += time.elapsed();
      std::cout << "Remeshing of " << poly_item->name().data()
                << " done in " << time.elapsed() << " ms" << std::endl;
    }
    std::cout << "Remeshing of all selected items done in "<< total_time << " ms" << std::endl;
    // change status
    for(Scene_facegraph_item* poly_item : selection)
    {
      //destroys the patch_id_map for the Surface_mesh_item to avoid assertions.
      poly_item->resetColors();
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
    // default cursor
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;

  Ui::Minimal_angle_remeshing_dialog
  remeshing_dialog(QDialog* dialog)
  {
    Ui::Minimal_angle_remeshing_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));
    return ui;
  }

  bool get_named_parameters(NamedParameters& np) {
    // Create dialog box
    QDialog dialog(this->mw);
    Ui::Minimal_angle_remeshing_dialog ui = remeshing_dialog(&dialog);
    // Get values
    int i = dialog.exec();
    if (i == QDialog::Rejected)
    {
      return false;
    }
    // general parameters
    np.max_error_threshold = ui.sb_max_error_threshold->value();
    np.min_angle_threshold = ui.sb_min_angle_threshold->value();
    np.max_mesh_complexity = ui.sb_max_mesh_complexity->value();
    np.smooth_angle_delta = ui.sb_smooth_angle_delta->value();

    np.apply_edge_flip = ui.cb_apply_edge_flip->isChecked();
    np.edge_flip_strategy = ui.cb_edge_flip_strategy->currentIndex() == 0 ?
      EdgeFlipStrategy::k_improve_valence : EdgeFlipStrategy::k_improve_angle;
    np.flip_after_split_and_collapse = ui.cb_flip_after_split_and_collapse->isChecked();
    np.relocate_after_local_operations = ui.cb_relocate_after_local_operations->isChecked();

    np.relocate_strategy = ui.cb_relocate_strategy->currentIndex() == 0 ?
      RelocateStrategy::k_barycenter : RelocateStrategy::k_cvt_barycenter;
    np.keep_vertex_in_one_ring = ui.cb_keep_vertex_in_one_ring->isChecked();
    np.use_local_aabb_tree = ui.cb_use_local_aabb_tree->isChecked();
    np.collapsed_list_size = ui.sb_collapsed_list_size->value();

    np.decrease_max_errors = ui.cb_decrease_max_errors->isChecked();
    np.verbose_progress = ui.cb_verbose_progress->isChecked();
    np.apply_initial_mesh_simplification = ui.cb_apply_initial_mesh_simplification->isChecked();
    np.apply_final_vertex_relocation = ui.cb_apply_final_vertex_relocation->isChecked();

    // sample parameters
    np.samples_per_face_in = ui.sb_samples_per_face_in->value();
    np.samples_per_face_out = ui.sb_samples_per_face_out->value();
    np.max_samples_per_area = ui.sb_max_samples_per_area->value();
    np.min_samples_per_triangle = ui.sb_min_samples_per_triangle->value();

    np.bvd_iteration_count = ui.sb_bvd_iteration_count->value();
    np.sample_number_strategy = ui.cb_sample_number_strategy->currentIndex() == 0 ?
      SampleNumberStrategy::k_fixed : SampleNumberStrategy::k_variable;
    np.sample_strategy = ui.cb_sample_strategy->currentIndex() == 0 ?
      SampleStrategy::k_uniform : SampleStrategy::k_adaptive;
    np.use_stratified_sampling = ui.cb_use_stratified_sampling->isChecked();

    // feature intensity parameters
    np.sum_theta = ui.sb_sum_theta->value();
    np.sum_delta = ui.sb_sum_delta->value();
    np.dihedral_theta = ui.sb_dihedral_theta->value();
    np.dihedral_delta = ui.sb_dihedral_delta->value();

    np.feature_difference_delta = ui.sb_feature_difference_delta->value();
    np.feature_control_delta = ui.sb_feature_control_delta->value();
    np.inherit_element_types = ui.cb_inherit_element_types->isChecked();
    np.use_feature_intensity_weights = ui.cb_use_feature_intensity_weights->isChecked();

    // vertex optimization parameters
    np.vertex_optimize_count = ui.sb_vertex_optimize_count->value();
    np.vertex_optimize_ratio = ui.sb_vertex_optimize_ratio->value();
    np.stencil_ring_size = ui.sb_stencil_ring_size->value();
    np.optimize_strategy = ui.cb_optimize_strategy->currentIndex() == 0 ?
      OptimizeStrategy::k_approximation : OptimizeStrategy::k_interpolation;

    np.face_optimize_type = (OptimizeType)(ui.cb_face_optimize_type->currentIndex());
    np.edge_optimize_type = (OptimizeType)(ui.cb_edge_optimize_type->currentIndex());
    np.vertex_optimize_type = (OptimizeType)(ui.cb_vertex_optimize_type->currentIndex());
    np.optimize_after_local_operations = ui.cb_optimize_after_local_operations->isChecked();

    return true;
  }

  void perform_minimal_angle_remeshing(FaceGraph& pmesh, const NamedParameters& np) {
    CGAL::Polygon_mesh_processing::minimal_angle_remeshing(pmesh,
      CGAL::Polygon_mesh_processing::parameters::
      max_error_threshold(np.max_error_threshold)       // general parameters
      .min_angle_threshold(np.min_angle_threshold)
      .max_mesh_complexity(np.max_mesh_complexity)
      .smooth_angle_delta(np.smooth_angle_delta)
      .apply_edge_flip(np.apply_edge_flip)
      .edge_flip_strategy(np.edge_flip_strategy)
      .flip_after_split_and_collapse(np.flip_after_split_and_collapse)
      .relocate_after_local_operations(np.relocate_after_local_operations)
      .relocate_strategy(np.relocate_strategy)
      .keep_vertex_in_one_ring(np.keep_vertex_in_one_ring)
      .use_local_aabb_tree(np.use_local_aabb_tree)
      .collapsed_list_size(np.collapsed_list_size)
      .decrease_max_errors(np.decrease_max_errors)
      .verbose_progress(np.verbose_progress)
      .apply_initial_mesh_simplification(np.apply_initial_mesh_simplification)
      .apply_final_vertex_relocation(np.apply_final_vertex_relocation)
      .samples_per_face_in(np.samples_per_face_in)       // sample parameters
      .samples_per_face_out(np.samples_per_face_out)
      .max_samples_per_area(np.max_samples_per_area)
      .min_samples_per_triangle(np.min_samples_per_triangle)
      .bvd_iteration_count(np.bvd_iteration_count)
      .sample_number_strategy(np.sample_number_strategy)
      .sample_strategy(np.sample_strategy)
      .use_stratified_sampling(np.use_stratified_sampling)
      .sum_theta(np.sum_theta)                           // feature paramemters
      .sum_delta(np.sum_delta)
      .dihedral_theta(np.dihedral_theta)
      .dihedral_delta(np.dihedral_delta)
      .feature_difference_delta(np.feature_difference_delta)
      .feature_control_delta(np.feature_control_delta)
      .inherit_element_types(np.inherit_element_types)
      .use_feature_intensity_weights(np.use_feature_intensity_weights)
      .vertex_optimize_count(np.vertex_optimize_count)   // vertex optimization parameters
      .vertex_optimize_ratio(np.vertex_optimize_ratio)
      .stencil_ring_size(np.stencil_ring_size)
      .optimize_strategy(np.optimize_strategy)
      .face_optimize_type(np.face_optimize_type)
      .edge_optimize_type(np.edge_optimize_type)
      .vertex_optimize_type(np.vertex_optimize_type)
      .optimize_after_local_operations(np.optimize_after_local_operations));
  }

private:
  QAction* actionMinimalAngleRemeshing_;

}; // end Polyhedron_demo_minimal_angle_remeshing_plugin

#include "Minimal_angle_remeshing_plugin.moc"
