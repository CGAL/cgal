#include "main_window.h"
// Qt
#include <QFileDialog>
#include <qsettings.h>
// local
#include "ui_main_window.h"
#include "ui_parameter_settings.h"
#include "parameter_settings.h"

MainWindow::MainWindow(QWidget *parent)
  : CGAL::Qt::DemosMainWindow(parent) {
  ui = new Ui::MainWindow;
  ui->setupUi(this);

  // saves some pointers from ui, for latter use.
  m_pViewer = ui->viewer;

  // does not save the state of the viewer 
  m_pViewer->setStateFileName(QString::null);

  // accepts drop events
  setAcceptDrops(true);
  // setups scene
  m_pScene = new Scene();
  m_pViewer->setScene(m_pScene);
  m_pViewer->setManipulatedFrame(m_pScene->manipulatedFrame());

  // connects actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionFile_quit, SIGNAL(triggered()), this, SLOT(quit()));
  this->addRecentFiles(ui->menuFile, ui->actionFile_quit);
  connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));

  readSettings();
  update_menu_items();
}

MainWindow::~MainWindow() {
  delete ui;
}

void MainWindow::updateViewerBBox() {
  //m_pScene->update_bbox();
  const Bbox &bbox = m_pScene->bbox();
  const double xmin = bbox.xmin();
  const double ymin = bbox.ymin();
  const double zmin = bbox.zmin();
  const double xmax = bbox.xmax();
  const double ymax = bbox.ymax();
  const double zmax = bbox.zmax();
  CGAL::qglviewer::Vec vec_min(xmin, ymin, zmin), vec_max(xmax, ymax, zmax);
  m_pViewer->setSceneBoundingBox(vec_min, vec_max);
  m_pViewer->camera()->showEntireScene();
}

void MainWindow::open(QString file_name) {
  open(file_name, OpenType::k_open_both);
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m) {
  m_pViewer->setAddKeyFrameKeyboardModifiers(m);
}

void MainWindow::quit() {
  writeSettings();
  close();
}

void MainWindow::readSettings() {
  this->readState("MainWindow", Size | State);
}

void MainWindow::writeSettings() {
  this->writeState("MainWindow");
  std::cerr << "Write setting... done.\n";
}

void MainWindow::dropEvent(QDropEvent *event) {
  Q_FOREACH(QUrl url, event->mimeData()->urls()) {
    QString filename = url.toLocalFile();
    if (!filename.isEmpty()) {
      QTextStream(stderr) << QString("dropEvent(\"%1\")\n").arg(filename);
      open(filename);
    }
  }
  event->acceptProposedAction();
}

void MainWindow::closeEvent(QCloseEvent *event) {
  writeSettings();
  event->accept();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event) {
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

// File menu:
void MainWindow::on_actionFile_open_triggered() {
  QSettings settings;
  QString directory = settings.value("OFF open directory",
    QDir::current().dirName()).toString();
  QStringList filenames =
    QFileDialog::getOpenFileNames(this,
    tr("Load polyhedron surface..."),
    directory,
    tr("OFF files (*.off)\n"
    "All files (*)"));
  if (!filenames.isEmpty()) {
    Q_FOREACH(QString filename, filenames) {
      open(filename, OpenType::k_open_both);
    }
  }
}

void MainWindow::on_actionFile_open_input_triggered() {
  QSettings settings;
  QString directory = settings.value("OFF open directory",
    QDir::current().dirName()).toString();
  QStringList filenames =
    QFileDialog::getOpenFileNames(this,
    tr("Load polyhedron surface..."),
    directory,
    tr("OFF files (*.off)\n"
    "All files (*)"));
  if (!filenames.isEmpty()) {
    Q_FOREACH(QString filename, filenames) {
      open(filename, OpenType::k_open_input);
    }
  }
}

void MainWindow::on_actionFile_open_remesh_triggered() {
  QSettings settings;
  QString directory = settings.value("OFF open directory",
    QDir::current().dirName()).toString();
  QStringList filenames =
    QFileDialog::getOpenFileNames(this,
    tr("Load polyhedron surface..."),
    directory,
    tr("OFF files (*.off)\n"
    "All files (*)"));
  if (!filenames.isEmpty()) {
    Q_FOREACH(QString filename, filenames) {
      open(filename, OpenType::k_open_remesh);
    }
  }
}

void MainWindow::on_actionFile_save_remesh_as_triggered() {
  QSettings settings;
  QString directory = settings.value("OFF open directory",
    QDir::current().dirName()).toString();

  QString filters("Off files (*.off);;Mesh files (*.mesh);;All files (*.*)");
  QString defaultFilter("Off files (*.off)");

  QString filename =
    QFileDialog::getSaveFileName(this,
    tr("Save polyhedral surface as..."),
    directory,
    filters, &defaultFilter);
  if (!filename.isEmpty()) {
    m_pScene->save_remesh_as(filename);
  }
}

// Edit menu:
void MainWindow::on_actionEdit_copy_snapshot_triggered() {
  // copy snapshot to clipboard
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QClipboard *qb = QApplication::clipboard();
  m_pViewer->makeCurrent();
  m_pViewer->raise();
  QImage snapshot = m_pViewer->grabFramebuffer();
  qb->setImage(snapshot);
  QApplication::restoreOverrideCursor();
}

void MainWindow::on_actionEdit_save_snapshot_triggered() {
  return;
  // save snapshot to file
  /*QApplication::setOverrideCursor(Qt::WaitCursor);
  QString filename = QFileDialog::getSaveFileName(this,
    tr("Save snapshot to file..."), "snapshot00.png", "*.png");
  m_pViewer->saveSnapshot(filename);
  QApplication::restoreOverrideCursor();*/
}

void MainWindow::on_actionEdit_parameter_settings_triggered() {
  ParameterSettings ps;
  // general parameters
  ps.set_max_error_threshold(m_pScene->get_max_error_threshold());        
  ps.set_min_angle_threshold(m_pScene->get_min_angle_threshold());
  ps.set_max_mesh_complexity(m_pScene->get_max_mesh_complexity());
  ps.set_smooth_angle_delta(m_pScene->get_smooth_angle_delta());
  ps.set_apply_edge_flip(m_pScene->get_apply_edge_flip());
  EdgeFlipStrategy efs = m_pScene->get_edge_flip_strategy();
  int current_index = efs == EdgeFlipStrategy::k_improve_valence ? 0 : 1;
  ps.set_edge_flip_strategy(current_index);
  ps.set_flip_after_split_and_collapse(
      m_pScene->get_flip_after_split_and_collapse());
  ps.set_relocate_after_local_operations(
      m_pScene->get_relocate_after_local_operations());
  RelocateStrategy rs = m_pScene->get_relocate_strategy();
  current_index = rs == RelocateStrategy::k_barycenter ? 0 : 1;
  ps.set_relocate_strategy(current_index);
  ps.set_keep_vertex_in_one_ring(m_pScene->get_keep_vertex_in_one_ring());
  ps.set_use_local_aabb_tree(m_pScene->get_use_local_aabb_tree());
  ps.set_collapsed_list_size(m_pScene->get_collapsed_list_size());
  ps.set_decrease_max_errors(m_pScene->get_decrease_max_errors());
  ps.set_track_information(m_pScene->get_track_information());
  ps.set_apply_initial_mesh_simplification(
      m_pScene->get_apply_initial_mesh_simplification());
  ps.set_apply_final_vertex_relocation(
      m_pScene->get_apply_final_vertex_relocation());
  // sample parameters
  ps.set_samples_per_facet_in(m_pScene->get_samples_per_facet_in());
  ps.set_samples_per_facet_out(m_pScene->get_samples_per_facet_out());
  ps.set_max_samples_per_area(m_pScene->get_max_samples_per_area());
  ps.set_min_samples_per_triangle(m_pScene->get_min_samples_per_triangle());
  ps.set_bvd_iteration_count(m_pScene->get_bvd_iteration_count());
  SampleNumberStrategy sns = m_pScene->get_sample_number_strategy();
  current_index = sns == SampleNumberStrategy::k_fixed ? 0 : 1;
  ps.set_sample_number_strategy(current_index);
  SampleStrategy ss = m_pScene->get_sample_strategy();
  current_index = ss == SampleStrategy::k_uniform ? 0 : 1;
  ps.set_sample_strategy(current_index);
  ps.set_use_stratified_sampling(m_pScene->get_use_stratified_sampling());
  // feature function parameters
  ps.set_sum_theta(m_pScene->get_sum_theta());
  ps.set_sum_delta(m_pScene->get_sum_delta());
  ps.set_dihedral_theta(m_pScene->get_dihedral_theta());
  ps.set_diheral_delta(m_pScene->get_dihedral_delta());
  ps.set_feature_difference_delta(m_pScene->get_feature_difference_delta());
  ps.set_feature_control_delta(m_pScene->get_feature_control_delta());
  ps.set_inherit_element_types(m_pScene->get_inherit_element_types());
  ps.set_use_feature_intensity_weights(
      m_pScene->get_use_feature_intensity_weights());
  // vertex optimization parameters
  ps.set_vertex_optimize_count(m_pScene->get_vertex_optimize_count());
  ps.set_vertex_optimize_ratio(m_pScene->get_vertex_optimize_ratio());
  ps.set_stencil_ring_size(m_pScene->get_stencil_ring_size());
  OptimizeStrategy os = m_pScene->get_optimize_strategy();
  current_index = os == OptimizeStrategy::k_approximation ? 0 : 1;
  ps.set_optimize_strategy(current_index);
  OptimizeType ot = m_pScene->get_facet_optimize_type(); 
  current_index = m_pScene->get_optimize_type_index(ot);
  ps.set_facet_optimize_type(current_index);
  ot = m_pScene->get_edge_optimize_type();                
  current_index = m_pScene->get_optimize_type_index(ot);
  ps.set_edge_optimize_type(current_index);
  ot = m_pScene->get_vertex_optimize_type();             
  current_index = m_pScene->get_optimize_type_index(ot);
  ps.set_vertex_optimize_type(current_index);
  ps.set_optimize_after_local_operations(
      m_pScene->get_optimize_after_local_operations());
  if (ps.exec() == QDialog::Accepted) {
    // general parameters
    m_pScene->set_max_error_threshold(ps.get_max_error_threshold());
    m_pScene->set_min_angle_threshold(ps.get_min_angle_threshold());
    m_pScene->set_max_mesh_complexity(ps.get_max_mesh_complexity());
    m_pScene->set_smooth_angle_delta(ps.get_smooth_angle_delta());
    m_pScene->set_apply_edge_flip(ps.get_apply_edge_flip());
    efs = ps.get_edge_flip_strategy() == 0 ? EdgeFlipStrategy::k_improve_valence : 
                                             EdgeFlipStrategy::k_improve_angle;
    m_pScene->set_edge_flip_strategy(efs);
    m_pScene->set_flip_after_split_and_collapse(
        ps.get_flip_after_split_and_collapse());
    m_pScene->set_relocate_after_local_operations(
        ps.get_relocate_after_local_operations());
    rs = ps.get_relocate_strategy() == 0 ? RelocateStrategy::k_barycenter : 
                                           RelocateStrategy::k_cvt_barycenter;
    m_pScene->set_relocate_strategy(rs);
    m_pScene->set_keep_vertex_in_one_ring(ps.get_keep_vertex_in_one_ring());
    m_pScene->set_use_local_aabb_tree(ps.get_use_local_aabb_tree());
    m_pScene->set_collapsed_list_size(ps.get_collapsed_list_size());
    m_pScene->set_decrease_max_errors(ps.get_decrease_max_errors());
    m_pScene->set_track_information(ps.get_track_information());
    m_pScene->set_apply_initial_mesh_simplification(
        ps.get_apply_initial_mesh_simplification());
    m_pScene->set_apply_final_vertex_relocation(
        ps.get_apply_final_vertex_relocation());
    // sample parameters
    sns = ps.get_sample_number_strategy() == 0 ? 
        SampleNumberStrategy::k_fixed : SampleNumberStrategy::k_variable;
    ss = ps.get_sample_strategy() == 0 ? 
        SampleStrategy::k_uniform : SampleStrategy::k_adaptive;
    bool samples_changed =
      (m_pScene->get_samples_per_facet_in() != ps.get_samples_per_facet_in()) ||
      (m_pScene->get_samples_per_facet_out() != ps.get_samples_per_facet_out()) ||
      (m_pScene->get_max_samples_per_area() != ps.get_max_samples_per_area()) ||
      (m_pScene->get_min_samples_per_triangle() != ps.get_min_samples_per_triangle()) ||
      (m_pScene->get_bvd_iteration_count() != ps.get_bvd_iteration_count()) ||
      (m_pScene->get_sample_number_strategy() != sns) ||
      (m_pScene->get_sample_strategy() != ss) ||
      (m_pScene->get_use_stratified_sampling() != ps.get_use_stratified_sampling());
    if (samples_changed) {
      m_pScene->set_samples_per_facet_in(ps.get_samples_per_facet_in());
      m_pScene->set_samples_per_facet_out(ps.get_samples_per_facet_out());
      m_pScene->set_max_samples_per_area(ps.get_max_samples_per_area());
      m_pScene->set_min_samples_per_triangle(ps.get_min_samples_per_triangle());
      m_pScene->set_bvd_iteration_count(ps.get_bvd_iteration_count());
      m_pScene->set_sample_number_strategy(sns);
      m_pScene->set_sample_strategy(ss);
      m_pScene->set_use_stratified_sampling(ps.get_use_stratified_sampling());
    }
    // feature function parameters
    bool feature_changed = 
      (m_pScene->get_sum_theta() != ps.get_sum_theta()) ||
      (m_pScene->get_sum_delta() != ps.get_sum_delta()) ||
      (m_pScene->get_dihedral_theta() != ps.get_dihedral_theta()) ||
      (m_pScene->get_dihedral_delta() != ps.get_dihedral_delta()) ||
      (m_pScene->get_feature_difference_delta() != ps.get_feature_difference_delta()) ||
      (m_pScene->get_feature_control_delta() != ps.get_feature_control_delta()) ||
      (m_pScene->get_inherit_element_types() != ps.get_inherit_element_types()) ||
      (m_pScene->get_use_feature_intensity_weights() != ps.get_use_feature_intensity_weights());
    if (feature_changed) {
      m_pScene->set_sum_theta(ps.get_sum_theta());
      m_pScene->set_sum_delta(ps.get_sum_delta());
      m_pScene->set_dihedral_theta(ps.get_dihedral_theta());
      m_pScene->set_dihedral_delta(ps.get_dihedral_delta());
      m_pScene->set_feature_difference_delta(ps.get_feature_difference_delta());
      m_pScene->set_feature_control_delta(ps.get_feature_control_delta());
      m_pScene->set_inherit_element_types(ps.get_inherit_element_types());
      m_pScene->set_use_feature_intensity_weights(ps.get_use_feature_intensity_weights());
    }
    if (samples_changed || feature_changed) {
      m_pScene->update_feature_intensities_and_clear_links();
      m_pScene->toggle_view_polyhedron_facets();
      m_pViewer->update();
      update_menu_items();
    }
    // vertex optimization parameters
    m_pScene->set_vertex_optimize_count(ps.get_vertex_optimize_count());
    m_pScene->set_vertex_optimize_ratio(ps.get_vertex_optimize_ratio());
    m_pScene->set_stencil_ring_size(ps.get_stencil_ring_size());
    os = ps.get_optimize_strategy() == 0 ? OptimizeStrategy::k_approximation :
                                           OptimizeStrategy::k_Interpolation;
    m_pScene->set_optimize_strategy(os);
    ot = m_pScene->get_optimize_type(ps.get_facet_optimize_type());
    m_pScene->set_facet_optimize_type(ot);
    ot = m_pScene->get_optimize_type(ps.get_edge_optimize_type());
    m_pScene->set_edge_optimize_type(ot);
    ot = m_pScene->get_optimize_type(ps.get_vertex_optimize_type());
    m_pScene->set_vertex_optimize_type(ot);
    m_pScene->set_optimize_after_local_operations(
        ps.get_optimize_after_local_operations());
  }
}

// Input
void MainWindow::on_actionInput_eliminate_degenerations_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->eliminate_degenerations();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionInput_split_long_edges_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->split_input_long_edges();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

void MainWindow::on_actionInput_properties_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->input_properties();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
}

// remesh
void MainWindow::on_actionRemesh_reset_from_input_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->reset_from_input();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionRemesh_generate_links_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->generate_links_and_types();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionRemesh_properties_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->remesh_properties();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

// isotropic
void MainWindow::on_actionIsotropic_remeshing_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->isotropic_remeshing();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionIsotropic_initial_mesh_simplification_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->initial_mesh_simplification();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionIsotropic_split_local_longest_edge_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->split_local_longest_edge();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionIsotropic_increase_minimal_angle_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->increase_minimal_angle();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionIsotropic_maximize_minimal_angle_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->maximize_minimal_angle();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

void MainWindow::on_actionIsotropic_final_Vertex_relocation_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  m_pScene->final_vertex_relocation();
  QApplication::restoreOverrideCursor();
  m_pViewer->update();
  update_menu_items();
}

// view
void MainWindow::on_actionView_input_triggered() {
  m_pScene->toggle_view_input();
  m_pViewer->update();
}

void MainWindow::on_actionView_remesh_triggered() {
  m_pScene->toggle_view_remesh();
  m_pViewer->update();
}

void MainWindow::on_actionView_toggle_input_remesh_triggered() {
  m_pScene->toggle_view_input_remesh();
  m_pViewer->update();
}

void MainWindow::on_actionView_polyhedron_edges_triggered() {
  m_pScene->toggle_view_polyhedron_edges();
  m_pViewer->update();
}

void MainWindow::on_actionView_minimal_angle_triggered() {
  m_pScene->toggle_view_minimal_angle();
  m_pViewer->update();
}

void MainWindow::on_actionView_polyhedron_facets_plain_triggered() {
  m_pScene->toggle_view_polyhedron_facets();
  m_pViewer->update();
}

void MainWindow::on_actionView_polyhedron_facets_ifi_triggered() {
  m_pScene->toggle_view_interpolated_feature_intensities();
  m_pViewer->update();
}

void MainWindow::on_actionView_polyhedron_facets_errors_triggered() {
  m_pScene->toggle_view_facet_errors();
  m_pViewer->update();
}

void MainWindow::on_actionView_all_sample_feature_intensities_triggered() {
  m_pScene->toggle_view_all_sample_feature_intensities();
  m_pViewer->update();
}

void MainWindow::on_actionView_all_sample_capacities_triggered() {
  m_pScene->toggle_view_all_sample_capacities();
  m_pViewer->update();
}

void MainWindow::on_actionView_all_sample_weights_triggered() {
  m_pScene->toggle_view_all_sample_weights();
  m_pViewer->update();
}

void MainWindow::on_actionView_element_classifications_triggered() {
  m_pScene->toggle_view_element_classifications();
  m_pViewer->update();
}

void MainWindow::on_actionView_gaussian_curvatures_triggered() {
  m_pScene->toggle_view_gaussian_curvatures();
  m_pViewer->update();
}

void MainWindow::on_actionView_maximal_normal_dihedrals_triggered() {
  m_pScene->toggle_view_maximal_normal_dihedrals();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_feature_intensities_triggered() {
  m_pScene->toggle_view_vertex_feature_intensities();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_capacities_triggered() {
  m_pScene->toggle_view_vertex_capacities();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_weights_triggered() {
  m_pScene->toggle_view_vertex_weights();
  m_pViewer->update();
}

void MainWindow::on_actionView_normal_dihedrals_triggered() {
  m_pScene->toggle_view_normal_dihedrals();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_feature_intensities_triggered() {
  m_pScene->toggle_view_edge_feature_intensities();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_capacities_triggered() {
  m_pScene->toggle_view_edge_capacities();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_weights_triggered() {
  m_pScene->toggle_view_edge_weights();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_feature_intensities_triggered() {
  m_pScene->toggle_view_facet_feature_intensities();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_capacities_triggered() {
  m_pScene->toggle_view_facet_capacities();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_weights_triggered() {
  m_pScene->toggle_view_facet_weights();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_in_start_points_triggered() {
  m_pScene->toggle_view_facet_in_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_in_end_points_triggered() {
  m_pScene->toggle_view_facet_in_end_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_in_links_triggered() {
  m_pScene->toggle_view_facet_in_links();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_out_start_points_triggered() {
  m_pScene->toggle_view_facet_out_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_facet_out_end_points_triggered() {
  m_pScene->toggle_view_facet_out_end_points();
  m_pViewer->update();
}
void MainWindow::on_actionView_facet_out_links_triggered() {
  m_pScene->toggle_view_facet_out_links();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_in_start_points_triggered() {
  m_pScene->toggle_view_edge_in_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_in_end_points_triggered() {
  m_pScene->toggle_view_edge_in_end_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_in_links_triggered() {
  m_pScene->toggle_view_edge_in_links();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_out_start_points_triggered() {
  m_pScene->toggle_view_edge_out_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_out_end_points_triggered() {
  m_pScene->toggle_view_edge_out_end_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_edge_out_links_triggered() {
  m_pScene->toggle_view_edge_out_links();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_in_start_points_triggered() {
  m_pScene->toggle_view_vertex_in_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_in_end_points_triggered() {
  m_pScene->toggle_view_vertex_in_end_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_in_links_triggered() {
  m_pScene->toggle_view_vertex_in_links();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_out_start_points_triggered() {
  m_pScene->toggle_view_vertex_out_start_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_out_end_points_triggered() {
  m_pScene->toggle_view_vertex_out_end_points();
  m_pViewer->update();
}

void MainWindow::on_actionView_vertex_out_links_triggered() {
  m_pScene->toggle_view_vertex_out_links();
  m_pViewer->update();
}

void MainWindow::open(QString file_name, OpenType open_type) {
  QFileInfo file_info(file_name);
  if (file_info.isFile() && file_info.isReadable()) {
    QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));
    bool suc = false;
    switch (open_type) {
    case OpenType::k_open_both:
      suc = m_pScene->open(file_name);
      break;
    case OpenType::k_open_input:
      suc = m_pScene->open_input(file_name);
      break;
    case OpenType::k_open_remesh:
      suc = m_pScene->open_remesh(file_name);
      break;
    default:
      break;
    }
    if (suc) {
      QSettings settings;
      settings.setValue("OFF open directory",
        file_info.absoluteDir().absolutePath());
      this->addToRecentFiles(file_name);
      updateViewerBBox();
      m_pViewer->update();
      update_menu_items();
    }
    QApplication::restoreOverrideCursor();
  }
}

void MainWindow::update_menu_items() {
  bool link_initialized = m_pScene->get_link_initialized();
  if (link_initialized) {
    // menu of "Polyhedron properties"
    ui->actionView_polyhedron_facets_errors->setEnabled(true);
    ui->menuView_samples_and_links->setEnabled(true);
    // menu of sample properties
    bool use_stratified_sampling = m_pScene->get_use_stratified_sampling();
    if (use_stratified_sampling) {
      ui->menuView_all_sample_properties->setEnabled(false);
      ui->menuView_vertex_sample_properties->setEnabled(true);
      ui->menuView_edge_sample_properties->setEnabled(true);
      ui->menuView_facet_sample_properties->setEnabled(true);
    }
    else {
      ui->menuView_all_sample_properties->setEnabled(true);
      ui->menuView_vertex_sample_properties->setEnabled(false);
      ui->menuView_edge_sample_properties->setEnabled(false);
      ui->menuView_facet_sample_properties->setEnabled(false);
    }
  }
  else {
    // menu of "Polyhedron properties"
    ui->actionView_polyhedron_facets_errors->setEnabled(false);
    ui->menuView_samples_and_links->setEnabled(false);
    // menu of sample properties
    ui->menuView_all_sample_properties->setEnabled(false);
    ui->menuView_vertex_sample_properties->setEnabled(false);
    ui->menuView_edge_sample_properties->setEnabled(false);
    ui->menuView_facet_sample_properties->setEnabled(false);
  }
}


