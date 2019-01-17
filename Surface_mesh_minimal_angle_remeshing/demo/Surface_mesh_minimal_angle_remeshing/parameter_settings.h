#ifndef PARAMETERSETTINGS_H_
#define PARAMETERSETTINGS_H_

#include "ui_parameter_settings.h"

class ParameterSettings : public QDialog, private Ui_ParameterSettingsDialog {
  Q_OBJECT
 public:
  explicit ParameterSettings(QWidget *parent = 0) {
    setupUi(this);
    // fill the items in cb_sample_number_strategy
    QStringList sample_number_strategy_types;
    sample_number_strategy_types.append("Fixed");
    sample_number_strategy_types.append("Variable");
    cb_sample_number_strategy->addItems(sample_number_strategy_types);
    cb_sample_number_strategy->setCurrentIndex(0);
    // fill the items in cb_sample_strategy
    QStringList sample_strategy_types;
    sample_strategy_types.append("Uniform");
    sample_strategy_types.append("Adaptive");
    cb_sample_strategy->addItems(sample_strategy_types);
    cb_sample_strategy->setCurrentIndex(1);
    // fill the items in cb_facet_relocate_type, 
    //                   cb_edge_relocate_type,
    //                   cb_vertex_relocate_type
    QStringList optimize_types;
    optimize_types.append("None");
    optimize_types.append("Input to Remesh");
    optimize_types.append("Remesh to Input");
    optimize_types.append("Both");
    cb_facet_optimize_type->addItems(optimize_types);
    cb_facet_optimize_type->setCurrentIndex(3);
    cb_edge_optimize_type->addItems(optimize_types);
    cb_edge_optimize_type->setCurrentIndex(3);
    cb_vertex_optimize_type->addItems(optimize_types);
    cb_vertex_optimize_type->setCurrentIndex(3);
    // fill the items in cb_relocate_strategy
    QStringList optimize_strategy_types;
    optimize_strategy_types.append("Approximation");
    optimize_strategy_types.append("Interpolation");
    cb_optimize_strategy->addItems(optimize_strategy_types);
    cb_optimize_strategy->setCurrentIndex(0);
    // fill the items in cb_edge_flip_strategy
    QStringList edge_flip_strategy_types;
    edge_flip_strategy_types.append("Improve valence");
    edge_flip_strategy_types.append("Improve angle");
    cb_edge_flip_strategy->addItems(edge_flip_strategy_types);
    cb_edge_flip_strategy->setCurrentIndex(1);
    // fill the items in cb_relocate_strategy
    QStringList relocate_strategy_types;
    relocate_strategy_types.append("Barycenter");
    relocate_strategy_types.append("Cvt barycenter");
    cb_relocate_strategy->addItems(relocate_strategy_types);
    cb_relocate_strategy->setCurrentIndex(1);
  }

  // disable copy/move construction
  ParameterSettings(const ParameterSettings &) = delete;
  ParameterSettings(const ParameterSettings &&) = delete;
  ParameterSettings &operator = (const ParameterSettings &) = delete;
  ParameterSettings &operator = (const ParameterSettings &&) = delete;

  // access functions

  // 1) general parameters
  double get_max_error_threshold() const { return sb_max_error_threshold->value(); }
  void set_max_error_threshold(double value) { sb_max_error_threshold->setValue(value); }

  double get_min_angle_threshold() const { return sb_min_angle_threshold->value(); }
  void set_min_angle_threshold(double value) { sb_min_angle_threshold->setValue(value); }

  int get_max_mesh_complexity() const { return sb_max_mesh_complexity->value(); }
  void set_max_mesh_complexity(int value) { sb_max_mesh_complexity->setValue(value); }

  double get_smooth_angle_delta() const { return sb_smooth_angle_delta->value(); }
  void set_smooth_angle_delta(double value) { sb_smooth_angle_delta->setValue(value); }

  bool get_apply_edge_flip() const { return cb_apply_edge_flip->isChecked(); }
  void set_apply_edge_flip(bool value) { cb_apply_edge_flip->setChecked(value); }

  int get_edge_flip_strategy() const { return cb_edge_flip_strategy->currentIndex(); }
  void set_edge_flip_strategy(int value) { cb_edge_flip_strategy->setCurrentIndex(value); }

  bool get_flip_after_split_and_collapse() const { return cb_flip_after_split_and_collapse->isChecked(); }
  void set_flip_after_split_and_collapse(bool value) { cb_flip_after_split_and_collapse->setChecked(value); }

  bool get_relocate_after_local_operations() const { return cb_relocate_after_local_operations->isChecked(); }
  void set_relocate_after_local_operations(bool value) { cb_relocate_after_local_operations->setChecked(value); }

  int get_relocate_strategy() const { return cb_relocate_strategy->currentIndex(); }
  void set_relocate_strategy(int value) { return cb_relocate_strategy->setCurrentIndex(value); }

  bool get_keep_vertex_in_one_ring() const { return cb_keep_vertex_in_one_ring->isChecked(); }
  void set_keep_vertex_in_one_ring(bool value) { cb_keep_vertex_in_one_ring->setChecked(value); }

  bool get_use_local_aabb_tree() const { return cb_use_local_aabb_tree->isChecked(); }
  void set_use_local_aabb_tree(bool value) { cb_use_local_aabb_tree->setChecked(value); }

  int get_collapsed_list_size() const { return sb_collapsed_list_size->value(); }
  void set_collapsed_list_size(int value) { sb_collapsed_list_size->setValue(value); }

  bool get_decrease_max_errors() const { return cb_decrease_max_errors->isChecked(); }
  void set_decrease_max_errors(bool value) { cb_decrease_max_errors->setChecked(value); }

  bool get_track_information() const { return cb_track_information->isChecked(); }
  void set_track_information(bool value) { cb_track_information->setChecked(value); }

  bool get_apply_initial_mesh_simplification() const { return cb_apply_initial_mesh_simplification->isChecked(); }
  void set_apply_initial_mesh_simplification(bool value) { cb_apply_initial_mesh_simplification->setChecked(value); }

  bool get_apply_final_vertex_relocation() const { return cb_apply_final_vertex_relocation->isChecked(); }
  void set_apply_final_vertex_relocation(bool value) { cb_apply_final_vertex_relocation->setChecked(value); }

  // 2) sample parameters
  int get_samples_per_facet_in() const { return sb_samples_per_facet_in->value(); }
  void set_samples_per_facet_in(int value) { sb_samples_per_facet_in->setValue(value); }

  int get_samples_per_facet_out() const { return sb_samples_per_facet_out->value(); }
  void set_samples_per_facet_out(int value) { sb_samples_per_facet_out->setValue(value); }

  int get_max_samples_per_area() const { return sb_max_samples_per_area->value(); }
  void set_max_samples_per_area(int value) { sb_max_samples_per_area->setValue(value); }

  int get_min_samples_per_triangle() const { return sb_min_samples_per_triangle->value(); }
  void set_min_samples_per_triangle(int value) { sb_min_samples_per_triangle->setValue(value); }

  int get_bvd_iteration_count() const { return sb_bvd_iteration_count->value(); }
  void set_bvd_iteration_count(int value) { sb_bvd_iteration_count->setValue(value); }

  int get_sample_number_strategy() const { return cb_sample_number_strategy->currentIndex(); }
  void set_sample_number_strategy(int value) { cb_sample_number_strategy->setCurrentIndex(value); }

  int get_sample_strategy() const { return cb_sample_strategy->currentIndex(); }
  void set_sample_strategy(int value) { cb_sample_strategy->setCurrentIndex(value); }

  bool get_use_stratified_sampling() const { return cb_use_stratified_sampling->isChecked(); }
  void set_use_stratified_sampling(bool value) { cb_use_stratified_sampling->setChecked(value); }

  // 3) feature function parameters
  double get_sum_theta() const { return sb_sum_theta->value(); }
  void set_sum_theta(double value) { sb_sum_theta->setValue(value); }

  double get_sum_delta() const { return sb_sum_delta->value(); }
  void set_sum_delta(double value) { sb_sum_delta->setValue(value); }

  double get_dihedral_theta() const { return sb_dihedral_theta->value(); }
  void set_dihedral_theta(double value) { sb_dihedral_theta->setValue(value); }

  double get_dihedral_delta() const { return sb_dihedral_delta->value(); }
  void set_diheral_delta(double value) { sb_dihedral_delta->setValue(value); }

  double get_feature_difference_delta() const { return sb_feature_difference_delta->value(); }
  void set_feature_difference_delta(double value) { sb_feature_difference_delta->setValue(value); }
 
  double get_feature_control_delta() const { return sb_feature_control_delta->value(); }
  void set_feature_control_delta(double value) { sb_feature_control_delta->setValue(value); }

  bool get_inherit_element_types() const { return cb_inherit_element_types->isChecked(); }
  void set_inherit_element_types(bool value) { cb_inherit_element_types->setChecked(value); }

  bool get_use_feature_intensity_weights() const { return cb_use_feature_intensity_weights->isChecked(); }
  void set_use_feature_intensity_weights(bool value) { cb_use_feature_intensity_weights->setChecked(value); }

  // 4) vertex optimization parameters
  int get_vertex_optimize_count() const { return sb_vertex_optimize_count->value(); }
  void set_vertex_optimize_count(int value) { sb_vertex_optimize_count->setValue(value); }

  double get_vertex_optimize_ratio() const { return sb_vertex_optimize_ratio->value(); }
  void set_vertex_optimize_ratio(double value) { sb_vertex_optimize_ratio->setValue(value); }

  int get_stencil_ring_size() const { return sb_stencil_ring_size->value(); }
  void set_stencil_ring_size(int value) { sb_stencil_ring_size->setValue(value); }

  int get_optimize_strategy() const { return cb_optimize_strategy->currentIndex(); }
  void set_optimize_strategy(int value) { cb_optimize_strategy->setCurrentIndex(value); }

  int get_facet_optimize_type() const { return cb_facet_optimize_type->currentIndex(); }
  void set_facet_optimize_type(int value) { cb_facet_optimize_type->setCurrentIndex(value); }

  int get_edge_optimize_type() const { return cb_edge_optimize_type->currentIndex(); }
  void set_edge_optimize_type(int value) { cb_edge_optimize_type->setCurrentIndex(value); }

  int get_vertex_optimize_type() const { return cb_vertex_optimize_type->currentIndex(); }
  void set_vertex_optimize_type(int value) { cb_vertex_optimize_type->setCurrentIndex(value); }

  bool get_optimize_after_local_operations() const { return cb_optimize_after_local_operations->isChecked(); }
  void set_optimize_after_local_operations(bool value) { cb_optimize_after_local_operations->setChecked(value); }
};

#endif // PARAMETERSETTINGS_H_
