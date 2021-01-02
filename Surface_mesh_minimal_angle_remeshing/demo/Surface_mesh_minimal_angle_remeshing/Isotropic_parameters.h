#ifndef CGAL_ISOTROPIC_PARAMETERS_H
#define CGAL_ISOTROPIC_PARAMETERS_H

#include "ui_Isotropic_parameters.h"

class IsotropicParameters : public QDialog, private Ui_IsotropicParametersDialog {
  Q_OBJECT
public:
  explicit IsotropicParameters(QWidget *parent = 0) {
    setupUi(this);
  }

  // disable copy/move construction
  IsotropicParameters(const IsotropicParameters &) = delete;
  IsotropicParameters(const IsotropicParameters &&) = delete;
  IsotropicParameters &operator = (const IsotropicParameters &) = delete;
  IsotropicParameters &operator = (const IsotropicParameters &&) = delete;

  // access functions

  // 1) general parameters
  double get_target_edge_length() const { return sb_target_edge_length->value(); }
  void set_target_edge_length(double value) { sb_target_edge_length->setValue(value); }

  int get_smooth_iteration_count() const { return sb_smooth_iteration_count->value(); }
  void set_smooth_iteration_count(double value) { sb_smooth_iteration_count->setValue(value); }
};

#endif // CGAL_ISOTROPIC_PARAMETERS_H
