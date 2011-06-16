#ifndef SHOW_POINT_DIALOG_H
#define SHOW_POINT_DIALOG_H
#include "config.h"

#include "Point_dialog_config.h"

#include <QDialog>

namespace Ui {
  class Show_point_dialog;
}

class POINT_DIALOG_EXPORT Show_point_dialog :
  public QDialog
{
  Q_OBJECT
public:
  Show_point_dialog(QWidget* parent = 0);
  ~Show_point_dialog();

  bool has_correct_coordinates() const;

  double get_x() const;
  double get_y() const;
  double get_z() const;

protected slots:
  void interprete_string(const QString&);

protected:
  Ui::Show_point_dialog* ui;
  bool m_has_correct_coordinates;
};

#endif
