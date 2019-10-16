#ifndef CGAL_DOUBLE_EDIT_H
#define CGAL_DOUBLE_EDIT_H



#include <QObject>
#include <QtCore/qglobal.h>
#include <QLineEdit>
#include "Scene_config.h"


class DoubleValidator;
class SCENE_EXPORT DoubleEdit : public QLineEdit {

  Q_OBJECT

public:
  DoubleEdit(QWidget* parent = nullptr);
  ~DoubleEdit();
  double value() const;
  void setValue(double d);
  void setMinimum(double d);
  void setMaximum(double d);
  void setRange(double min, double max);
private:
  DoubleValidator* validator;
};
#endif // CGAL_DOUBLE_EDIT_H
