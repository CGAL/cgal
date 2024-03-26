#ifndef CGAL_DOUBLE_EDIT_H
#define CGAL_DOUBLE_EDIT_H



#include <QObject>
#include <QtCore/qglobal.h>
#include <QLineEdit>
#include "Scene_config.h"


class DoubleValidator;
class SCENE_EXPORT DoubleEdit : public QLineEdit {

  Q_OBJECT
  Q_PROPERTY(double value READ getValue WRITE setValue)
  Q_PROPERTY(double minimum READ getMinimum WRITE setMinimum)
  Q_PROPERTY(double maximum READ getMaximum WRITE setMaximum)
public:
  DoubleEdit(QWidget* parent = nullptr);
  ~DoubleEdit();
  double value() const;
  void setValue(double d);
  void setMinimum(double d);
  void setMaximum(double d);
  void setRange(double rmin, double rmax);
  double getValue();
  double getMinimum();
  double getMaximum();
private:
  DoubleValidator* validator;
};
#endif // CGAL_DOUBLE_EDIT_H
