#ifndef CGAL_DOUBLE_EDIT_H
#define CGAL_DOUBLE_EDIT_H



#include <QObject>
#include <QtCore/qglobal.h>
#include <QLineEdit>
#ifdef cgal_double_edit_EXPORTS
#  define CGAL_DOUBLE_EDIT_EXPORT Q_DECL_EXPORT
#else
#  define CGAL_DOUBLE_EDIT_EXPORT Q_DECL_IMPORT
#endif

class QDoubleValidator;
class CGAL_DOUBLE_EDIT_EXPORT DoubleEdit : public QLineEdit {

  Q_OBJECT

public:
  DoubleEdit(QWidget* parent = nullptr);
  ~DoubleEdit();
  double value() const;
  void setValue(double d);
  void setMinimum(double d);
  void setMaximum(double d);
private:
  QDoubleValidator* validator;
};
#endif // CGAL_DOUBLE_EDIT_H
