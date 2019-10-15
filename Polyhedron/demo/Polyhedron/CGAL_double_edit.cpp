#include "CGAL_double_edit.h"

#include <QDoubleValidator>

  DoubleEdit::DoubleEdit(QWidget* parent = nullptr)
    : QLineEdit(parent)
  {
    QDoubleValidator* validator = new QDoubleValidator(parent);
    validator->setLocale(QLocale::C);
    this->setValidator(validator);
  }

  #include "CGAL_double_edit.moc"

