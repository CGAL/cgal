#include "CGAL_double_edit.h"

#include <QDoubleValidator>

  DoubleEdit::DoubleEdit(QWidget *parent)
    : QLineEdit()
  {
    validator = new QDoubleValidator(this);
    validator->setLocale(QLocale::C);
    this->setValidator(validator);
  }

  DoubleEdit::~DoubleEdit()
  {
    delete validator;
  }
  double DoubleEdit::value() const
  {
    return this->text().toDouble();
  }

  void DoubleEdit::setValue(double d)
  {
    this->setText(tr("%1").arg(d));
  }

  void DoubleEdit::setMinimum(double d)
  {
    this->validator->setBottom(d);
  }

  void DoubleEdit::setMaximum(double d)
  {
    this->validator->setTop(d);
  }
  #include "CGAL_double_edit.moc"

