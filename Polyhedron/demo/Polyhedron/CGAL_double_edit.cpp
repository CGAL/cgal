#include "CGAL_double_edit.h"
#include <iostream>
#include <QDoubleValidator>

class DoubleValidator : public QDoubleValidator
{
public:
  DoubleValidator(QObject* parent = nullptr)
    : QDoubleValidator(parent)
  {
    setLocale(QLocale::C);
  }

  void fixup ( QString & input ) const
  {
    input.replace(".", locale().decimalPoint());
    input.replace(",", locale().decimalPoint());
    QDoubleValidator::fixup(input);
  }
  QValidator::State validate ( QString & input, int & pos ) const
  {
    fixup(input);
    return QDoubleValidator::validate(input, pos);
  }
};

  DoubleEdit::DoubleEdit(QWidget *parent)
    : QLineEdit(parent)
  {
    validator = new DoubleValidator(this);
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

  void DoubleEdit::setRange(double rmin, double rmax)
  {
    this->validator->setRange(rmin, rmax, this->validator->decimals());
  }

  double DoubleEdit::getValue()
  {
    return this->value();
  }

  double DoubleEdit::getMinimum()
  {
    return this->validator->bottom();
  }

  double DoubleEdit::getMaximum()
  {
    return this->validator->top();
  }
