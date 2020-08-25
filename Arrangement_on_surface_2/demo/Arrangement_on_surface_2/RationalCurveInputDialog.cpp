#include "RationalCurveInputDialog.h"
#include "ui_RationalCurveInputDialog.h"

RationalCurveInputDialog::RationalCurveInputDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RationalCurveInputDialog)
{
    ui->setupUi(this);
}

RationalCurveInputDialog::~RationalCurveInputDialog()
{
    delete ui;
}

std::string RationalCurveInputDialog::getNumeratorText()
{
  return ui->numeratorLineEdit->text().toStdString();
}

std::string RationalCurveInputDialog::getDenominatorText()
{
  return ui->denominatorLineEdit->text().toStdString();
}
