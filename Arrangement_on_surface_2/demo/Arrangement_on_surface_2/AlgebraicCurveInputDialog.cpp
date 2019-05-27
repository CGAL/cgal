#include "AlgebraicCurveInputDialog.h"
#include "ui_AlgebraicCurveInputDialog.h"

AlgebraicCurveInputDialog::AlgebraicCurveInputDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AlgebraicCurveInputDialog)
{
    ui->setupUi(this);
}

AlgebraicCurveInputDialog::~AlgebraicCurveInputDialog()
{
    delete ui;
}

//! member function to get the expression entered in the dialogue box
/*!
  \return string value of the polynomial expression entered
*/
std::string AlgebraicCurveInputDialog::getLineEditText()
{
    QString lineEditText = ui->lineEdit->text();
    return lineEditText.toStdString();
}
