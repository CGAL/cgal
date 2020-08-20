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

std::string AlgebraicCurveInputDialog::getLineEditText()
{
    QString lineEditText = ui->lineEdit->text();
    return lineEditText.toStdString();
}
