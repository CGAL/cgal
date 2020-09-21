// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#include "RationalCurveInputDialog.h"
#include "ui_RationalCurveInputDialog.h"

RationalCurveInputDialog::RationalCurveInputDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RationalCurveInputDialog)
{
    ui->setupUi(this);
    setTabOrder(ui->numeratorLineEdit, ui->denominatorLineEdit);
    ui->numeratorLineEdit->setFocus();
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
