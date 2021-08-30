// Copyright (c) 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include "AlgebraicCurveInputDialog.h"
#include "ui_AlgebraicCurveInputDialog.h"

AlgebraicCurveInputDialog::AlgebraicCurveInputDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AlgebraicCurveInputDialog)
{
    ui->setupUi(this);
    ui->lineEdit->setFocus();
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
