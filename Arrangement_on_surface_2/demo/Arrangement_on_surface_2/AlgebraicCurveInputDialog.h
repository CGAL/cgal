// Copyright (c) 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#ifndef ALGEBRAICCURVEINPUTDIALOG_H
#define ALGEBRAICCURVEINPUTDIALOG_H

#include <QDialog>

namespace Ui
{
class AlgebraicCurveInputDialog;
}

class AlgebraicCurveInputDialog : public QDialog
{
    Q_OBJECT

public:
    explicit AlgebraicCurveInputDialog(QWidget *parent = nullptr);
    ~AlgebraicCurveInputDialog();
    std::string getLineEditText();
    Ui::AlgebraicCurveInputDialog* getUi(){return this->ui;}

private:
    Ui::AlgebraicCurveInputDialog *ui;
};

#endif // ALGEBRAICCURVEINPUTDIALOG_H
