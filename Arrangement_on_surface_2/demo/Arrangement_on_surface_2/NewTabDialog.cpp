// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "NewTabDialog.h"
#include "ArrangementDemoWindow.h"
#include "ArrangementTypesUtils.h"
#include "ui_NewTabDialog.h"
#include <QButtonGroup>


NewTabDialog::NewTabDialog( QWidget* parent ) :
  QDialog( parent ),
  ui( new Ui::NewTabDialog ),
  buttonGroup( new QButtonGroup )
{
  using TraitsType = demo_types::TraitsType;

  this->ui->setupUi(this);

  this->buttonGroup->addButton(
    this->ui->segmentRadioButton, static_cast<int>(TraitsType::SEGMENT_TRAITS));
  this->buttonGroup->addButton(
    this->ui->polylineRadioButton,
    static_cast<int>(TraitsType::POLYLINE_TRAITS));
  this->buttonGroup->addButton(
    this->ui->linearRadioButton, static_cast<int>(TraitsType::LINEAR_TRAITS));
#ifdef CGAL_USE_CORE
  this->buttonGroup->addButton(
    this->ui->conicRadioButton, static_cast<int>(TraitsType::CONIC_TRAITS));
  this->buttonGroup->addButton(
    this->ui->algebraicRadioButton,
    static_cast<int>(TraitsType::ALGEBRAIC_TRAITS));
  this->buttonGroup->addButton(
    this->ui->bezierRadioButton, static_cast<int>(TraitsType::BEZIER_TRAITS));
  this->buttonGroup->addButton(
    this->ui->rationalFunctionRadioButton,
    static_cast<int>(TraitsType::RATIONAL_FUNCTION_TRAITS));
#else
  this->ui->conicRadioButton->hide();
  this->ui->algebraicRadioButton->hide();
  this->ui->bezierRadioButton->hide();
  this->ui->rationalFunctionRadioButton->hide();
#endif
}

int NewTabDialog::checkedId( ) const
{
  return this->buttonGroup->checkedId( );
}
