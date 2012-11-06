// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "NewTabDialog.h"
#include "ArrangementDemoWindow.h"
#include "ui_NewTabDialog.h"

NewTabDialog::NewTabDialog( QWidget* parent, Qt::WindowFlags f ) :
  QDialog( parent, f ),
  ui( new Ui::NewTabDialog ),
  buttonGroup( new QButtonGroup )
{
  this->ui->setupUi( this );
    
  this->buttonGroup->addButton( this->ui->segmentRadioButton,
                                ArrangementDemoWindow::SEGMENT_TRAITS );
  this->buttonGroup->addButton( this->ui->polylineRadioButton,
                                ArrangementDemoWindow::POLYLINE_TRAITS );
  this->buttonGroup->addButton( this->ui->conicRadioButton,
                                ArrangementDemoWindow::CONIC_TRAITS );
  this->buttonGroup->addButton( this->ui->linearRadioButton,
                                ArrangementDemoWindow::LINEAR_TRAITS );
  this->buttonGroup->addButton( this->ui->circularArcRadioButton,
                                ArrangementDemoWindow::CIRCULAR_ARC_TRAITS );
  // this->buttonGroup->addButton( this->ui->algebraicRadioButton,
  //                               ArrangementDemoWindow::ALGEBRAIC_TRAITS );
}

int NewTabDialog::checkedId( ) const
{
  return this->buttonGroup->checkedId( );
}
