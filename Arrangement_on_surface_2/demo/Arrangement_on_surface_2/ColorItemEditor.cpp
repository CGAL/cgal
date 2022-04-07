/****************************************************************************
 **
 ** Copyright (C) 2009 Nokia Corporation and/or its subsidiary(-ies).
 ** All rights reserved.
 ** Contact: Nokia Corporation (qt-info@nokia.com)
 **
 ** This file is part of the examples of the Qt Toolkit.
 **
 ** $QT_BEGIN_LICENSE:LGPL$
 ** Commercial Usage
 ** Licensees holding valid Qt Commercial licenses may use this file in
 ** accordance with the Qt Commercial License Agreement provided with the
 ** Software or, alternatively, in accordance with the terms contained in
 ** a written agreement between you and Nokia.
 **
 ** GNU Lesser General Public License Usage
 ** Alternatively, this file may be used under the terms of the GNU Lesser
 ** General Public License version 2.1 as published by the Free Software
 ** Foundation and appearing in the file LICENSE.LGPL included in the
 ** packaging of this file.  Please review the following information to
 ** ensure the GNU Lesser General Public License version 2.1 requirements
 ** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
 **
 ** In addition, as a special exception, Nokia gives you certain additional
 ** rights.  These rights are described in the Nokia Qt LGPL Exception
 ** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
 **
 ** GNU General Public License Usage
 ** Alternatively, this file may be used under the terms of the GNU
 ** General Public License version 3.0 as published by the Free Software
 ** Foundation and appearing in the file LICENSE.GPL included in the
 ** packaging of this file.  Please review the following information to
 ** ensure the GNU General Public License version 3.0 requirements will be
 ** met: http://www.gnu.org/copyleft/gpl.html.
 **
 ** If you have questions regarding the use of this file, please contact
 ** Nokia at qt-info@nokia.com.
 ** $QT_END_LICENSE$
 **
 ** SPDX-License-Identifier: LGPL-2.1-only
 **
 ****************************************************************************/

#include <QtGui>
#include <iostream>
#include <QColorDialog>


#include "ColorItemEditor.h"

//! a pick color option
/*!
  \param widget A QWidget pointer
  \return a starting point to ask the user for color
*/
ColorItemEditor::ColorItemEditor( QWidget* widget ) : QPushButton( widget )
{
  this->setText( tr("Select a color") );
}

//! get the color reference
/*!
  \return the Qcolor object of the selected color
*/
QColor ColorItemEditor::color( ) const
{
  return this->m_color;
}

//! a pick color option
/*!
  \param color A QColor that the user selected to be changed
*/
void ColorItemEditor::setColor( QColor color )
{
  this->m_color = color;
}

//! checking if the color reference that user wanted is valid or not
/*!
*/
void ColorItemEditor::mousePressEvent(QMouseEvent* /* e */)
{
  QColor selectedColor = QColorDialog::getColor(this->m_color);
  if (selectedColor.isValid()) {
    // std::cout << selectedColor.name().toStdString() << std::endl;
    this->setColor( selectedColor );
  }

  Q_EMIT confirmed();
}
