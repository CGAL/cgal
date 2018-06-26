/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.
 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 version 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0

#ifndef CGAL_IMAGE_INTERFACE_H
#define CGAL_IMAGE_INTERFACE_H
#include <CGAL/license/GraphicsView.h>

#include <QDialog>
#include <QWidget>
#include <QApplication>

#include "ui_ImageInterface.h"
class ImageInterface: public QDialog, public Ui::ImageInterface
{
  Q_OBJECT
  qreal ratio;
  QWidget *currentlyFocused;
public:
  ImageInterface(QWidget *parent, qreal ratio)
    : QDialog(parent), ratio(ratio)
  {
    currentlyFocused = NULL;
    setupUi(this);
    connect(imgHeight, SIGNAL(valueChanged(int)),
            this, SLOT(imgHeightValueChanged(int)));

    connect(imgWidth, SIGNAL(valueChanged(int)),
            this, SLOT(imgWidthValueChanged(int)));

    connect(qApp, SIGNAL(focusChanged(QWidget*, QWidget*)),
            this, SLOT(onFocusChanged(QWidget*, QWidget*)));
  }
private Q_SLOTS:
  void imgHeightValueChanged(int i)
  {
    if(currentlyFocused == imgHeight
       && ratioCheckBox->isChecked())
    {imgWidth->setValue(int(i*ratio));}
  }

  void imgWidthValueChanged(int i)
  {
    if(currentlyFocused == imgWidth
       && ratioCheckBox->isChecked())
    {imgHeight->setValue(int(i/ratio));}
  }

  void onFocusChanged(QWidget*, QWidget* now)
  {
    currentlyFocused = now;
  }
};
#endif // CGAL_IMAGE_INTERFACE_H
