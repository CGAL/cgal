// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_WIDGET_STYLE_EDITOR_AUX_H
#define CGAL_QT_WIDGET_STYLE_EDITOR_AUX_H

#include <qcolor.h>
#include <qcolordialog.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qcombobox.h>
#include <qscrollview.h>
#include <qpixmap.h>

class Color_selector : public QPushButton
{
  Q_OBJECT
public:
  Color_selector(QColor c = Qt::black,
		 QWidget* parent = 0, const char* name = 0)
    : QPushButton(parent, name)
  {
    setColor(c);
    connect(this, SIGNAL(clicked()),
	    this, SLOT(color_dialog()) );
  }

  virtual ~Color_selector() {};

  QColor value() const
  {
    return color;
  }

public slots:
  void setColor(QColor c)
  {
    color = c;

    QPixmap pix(24,20);
    pix.fill(c);
    setPixmap(pix);

    emit newColor(c);
  }

signals:
  void newColor(QColor);

private slots:
  void color_dialog()
  {
    QColor c = QColorDialog::getColor(value());
    if( c.isValid() )
      setColor(c);
  }

private:
  QColor color;
};

class Int_selector : public QSpinBox
{
  Q_OBJECT
public:
  Int_selector(int i, QWidget *parent = 0, const char *name = 0)
    : QSpinBox(-INT_MAX, INT_MAX, 1, parent, name)
  {
    setValue(i);
  }

  virtual ~Int_selector() {};
};

class Bool_selector : public QComboBox
{
  Q_OBJECT
public:
  Bool_selector(bool b_, QWidget *parent = 0, const char *name = 0)
    : QComboBox(false, parent, name)
  {
    insertItem("False");
    insertItem("True");

    if(b_)
      setCurrentItem(1);
    else
      setCurrentItem(0);
  }

  virtual ~Bool_selector() {};

  bool value() const
  {
    return currentItem() == 1;
  }
};

class Point_style_selector : public QComboBox
{
  Q_OBJECT
public:
  typedef ::CGAL::PointStyle PointStyle;
  Point_style_selector(PointStyle s,
		       QWidget *parent = 0, const char *name = 0)
    : QComboBox(false, parent, name)
  {
    insertItem("Pixel");
    insertItem("Cross");
    insertItem("Plus");
    insertItem("Circle");
    insertItem("Disc");
    insertItem("Rect");
    insertItem("Box");

    setCurrentItem(static_cast<int>(s));
  }

  virtual ~Point_style_selector() {};

  PointStyle value() const
  {
    return PointStyle(currentItem());
  }
};

#endif // CGAL_QT_WIDGET_STYLE_EDITOR_AUX_H
