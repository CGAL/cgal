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

#include <CGAL/basic.h>


#include "Qt_widget_style_editor.h"
#include "Qt_widget_style_editor-aux.h"

#include <qcolor.h>
#include <qlabel.h>
#include <qvbox.h>
#include <qhbox.h>
#include <qlayout.h>
#include <qgrid.h>
#include <qvariant.h>

namespace CGAL {

Qt_widget_style_editor::Qt_widget_style_editor(Style* style,
					       QWidget *parent,
					       const char *name)
  : QFrame(parent, name), style(style)
{
  typedef Style::const_iterator iterator;

  QGridLayout* layout = new QGridLayout(this);
  layout->addColSpacing(1,5);

  const int labels_col = 0; // column number of labels
  const int selectors_col = 2; // column number of selectors

  int row = 0;
  for(iterator it=style->begin();
      it != style->end();
      ++it)
    {
      QLabel* label = new QLabel( it.key(), this);
      layout->addWidget(label, row, labels_col);

      QWidget* selector = 0;
      switch( it.data().type() ) {
      case QVariant::Color:
	selector = new Color_selector(it.data().toColor(), this);
	connect(selector, SIGNAL(newColor(QColor)),
		this, SLOT(map(QColor)));
	break;
      case QVariant::Int:
	selector = new Int_selector(it.data().toInt(), this);
	connect(selector, SIGNAL(valueChanged(int)),
		this, SLOT(map(int)));
	break;
      case QVariant::Bool:
	selector = new Bool_selector(it.data().toBool(),
				     this);
	connect(selector, SIGNAL(toggled(bool)),
		this, SLOT(map(bool)));
	break;
      case QVariant::UInt:
	selector =
	  new Point_style_selector(PointStyle(it.data().toUInt()),
				   this);
	connect(selector, SIGNAL(activated(int)),
		this, SLOT(pointstyle(int)));
	break;
      default:
	CGAL_error();
	break;
      }

      mapper[selector]=it.key();

      layout->addWidget(selector, row, selectors_col);

      ++row;
    }
}

void Qt_widget_style_editor::map(QColor c)
{
  const QObject* s = sender();
  if( mapper.contains(s) )
    style->setColor(mapper[s], c);
  emit styleChanged();
}

void Qt_widget_style_editor::map(int i)
{
  const QObject* s = sender();
  if( mapper.contains(s) )
    style->setInt(mapper[s], i);
  emit styleChanged();
}

void Qt_widget_style_editor::map(bool b)
{
  const QObject* s = sender();
  if( mapper.contains(s) )
    style->setBool(mapper[s], b);
  emit styleChanged();
}

} // end namespace CGAL

// moc_source_file: Qt_widget_style_editor.h
#include "Qt_widget_style_editor.moc"

// moc_source_file: Qt_widget_style_editor-aux.h
#include "Qt_widget_style_editor-aux.moc"

