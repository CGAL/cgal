/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// include files for Qt
#include <qprinter.h>
#include <qpainter.h>
#include <stdio.h>

// application specific includes
#include "toolview.h"
#include "tooldoc.h"


ToolView::ToolView(ToolDoc* pDoc,
                   const QGLFormat& format,
                   QWidget *parent,
                   const char* name,
                   int wflags)
 : QGLView(format, parent, name, wflags)
{
    doc=pDoc;
    mousePressed = false;
    m_LeftButtonDown = false;
    m_RightButtonDown = false;
    m_Move = false;
}

ToolView::~ToolView()
{
}

ToolDoc *ToolView::getDocument() const
{
	return doc;
}

void ToolView::update(ToolView* pSender){
	if(pSender != this)
		repaint();
}

void ToolView::print(QPrinter *pPrinter)
{
  if (pPrinter->setup(this))
  {
		QPainter p;
		p.begin(pPrinter);
		
		p.end();
  }
}

void ToolView::focusInEvent(QFocusEvent *)
{
  //This is called when the widget gets the focus.
	emit(invalidate_view(this));
}

void ToolView::closeEvent(QCloseEvent*)
{
  // LEAVE THIS EMPTY: THE EVENT FILTER IN THE ToolApp CLASS TAKES CARE FOR CLOSING
  // QWidget closeEvent must be prevented.
}


/////////////////////////////////////////////////////
// MOUSE EVENTS
/////////////////////////////////////////////////////

void ToolView::HandleMouseButton(int x, int y)
{
	CVector3d theVec = theArcball.Intersect(x,theViewport.yRes()-y,
		theCamera,theViewport);
	theArcball.EndDrag(theVec);
	theArcball.SetMode(m_LeftButtonDown+2*m_RightButtonDown);
	theVec = theArcball.Intersect(x,theViewport.yRes()-y,
		theCamera,theViewport);
	theArcball.BeginDrag(theVec);
}

void ToolView::mousePressEvent(QMouseEvent *event)
{
  if(event->button() & LeftButton)
    m_LeftButtonDown = true;
  if(event->button() & RightButton)
    m_RightButtonDown = true;

  if(m_LeftButtonDown || m_RightButtonDown)
  {
    QPoint point = event->pos();
  	HandleMouseButton(point.x(),point.y());
	}
}

void ToolView::mouseMoveEvent(QMouseEvent *event)
{
  if(m_LeftButtonDown || m_RightButtonDown)
  {
    QPoint point = event->pos();
		CVector3d theVec = theArcball.Intersect(point.x(),
		                                        theViewport.yRes()-point.y(),
										                        theCamera,
										                        theViewport);
		theArcball.Motion(theVec);
  	m_Move = true;
	  updateGL();
  }
}

void ToolView::mouseReleaseEvent(QMouseEvent *event)
{
  m_LeftButtonDown = false;
  m_RightButtonDown = false;
	
	m_Move = false;
  QPoint point = event->pos();
	HandleMouseButton(point.x(),point.y());
	updateGL();
}

void ToolView::drawScene(bool superimposededges, bool superimposedvertices,
                         bool first)
{
  doc->drawScene(m_Move, superimposededges, superimposedvertices, first, m_Smooth, m_Lighting,
    m_fore_color[0], m_fore_color[1], m_fore_color[2]);
}


#include "toolview.moc"

