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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/mouseGrabber.h>

namespace CGAL{
namespace qglviewer{


// Static private variable
CGAL_INLINE_FUNCTION
 QList<MouseGrabber *> &MouseGrabber::MouseGrabberPool() {
  static QList<MouseGrabber*> MouseGrabberPool_;
  void* p = qApp->property("qglviewer mouse grabber pool").value<void*>();
  if(p == 0) {
    p = (void*)(&MouseGrabberPool_);
    qApp->setProperty("qglviewer mouse grabber pool", QVariant::fromValue(p));
  }
  return *static_cast<QList<MouseGrabber *> * >(p);
}

/*! Default constructor.

Adds the created MouseGrabber in the MouseGrabberPool(). grabsMouse() is set to
\c false. */
CGAL_INLINE_FUNCTION
MouseGrabber::MouseGrabber() : grabsMouse_(false) { addInMouseGrabberPool(); }

/*! Adds the MouseGrabber in the MouseGrabberPool().

All created MouseGrabber are automatically added in the MouseGrabberPool() by
the constructor. Trying to add a MouseGrabber that already
isInMouseGrabberPool() has no effect.

Use removeFromMouseGrabberPool() to remove the MouseGrabber from the list, so
that it is no longer tested with checkIfGrabsMouse() by the CGAL::QGLViewer, and hence
can no longer grab mouse focus. Use isInMouseGrabberPool() to know the current
state of the MouseGrabber. */
CGAL_INLINE_FUNCTION
void MouseGrabber::addInMouseGrabberPool() {
  if (!isInMouseGrabberPool())
    MouseGrabber::MouseGrabberPool().append(this);
}

/*! Removes the MouseGrabber from the MouseGrabberPool().

See addInMouseGrabberPool() for details. Removing a MouseGrabber that is not in
MouseGrabberPool() has no effect. */
CGAL_INLINE_FUNCTION
void MouseGrabber::removeFromMouseGrabberPool() {
  if (isInMouseGrabberPool())
    MouseGrabber::MouseGrabberPool().removeAll(const_cast<MouseGrabber *>(this));
}

/*! Clears the MouseGrabberPool().

 Use this method only if it is faster to clear the MouseGrabberPool() and then
 to add back a few MouseGrabbers than to remove each one independently. Use
 CGAL::QGLViewer::setMouseTracking(false) instead if you want to disable mouse
 grabbing.

 When \p autoDelete is \c true, the MouseGrabbers of the MouseGrabberPool() are
 actually deleted (use this only if you're sure of what you do). */
CGAL_INLINE_FUNCTION
void MouseGrabber::clearMouseGrabberPool(bool autoDelete) {
  if (autoDelete)
    qDeleteAll(MouseGrabber::MouseGrabberPool());
  MouseGrabber::MouseGrabberPool().clear();
}
}}
