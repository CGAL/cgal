// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//

#include "Viewer.h"
#include <CGAL/Qt/vec.h>

Viewer::Viewer(QWidget* parent) :
    Base(parent, NULL, ""),
    m_previous_scene_empty(true)
{}

void Viewer::setScene(Scene* scene_, bool doredraw)
{
  scene = scene_;
  set_lcc(scene->lcc, doredraw);
}

void Viewer::sceneChanged()
{
  Base::compute_elements();
  this->camera()->
      setSceneBoundingBox(CGAL::qglviewer::Vec(m_bounding_box.xmin(),
                                               m_bounding_box.ymin(),
                                               m_bounding_box.zmin()),
                          CGAL::qglviewer::Vec(m_bounding_box.xmax(),
                                               m_bounding_box.ymax(),
                                               m_bounding_box.zmax()));
  Base::redraw();
  if (m_previous_scene_empty)
  { this->showEntireScene(); }

  m_previous_scene_empty = scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
  // const Qt::KeyboardModifiers modifiers = e->modifiers();
  Base::keyPressEvent(e);
}

QString Viewer::helpString() const
{ return Base::helpString("LCC Demo"); }

