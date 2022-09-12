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
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>
//

#include "Viewer.h"
#include <CGAL/Qt/vec.h>

Viewer::Viewer(QWidget *parent)
  : Base(parent, m_graphic_buffer, ""), 
    m_nofaces_(false), m_random_face_color_(false),
    m_previous_scene_empty(true) {}

void Viewer::setScene(Scene *scene_, bool doredraw) {
  scene = scene_;

  if (scene->lcc != nullptr) {
    CGAL::add_in_graphic_buffer_lcc(gBuffer, m_drawing_functor_, scene->lcc, m_nofaces_,
                              m_random_face_color_);
  }

  if (doredraw) {
    Base::redraw();
  }
}

void Viewer::sceneChanged() {
  gBuffer.clear();
  CGAL::add_in_graphic_buffer_lcc(gBuffer, m_drawing_functor_, scene->lcc,
                                  m_nofaces_, m_random_face_color_);

  this->camera()->setSceneBoundingBox(
      CGAL::qglviewer::Vec(gBuffer.get_bounding_box().xmin(),
                           gBuffer.get_bounding_box().ymin(),
                           gBuffer.get_bounding_box().zmin()),
      CGAL::qglviewer::Vec(gBuffer.get_bounding_box().xmax(),
                           gBuffer.get_bounding_box().ymax(),
                           gBuffer.get_bounding_box().zmax()));
  Base::redraw();
  if (m_previous_scene_empty) {
    this->showEntireScene();
  }

  m_previous_scene_empty =
      scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::keyPressEvent(QKeyEvent *e) {
  // const Qt::KeyboardModifiers modifiers = e->modifiers();
  Base::keyPressEvent(e);
}

QString Viewer::helpString() const { return Base::helpString("LCC Demo"); }
