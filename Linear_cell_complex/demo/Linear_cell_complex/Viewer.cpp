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
    // TODO: add a new constructor that does not take graphic buffer.
    : Base(parent, ""), m_drawing_functor(MyDrawingFunctorLCC()),
      m_nofaces(false), m_random_face_color(false),
      m_previous_scene_empty(true) {}

void Viewer::setScene(Scene *scene_, bool doredraw) {
  scene = scene_;

  if (scene->lcc != nullptr) {
    compute_elements(gBuffer, scene->lcc, m_drawing_functor, m_nofaces,
                     m_random_face_color);
  }

  if (doredraw) {
    Base::redraw();
  }
}

void Viewer::sceneChanged() {
  compute_elements(gBuffer, scene->lcc, m_drawing_functor, m_nofaces,
                   m_random_face_color);

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
