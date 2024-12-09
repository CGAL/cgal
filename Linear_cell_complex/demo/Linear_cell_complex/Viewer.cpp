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
    m_previous_scene_empty(true)
{
  m_gs_options.face_color=[](const LCC & alcc,
                             Dart_const_descriptor dh)->CGAL::IO::Color
  {
    if(alcc.template is_free<3>(dh))
    { return alcc.template info<3>(dh).color(); }

    if(!alcc.template info<3>(dh).is_visible() ||
       !alcc.template info<3>(dh).is_filled())
    { return alcc.template info<3>(alcc.template beta<3>(dh)).color(); }

    if(!alcc.template info<3>(alcc.template beta<3>(dh)).is_visible() ||
       !alcc.template info<3>(alcc.template beta<3>(dh)).is_filled())
    { return alcc.template info<3>(dh).color(); }

    const CGAL::IO::Color& c1=alcc.template info<3>(dh).color();
    const CGAL::IO::Color& c2=alcc.template info<3>(alcc.template beta<3>(dh)).color();
    return CGAL::IO::Color((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2);
  };

  m_gs_options.colored_face=[](const LCC &, Dart_const_descriptor)->bool
  { return true; };

  m_gs_options.draw_volume=[](const LCC & alcc, Dart_const_descriptor dh)->bool
  { return alcc.template info<3>(dh).is_visible(); };

  m_gs_options.volume_wireframe=[](const LCC& alcc, Dart_const_descriptor dh)->bool
  { return !(alcc.template info<3>(dh).is_filled()); };
}

void Viewer::setScene(Scene *scene_, bool doredraw)
{
  scene = scene_;

  if (scene->lcc!=nullptr)
  { CGAL::add_to_graphics_scene(*scene->lcc, m_graphic_buffer, m_gs_options); }

  if (doredraw)
  { Base::redraw(); }
}

void Viewer::sceneChanged()
{
  m_graphic_buffer.clear();
  CGAL::add_to_graphics_scene(*scene->lcc, m_graphic_buffer, m_gs_options);

  this->camera()->setSceneBoundingBox(
      CGAL::qglviewer::Vec(gBuffer.bounding_box().xmin(),
                           gBuffer.bounding_box().ymin(),
                           gBuffer.bounding_box().zmin()),
      CGAL::qglviewer::Vec(gBuffer.bounding_box().xmax(),
                           gBuffer.bounding_box().ymax(),
                           gBuffer.bounding_box().zmax()));
  Base::redraw();
  if (m_previous_scene_empty)
  { this->showEntireScene(); }

  m_previous_scene_empty=scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
  // const Qt::KeyboardModifiers modifiers = e->modifiers();
  Base::keyPressEvent(e);
}

QString Viewer::helpString() const
{ return Base::helpString("LCC Demo"); }
