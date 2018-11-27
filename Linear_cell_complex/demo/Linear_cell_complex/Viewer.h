// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//
#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <CGAL/draw_linear_cell_complex.h>

// Functor used by SimpleLCCViewerQt to colorize of not elements.
struct MyDrawingFunctorLCC
{
  /// @return true iff the volume containing dh is drawn.
  template<typename LCC>
  bool draw_volume(const LCC& alcc,
                   typename LCC::Dart_const_handle dh) const
  { return alcc.template info<3>(dh).is_visible(); }
  /// @return true iff the face containing dh is drawn.
  template<typename LCC>
  bool draw_face(const LCC& alcc,
                 typename LCC::Dart_const_handle dh) const
  { return true; }
  /// @return true iff the edge containing dh is drawn.
  template<typename LCC>
  bool draw_edge(const LCC& alcc,
                 typename LCC::Dart_const_handle dh) const
  { return true; }
  /// @return true iff the vertex containing dh is drawn.
  template<typename LCC>
  bool draw_vertex(const LCC& alcc,
                   typename LCC::Dart_const_handle dh) const
  { return true; }

  /// @return true iff the volume containing dh is drawn in wireframe.
  template<typename LCC>
  bool volume_wireframe(const LCC& alcc,
                        typename LCC::Dart_const_handle dh) const
  { return false; }
  /// @return true iff the face containing dh is drawn in wireframe.
  template<typename LCC>
  bool face_wireframe(const LCC& alcc,
                        typename LCC::Dart_const_handle dh) const
  { return false; }

  /// @return true iff the volume containing dh is colored.
  template<typename LCC>
  bool colored_volume(const LCC& alcc,
                      typename LCC::Dart_const_handle dh) const
  { return true; }
  /// @return true iff the face containing dh is colored.
  ///  if we have also colored_volume(alcc, dh), the volume color is
  ///  ignored and only the face color is considered.
  template<typename LCC>
  bool colored_face(const LCC& alcc,
                    typename LCC::Dart_const_handle dh) const
  { return false; }
  /// @return true iff the edge containing dh is colored.
  template<typename LCC>
  bool colored_edge(const LCC& alcc,
                    typename LCC::Dart_const_handle dh) const
  { return false; }
  /// @return true iff the vertex containing dh is colored.
  template<typename LCC>
  bool colored_vertex(const LCC& alcc,
                      typename LCC::Dart_const_handle dh) const
  { return false; }

  /// @return the color of the volume containing dh
  ///  used only if colored_volume(alcc, dh) and !colored_face(alcc, dh)
  template<typename LCC>
  CGAL::Color volume_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  { return alcc.template info<3>(dh).color(); }
  /// @return the color of the face containing dh
  ///  used only if colored_face(alcc, dh)
  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the edge containing dh
  ///  used only if colored_edge(alcc, dh)
  template<typename LCC>
  CGAL::Color edge_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  { return CGAL::Color(0, 0, 0); }
  /// @return the color of the vertex containing dh
  ///  used only if colored_vertex(alcc, dh)
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  { return CGAL::Color(0, 0, 0); }
};


class Viewer : public CGAL::SimpleLCCViewerQt<LCC, MyDrawingFunctorLCC>
{
  Q_OBJECT

  typedef CGAL::SimpleLCCViewerQt<LCC, MyDrawingFunctorLCC> Base;

public:
  Viewer(QWidget* parent);
  void setScene(Scene* scene_, bool doredraw=true);
  void keyPressEvent(QKeyEvent *e);
  virtual QString helpString() const;

public Q_SLOTS:
  void sceneChanged();

private:
  Scene* scene;
  bool m_previous_scene_empty;
};

#endif
