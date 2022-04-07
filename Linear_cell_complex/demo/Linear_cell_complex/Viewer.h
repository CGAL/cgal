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
  bool draw_face(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the edge containing dh is drawn.
  template<typename LCC>
  bool draw_edge(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the vertex containing dh is drawn.
  template<typename LCC>
  bool draw_vertex(const LCC&,
                   typename LCC::Dart_const_handle) const
  { return true; }

  /// @return true iff the volume containing dh is drawn in wireframe.
  template<typename LCC>
  bool volume_wireframe(const LCC& alcc,
                        typename LCC::Dart_const_handle dh) const
  { return !(alcc.template info<3>(dh).is_filled()); }
  /// @return true iff the face containing dh is drawn in wireframe.
  template<typename LCC>
  bool face_wireframe(const LCC&,
                        typename LCC::Dart_const_handle) const
  { return false; }

  /// @return true iff the volume containing dh is colored.
  template<typename LCC>
  bool colored_volume(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the face containing dh is colored.
  ///  if we have also colored_volume(alcc, dh), the volume color is
  ///  ignored and only the face color is considered.
  template<typename LCC>
  bool colored_face(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the edge containing dh is colored.
  template<typename LCC>
  bool colored_edge(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the vertex containing dh is colored.
  template<typename LCC>
  bool colored_vertex(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return false; }

  /// @return the color of the volume containing dh
  ///  used only if colored_volume(alcc, dh) and !colored_face(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color volume_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  { return alcc.template info<3>(dh).color(); }
  /// @return the color of the face containing dh
  ///  used only if colored_face(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color face_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the edge containing dh
  ///  used only if colored_edge(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color edge_color(const LCC&,
                         typename LCC::Dart_const_handle) const
  { return CGAL::IO::Color(0, 0, 0); }
  /// @return the color of the vertex containing dh
  ///  used only if colored_vertex(alcc, dh)
  template<typename LCC>
  CGAL::IO::Color vertex_color(const LCC&,
                           typename LCC::Dart_const_handle) const
  { return CGAL::IO::Color(0, 0, 0); }
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
