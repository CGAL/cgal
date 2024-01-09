// Copyright (c) 2018-2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawSurfaceMesh

Open a new window and draw `asm`, an instance of the `CGAL::Surface_mesh` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam SM an instance of the `CGAL::Surface_mesh` class.
\param asm the surface mesh to draw.

*/
template<class SM>
void draw(const SM& asm);

#else // DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/use.h>

namespace CGAL {

// Check if there are any color maps that could be used, random otherwise
template <typename K>
struct Surface_mesh_basic_viewer_color_map
  : DefaultColorFunctorFaceGraph
{
  using Base = DefaultColorFunctorFaceGraph;

  using SM = ::CGAL::Surface_mesh<K>;
  using vertex_descriptor = typename boost::graph_traits<SM>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<SM>::edge_descriptor;
  using face_descriptor = typename boost::graph_traits<SM>::face_descriptor;

  using vertex_colors = typename SM::template Property_map<vertex_descriptor, CGAL::IO::Color>;
  using edge_colors = typename SM::template Property_map<edge_descriptor, CGAL::IO::Color>;
  using face_colors = typename SM::template Property_map<face_descriptor, CGAL::IO::Color>;

  Surface_mesh_basic_viewer_color_map(const SM& amesh)
  {
    bool found = false;
    std::tie(vcolors, found) = amesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");
    std::tie(ecolors, found) = amesh.template property_map<edge_descriptor, CGAL::IO::Color>("e:color");
    std::tie(fcolors, found) = amesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");
    CGAL_USE(found);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const vertex_descriptor v) const
  {
    return vcolors ? get(vcolors, v) : Base::operator()(amesh, v);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const edge_descriptor e) const
  {
    return ecolors ? get(ecolors, e) : Base::operator()(amesh, e);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const face_descriptor f) const
  {
    return fcolors ? get(fcolors, f) : Base::operator()(amesh, f);
  }

private:
  vertex_colors vcolors;
  edge_colors ecolors;
  face_colors fcolors;
};

// Specialization of draw function.
template<class K>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer",
          bool nofill=false)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    CGAL::Qt::init_ogl_context(4,3);
    int argc=1;
    const char* argv[2]={"surface_mesh_viewer", nullptr};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleFaceGraphViewerQt mainwindow(app.activeWindow(), amesh, title, nofill,
                                       Surface_mesh_basic_viewer_color_map<K>(amesh));
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
