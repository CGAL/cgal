// Copyright (c) 2018 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H
#define CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H

#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Surface_mesh.h>
#include <CGAL/Random.h>

template<class SM>
class SimpleSurfaceMeshWithSmallFacesViewerQt : public CGAL::Basic_viewer_qt
{
  typedef CGAL::Basic_viewer_qt Base;
  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  typedef typename SM::Vertex_index vertex_descriptor;
  typedef typename SM::Face_index face_descriptor;
  typedef typename SM::Edge_index edge_descriptor;
  typedef typename SM::Halfedge_index halfedge_descriptor;
  typedef typename Kernel::FT FT;

public:
  /// Construct the viewer.
  /// @param amesh the surface mesh to view
  SimpleSurfaceMeshWithSmallFacesViewerQt(QWidget* parent,
                                          SM& amesh) :
    // First draw: no vertex; no edge, faces; multi-color; inverse normal
    Base(parent, "Surface mesh viewer with small faces", false, false, true, false, false),
    sm(amesh),
    m_threshold(85),
    m_draw_small_faces(true),
    m_draw_big_faces(true)
  {
    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_I, "Increment threshold for small faces");
    setKeyDescription(Qt::Key_D, "Decrement threshold for small faces");
    setKeyDescription(Qt::Key_S, "Draw small faces only , big faces only, both");

    if (sm.faces().begin()!=sm.faces().end())
    {
      bool exist;
      typename SM::template Property_map<face_descriptor, FT> faces_size;
      boost::tie(faces_size, exist)=sm.template property_map<face_descriptor, FT>("f:size");
      CGAL_assertion(exist);

      m_min_size=faces_size[*(sm.faces().begin())];
      m_max_size=m_min_size;
      FT cur_size;

      for (typename SM::Face_range::iterator f=sm.faces().begin(); f!=sm.faces().end(); ++f)
      {
        cur_size=faces_size[*f];
        if (cur_size<m_min_size) m_min_size=cur_size;
        if (cur_size>m_max_size) m_max_size=cur_size;
      }
    }

    compute_elements();
  }

protected:
  void compute_face(face_descriptor fh)
  {
    // [Face creation]
    bool issmall=false;

    // Default color of faces
    CGAL::IO::Color c(75,160,255);

    // Compare the size of the face with the % m_threshold
    bool exist;
    typename SM::template Property_map<face_descriptor, FT> faces_size;
    boost::tie(faces_size, exist)=sm.template property_map<face_descriptor, FT>("f:size");
    CGAL_assertion(exist);

    // It it is smaller, color the face in red.
    if (get(faces_size, fh)<m_min_size+((m_max_size-m_min_size)/(100-m_threshold)))
    {
      c=CGAL::IO::Color(255,20,20);
      issmall=true;
    }

    if ((issmall && !m_draw_small_faces) || (!issmall && !m_draw_big_faces))
    { return; }

    // Add the color of the face, then all its points.
    face_begin(c);
    halfedge_descriptor hd=sm.halfedge(fh);
    do
    {
      add_point_in_face(sm.point(sm.source(hd)), get_vertex_normal(hd));
      hd=sm.next(hd);
    }
    while(hd!=sm.halfedge(fh));
    face_end();
    /// [Face creation]
  }

  // Copy from draw_surface_mesh.h
  void compute_edge(edge_descriptor e)
  {
    /// [Edge creation]
    add_segment(sm.point(sm.source(sm.halfedge(e))),
                sm.point(sm.target(sm.halfedge(e))));
    /// [Edge creation]
  }

  void compute_vertex(vertex_descriptor vh)
  {
    /// [Vertex creation]
    add_point(sm.point(vh));
    /// [Vertex creation]
  }

  void compute_elements()
  {
    clear();

    for (typename SM::Face_range::iterator f=sm.faces().begin();
         f!=sm.faces().end(); ++f)
    {
      if (*f!=boost::graph_traits<SM>::null_face())
      { compute_face(*f); }
    }

    for (typename SM::Edge_range::iterator e=sm.edges().begin();
         e!=sm.edges().end(); ++e)
    { compute_edge(*e); }

    for (typename SM::Vertex_range::iterator v=sm.vertices().begin();
         v!=sm.vertices().end(); ++v)
    { compute_vertex(*v); }
  }

  // Call: * compute_elements() if the model changed, followed by
  //       * redraw() if some viewing parameters changed that implies some
  //                  modifications of the buffers
  //                  (eg. type of normal, color/mono)
  //       * update() just to update the drawing
  virtual void keyPressEvent(QKeyEvent *e)
  {
    /// [Keypress]
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    if ((e->key()==Qt::Key_I) && (modifiers==Qt::NoButton))
    {
      if (m_threshold<100) { ++m_threshold; }
      displayMessage(QString("Threshold percent=%1%.").arg(m_threshold));
      compute_elements();
      redraw();
    }
    else if ((e->key()==Qt::Key_D) && (modifiers==Qt::NoButton))
    {
      if (m_threshold>0) { --m_threshold; }
      displayMessage(QString("Threshold percent=%1%.").arg(m_threshold));
      compute_elements();
      redraw();
    }
    else if ((e->key()==Qt::Key_S) && (modifiers==Qt::NoButton))
    {
      QString msg;
      if (m_draw_small_faces)
      {
        if (m_draw_big_faces)
        {
          m_draw_big_faces=false;
          msg=QString("Draw small faces only.");
        }
        else
        {
          m_draw_big_faces=true; m_draw_small_faces=false;
          msg=QString("Draw big faces only.");
        }
      }
      else
      {
        assert(m_draw_big_faces);
        m_draw_small_faces=true;
        msg=QString("Draw small and big faces.");
      }

      displayMessage(msg);
      compute_elements();
      redraw();
    }
    else
    {
      // Call the base method to process others/classicals key
      Base::keyPressEvent(e);
    }
    /// [Keypress]
  }

protected:
  typename Kernel::Vector_3 get_face_normal(halfedge_descriptor he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    unsigned int nb=0;
    do
    {
      CGAL::internal::newell_single_step_3(sm.point(sm.source(he)),
                                           sm.point(sm.target(he)), normal);
      ++nb;
      he=sm.next(he);
    }
    while (he!=end);
    assert(nb>0);
    return (typename Kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }

  typename Kernel::Vector_3 get_vertex_normal(halfedge_descriptor he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    do
    {
      if (!sm.is_border(he))
      {
        typename Kernel::Vector_3 n=get_face_normal(he);
        normal=typename Kernel::Construct_sum_of_vectors_3()(normal, n);
      }
      he=sm.next(sm.opposite(he));
    }
    while (he!=end);

    if (!typename Kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
    { normal=(typename Kernel::Construct_scaled_vector_3()
              (normal, 1.0/CGAL::sqrt(normal.squared_length()))); }

    return normal;
  }

protected:
  SM& sm;
  unsigned int m_threshold;
  FT m_min_size, m_max_size;
  bool m_draw_small_faces;
  bool m_draw_big_faces;
};

template<class K>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>& amesh)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"surface_mesh_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleSurfaceMeshWithSmallFacesViewerQt<CGAL::Surface_mesh<K>>
      mainwindow(app.activeWindow(), amesh);
    mainwindow.show();
    app.exec();
  }
}

#else // CGAL_USE_BASIC_VIEWER

template<class K>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>&)
{
  std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
}

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H
