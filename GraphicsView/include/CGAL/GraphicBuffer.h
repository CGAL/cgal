// Copyright (c) 2022 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_GRAPHIC_BUFFER_H
#define CGAL_GRAPHIC_BUFFER_H

#include <CGAL/license/GraphicsView.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>

#include <cstdlib>
#include <map>
#include <queue>

#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

namespace CGAL {

// This class is responsible to deal with CGAL data structures,
// filling mesh, and handling buffers.
template <typename BufferType = float> class GraphicBuffer {

public:
  GraphicBuffer(std::vector<BufferType> (&pos)[20], CGAL::Bbox_3 &bbox)
      : m_buffer_for_mono_points(&pos[POS_MONO_POINTS], nullptr, &bbox, nullptr,
                                 nullptr, nullptr),
        m_buffer_for_colored_points(&pos[POS_COLORED_POINTS], nullptr, &bbox,
                                    &pos[COLOR_POINTS], nullptr, nullptr),
        m_buffer_for_mono_segments(&pos[POS_MONO_SEGMENTS], nullptr, &bbox,
                                   nullptr, nullptr, nullptr),
        m_buffer_for_colored_segments(&pos[POS_COLORED_SEGMENTS], nullptr,
                                      &bbox, &pos[COLOR_SEGMENTS], nullptr,
                                      nullptr),
        m_buffer_for_mono_rays(&pos[POS_MONO_RAYS], nullptr, &bbox, nullptr,
                               nullptr),
        m_buffer_for_colored_rays(&pos[POS_COLORED_RAYS], nullptr, &bbox,
                                  &pos[COLOR_RAYS], nullptr, nullptr),
        m_buffer_for_mono_lines(&pos[POS_MONO_RAYS], nullptr, &bbox, nullptr,
                                nullptr),
        m_buffer_for_colored_lines(&pos[POS_COLORED_LINES], nullptr, &bbox,
                                   &pos[COLOR_LINES], nullptr, nullptr),
        m_buffer_for_mono_faces(&pos[POS_MONO_FACES], nullptr, &bbox, nullptr,
                                &pos[FLAT_NORMAL_MONO_FACES],
                                &pos[SMOOTH_NORMAL_MONO_FACES]),
        m_buffer_for_colored_faces(&pos[POS_COLORED_FACES], nullptr, &bbox,
                                   &pos[COLOR_FACES],
                                   &pos[FLAT_NORMAL_COLORED_FACES],
                                   &pos[SMOOTH_NORMAL_COLORED_FACES]) {}

  const Buffer_for_vao<BufferType> &get_buffer_for_mono_points() const {
    return m_buffer_for_mono_points;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_colored_points() const {
    return m_buffer_for_colored_points;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_mono_segments() const {
    return m_buffer_for_mono_segments;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_colored_segments() const {
    return m_buffer_for_colored_segments;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_mono_rays() const {
    return m_buffer_for_mono_rays;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_colored_rays() const {
    return m_buffer_for_colored_rays;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_mono_lines() const {
    return m_buffer_for_mono_lines;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_colored_lines() const {
    return m_buffer_for_colored_lines;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_mono_faces() const {
    return m_buffer_for_mono_faces;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_colored_faces() const {
    return m_buffer_for_colored_faces;
  }

  const Buffer_for_vao<BufferType> &get_buffer_for_clipping_plane() const {
    return m_buffer_for_clipping_plane;
  }

  template <typename KPoint> void add_point(const KPoint &p) {
    m_buffer_for_mono_points.add_point(p);
  }

  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor) {
    m_buffer_for_colored_points.add_point(p, acolor);
  }

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2) {
    m_buffer_for_mono_segments.add_segment(p1, p2);
  }

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor) {
    m_buffer_for_colored_segments.add_segment(p1, p2, acolor);
  }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v) {
    double bigNumber = 1e30;
    m_buffer_for_mono_rays.add_ray_segment(p, (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor) {
    double bigNumber = 1e30;
    m_buffer_for_colored_rays.add_ray_segment(p, (p + (bigNumber)*v), acolor);
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v) {
    double bigNumber = 1e30;
    m_buffer_for_mono_lines.add_line_segment((p - (bigNumber)*v),
                                             (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor) {
    double bigNumber = 1e30;
    m_buffer_for_colored_lines.add_line_segment((p - (bigNumber)*v),
                                                (p + (bigNumber)*v), acolor);
  }

  template <typename KPoint> bool add_point_in_face(const KPoint &kp) {
    if (m_buffer_for_mono_faces.is_a_face_started()) {
      return m_buffer_for_mono_faces.add_point_in_face(kp);
    } else if (m_buffer_for_colored_faces.is_a_face_started()) {
      return m_buffer_for_colored_faces.add_point_in_face(kp);
    }
    return false;
  }

  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal) {
    if (m_buffer_for_mono_faces.is_a_face_started()) {
      return m_buffer_for_mono_faces.add_point_in_face(kp, p_normal);
    } else if (m_buffer_for_colored_faces.is_a_face_started()) {
      return m_buffer_for_colored_faces.add_point_in_face(kp, p_normal);
    }
    return false;
  }

protected:
  // The following enum gives the indices of different elements of arrays
  // vectors.
  enum {
    BEGIN_POS = 0,
    POS_MONO_POINTS = BEGIN_POS,
    POS_COLORED_POINTS,
    POS_MONO_SEGMENTS,
    POS_COLORED_SEGMENTS,
    POS_MONO_RAYS,
    POS_COLORED_RAYS,
    POS_MONO_LINES,
    POS_COLORED_LINES,
    POS_MONO_FACES,
    POS_COLORED_FACES,
    POS_CLIPPING_PLANE,
    END_POS,
    BEGIN_COLOR = END_POS,
    COLOR_POINTS = BEGIN_COLOR,
    COLOR_SEGMENTS,
    COLOR_RAYS,
    COLOR_LINES,
    COLOR_FACES,
    END_COLOR,
    BEGIN_NORMAL = END_COLOR,
    SMOOTH_NORMAL_MONO_FACES = BEGIN_NORMAL,
    FLAT_NORMAL_MONO_FACES,
    SMOOTH_NORMAL_COLORED_FACES,
    FLAT_NORMAL_COLORED_FACES,
    END_NORMAL,
    LAST_INDEX = END_NORMAL
  };

  Buffer_for_vao<BufferType> m_buffer_for_mono_points;
  Buffer_for_vao<BufferType> m_buffer_for_colored_points;
  Buffer_for_vao<BufferType> m_buffer_for_mono_segments;
  Buffer_for_vao<BufferType> m_buffer_for_colored_segments;
  Buffer_for_vao<BufferType> m_buffer_for_mono_rays;
  Buffer_for_vao<BufferType> m_buffer_for_colored_rays;
  Buffer_for_vao<BufferType> m_buffer_for_mono_lines;
  Buffer_for_vao<BufferType> m_buffer_for_colored_lines;
  Buffer_for_vao<BufferType> m_buffer_for_mono_faces;
  Buffer_for_vao<BufferType> m_buffer_for_colored_faces;
  Buffer_for_vao<BufferType> m_buffer_for_clipping_plane;
};

} // namespace CGAL

#endif // CGAL_GRAPHIC_BUFFER_H
