// Copyright (c) 2022 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//            Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_GRAPHICS_SCENE_H
#define CGAL_GRAPHICS_SCENE_H

// TODO #include <CGAL/license/GraphicsView.h>

#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <CGAL/Random.h>

#include <cstdlib>

#include <CGAL/Buffer_for_vao.h>

namespace CGAL {

//------------------------------------------------------------------------------
inline CGAL::IO::Color get_random_color(CGAL::Random& random)
{
  CGAL::IO::Color res;
  do
  {
    res=CGAL::IO::Color(random.get_int(0,256),
                        random.get_int(0,256),
                        random.get_int(0,256));
  }
  while(res.red()==255 && res.green()==255 && res.blue()==255);
  return res;
}
//------------------------------------------------------------------------------
// This class is responsible for dealing with available CGAL data structures and
// handling buffers.
class Graphics_scene
{
public:
  using BufferType=float;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef Local_kernel::Point_3 Local_point;
  typedef Local_kernel::Vector_3 Local_vector;

  Graphics_scene()
      : m_buffer_for_mono_points(&arrays[POS_MONO_POINTS], nullptr,
                                 &m_bounding_box, nullptr, nullptr, nullptr),
        m_buffer_for_colored_points(&arrays[POS_COLORED_POINTS], nullptr,
                                    &m_bounding_box, &arrays[COLOR_POINTS],
                                    nullptr, nullptr),
        m_buffer_for_mono_segments(&arrays[POS_MONO_SEGMENTS], nullptr,
                                   &m_bounding_box, nullptr, nullptr, nullptr),
        m_buffer_for_colored_segments(&arrays[POS_COLORED_SEGMENTS], nullptr,
                                      &m_bounding_box, &arrays[COLOR_SEGMENTS],
                                      nullptr, nullptr),
        m_buffer_for_mono_rays(&arrays[POS_MONO_RAYS], nullptr, &m_bounding_box,
                               nullptr, nullptr),
        m_buffer_for_colored_rays(&arrays[POS_COLORED_RAYS], nullptr,
                                  &m_bounding_box, &arrays[COLOR_RAYS], nullptr,
                                  nullptr),
        m_buffer_for_mono_lines(&arrays[POS_MONO_RAYS], nullptr,
                                &m_bounding_box, nullptr, nullptr),
        m_buffer_for_colored_lines(&arrays[POS_COLORED_LINES], nullptr,
                                   &m_bounding_box, &arrays[COLOR_LINES],
                                   nullptr, nullptr),
        m_buffer_for_mono_faces(
            &arrays[POS_MONO_FACES], nullptr, &m_bounding_box, nullptr,
            &arrays[FLAT_NORMAL_MONO_FACES], &arrays[SMOOTH_NORMAL_MONO_FACES]),
        m_buffer_for_colored_faces(&arrays[POS_COLORED_FACES], nullptr,
                                   &m_bounding_box, &arrays[COLOR_FACES],
                                   &arrays[FLAT_NORMAL_COLORED_FACES],
                                   &arrays[SMOOTH_NORMAL_COLORED_FACES])
  {}

  const Buffer_for_vao &get_buffer_for_mono_points() const
  { return m_buffer_for_mono_points; }

  const Buffer_for_vao &get_buffer_for_colored_points() const
  { return m_buffer_for_colored_points; }

  const Buffer_for_vao &get_buffer_for_mono_segments() const
  { return m_buffer_for_mono_segments; }

  const Buffer_for_vao &get_buffer_for_colored_segments() const
  { return m_buffer_for_colored_segments; }

  const Buffer_for_vao &get_buffer_for_mono_rays() const
  { return m_buffer_for_mono_rays; }

  const Buffer_for_vao &get_buffer_for_colored_rays() const
  { return m_buffer_for_colored_rays; }

  const Buffer_for_vao &get_buffer_for_mono_lines() const
  { return m_buffer_for_mono_lines; }

  const Buffer_for_vao &get_buffer_for_colored_lines() const
  { return m_buffer_for_colored_lines; }

  const Buffer_for_vao &get_buffer_for_mono_faces() const
  { return m_buffer_for_mono_faces; }

  const Buffer_for_vao &get_buffer_for_colored_faces() const
  { return m_buffer_for_colored_faces; }

  const CGAL::Bbox_3 &bounding_box() const { return m_bounding_box; }

  const std::vector<BufferType> &get_array_of_index(int index) const
  { assert(index<LAST_INDEX); return arrays[index]; }

  int get_size_of_index(int index) const
  { return static_cast<int>(arrays[index].size()*sizeof(BufferType)); }

  unsigned int number_of_elements(int index) const
  { return static_cast<unsigned int>(arrays[index].size()/3); }

  void initiate_bounding_box(const CGAL::Bbox_3& new_bounding_box)
  { m_bounding_box = new_bounding_box; }

  void update_bounding_box(const CGAL::Bbox_3 &box) { m_bounding_box+=box; }

  template <typename KPoint, typename KVector>
  void update_bounding_box_for_ray(const KPoint &p, const KVector &v)
  {
    Local_point lp = get_local_point(p);
    Local_vector lv = get_local_vector(v);
    update_bounding_box((lp + lv).bbox());
  }

  template <typename KPoint, typename KVector>
  void update_bounding_box_for_line(const KPoint &p, const KVector &v,
                                    const KVector &pv)
  {
    Local_point lp = get_local_point(p);
    Local_vector lv = get_local_vector(v);
    Local_vector lpv = get_local_vector(pv);
    update_bounding_box(lp.bbox() + (lp + lv).bbox() + (lp + lpv).bbox());
  }

  void reverse_all_normals() const
  {
    m_buffer_for_mono_faces.negate_normals();
    m_buffer_for_colored_faces.negate_normals();
  }

  template <typename KPoint> void add_point(const KPoint &p)
  { m_buffer_for_mono_points.add_point(p); }

  template <typename KPoint>
  void add_point(const KPoint &p, const CGAL::IO::Color &acolor)
  { m_buffer_for_colored_points.add_point(p, acolor); }

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2)
  { m_buffer_for_mono_segments.add_segment(p1, p2); }

  template <typename KPoint>
  void add_segment(const KPoint &p1, const KPoint &p2,
                   const CGAL::IO::Color &acolor)
  { m_buffer_for_colored_segments.add_segment(p1, p2, acolor); }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v)
  {
    double bigNumber = 1e30;
    m_buffer_for_mono_rays.add_ray_segment(p, (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v,
               const CGAL::IO::Color &acolor)
  {
    double bigNumber = 1e30;
    m_buffer_for_colored_rays.add_ray_segment(p, (p + (bigNumber)*v), acolor);
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v)
  {
    double bigNumber = 1e30;
    m_buffer_for_mono_lines.add_line_segment((p - (bigNumber)*v),
                                             (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v,
                const CGAL::IO::Color &acolor)
  {
    double bigNumber = 1e30;
    m_buffer_for_colored_lines.add_line_segment((p - (bigNumber)*v),
                                                (p + (bigNumber)*v), acolor);
  }

  template <typename KPoint> bool add_point_in_face(const KPoint &kp)
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { return m_buffer_for_mono_faces.add_point_in_face(kp); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { return m_buffer_for_colored_faces.add_point_in_face(kp); }
    return false;
  }

  template <typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint &kp, const KVector &p_normal)
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { return m_buffer_for_mono_faces.add_point_in_face(kp, p_normal); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { return m_buffer_for_colored_faces.add_point_in_face(kp, p_normal); }
    return false;
  }

  bool a_face_started() const
  {
    return m_buffer_for_mono_faces.is_a_face_started() ||
           m_buffer_for_colored_faces.is_a_face_started();
  }

  void face_begin()
  {
    if (a_face_started())
    {
      std::cerr
          << "You cannot start a new face before to finish the previous one."
          << std::endl;
    }
    else
    { m_buffer_for_mono_faces.face_begin(); }
  }

  void face_begin(const CGAL::IO::Color &acolor)
  {
    if (a_face_started())
    {
      std::cerr
          << "You cannot start a new face before to finish the previous one."
          << std::endl;
    }
    else
    { m_buffer_for_colored_faces.face_begin(acolor); }
  }

  void face_end()
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { m_buffer_for_mono_faces.face_end(); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { m_buffer_for_colored_faces.face_end(); }
  }

  template <typename KPoint>
  void add_text(const KPoint &kp, const std::string &txt)
  {
    Local_point p = get_local_point(kp);
    m_texts.push_back(std::make_tuple(p, txt));
  }

  template <typename KPoint>
  void add_text(const KPoint &kp, const char *txt)
  { add_text(kp, std::string(txt)); }

  bool empty() const
  {
    return (m_buffer_for_mono_points.is_empty() &&
            m_buffer_for_colored_points.is_empty() &&
            m_buffer_for_mono_segments.is_empty() &&
            m_buffer_for_colored_segments.is_empty() &&
            m_buffer_for_mono_rays.is_empty() &&
            m_buffer_for_colored_rays.is_empty() &&
            m_buffer_for_mono_lines.is_empty() &&
            m_buffer_for_colored_lines.is_empty() &&
            m_buffer_for_mono_faces.is_empty() &&
            m_buffer_for_colored_faces.is_empty());
  }

  bool has_zero_x() const
  {
    return m_buffer_for_mono_points.has_zero_x() &&
           m_buffer_for_colored_points.has_zero_x() &&
           m_buffer_for_mono_segments.has_zero_x() &&
           m_buffer_for_colored_segments.has_zero_x() &&
           m_buffer_for_mono_faces.has_zero_x() &&
           m_buffer_for_colored_faces.has_zero_x() &&
           m_buffer_for_mono_rays.has_zero_x() &&
           m_buffer_for_colored_rays.has_zero_x() &&
           m_buffer_for_mono_lines.has_zero_x() &&
           m_buffer_for_colored_lines.has_zero_x();
  }

  bool has_zero_y() const
  {
    return m_buffer_for_mono_points.has_zero_y() &&
           m_buffer_for_colored_points.has_zero_y() &&
           m_buffer_for_mono_segments.has_zero_y() &&
           m_buffer_for_colored_segments.has_zero_y() &&
           m_buffer_for_mono_faces.has_zero_y() &&
           m_buffer_for_colored_faces.has_zero_y() &&
           m_buffer_for_mono_rays.has_zero_y() &&
           m_buffer_for_colored_rays.has_zero_y() &&
           m_buffer_for_mono_lines.has_zero_y() &&
           m_buffer_for_colored_lines.has_zero_y();
  }

  bool has_zero_z() const
  {
    return m_buffer_for_mono_points.has_zero_z() &&
           m_buffer_for_colored_points.has_zero_z() &&
           m_buffer_for_mono_segments.has_zero_z() &&
           m_buffer_for_colored_segments.has_zero_z() &&
           m_buffer_for_mono_faces.has_zero_z() &&
           m_buffer_for_colored_faces.has_zero_z() &&
           m_buffer_for_mono_rays.has_zero_z() &&
           m_buffer_for_colored_rays.has_zero_z() &&
           m_buffer_for_mono_lines.has_zero_z() &&
           m_buffer_for_colored_lines.has_zero_z();
  }

  // Returns true if the data structure lies on a XY or XZ or YZ plane
  bool is_two_dimensional() const
  {
    return (!empty() && (has_zero_x() || has_zero_y() || has_zero_z()));
  }

  void clear()
  {
    m_buffer_for_mono_points.clear();
    m_buffer_for_colored_points.clear();
    m_buffer_for_mono_segments.clear();
    m_buffer_for_colored_segments.clear();
    m_buffer_for_mono_rays.clear();
    m_buffer_for_colored_rays.clear();
    m_buffer_for_mono_lines.clear();
    m_buffer_for_colored_lines.clear();
    m_buffer_for_mono_faces.clear();
    m_buffer_for_colored_faces.clear();
    m_texts.clear();
    m_bounding_box=CGAL::Bbox_3();
  }

  template <typename KPoint>
  static Local_point get_local_point(const KPoint &p)
  {
    return internal::Geom_utils<typename CGAL::Kernel_traits<KPoint>::Kernel,
                                Local_kernel>::get_local_point(p);
  }

  template <typename KVector>
  static Local_vector get_local_vector(const KVector &v)
  {
    return internal::Geom_utils<typename CGAL::Kernel_traits<KVector>::Kernel,
                                Local_kernel>::get_local_vector(v);
  }

  void m_texts_clear()
  { m_texts.clear(); }

  std::size_t m_texts_size() const
  { return m_texts.size(); }

  const std::vector<std::tuple<Local_point, std::string>>& get_m_texts() const
  { return m_texts; }

public:
  // The following enum gives the indices of different elements of arrays
  // vectors.
  enum Buffers {
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

protected:
  Buffer_for_vao m_buffer_for_mono_points;
  Buffer_for_vao m_buffer_for_colored_points;
  Buffer_for_vao m_buffer_for_mono_segments;
  Buffer_for_vao m_buffer_for_colored_segments;
  Buffer_for_vao m_buffer_for_mono_rays;
  Buffer_for_vao m_buffer_for_colored_rays;
  Buffer_for_vao m_buffer_for_mono_lines;
  Buffer_for_vao m_buffer_for_colored_lines;
  Buffer_for_vao m_buffer_for_mono_faces;
  Buffer_for_vao m_buffer_for_colored_faces;

  std::vector<std::tuple<Local_point, std::string>> m_texts;

  std::vector<BufferType> arrays[LAST_INDEX];

  CGAL::Bbox_3 m_bounding_box;
};

} // namespace CGAL

#endif // CGAL_GRAPHICS_SCENE_H
