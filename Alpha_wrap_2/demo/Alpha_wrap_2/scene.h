// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_DEMO_SCENE_H
#define CGAL_ALPHA_WRAP_2_DEMO_SCENE_H

#include "types.h"

#include <QtOpenGL>
#include <QString>

#include <fstream>
#include <algorithm>

struct Stats_data_structure {
  double nb_proj = 0.;
  double nb_steiner = 0.;
  int nb_gate_traversed = 0;
};

class Scene
{
private:
  Pslg m_input;
  Pslg m_wrap;

  double m_alpha = 20.;
  double m_offset = 600.;
  bool is_alpha_relative = true;
  bool is_offset_relative = true;
  bool is_step_by_step = false;

  Alpha_wrapper m_wrapper;
  Stats_data_structure m_stats;

  Point_2 m_mouse_pos;

  bool view_input = true;
  bool view_alpha_wrap = true;
  bool view_dt2_inside_outside = true;
  bool view_dt2_edge = true;
  bool view_voronoi = true;
  bool view_empty_alpha_pencils = false;
  bool view_steiner_point = true;
  bool view_next_gate = true;
  bool view_next_gate_pencil = false;

  QString screenshot_folder = "";
  QString screenshot_filename = "";
  int screenshot_number = 0;

public:
  Scene() { }

  ~Scene() { clear(); }

  void toggle_view_input() { view_input = !view_input; }

  void toggle_view_alpha_wrap() { view_alpha_wrap = !view_alpha_wrap; }

  void toggle_view_dt2_inside_outside() { view_dt2_inside_outside = !view_dt2_inside_outside; }

  void toggle_view_dt2_edge() { view_dt2_edge = !view_dt2_edge; }

  void toggle_view_voronoi() { view_voronoi = !view_voronoi; }

  void toggle_view_empty_alpha_pencils() { view_empty_alpha_pencils = !view_empty_alpha_pencils; }

  void toggle_view_steiner_point() { view_steiner_point = !view_steiner_point; }

  void toggle_view_next_gate() { view_next_gate = !view_next_gate; }

  void toggle_view_next_gate_pencil() { view_next_gate_pencil = !view_next_gate_pencil; }

  void clear()
  {
    m_input.clear();
    m_wrap.clear();
    m_wrapper.clear();
    m_stats = {};
  }

  const Pslg& get_input() const {
    return m_input;
  }

  Pslg& get_alpha_wrap() {
    return m_wrap;
  }

  const Alpha_wrapper& get_wrapper() const {
    return m_wrapper;
  }

  Alpha_wrapper& get_wrapper() {
    return m_wrapper;
  }

  double get_alpha() const {
    return m_alpha;
  }

  void set_alpha(const double a) {
    m_alpha = a;
  }

  double get_offset() const {
    return m_offset;
  }

  void set_offset(const double o) {
    m_offset = o;
  }

  bool get_is_alpha_relative() const {
    return is_alpha_relative;
  }

  void set_is_alpha_relative(const bool b) {
    is_alpha_relative = b;
  }

  bool get_is_offset_relative() const {
    return is_offset_relative;
  }

  void set_is_offset_relative(const bool b) {
    is_offset_relative = b;
  }

  bool get_is_step_by_step() const {
    return is_step_by_step;
  }

  void set_is_step_by_step(const bool b) {
    is_step_by_step = b;
  }

  QString get_screenshot_folder() const {
    return screenshot_folder;
  }

  void set_screenshot_folder(QString s) {
    screenshot_folder = s;
  }

  QString get_screenshot_filename() const {
    return screenshot_filename;
  }

  void set_screenshot_filename(const QString& s) {
    screenshot_filename = s;
  }

  int get_screenshot_number() const {
    return screenshot_number;
  }

  void set_screenshot_number(const int i) {
    screenshot_number = i;
  }

  double get_nb_proj() const {
    return m_stats.nb_proj;
  }

  double get_nb_steiner() const {
    return m_stats.nb_steiner;
  }

  int get_nb_gate_traversed() const {
    return m_stats.nb_gate_traversed;
  }

  void reset_stats() {
    m_stats = Stats_data_structure();
  }

  void set_mouse_pos(const Point_2& pos) {  m_mouse_pos = pos; }

  bool load(const QString& filename)
  {
    Pslg pslg;
    AW2::IO::IO_exit_code success_read = AW2::IO::read_input_polylines_file(filename.toUtf8().constData(), pslg);
    if(success_read == AW2::IO::UNREADABLE_INPUT) {
      std::cerr << "Warning: the input is unreadable \n";
      return false;
    }

    if(success_read == AW2::IO::INPUT_IS_EMPTY) {
      std::cout << "The input is empty \n";
      return false;
    }

    clear();
    m_input = std::move(pslg);

    return true;
  }

  void render()
  {
    if(view_input)
      gl_draw_pslg(m_input, true, 250, 0, 0);
    if(view_alpha_wrap)
      gl_draw_pslg(m_wrap, false, 0, 0, 250);
    if(view_dt2_inside_outside)
      gl_draw_dt2_inside_outside();
    if(view_dt2_edge)
      gl_draw_dt2_edge();
    if(view_voronoi)
      gl_draw_voronoi();
    if(view_empty_alpha_pencils)
      gl_draw_all_pencils();
    if(view_steiner_point)
      gl_draw_steiner_point();
    if(view_next_gate)
      gl_draw_next_gate();
    if(view_next_gate_pencil) {
      gl_draw_next_gate_pencil();
    }
  }

  void gl_draw_pslg(const Pslg& pslg,
                    bool /*force_close*/,
                    const unsigned char r,
                    const unsigned char g,
                    const unsigned char b)
  {
    if(pslg.empty())
      return;

    ::glLineWidth(3.0f);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);
    for(const auto& cmp : pslg) {
      if(cmp.size() < 2)
        continue;
      for(size_t i = 1; i < cmp.size(); ++i) {
        const Point_2& start = cmp[i-1];
        const Point_2& end = cmp[i];
        ::glVertex2d(start.x(),start.y());
        ::glVertex2d(end.x(),end.y());
      }
    }
    ::glEnd();
  }

  void gl_draw_dt2_inside_outside() const
  {
    const Triangulation& tr = m_wrapper.triangulation();

    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glColor4f(0.25,0.25,0.25,0.8);
    ::glBegin(GL_TRIANGLES);
    for(const Face_handle& fh2 : tr.finite_face_handles()) {
      if(fh2->is_outside())
        continue;
      const Point_2& p1 = fh2->vertex(0)->point();
      const Point_2& p2 = fh2->vertex(1)->point();
      const Point_2& p3 = fh2->vertex(2)->point();
      ::glVertex2d(p1.x(),p1.y());
      ::glVertex2d(p2.x(),p2.y());
      ::glVertex2d(p3.x(),p3.y());
    }
    ::glEnd();
    ::glDisable(GL_BLEND);
  }

  void gl_draw_dt2_edge() const
  {
    const Triangulation& tr = m_wrapper.triangulation();

    ::glColor3ub(127,127,127);
    ::glLineWidth(1.0f);
    ::glBegin(GL_LINES);
    for(const Face_handle& fh2 : tr.finite_face_handles()) {
      const Point_2& p1 = fh2->vertex(0)->point();
      const Point_2& p2 = fh2->vertex(1)->point();
      const Point_2& p3 = fh2->vertex(2)->point();
      ::glVertex2d(p1.x(),p1.y());
      ::glVertex2d(p2.x(),p2.y());
      ::glVertex2d(p2.x(),p2.y());
      ::glVertex2d(p3.x(),p3.y());
      ::glVertex2d(p3.x(),p3.y());
      ::glVertex2d(p1.x(),p1.y());
    }
    ::glEnd();
  }

  void gl_draw_voronoi() const
  {
    const Triangulation& tr = m_wrapper.triangulation();

    ::glColor3ub(133, 193, 233);
    ::glLineWidth(2.5f);

    ::glBegin(GL_LINES);
    for(const Edge& eh2 : tr.finite_edges()) {
      const Face_handle& f = eh2.first;
      const Face_handle& n = f->neighbor(eh2.second);
      Point_2 cc1 = f->circumcenter();
      Point_2 cc2 = n->circumcenter();
      if(tr.is_infinite(f)) {
        for(int i=0; i<3; ++i) {
          if(!tr.is_infinite(f->vertex(i)))
            continue;
          Edge boundary_edge(f, i);
          Segment_2 boundary_seg = AW2i::delaunay_edge_to_segment(boundary_edge);
          Line_2 boundary_edge_line(boundary_seg);
          Line_2 perpendicular = boundary_edge_line.perpendicular(cc1);

          // look for the opposite vertex to get the direction where the triangulation is finite
          Edge mirror_boundary_edge = tr.mirror_edge(boundary_edge);
          Vertex_handle oposite_vertex = n->vertex(mirror_boundary_edge.second);
          Point_2 oposite_vertex_projection = perpendicular.projection(tr.point(oposite_vertex));
          Vector_2 infinite_cc_unit_vec = (cc1 - oposite_vertex_projection);
          infinite_cc_unit_vec /= CGAL::sqrt(CGAL::to_double(infinite_cc_unit_vec.squared_length()));
          cc1 = (cc1 + (infinite_cc_unit_vec * 1e5));
        }
      }
      if(tr.is_infinite(n)) {
        for(int i=0; i<3; ++i) {
          if(!tr.is_infinite(n->vertex(i)))
            continue;
          Edge boundary_edge(n, i);
          Segment_2 boundary_seg = AW2i::delaunay_edge_to_segment(boundary_edge);
          Line_2 boundary_edge_line(boundary_seg);
          Line_2 perpendicular = boundary_edge_line.perpendicular(cc2);

          // look for the opposite vertex to get the direction where the triangulation is finite
          Edge mirror_boundary_edge = tr.mirror_edge(boundary_edge);
          Vertex_handle oposite_vertex = f->vertex(mirror_boundary_edge.second);
          Point_2 oposite_vertex_projection = perpendicular.projection(tr.point(oposite_vertex));
          Vector_2 infinite_cc_unit_vec = (cc2-oposite_vertex_projection);
          infinite_cc_unit_vec /= CGAL::sqrt(CGAL::to_double(infinite_cc_unit_vec.squared_length()));
          cc2 = (cc2 + (infinite_cc_unit_vec * 1e5));
        }
      }
      ::glVertex2d(cc1.x(),cc1.y());
      ::glVertex2d(cc2.x(),cc2.y());
    }
    ::glEnd();
  }

  void gl_draw_circle(const Point_2& center, const FT& sq_radius) const
  {
    double radius = std::sqrt(CGAL::to_double(sq_radius));
    ::glBegin(GL_LINE_LOOP);
    const double deg_to_rad = 3.14159/180;
    for(int i=0; i < 360; ++i) {
       double degree_radian = i*deg_to_rad;
       ::glVertex2d(center.x() + cos(degree_radian)*radius,
                    center.y() + sin(degree_radian)*radius);
    }
    ::glEnd();
  }

  void gl_draw_circle_cut_by_edge(const Point_2& center, FT sq_radius,
                                  const Edge& eh2) const
  {
    const Point_2& p1 = eh2.first->vertex((eh2.second+1)%3)->point();
    const Point_2& p2 = eh2.first->vertex((eh2.second+2)%3)->point();
    Point_2 p3;
    if(m_wrapper.triangulation().is_infinite(eh2.first->vertex(eh2.second))) {
      p3 = m_wrapper.triangulation().mirror_edge(eh2).first->vertex(
             m_wrapper.triangulation().mirror_edge(eh2).second)->point();
    } else {
      p3 = eh2.first->vertex(eh2.second)->point();
    }

    Line_2 eh2_line(p1,p2);
    bool is_negative_side_interior = false;
    if((eh2_line.oriented_side(p3) == CGAL::ON_NEGATIVE_SIDE) && !eh2.first->is_outside()) {
      is_negative_side_interior = true;
    } else if((eh2_line.oriented_side(p3) == CGAL::ON_POSITIVE_SIDE) && eh2.first->is_outside()) {
      is_negative_side_interior = true;
    }
    double radius = std::sqrt(CGAL::to_double(sq_radius));
    ::glBegin(GL_LINE_LOOP);
    const double deg_to_rad = 3.14159/180;
    for(int i=0; i<360; ++i)
    {
      double degree_radian = i*deg_to_rad;
      Point_2 circle_pt(center.x() + cos(degree_radian)*radius,
                        center.y() + sin(degree_radian)*radius);
      if((eh2_line.oriented_side(circle_pt) == CGAL::ON_NEGATIVE_SIDE) && is_negative_side_interior)
        continue;
      if((eh2_line.oriented_side(circle_pt) == CGAL::ON_POSITIVE_SIDE) && !is_negative_side_interior)
        continue;
      ::glVertex2d(circle_pt.x(),circle_pt.y());
    }
    ::glEnd();
  }

  void gl_draw_first_empty_circle_on_the_voronoi_edge(FT min_empty_circle_sq_radius,
                                                      const Edge& delaunay_edge,
                                                      const Triangulation& tr,
                                                      bool is_cut_by_delaunay_edge) const
  {
    const Point_2& p1 = delaunay_edge.first->vertex((delaunay_edge.second+1)%3)->point();
    const Point_2& p2 = delaunay_edge.first->vertex((delaunay_edge.second+2)%3)->point();
    const Point_2& m = CGAL::midpoint(p1,p2);
    FT edge_half_length_squared = CGAL::squared_distance(p1,p2) / 4.0;
    FT circle_center_to_edge_midpoint_sq_distance = min_empty_circle_sq_radius - edge_half_length_squared;

    const Face_handle& f = delaunay_edge.first;
    const Face_handle& n = f->neighbor(delaunay_edge.second);

    Point_2 cc1 = f->circumcenter();
    Point_2 cc2 = n->circumcenter();
    if(tr.is_infinite(f))
      cc1 = CGAL::midpoint(CGAL::midpoint(p1,p2),cc1);
    if(tr.is_infinite(n))
      cc2 = CGAL::midpoint(CGAL::midpoint(p1,p2),cc2);

    Vector_2 unit;
    if(CGAL::squared_distance(cc1,m) < circle_center_to_edge_midpoint_sq_distance) {
      unit = cc2 - cc1;
    } else {
      unit = cc1 - cc2;
    }
    FT cc1_cc2_dist = CGAL::sqrt(CGAL::to_double(unit.squared_length()));
    unit /= cc1_cc2_dist;

    Point_2 circle_center;
    if(circle_center_to_edge_midpoint_sq_distance > 0.) {
      circle_center = m + (unit * CGAL::sqrt(CGAL::to_double(circle_center_to_edge_midpoint_sq_distance)));
    } else {
      circle_center = m;
    }

    if(is_cut_by_delaunay_edge) {
      gl_draw_circle_cut_by_edge(circle_center, min_empty_circle_sq_radius, delaunay_edge);
    } else {
      gl_draw_circle(circle_center, min_empty_circle_sq_radius);
    }
  }

  void gl_draw_all_pencils()
  {
    const Triangulation& tr = m_wrapper.triangulation();
    Alpha_PQ& queue = m_wrapper.queue();
    FT sq_alpha = CGAL::square(m_wrapper  .alpha());

    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glColor4f(0.25,0.25,0.25,0.25);
    ::glLineWidth(0.5f);

    std::vector<Gate> gates;
    gates.reserve(queue.size());
    while(!queue.empty()) {
      gates.push_back(std::move(const_cast<Gate&>(queue.top())));
      const Gate& gate = gates.back();
      const Edge& eh2 = gate.edge();
      if(gate.has_steiner_point()) {
        gl_draw_first_empty_circle_on_the_voronoi_edge(
          CGAL::Alpha_wraps_2::internal::smallest_squared_radius_2(eh2, tr), eh2, tr, true);
      } else {
        gl_draw_first_empty_circle_on_the_voronoi_edge(
          CGAL::Alpha_wraps_2::internal::smallest_squared_radius_2(eh2, tr), eh2, tr, false);
      }
      queue.pop();
    }

    // restore queue
    for(auto it = gates.rbegin(); it != gates.rend(); ++it) {
      queue.push(std::move(*it));
    }

    for(const Edge& eh2 : tr.finite_edges()) {
      const Face_handle& f = eh2.first;
      const Face_handle& n = f->neighbor(eh2.second);
      if(f->is_outside() == n->is_outside())
        continue;
      if(sq_alpha < CGAL::Alpha_wraps_2::internal::smallest_squared_radius_2(eh2, tr))
        continue;

      const Point_2& cc1 = f->circumcenter();
      const Point_2& cc2 = n->circumcenter();
      const Point_2& p1 = eh2.first->vertex((eh2.second+1)%3)->point();
      FT empty_circle_cc1_sq_radius = CGAL::squared_distance(cc1,p1);
      FT empty_circle_cc2_sq_radius = CGAL::squared_distance(cc2,p1);
      FT max_empty_circle_sq_radius = CGAL::max(empty_circle_cc1_sq_radius,
                                                empty_circle_cc2_sq_radius);
      FT alpha_or_max_empty_circle_sq_radius = CGAL::min(max_empty_circle_sq_radius, sq_alpha);
      ::glColor3ub(255, 166, 216);
      ::glLineWidth(1.5f);
      gl_draw_first_empty_circle_on_the_voronoi_edge(alpha_or_max_empty_circle_sq_radius,
                                                     eh2, tr, true);
      ::glColor4f(0.25,0.25,0.25,0.25);
      ::glLineWidth(0.5f);
    }
  }

  void gl_draw_next_gate_pencil() const
  {
    const Triangulation& tr = m_wrapper.triangulation();
    const Alpha_PQ& queue = m_wrapper.queue();
    if(queue.empty())
      return;

    ::glLineWidth(1.9f);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glColor4f(0.25,0.25,0.25,0.8);
    const Gate& gate = queue.top();
    const Edge& eh2 = gate.edge();
    if(gate.has_steiner_point()) {
      gl_draw_first_empty_circle_on_the_voronoi_edge(
        CGAL::Alpha_wraps_2::internal::smallest_squared_radius_2(eh2, tr), eh2, tr, true);
    } else {
      gl_draw_first_empty_circle_on_the_voronoi_edge(
        CGAL::Alpha_wraps_2::internal::smallest_squared_radius_2(eh2, tr), eh2, tr, false);
    }
  }

  void gl_draw_steiner_point()
  {
    Alpha_PQ& queue = m_wrapper.queue();
    if(queue.empty())
      return;

    std::vector<Gate> gates;
    gates.reserve(queue.size());
    while(!queue.empty()) {
      gates.push_back(std::move(const_cast<Gate&>(queue.top())));
      queue.pop();
    }

    if(gates.empty())
      return;

    Gate top_gate = gates.front();

    // projection points
    ::glColor3ub(35, 155, 86 );
    ::glPointSize(10.f);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glEnable(GL_POINT_SMOOTH);
    ::glBegin(GL_POINTS);

    for(const Gate& gate : gates) {
      if(gate.has_steiner_from_projection()) {
        const Point_2& steiner_point = gate.steiner_point();
        ::glVertex2d(steiner_point.x(), steiner_point.y());
      }
    }
    ::glEnd();

    // intersection points
    ::glColor3ub(52, 152, 219);
    ::glPointSize(10.f);
    ::glBegin(GL_POINTS);
    for(const Gate& gate : gates) {
      if(gate.has_steiner_from_intersection()) {
        const Point_2& steiner_point = gate.steiner_point();
        ::glVertex2d(steiner_point.x(), steiner_point.y());
      }
    }
    ::glEnd();

    // next point
    const Gate& gate = top_gate;
    if(gate.has_steiner_point()) {
      ::glPointSize(15.f);
      ::glBegin(GL_POINTS);
      bool is_projection = gate.has_steiner_from_projection();
      is_projection ? ::glColor3ub(35, 155, 86 ) : ::glColor3ub(52, 152, 219);
      const Point_2& steiner_point = gate.steiner_point();
      ::glVertex2d(steiner_point.x(),steiner_point.y());
      ::glEnd();
    }

    // restore queue
    for(auto it = gates.rbegin(); it != gates.rend(); ++it) {
      queue.push(std::move(*it));
    }
  }

  void gl_draw_next_gate() const
  {
    const Alpha_PQ& queue = m_wrapper.queue();
    if(queue.empty())
      return;
    const Gate& gate = queue.top();
    const Edge& e = gate.edge();
    Segment_2 seg = AW2i::delaunay_edge_to_segment(e);
    const Point_2& p1 = seg.source();
    const Point_2& p2 = seg.target();
    ::glColor3ub(35, 155, 86 );
    ::glLineWidth(7.f);
    ::glBegin(GL_LINES);
    ::glVertex2d(p1.x(),p1.y());
    ::glVertex2d(p2.x(),p2.y());
    ::glEnd();
  }

  void add_vertex(const Point_2& p, bool new_cmp, bool is_closed)
  {
    if(m_input.empty() || new_cmp) {
      Pslg_component cmp;
      cmp.push_back(p);

      cmp.set_is_closed(is_closed);
      if(is_closed) {
        cmp.push_back(p);
      }
      m_input.push_back(cmp);
    } else {
      m_input[m_input.size()-1].set_is_closed(is_closed);
      if(is_closed) {
        auto pos_before_last = m_input[m_input.size()-1].end()-1;
        m_input[m_input.size()-1].insert(pos_before_last,p);
      } else {
        m_input[m_input.size()-1].push_back(p);
      }
    }
  }

  void close_input()
  {
    m_input[m_input.size()-1].push_back(m_input[m_input.size()-1][0]);
    m_input[m_input.size()-1].set_is_closed(true);
  }

  void open_input()
  {
    m_input[m_input.size()-1].pop_back();
    m_input[m_input.size()-1].set_is_closed(false);
  }

  CGAL::Bbox_2 get_bbox() const
  {
    CGAL::Bbox_2 bbox = m_input.bbox_2();
    if(!m_wrap.empty())
      bbox += m_wrap.bbox_2();
    return bbox;
  }

  void explode_input()
  {
    Pslg exploded_input;
    for(const Pslg_component& cmp : m_input) {
      for(size_t i = 1; i < cmp.size(); ++i) {
        Pslg_component exploded_cmp;
        exploded_cmp.push_back(random_point(cmp[i-1]));
        exploded_cmp.push_back(random_point(cmp[i]));
        exploded_input.push_back(exploded_cmp);
      }
    }
    m_input = exploded_input;
  }

  void smooth_input()
  {
    for(auto& cmp : m_input)
      cmp.resample(cmp.length()/100.);
    m_input.smooth(1);
  }

  double get_diagonal_bbox() const
  {
    CGAL::Bbox_2 bbox = m_input.bbox_2();
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()));
    return diag_length;
  }

  Point_2 random_point(const Point_2& in) const
  {
    double random_offset = get_diagonal_bbox() / 100;
    FT x = in.x() + ((static_cast<double>(rand()) / RAND_MAX) * random_offset);
    FT y = in.y() + ((static_cast<double>(rand()) / RAND_MAX) * random_offset);
    return Point_2(x,y);
  }

  void init_alpha_data_structure(const double alpha_value,
                                 const double offset_value)
  {
    Oracle oracle = m_wrapper.oracle();
    oracle.clear();

    std::vector<Segment_2> segments;
    for(const auto& cmp : m_input) {
      if(cmp.size() < 2)
        continue;
      for(size_t i=1; i<cmp.size(); ++i) {
        const Point_2& start = cmp[i-1];
        const Point_2& end = cmp[i];
        Segment_2 seg(start, end);
        if(seg.is_degenerate())
          continue; // ignore input segments that are problematic for aabb tree
        segments.push_back(seg);
      }
    }

    oracle.add_segments(segments);

    m_wrapper.initialize(alpha_value, offset_value, false /*refining*/);
  }

  void alpha_flood_fill(int max_iter = -1)
  {
    struct Demo_visitor
      : CGAL::Alpha_wraps_2::internal::Wrapping_default_visitor
    {
      Demo_visitor(Stats_data_structure& stats, int max_iter)
        : m_stats(stats), m_max_iter(max_iter), m_iter(0)
      { }

    public:
      // Whether the flood filling process should continue
      constexpr bool go_further(const Alpha_wrapper&) {
        return (m_iter < m_max_iter);
      }

      void before_edge_treatment(const Alpha_wrapper&, const Gate& gate)
      {
        ++m_stats.nb_gate_traversed;
        if(gate.has_steiner_point()) {
          ++m_iter;
          if(gate.has_steiner_from_projection())
            ++m_stats.nb_proj;
          ++m_stats.nb_steiner;
        }
      }

    private:
      Stats_data_structure& m_stats;
      std::size_t m_max_iter, m_iter;
    };

    Demo_visitor visitor(m_stats, max_iter);
    m_wrapper.alpha_flood_fill(visitor);
  }

  void extract_pslg_2_soup()
  {
    m_wrap = AW2i::extract_pslg_2_soup<EPICK>(m_wrapper.triangulation());
  }
};

#endif // CGAL_ALPHA_WRAP_2_DEMO_SCENE_H
