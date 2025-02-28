// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "GUI_country_pick_handler.h"

//#include <qvector3d.h>

//! \brief
GUI_country_pick_handler::GUI_country_pick_handler(Main_widget& main_widget) :
  m_main_widget(main_widget),
  m_camera(main_widget.get_camera())
{}

//! \brief
void GUI_country_pick_handler::mouse_press_event(QMouseEvent* e) {
  // handle country selection
  if (e->button() == Qt::RightButton) {
    auto p = e->pos();
    QVector3D sp0(p.x(), m_vp_height - p.y(), 0);
    QVector3D sp1(p.x(), m_vp_height - p.y(), 1);

    auto proj = m_camera.get_projection_matrix();
    auto view = m_camera.get_view_matrix();
    auto model_view = view * m_main_widget.get_model_matrix();
    QRect viewport(0, 0, m_vp_width, m_vp_height);
    auto wp0 = sp0.unproject(model_view, proj, viewport);
    auto wp1 = sp1.unproject(model_view, proj, viewport);

    // ASSERTION!!!
    m_main_widget.set_mouse_pos(wp0);

    // define a ray from the camera pos to the world-point
    //auto o = m_camera.get_pos();
    //auto u = wp - o;
    auto o = wp0;
    auto u = wp1 - wp0;

    // solve the quadratic equation to check for intersection of ray with sphere
    auto a = QVector3D::dotProduct(u, u);
    auto b = 2 * QVector3D::dotProduct(u, o);
    auto c = QVector3D::dotProduct(o, o) - 1;
    auto d = b * b - 4 * a * c;

    float ti = -1;
    if (abs(d) < std::numeric_limits<float>::epsilon()) {
      // single intersection
      ti = -b / (2 * a);
    }
    else {
      if (d < 0) {
        // no intersection
        return;
      }
      else {
        // two intersections
        auto sd = sqrt(d);
        auto t1 = (-b - sd) / (2 * a);
        auto t2 = (-b + sd) / (2 * a);
        if (t1 > 0 && t2 > 0) ti = (std::min)(t1, t2);
        else if (t1 > 0) ti = t1;
        else ti = t2;
      }
    }

    //m_mouse_pos = o + ti * u;
    auto pos = o + ti * u;
    m_main_widget.set_mouse_pos(pos);
    static std::string prev_picked_country;
    auto& arrh = m_main_widget.get_arr_handle();
    auto picked_country = Aos::locate_country(arrh, pos);

    m_main_widget.hightlight_country(picked_country);
    //  if (!prev_picked_country.empty())
  //  {
  //    // dim the previous country color
  //    auto& prev_country = m_country_triangles[prev_picked_country];
  //    auto color = prev_country->get_color();
  //    color *= s_dimming_factor;
  //    color.setW(1);
  //    prev_country->set_color(color);
  //  }

  //  if (!picked_country.empty())
  //  {
  //    // highlight the current country color
  //    auto& curr_country = m_country_triangles[picked_country];
  //    auto color = curr_country->get_color();
  //    color /= s_dimming_factor;
  //    color.setW(1);
  //    curr_country->set_color(color);
  //    qDebug() << "SELECTED COUNTRY: " << picked_country;
  //  }

  //  prev_picked_country = picked_country;
  }
}

//! \brief
void GUI_country_pick_handler::resize(int w, int h) {
  m_vp_width = w;
  m_vp_height = h;
}
