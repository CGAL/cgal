// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Verification.h"

#include <iostream>

#include "Kml_reader.h"

void Verification::find_minimum_projected_error_on_sphere(float we,
                                                          Camera& cam,
                                                          int vp_width,
                                                          int vp_height) {
  QRect vp(0, 0, vp_width, vp_height);
  auto proj = cam.get_projection_matrix();
  auto view = cam.get_view_matrix();
  QMatrix4x4 model;
  auto model_view = view * model;

  float max_err = 0;
  float max_theta = -1;
  float max_phi = -1;

  int num_divs = 200;
  const float dtheta = M_PI_2 / num_divs;
  const float dphi = M_PI_2 / num_divs;

  const float r1 = 1.f;
  const float r2 = r1 - we;
  for (int i = 0; i <= num_divs; ++i) {
    const float theta = dtheta * i;
    const float cos_theta = std::cos(theta);
    const float sin_theta = std::sin(theta);

    for (int j = 0; j <= num_divs; ++j) {
      QVector3D p1, p2;
      const float phi = dphi * j;
      const float cos_phi = std::cos(phi);
      const float sin_phi = std::sin(phi);

      // p1
      const float r1xz = r1 * sin_phi;
      p1.setY(r1 * cos_phi);
      p1.setX(r1xz * cos_theta);
      p1.setZ(r1xz * sin_theta);

      // p2
      const float r2xz = r2 * sin_phi;
      p2.setY(r2 * cos_phi);
      p2.setX(r2xz * cos_theta);
      p2.setZ(r2xz * sin_theta);

      auto wp1 = p1.project(model_view, proj, vp);
      auto wp2 = p2.project(model_view, proj, vp);

      const auto pe = wp1.distanceToPoint(wp2);
      if (max_err < pe) {
        max_err = pe;
        max_theta = theta;
        max_phi = phi;
      }
    }
  }

  std::cout << "max err = " << max_err << std::endl;
  std::cout << "max phi = " << max_phi * 180 / M_PI << std::endl;
  std::cout << "max theta = " << max_theta * 180 / M_PI << std::endl;

  auto wp1 = QVector3D(0, r1, 0).project(model_view, proj, vp);
  auto wp2 = QVector3D(0, r2, 0).project(model_view, proj, vp);
  auto pe = wp1.distanceToPoint(wp2);
  std::cout << "polar err = " << pe << std::endl;

  wp1 = QVector3D(r1, 0, 0).project(model_view, proj, vp);
  wp2 = QVector3D(r2, 0, 0).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "x-axis err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, 1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, 1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "nearest proj err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, -1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, -1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "farthest proj err = " << pe << std::endl;

  // project the origin on the screen (to check if it projects to the mid-vp)
  //std::cout << QVector3D(0, 0, 0).project(model_view, proj, vp) << std::endl;
}

//! \brief
void Verification::verify_antarctica_node_is_redundant() {
  Kml::Node n1(178.277211542064, -84.4725179992025),
            n2(180.0, -84.71338),
            n3(-179.942499356179, -84.7214433735525);

  // 1) check if it is collinear with its neighboring nodes:
  // all of the vectors in 3D must lie in the same plane
  auto v1 = n1.get_coords_3f();
  auto v2 = n2.get_coords_3f();
  auto v3 = n3.get_coords_3f();
  auto n = QVector3D::crossProduct(v1, v3);
  n.normalize();
  std::cout << "** DOT PRODUCT = " << QVector3D::dotProduct(n, v2) << std::endl;

  // 2) check if it is between its neighbors (check if r,s > 0)
  auto det = [](float ax, float ay, float bx, float by)
             { return ax * by - ay * bx; };
  auto D = det(v1.x(), v1.y(), v3.x(), v3.y());
  auto Dr = det(v2.x(), v2.y(), v3.x(), v3.y());
  auto Ds = det(v1.x(), v1.y(), v2.x(), v2.y());
  auto r = Dr / D;
  auto s = Ds / D;
  std::cout << "r = " << r << "\ns=" << s << std::endl;
}
