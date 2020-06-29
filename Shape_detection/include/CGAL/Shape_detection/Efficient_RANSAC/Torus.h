// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_TORUS_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_TORUS_H

#include <CGAL/license/Shape_detection.h>

#include <cmath>
#include <CGAL/Circle_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>

namespace CGAL {
  namespace Shape_detection {

    /*!
      \ingroup PkgShapeDetectionRANSACShapes
      \brief Torus implements Shape_base. The torus is represented by the
      symmetry axis, its center on the axis, and the major and minor radii.
     \tparam Traits must be a model of `EfficientRANSACTraits` with the additional
             requirement for tori (see `EfficientRANSACTraits` documentation).
     */
  template <class Traits>
  class Torus : public Shape_base<Traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Point_map Point_map;
     ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
     ///< property map to access the unoriented normal of an input point.
    typedef typename Traits::FT FT;
    ///< number type.
    typedef typename Traits::Point_3 Point_3;
    ///< point type.
    typedef typename Traits::Vector_3 Vector_3;
    ///< vector type.
    typedef typename Traits::Vector_2 Vector_2;
    ///< 2D vector type.
    typedef typename Traits::Point_2 Point_2;
    ///< 2D point type used during construction.
    typedef typename Traits::Circle_2 Circle_2;
     ///< circle type used during construction.
    /// \endcond

    Torus() : Shape_base<Traits>() {}

    /*!
      Direction of symmetry axis.
     */
    Vector_3 axis() const {
      return m_axis;
    }

    /*!
      Center point on symmetry axis.
     */
    Point_3 center() const {
      return m_center;
    }

    /*!
      Major radius of the torus.
     */
    FT major_radius() const {
      return m_majorRad;
    }

    /*!
      Minor radius of the torus.
      */
    FT minor_radius() const {
      return m_minorRad;
    }

    /// \cond SKIP_IN_MANUAL

    /*!
      Helper function to write center point, symmetry axis
      and the two radii into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      sstr << "Type: torus center(" << this->get_x(m_center) << ", " << this->get_y(m_center);
      sstr << ", " << this->get_z(m_center) << ") axis(" << this->get_x(m_axis) << ", ";
      sstr << this->get_y(m_axis) << ", " << this->get_z(m_axis) << ") major radius = ";
      sstr << m_majorRad << " minor radius = " << m_minorRad << " #Pts: ";
      sstr << this->m_indices.size();

      return sstr.str();
    }

    /*!
      Computes squared Euclidean distance from query point to the shape.
      */
    FT squared_distance(const Point_3 &p) const {
      const Vector_3 d = this->constr_vec(m_center, p);

      // height over symmetry plane
      const FT height = this->scalar_pdct(d, m_axis);

      // distance from axis in plane
      const FT l = CGAL::sqrt(this->scalar_pdct(d, d) - height * height);

      // inPlane distance from circle
      const FT l2 = m_majorRad - l;

      // distance from torus
      const FT squared_dist = CGAL::sqrt(height * height + l2 * l2) - m_minorRad;

      return squared_dist * squared_dist;
    }
    /// \endcond

  protected:
    /// \cond SKIP_IN_MANUAL

    // ------------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------------
    FT get_x_2(const Vector_2& v) const { return this->m_traits.compute_x_2_object()(v); }
    FT get_y_2(const Vector_2& v) const { return this->m_traits.compute_y_2_object()(v); }
    FT get_x_2(const Point_2& p) const { return this->m_traits.compute_x_2_object()(p); }
    FT get_y_2(const Point_2& p) const { return this->m_traits.compute_y_2_object()(p); }

    Circle_2 constr_circle(const Point_2& a, const Point_2& b,const Point_2& c) const
    { return this->m_traits.construct_circle_2_object()(a, b, c); }
    Point_2 circle_center(const Circle_2& s) const
    { return this->m_traits.construct_center_2_object()(s); }
    FT sqradius(const Circle_2& s) const
    { return this->m_traits.compute_squared_radius_2_object()(s); }
    Point_2 constr_point_2(FT a, FT b) const
    { return this->m_traits.construct_point_2_object()(a, b); }
    Vector_2 constr_vec_2(const Point_2 &a, const Point_2 &b) const
    { return this->m_traits.construct_vector_2_object()(a, b); }
    FT sqlen_2(const Vector_2& v) const
    { return this->m_traits.compute_squared_length_2_object()(v); }
    bool collinear_2(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return this->m_traits.collinear_2_object()(p, q, r); }

    void create_shape(const std::vector<std::size_t> &indices) {
      std::vector<Point_3> p;
      std::vector<Vector_3> n;

      p.resize(indices.size());
      n.resize(indices.size());

      for (std::size_t i = 0;i<indices.size();i++) {
        p[i] = this->point(indices[i]);
        n[i] = this->normal(indices[i]);
      }

      // Implemented method from 'Geometric least-squares fitting of spheres, cylinders, cones and tori' by G. Lukacs,A.D. Marshall, R. R. Martin
      const FT a01 = this->scalar_pdct(this->cross_pdct(n[0], n[1]), n[2]);
      const FT b01 = this->scalar_pdct(this->cross_pdct(n[0], n[1]), n[3]);
      const FT a0 = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[1], p[2]), n[0]), n[2]);
      const FT b0 = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[1], p[3]), n[0]), n[3]);
      const FT a1 = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[2], p[0]), n[1]), n[2]);
      const FT b1 = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[3], p[0]), n[1]), n[3]);
      const FT a = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[2], p[0]), this->constr_vec(p[0], p[1])), n[2]);
      const FT b = this->scalar_pdct(this->cross_pdct(this->constr_vec(p[3], p[0]), this->constr_vec(p[0], p[1])), n[3]);

      FT div = (b1 * a01 - b01 * a1);
      if (div == 0)
        return;

      div = (FT)1.0 / div;
      const FT r = ((a01 * b + b1 * a0 - b0 * a1 - b01 * a)) * div * (FT)0.5;
      const FT q = (b * a0 - b0 * a) * div;

      FT root = r * r - q;
      if (r * r - q < 0)
        root = 0;

      const FT y1 = -r - CGAL::sqrt(root);
      const FT y2 = -r + CGAL::sqrt(root);
      FT x1 = (a01 * y1 + a0);
      FT x2 =  (a01 * y2 + a0);

      if (x1 == 0 || x2 == 0)
        return;

      x1 = -(a1 * y1 + a) / x1;
      x2 = -(a1 * y2 + a) / x2;

      // 1. center + axis
      FT majorRad1 = FLT_MAX, minorRad1 = FLT_MAX, dist1 = FLT_MAX;
      Point_3 c1 = this->constr_pt();
      Vector_3 axis1 = this->constr_vec();
      if (is_finite(x1) && is_finite(y1)) {
        c1 = this->transl(p[0], this->scale(n[0], x1));
        axis1 = this->constr_vec(this->transl(p[1], this->scale(n[1], y1)), c1);

        FT l = this->sqlen(axis1);
        if (l > (FT)0.00001 && l == l) {
          axis1 = this->scale(axis1, FT(1.0) / CGAL::sqrt(l));
          dist1 = getCircle(c1, axis1, p, majorRad1, minorRad1);
        }
      }

      // 2. center + axis
      FT majorRad2 = 0, minorRad2 = 0, dist2 = FLT_MAX;
      Point_3 c2 = this->constr_pt();
      Vector_3 axis2 = this->constr_vec();
      if (is_finite(x2) && is_finite(y2)) {
        c2 = this->transl(p[0], this->scale(n[0], x2));
        axis2 = this->constr_vec(this->transl(p[1], this->scale(n[1], y2)), c2);

        FT l = this->sqlen(axis2);
        if (l > (FT)0.00001 && l == l) {
          axis2 = this->scale(axis2, FT(1.0) / CGAL::sqrt(l));
          dist2 = getCircle(c2, axis2, p, majorRad2, minorRad2);
        }
      }

      if (dist1 < dist2) {
        m_center = c1;
        m_axis = axis1;
        m_majorRad = majorRad1;
        m_minorRad = CGAL::sqrt(minorRad1);
      }
      else {
        m_center = c2;
        m_axis = axis2;
        m_majorRad = majorRad2;
        m_minorRad = CGAL::sqrt(minorRad2);
      }

      // Drop if shape is probably sphere
      if (m_majorRad < this->m_epsilon) {
        this->m_is_valid = false;

        return;
      }

      //validate points and normals
      for (std::size_t i = 0;i<indices.size();i++) {
        // check distance
        if (squared_distance(p[i]) > this->m_epsilon) {
          this->m_is_valid = false;
          return;
        }

        // check normal deviation
        Vector_3 d = this->constr_vec(m_center, p[i]);

        Vector_3 in_plane = this->cross_pdct(m_axis,
          this->cross_pdct(m_axis, d));
        if (this->scalar_pdct(in_plane, d) < 0)
          in_plane = this->scale(in_plane, FT(-1.0));

        FT length = CGAL::sqrt(this->sqlen(in_plane));
        if (length == 0)
          return;

        in_plane = this->scale(in_plane, FT(1.0) / length);

        d = this->constr_vec((this->transl(m_center, this->scale(in_plane, m_majorRad))), p[i]);

        length = CGAL::sqrt(this->sqlen(d));
        if (length == 0)
          return;

        d = this->scale(d, FT(1.0) / length);
        if (CGAL::abs(this->scalar_pdct(d, n[i])) < this->m_normal_threshold) {
          this->m_is_valid = false;
          return;
        }
      }

      this->m_is_valid = true;
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Point_3 po = this->point(indices[i]);
        Vector_3 d = this->constr_vec(m_center, po);
        // height over symmetry plane
        const FT p = this->scalar_pdct(d, m_axis);
        // distance from axis in plane
        FT l = CGAL::sqrt(this->scalar_pdct(d, d) - p * p);

        // inPlane distance from circle
        const FT l2 = m_majorRad - l;

        // distance from torus
        l = CGAL::sqrt(p * p + l2 * l2) - m_minorRad;
        dists[i] = l * l;
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices,
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 d = this->constr_vec(m_center, this->point(indices[i]));

        Vector_3 in_plane = this->cross_pdct(m_axis,
                                              this->cross_pdct(m_axis, d));
        if (this->scalar_pdct(in_plane, d) < 0)
          in_plane = this->scale(in_plane, FT(-1.0));

        FT length = (FT) CGAL::sqrt(this->sqlen(in_plane));

        // If length is 0 the point is on the axis, maybe in the apex. We
        // accept any normal for that position.
        if (length == 0) {
          angles[i] = (FT)1.0;
          continue;
        }

        in_plane = this->scale(in_plane,FT(1.0) / CGAL::sqrt(this->sqlen(in_plane)));

        d = this->constr_vec((this->transl(m_center, this->scale(in_plane, m_majorRad))), this->point(indices[i]));
        d = this->scale(d, FT(1.0) / CGAL::sqrt(this->sqlen(d)));
        angles[i] = CGAL::abs(this->scalar_pdct(d, this->normal(indices[i])));
      }
    }

    FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 d = this->constr_vec(m_center, p);

      Vector_3 in_plane = this->cross_pdct(m_axis,
                                           this->cross_pdct(m_axis, d));
      if (this->scalar_pdct(in_plane, d) < 0)
        in_plane = -in_plane;

      float length = CGAL::sqrt(this->sqlen(in_plane));

      // If length is 0 the point is on the axis, maybe in the apex. We
      // accept any normal for that position.
      if (length == 0) {
        return (FT)1.0;
      }

      in_plane = in_plane / CGAL::sqrt(this->sqlen(in_plane));

      d = p - (m_center + this->scalar_pdct(in_plane, m_majorRad));
      d = d / CGAL::sqrt(this->sqlen(d));

      return CGAL::abs(d * n);
    }

    virtual std::size_t minimum_sample_size() const {
        return 4;
    }

    virtual bool supports_connected_component() const {
      return false;
    }

  private:
    FT getCircle(Point_3 &center, const Vector_3 &axis, std::vector<Point_3> p, FT &majorRad, FT &minorRad) const {
      // create spin image
      std::vector<Point_2> pts;
      pts.resize(p.size());
      for (unsigned int i = 0;i<p.size();i++) {
        Vector_3 d = this->constr_vec(center, p[i]);
        FT e = this->scalar_pdct(d, axis);
        FT f = this->scalar_pdct(d, d) - e * e;
        if (f <= 0)
          pts[i] = this->constr_point_2(e, (FT) 0);
        else
          pts[i] = this->constr_point_2(e, CGAL::sqrt(this->scalar_pdct(d, d) - e * e));
      }

      if (this->collinear_2(pts[0], pts[1], pts[2])) {
        return (std::numeric_limits<FT>::max)();
      }

      Circle_2 c = this->constr_circle(pts[0], pts[1], pts[2]);
      minorRad = this->sqradius(c);
      majorRad = this->get_y_2(this->circle_center(c));
      center = this->transl(center, this->scale(axis, this->get_x_2(this->circle_center(c))));

      return CGAL::abs(
        this->sqlen_2(this->constr_vec_2(this->circle_center(c), pts[3])) - this->sqradius(c));
    }

    Point_3 m_center;
    Vector_3 m_axis;
    FT m_majorRad;
    FT m_minorRad;
    /// \endcond
  };
}
}
#endif
