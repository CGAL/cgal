// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_CANVAS_EUCLIDEAN_METRIC_FIELD
#define CGAL_CANVAS_EUCLIDEAN_METRIC_FIELD

#include <CGAL/Canvas/Metric_field.h>

namespace CGAL {
namespace Canvas {
namespace Metric_fields {

template<typename K>
class Euclidean_metric_field
  : public Metric_field<K>
{
public:
  using Base = Metric_field<K>;
  using FT = typename Base::FT;
  using Metric = typename Base::Metric;
  using Point_3 = typename Base::Point_3;
  using Vector_3 = typename Base::Vector_3;

private:
  double m_a, m_b, m_c;

public:
  double& a() { return m_a; }
  double& b() { return m_b; }
  double& c() { return m_c; }

public:
  virtual Metric compute_metric(const Point_3&) const
  {
    return this->build_metric(Vector_3(1, 0, 0),
                              Vector_3(0, 1, 0),
                              Vector_3(0, 0, 1),
                              m_a, m_b, m_c);
  }

public:
  Euclidean_metric_field(const double a = 1.,
                         const double b = 1.,
                         const double c = 1.,
                         FT epsilon_ = 1e-6,
                         const double en_factor = 1.)
  : Metric_field<K>(epsilon_, en_factor), m_a(a), m_b(b), m_c(c)
  { }
};

} // namespace Metric_fields
} // namespace Canvas
} // namespace CGAL

#endif
