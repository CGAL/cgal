/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-only

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/vec.h>

// Most of the methods are declared inline in vec.h

namespace CGAL{
namespace qglviewer{


/*! Projects the Vec on the axis of direction \p direction that passes through
the origin.

\p direction does not need to be normalized (but must be non null). */
CGAL_INLINE_FUNCTION
void Vec::projectOnAxis(const Vec &direction) {
#ifndef QT_NO_DEBUG
  if (direction.squaredNorm() < 1.0E-10)
    qWarning("Vec::projectOnAxis: axis direction is not normalized (norm=%f).",
             direction.norm());
#endif

  *this = (((*this) * direction) / direction.squaredNorm()) * direction;
}

/*! Projects the Vec on the plane whose normal is \p normal that passes through
the origin.

\p normal does not need to be normalized (but must be non null). */
CGAL_INLINE_FUNCTION
void Vec::projectOnPlane(const Vec &normal) {
#ifndef QT_NO_DEBUG
  if (normal.squaredNorm() < 1.0E-10)
    qWarning("Vec::projectOnPlane: plane normal is not normalized (norm=%f).",
             normal.norm());
#endif

  *this -= (((*this) * normal) / normal.squaredNorm()) * normal;
}

/*! Returns a Vec orthogonal to the Vec. Its norm() depends on the Vec, but is
 zero only for a null Vec. Note that the function that associates an
 orthogonalVec() to a Vec is not continous. */
CGAL_INLINE_FUNCTION
Vec Vec::orthogonalVec() const {
  // Find smallest component. Keep equal case for null values.
  if ((fabs(y) >= 0.9 * fabs(x)) && (fabs(z) >= 0.9 * fabs(x)))
    return Vec(0.0, -z, y);
  else if ((fabs(x) >= 0.9 * fabs(y)) && (fabs(z) >= 0.9 * fabs(y)))
    return Vec(-z, 0.0, x);
  else
    return Vec(-y, x, 0.0);
}

CGAL_INLINE_FUNCTION
std::ostream &operator<<(std::ostream &o, const Vec &v) {
  return o << v.x << '\t' << v.y << '\t' << v.z;
}

}} // namespace CGAL::qglviewer
