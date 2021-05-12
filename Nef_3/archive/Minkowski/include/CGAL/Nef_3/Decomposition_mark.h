// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_DECOMPOSITION_MARK_H
#define CGAL_DECOMPOSITION_MARK_H

namespace CGAL {

class Decomposition_mark {

  bool mark;

  friend std::ostream& operator<<(std::ostream& out, const Decomposition_mark& m);
  friend std::istream& operator>>(std::istream& in, Decomposition_mark& m);

 public:
  Decomposition_mark() {}
  Decomposition_mark(bool m) : mark(m) {}

  bool operator==(const Decomposition_mark& m) const {
    return mark == m.mark;
  }

  bool operator!=(const Decomposition_mark& m) const {
    return !*this==m;
  }

  Decomposition_mark operator!() const {
    return Decomposition_mark(!mark);
  }

  Decomposition_mark operator&&(const Decomposition_mark& m) const {
    return Decomposition_mark(mark && m.mark);
  }

  Decomposition_mark operator||(const Decomposition_mark& m) const {
    return Decomposition_mark(mark || m.mark);
  }
};

std::ostream& operator<<(std::ostream& out, const Decomposition_mark& m) {
  out << m.mark;
  return out;
}

std::istream& operator>>(std::istream& in, Decomposition_mark& m) {
  in >> m.mark;
  return in;
}

} //namespace CGAL
#endif // CGAL_DECOMPOSITION_MARK_H
