// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_MORE_FUNCTIONS_ON_SIGNS_H
#define CGAL_MORE_FUNCTIONS_ON_SIGNS_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/enum.h>

namespace CGAL {
inline
Sign
sign_of_first_root(const Sign &sb, const Sign &sc)
{
  // this method computes the sign of the second root of the quadratic
  // equation: a x^2 - b x + c == 0. It is assumed that the sign of a
  // is +1.

  if ( sc == NEGATIVE )
    return NEGATIVE;

  if ( sb == POSITIVE )  return sc;
  if ( sb == NEGATIVE )  return NEGATIVE;
  return opposite(sc);
}


inline
Sign
sign_of_second_root(const Sign &sb, const Sign &sc)
{
  // this method computes the sign of the first root of the quadratic
  // equation: a x^2 - b x + c == 0. It is assumed that the sign of a
  // is +1.

  if ( sc == NEGATIVE )
    return POSITIVE;

  if ( sb == POSITIVE )  return POSITIVE;
  if ( sb == NEGATIVE )  return opposite(sc);
  return sc;
}

} //namespace CGAL

#endif // CGAL_MORE_FUNCTIONS_ON_SIGNS_H
