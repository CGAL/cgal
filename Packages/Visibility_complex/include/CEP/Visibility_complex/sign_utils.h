// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef SIGN_UTILS_H
#define SIGN_UTILS_H

#include <utility>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

Sign
operator *(Sign s, Sign t) 
{
    if ((s == POSITIVE && t == POSITIVE) || 
	(s == NEGATIVE && t == NEGATIVE)) return POSITIVE;
    else if ((s == POSITIVE && t == NEGATIVE) ||
	     (s == NEGATIVE && t == POSITIVE)) return NEGATIVE;
    return ZERO;
}

// is s and t have the same sign (or one or both are ZERO) the first element of
// the pair is true and the sign is the sign of the sum.
// Otherwise the first element of the pair if false and the second is the sign
// of s.
std::pair<bool,Sign>
sign_of_sum(Sign s, Sign t)
{
    typedef std::pair<bool,Sign> p;
    if (s == POSITIVE) {
	if (t == POSITIVE || t == ZERO) return p(true,POSITIVE);
	else return p(false,POSITIVE);
    }
    else if (s == NEGATIVE) {
	if (t == NEGATIVE || t == ZERO) return p(true,NEGATIVE);
	else return p(false,NEGATIVE);
    }
    return p(true,t);
}

CGAL_END_NAMESPACE

#endif // SIGN_UTILS_H
