// Copyright (c) 1999,2007  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stefan Schirra, Michael Hemmer


#ifndef CGAL_NUMBER_TYPE_BASIC_H
#define CGAL_NUMBER_TYPE_BASIC_H

#include <CGAL/number_type_config.h>

#include <CGAL/basic.h>

// basic tools needed in several files
#include <boost/type_traits/is_same.hpp>
#include <functional>

#include <CGAL/Quotient_fwd.h>

#include <CGAL/Kernel/mpl.h>      // First_if_different
#include <CGAL/enum.h>            // CGAL::Sign etc.
#include <CGAL/tags.h>            // Tag_true / Tag_false

#include <CGAL/Coercion_traits.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Rational_traits.h>

#include <CGAL/Scalar_factor_traits.h>       // not part of CGAL 3.3
#include <CGAL/Algebraic_extension_traits.h> // not part of CGAL 3.3 

#include <CGAL/Needs_parens_as_product.h>

#include <CGAL/utils_classes.h>
#include <CGAL/utils.h>
#include <CGAL/FPU.h>

#include <CGAL/float.h>
#include <CGAL/double.h>
#include <CGAL/long_double.h>

#include <CGAL/Interval_nt.h> // needed by To_interval(long double), To_interval(long), To_interval(long long)

#include <CGAL/int.h>
#ifdef CGAL_USE_LONG_LONG
#include <CGAL/long_long.h>
#endif



#ifdef CGAL_USE_GMP
#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif // CGAL_USE_GMPXX
#endif // CGAL_USE_GMP

#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>

#endif // CGAL_NUMBER_TYPE_BASIC_H
