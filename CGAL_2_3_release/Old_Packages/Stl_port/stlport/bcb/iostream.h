/*
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/*
 *
 * This wrapper is needed for Borland C++ Builder to get 
 * rid of ambiguity introduced by "using namespace std;" clause
 * in its native <iostream.h>
 *
 */

#ifndef __STLPORT_IOSTREAM_H
#  define  __STLPORT_IOSTREAM_H
#  include  <..\iostream.>

using __STL_VENDOR_STD::istream;
using __STL_VENDOR_STD::ostream;
using __STL_VENDOR_STD::cin;
using __STL_VENDOR_STD::cout;
using __STL_VENDOR_STD::cerr;
using __STL_VENDOR_STD::clog;
using __STL_VENDOR_STD::endl;
using __STL_VENDOR_STD::ends;

using __STL_VENDOR_STD::ios;
using __STL_VENDOR_STD::flush;

#endif

// Local Variables:
// mode:C++
// End:
