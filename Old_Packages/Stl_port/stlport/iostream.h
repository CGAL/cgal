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

#ifndef __STLPORT_OLDSTD_iostream
# define __STLPORT_OLDSTD_iostream

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

# if ! defined (__STL_NO_IOSTREAMS)

# if defined ( __STL_REDEFINE_STD ) && defined (std) 
#    undef std
#    define __STL_RESUME_STD_FOR_iostream_H
#    define __STLPORT_NATIVE_PASS
# endif

# include __STL_NATIVE_HEADER(iostream.h)

# if defined ( __STL_RESUME_STD_FOR_iostream_H )
#    undef __STL_RESUME_STD_FOR_iostream_H
#    define std __STLPORT_NAMESPACE
#    undef __STLPORT_NATIVE_PASS
# endif


# if defined (__STL_USE_NAMESPACES) && !defined (__STL_BROKEN_USING_DIRECTIVE)

__STL_BEGIN_NAMESPACE

using ::istream;
using ::ostream;
using ::cin;
using ::cout;
using ::cerr;
using ::clog;
using ::endl;
using ::ends;

using ::ios;
using ::flush;

// using ::ws;

__STL_END_NAMESPACE

# endif /* __STL_USE_OWN_NAMESPACE */

# endif /* WINCE */

#endif /* __STLPORT_OLDSTD_iostream */

// Local Variables:
// mode:C++
// End:
