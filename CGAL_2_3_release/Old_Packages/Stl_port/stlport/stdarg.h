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

#ifndef __STLPORT_CSTD_stdarg
# define __STLPORT_CSTD_stdarg

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif


// GCC has some problems if this is being included 
// before <stddef.h>
# ifdef __GNUC__
#  include <stddef.h>
# endif

# if defined ( __STL_REDEFINE_STD ) && defined (std) 
#    undef std
#    define __STL_RESUME_STD_FOR_stdarg
#    define __STLPORT_NATIVE_PASS
# endif

# include __STL_NATIVE_C_HEADER(stdarg.h)

# if defined ( __STL_RESUME_STD_FOR_stdarg )
#    undef __STL_RESUME_STD_FOR_stdarg
#    define std __STLPORT_NAMESPACE
#    undef __STLPORT_NATIVE_PASS
# endif

#endif /* __STLPORT_stdarg */

// Local Variables:
// mode:C++
// End:
