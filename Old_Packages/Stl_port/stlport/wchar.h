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

#ifndef __STLPORT_CSTD_wchar
# define __STLPORT_CSTD_wchar

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

# if ! defined (__STL_WINCE)

# if defined ( __STL_REDEFINE_STD ) && defined (std) 
#    undef std
#    define __STL_RESUME_STD_FOR_wchar
#    define __STLPORT_NATIVE_PASS
# endif

#if defined (__SUNPRO_CC) && (__SUNPRO_CC >= 0x500) && !defined (__STLPORT_NEW_IOSTREAMS)
#if !defined(__SunOS_5_5_1) && !defined(__SunOS_5_6)
#define _MBSTATE_T
#define _STD_MBSTATE_T
#define _MBSTATET_H
#endif
namespace std { class mbstate_t; }
#endif

# include __STL_NATIVE_C_HEADER(wchar.h)

# if defined ( __STL_RESUME_STD_FOR_wchar )
#    undef __STL_RESUME_STD_FOR_wchar
#    define std __STLPORT_NAMESPACE
#    undef __STLPORT_NATIVE_PASS
# endif

# endif /* WINCE */

#endif /* __STLPORT_wchar */

// Local Variables:
// mode:C++
// End:
