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

#ifndef __STLPORT_CSTD_string
# define __STLPORT_CSTD_string 1

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

# if defined ( __STL_REDEFINE_STD ) && defined (std) 
#    undef std
#    define __STL_RESUME_STD_FOR_string
#    define __STLPORT_NATIVE_PASS
# endif

#   include __STL_NATIVE_C_HEADER(string.h)

# if defined ( __STL_RESUME_STD_FOR_string )
#    undef __STL_RESUME_STD_FOR_string
#    define std __STLPORT_NAMESPACE
#    undef __STLPORT_NATIVE_PASS
# endif


# if defined (__BORLANDC__) && defined (__STL_IMPORT_VENDOR_CSTD)

using __STL_VENDOR_CSTD::size_t;
using __STL_VENDOR_CSTD::memcpy;
using __STL_VENDOR_CSTD::memchr;
using __STL_VENDOR_CSTD::memmove;
using __STL_VENDOR_CSTD::memcmp;
using __STL_VENDOR_CSTD::memset;

using __STL_VENDOR_CSTD::strcat;
using __STL_VENDOR_CSTD::strchr;
using __STL_VENDOR_CSTD::strcmp;
using __STL_VENDOR_CSTD::strcoll;
using __STL_VENDOR_CSTD::strcpy;

using __STL_VENDOR_CSTD::strcspn;
using __STL_VENDOR_CSTD::strerror;
using __STL_VENDOR_CSTD::strlen;
using __STL_VENDOR_CSTD::strncat;
using __STL_VENDOR_CSTD::strncmp;

using __STL_VENDOR_CSTD::strncpy;
using __STL_VENDOR_CSTD::strpbrk;
using __STL_VENDOR_CSTD::strrchr;
using __STL_VENDOR_CSTD::strspn;
using __STL_VENDOR_CSTD::strstr;

using __STL_VENDOR_CSTD::strtok;
using __STL_VENDOR_CSTD::strxfrm;

# endif /* BORLAND */

#endif /* __STLPORT_string */
// #endif /* NATIVE */
// Local Variables:
// mode:C++
// End:
