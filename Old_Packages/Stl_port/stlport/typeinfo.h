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

#ifndef __STLPORT_OLDSTD_typeinfo
# define __STLPORT_OLDSTD_typeinfo

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

# if defined ( __STL_REDEFINE_STD ) && defined (std) 
#    undef std
#    define __STL_RESUME_STD_FOR_typeinfo_H
#    define __STLPORT_NATIVE_PASS
# endif

# if defined (__GNUC__) && (__GNUC_MINOR__ >= 8 )
#  include <../include/typeinfo.h>
# else
#  include __STL_NATIVE_HEADER(typeinfo.h)
# endif

# if defined  (__STL_USE_OWN_NAMESPACE)

__STL_BEGIN_NAMESPACE

using __STL_VENDOR_EXCEPT_STD::type_info;
using __STL_VENDOR_EXCEPT_STD::bad_typeid;
using __STL_VENDOR_EXCEPT_STD::bad_cast;

__STL_END_NAMESPACE

#endif /* __STL_OWN_NAMESPACE */

# if defined ( __STL_RESUME_STD_FOR_typeinfo_H )
#    undef __STL_RESUME_STD_FOR_typeinfo_H
#    define std __STLPORT_NAMESPACE
#    undef __STLPORT_NATIVE_PASS
# endif

#endif /* __STLPORT_OLDSTD_typeinfo */

// Local Variables:
// mode:C++
// End:
