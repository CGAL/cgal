/*
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
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
#ifndef __STL_STRING_FWD_C
#define __STL_STRING_FWD_C

__STL_BEGIN_NAMESPACE

template <class _CharT, class _Traits, class _Alloc>
const char* 
__get_c_string(const basic_string<_CharT,_Traits,_Alloc>& __str) { 
  return __str.c_str(); 
}

__STL_END_NAMESPACE

#endif /*  __STL_STRING_FWD_C */

// Local Variables:
// mode:C++
// End:
