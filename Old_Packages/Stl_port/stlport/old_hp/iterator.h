/*
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

#ifndef __SGI_STL_ITERATOR_H
#define __SGI_STL_ITERATOR_H

#ifndef __STL_CONFIG_H
#include <stl_config.h>
#endif

#if defined  (__STL_DEBUG) || defined (__STL_ASSERTIONS) && !defined (__STLPORT_DEBUG_H)
# include <stldebug.h>
#endif

#if defined (__STL_USE_NEW_STYLE_HEADERS)
# include <cstddef>
#else
# include <stddef.h>
#endif

# ifndef __STLPORT_NEW
#  include <new>
# endif

# ifndef __TYPE_TRAITS_H
#  include <type_traits.h>
# endif

#ifndef __SGI_STL_FUNCTION_H
#include <function.h>
#endif

# ifndef __STLPORT_IOSFWD
#  include <iosfwd>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_BASE_H
#  include <stl_iterator_base.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

#ifndef __SGI_STL_INTERNAL_CONSTRUCT_H
#include <stl_construct.h>
#endif

#ifndef __SGI_STL_INTERNAL_RAW_STORAGE_ITERATOR_H
#include <stl_raw_storage_iter.h>
#endif

# ifndef __SGI_STL_INTERNAL_STREAM_ITERATOR_H
#  include <stl_stream_iterator.h>
# endif

#ifdef __STL_USE_NAMESPACES

// Names from stl_iterator.h

# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else

using __STLPORT_STD::input_iterator_tag;
using __STLPORT_STD::output_iterator_tag;
using __STLPORT_STD::forward_iterator_tag;
using __STLPORT_STD::bidirectional_iterator_tag;
using __STLPORT_STD::random_access_iterator_tag;

using __STLPORT_STD::input_iterator;
using __STLPORT_STD::output_iterator;
using __STLPORT_STD::forward_iterator;
using __STLPORT_STD::bidirectional_iterator;
using __STLPORT_STD::random_access_iterator;

#ifdef __STL_CLASS_PARTIAL_SPECIALIZATION
using __STLPORT_STD::iterator_traits;
#endif

using __STLPORT_STD::iterator_category;
using __STLPORT_STD::distance_type;
using __STLPORT_STD::value_type;

using __STLPORT_STD::distance; 
using __STLPORT_STD::advance; 

using __STLPORT_STD::insert_iterator;
using __STLPORT_STD::front_insert_iterator;
using __STLPORT_STD::back_insert_iterator;
using __STLPORT_STD::inserter;
using __STLPORT_STD::front_inserter;
using __STLPORT_STD::back_inserter;

using __STLPORT_STD::reverse_iterator;
using __STLPORT_STD::reverse_bidirectional_iterator;

using __STLPORT_STD::istream_iterator;
using __STLPORT_STD::ostream_iterator;

// Names from stl_construct.h
using __STLPORT_STD::construct;
using __STLPORT_STD::destroy;

// Names from stl_raw_storage_iter.h
using __STLPORT_STD::raw_storage_iterator;
# endif

#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_ITERATOR_H */

// Local Variables:
// mode:C++
// End:
