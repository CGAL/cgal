//  (C) Copyright David Abrahams 2000. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//
//  The author gratefully acknowleges the support of Dragon Systems, Inc., in
//  producing this work.

//  Revision History:
//  04 Mar 01  Some fixes so it will compile with Intel C++ (Dave Abrahams)

#ifndef CONFIG_DWA052200_H_
# define CONFIG_DWA052200_H_

# include <boost/config.hpp>

# ifdef BOOST_NO_OPERATORS_IN_NAMESPACE
   // A gcc bug forces some symbols into the global namespace
#  define BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE
#  define BOOST_PYTHON_END_CONVERSION_NAMESPACE
#  define BOOST_PYTHON_CONVERSION
#  define BOOST_PYTHON_IMPORT_CONVERSION(x) using ::x
# else
#  define BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE namespace boost { namespace python {
#  define BOOST_PYTHON_END_CONVERSION_NAMESPACE }} // namespace boost::python
#  define BOOST_PYTHON_CONVERSION boost::python
#  define BOOST_PYTHON_IMPORT_CONVERSION(x) void never_defined() // so we can follow the macro with a ';'
# endif

# if defined(BOOST_MSVC)
#  if _MSC_VER <= 1200
#   define BOOST_MSVC6_OR_EARLIER 1
#  endif

# pragma warning (disable : 4786) // disable truncated debug symbols
# pragma warning (disable : 4251) // disable exported dll function
# pragma warning (disable : 4800) //'int' : forcing value to bool 'true' or 'false'
# pragma warning (disable : 4275) // non dll-interface class

# elif defined(__ICL) && __ICL < 600 // Intel C++ 5

#  pragma warning(disable: 985) // identifier was truncated in debug information

# endif

// The STLport puts all of the standard 'C' library names in std (as far as the
// user is concerned), but without it you need a fix if you're using MSVC or
// Intel C++
# if defined(BOOST_MSVC_STD_ITERATOR)
#  define BOOST_CSTD_
# else
#  define BOOST_CSTD_ std
# endif

/*****************************************************************************
 *
 *  Set up dll import/export options:
 *
 ****************************************************************************/

// backwards compatibility:
#ifdef BOOST_PYTHON_STATIC_LIB
#  define BOOST_PYTHON_STATIC_LINK
# elif !defined(BOOST_PYTHON_DYNAMIC_LIB)
#  define BOOST_PYTHON_DYNAMIC_LIB
#endif

#if defined(__MWERKS__) \
  || (defined(__DECCXX_VER) && __DECCXX_VER <= 60590002) \
  || (defined(__sgi) && defined(_COMPILER_VERSION) && _COMPILER_VERSION <= 730)
# define BOOST_PYTHON_NO_TEMPLATE_EXPORT
#endif

#if defined(BOOST_PYTHON_DYNAMIC_LIB) && (defined(_WIN32) || defined(__CYGWIN__))
#  if defined(BOOST_PYTHON_SOURCE)
#     define BOOST_PYTHON_DECL __declspec(dllexport)
#     define BOOST_PYTHON_BUILD_DLL
#  else
#     define BOOST_PYTHON_DECL __declspec(dllimport)
#  endif

// MinGW, at least, has some problems exporting template instantiations
#  if defined(__GNUC__) && __GNUC__ < 3 && !defined(__CYGWIN__)
#   define BOOST_PYTHON_NO_TEMPLATE_EXPORT
#  endif

#endif

#ifndef BOOST_PYTHON_DECL
#  define BOOST_PYTHON_DECL
#endif

#ifndef BOOST_PYTHON_EXPORT
# define BOOST_PYTHON_EXPORT extern
#endif 

#if !defined(BOOST_PYTHON_NO_TEMPLATE_EXPORT)
# define BOOST_PYTHON_EXPORT_CLASS_TEMPLATE(instantiation) BOOST_PYTHON_EXPORT template class BOOST_PYTHON_DECL instantiation
#else
# define BOOST_PYTHON_EXPORT_CLASS_TEMPLATE(instantiation) struct ThIsTyPeNeVeRuSeD
#endif

#if (defined(__DECCXX_VER) && __DECCXX_VER <= 60590031)
// Replace broken Tru64/cxx offsetof macro
# define BOOST_PYTHON_OFFSETOF(s_name, s_member) \
        ((size_t)__INTADDR__(&(((s_name *)0)->s_member)))
#else
# define BOOST_PYTHON_OFFSETOF offsetof
#endif

#endif // CONFIG_DWA052200_H_
