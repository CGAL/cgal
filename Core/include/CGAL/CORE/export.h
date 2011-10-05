#ifndef CGAL_CORE_EXPORT_H
#define CGAL_CORE_EXPORT_H

#include <boost/config.hpp>

#if defined(BOOST_MSVC) && ( ! defined(CGAL_EXPORTS) )

#if defined(CGAL_Core_EXPORTS) // add by CMake or in cpp files of the dll
#define	CGAL_CORE_EXPORT __declspec (dllexport)
#define CGAL_CORE_EXPIMP_TEMPLATE
#else
#define CGAL_CORE_EXPORT __declspec (dllimport)
#define CGAL_CORE_EXPIMP_TEMPLATE extern
#endif
 
#else 

#define  CGAL_CORE_EXPORT
#define CGAL_CORE_EXPIMP_TEMPLATE 
#endif

#endif //  CGAL_CORE_EXPORT_H


