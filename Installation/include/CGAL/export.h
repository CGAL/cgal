#ifndef CGAL_EXPORT_H
#define CGAL_EXPORT_H

#include <boost/config.hpp>

#if defined(BOOST_MSVC)

#if defined(CGAL_EXPORTS) // add by CMake or in cpp files of the dll
#define	CGAL_EXPORT __declspec (dllexport)
#define CGAL_EXPIMP_TEMPLATE
#else
#define CGAL_EXPORT __declspec (dllimport)
#define CGAL_EXPIMP_TEMPLATE extern
#endif
 
#else 

#define  CGAL_EXPORT
#define CGAL_EXPIMP_TEMPLATE 
#endif

#endif //  CGAL_EXPORT_H


