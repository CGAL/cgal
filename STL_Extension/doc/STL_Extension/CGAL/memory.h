#include <memory>

/*!
  \file memory.h
*/

/** 
 * \ingroup PkgStlExtension
 *
 * A define for the allocator used by %CGAL. This is only defined if
 * there is no user defined version before \ref memory.h is included
 * the first time.
 *
 */
#define CGAL_ALLOCATOR(T) std::allocator< T >
