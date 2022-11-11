#include <memory>

/**
 * \ingroup PkgSTLExtensionRef
 *
 * A define for the allocator used by \cgal. This is only defined if
 * there is no user defined version before `memory.h` is included
 * the first time.
 *
 */
#define CGAL_ALLOCATOR(T) std::allocator< T >
