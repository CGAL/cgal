// A macro which helps removing "unused variable" warnings.

#ifndef CGAL_USE_H
#define CGAL_USE_H

namespace CGAL {

template < typename T >
inline void use(const T&) {}

} // namespace CGAL

#define CGAL_USE(x) CGAL::use(x);

#endif
