#ifndef CGAL_MINMAX_FIX_H
#define CGAL_MINMAX_FIX_H

#include <cmath>
#include <cstdlib>
#include <cstddef>

#ifdef atoi
#undef atoi
#endif

namespace std {

using ::fabs;
using ::atoi;
using ::size_t;
using ::sqrt;

template <class T>
inline const T& min(const T& a, const T&b)
{ return (a<b) ? a : b;}

template <class T, class Cmp>
inline const T& min(const T& a, const T&b, Cmp cmp)
{ return (cmp(b,a)) ? a : b;}

template <class T>
inline const T& max(const T& a, const T&b)
{ return (a<b) ? b : a;}

}

#endif

