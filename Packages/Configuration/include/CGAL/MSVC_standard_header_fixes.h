#ifndef CGAL_MSVC_STANDARD_HEADER_FIXES_H
#define CGAL_MSVC_STANDARD_HEADER_FIXES_H

#pragma warning(once: 4291)

#include <cmath>
namespace std {
	using ::fabs;
	using ::sqrt;
	using ::log;
}

namespace std{
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


#include <cstddef>
namespace std{
using ::size_t;
using ::ptrdiff_t;
}

#include <cstdlib>
namespace std{
using ::atoi;
}

#include <iterator>
namespace std{
template < class T, class Dist> struct input_iterator: public iterator<input_iterator_tag,T,Dist> {};
template <class T, class Dist> struct forward_iterator : public iterator<forward_iterator_tag,T,Dist> {};
template <class T, class Dist> struct bidirectional_iterator : public iterator<bidirectional_iterator_tag,T,Dist> {};
template <class T, class Dist> struct random_access_iterator : public iterator<random_access_iterator_tag,T,Dist> {};
struct output_iterator : public iterator<output_iterator_tag,void,void> {};
}

#include <cstring>
namespace std{
using ::strcat;
using ::strcpy;
}

#endif // CGAL_MSVC_STANDARD_HEADER_FIXES_H
