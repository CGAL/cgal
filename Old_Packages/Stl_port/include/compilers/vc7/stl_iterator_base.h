#ifndef CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC
#include <iterator>
// #include <utility>

#define CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(T)                    \
namespace std {                                                        \
    template <>                                                        \
    struct iterator_traits<const T*> {                                 \
	typedef random_access_iterator_tag iterator_category;          \
	typedef T                          value_type;                 \
	typedef ptrdiff_t                  difference_type;            \
	typedef const T*                   pointer;                    \
	typedef const T&                   reference;                  \
    };                                                                 \
    template <>                                                        \
    struct iterator_traits<T*> {                                       \
	typedef random_access_iterator_tag iterator_category;          \
	typedef T                          value_type;                 \
	typedef ptrdiff_t                  difference_type;            \
	typedef T*                         pointer;                    \
	typedef T&                         reference;                  \
    };                                                                 \
}

// add more stuff accoring to taste...
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int)
namespace std {
    template <>                                                        \
    struct iterator_traits<const void*> {                                 \
	typedef random_access_iterator_tag iterator_category;          \
	typedef ptrdiff_t                  difference_type;            \
	typedef const void*                   pointer;                    \
    };                                                                 \
    template <>                                                        \
    struct iterator_traits<void*> {                                       \
	typedef random_access_iterator_tag iterator_category;          \
	typedef ptrdiff_t                  difference_type;            \
	typedef void*                         pointer;                    \
    };                                                                 \
}

  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(void*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(void**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short**)

#endif

