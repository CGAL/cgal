/*
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

#ifndef __TYPE_TRAITS_H
#define __TYPE_TRAITS_H

/*
This header file provides a framework for allowing compile time dispatch
based on type attributes. This is useful when writing template code.
For example, when making a copy of an array of an unknown type, it helps
to know if the type has a trivial copy constructor or not, to help decide
if a memcpy can be used.

The class template __type_traits provides a series of typedefs each of
which is either __true_type or __false_type. The argument to
__type_traits can be any type. The typedefs within this template will
attain their correct values by one of these means:
    1. The general instantiation contain conservative values which work
       for all types.
    2. Specializations may be declared to make distinctions between types.
    3. Some compilers (such as the Silicon Graphics N32 and N64 compilers)
       will automatically provide the appropriate specializations for all
       types.

EXAMPLE:

//Copy an array of elements which have non-trivial copy constructors
template <class T> void copy(T* source, T* destination, int n, __false_type);
//Copy an array of elements which have trivial copy constructors. Use memcpy.
template <class T> void copy(T* source, T* destination, int n, __true_type);

//Copy an array of any type by using the most efficient copy mechanism
template <class T> inline void copy(T* source,T* destination,int n) {
   copy(source, destination, n,
        typename __type_traits<T>::has_trivial_copy_constructor());
}
*/


struct __true_type {
};

struct __false_type {
};

template <class _Tp>
struct __type_traits { 
   typedef __true_type     this_dummy_member_must_be_first;
                   /* Do not remove this member. It informs a compiler which
                      automatically specializes __type_traits that this
                      __type_traits template is special. It just makes sure that
                      things work if an implementation is using a template
                      called __type_traits for something unrelated. */

   /* The following restrictions should be observed for the sake of
      compilers which automatically produce type specific specializations 
      of this class:
          - You may reorder the members below if you wish
          - You may remove any of the members below if you wish
          - You must not rename members without making the corresponding
            name change in the compiler
          - Members you add will be treated like regular members unless
            you add the appropriate support in the compiler. */
 

   typedef __false_type    has_trivial_default_constructor;
   typedef __false_type    has_trivial_copy_constructor;
   typedef __false_type    has_trivial_assignment_operator;
   typedef __false_type    has_trivial_destructor;
   typedef __false_type    is_POD_type;
};



// Provide some specializations.  This is harmless for compilers that
//  have built-in __types_traits support, and essential for compilers
//  that don't.

#ifndef __STL_NO_BOOL

__STL_TEMPLATE_NULL struct __type_traits<bool> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#endif /* __STL_NO_BOOL */

__STL_TEMPLATE_NULL struct __type_traits<char> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#ifndef __STL_NO_SIGNED_BUILTINS

__STL_TEMPLATE_NULL struct __type_traits<signed char> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

# endif

__STL_TEMPLATE_NULL struct __type_traits<unsigned char> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#if defined ( __STL_HAS_WCHAR_T ) && ! defined (__STL_WCHAR_T_IS_USHORT)

__STL_TEMPLATE_NULL struct __type_traits<wchar_t> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#endif /* __STL_HAS_WCHAR_T */

__STL_TEMPLATE_NULL struct __type_traits<short> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<unsigned short> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<int> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<unsigned int> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<long> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<unsigned long> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#ifdef __STL_LONG_LONG

__STL_TEMPLATE_NULL struct __type_traits<long long> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<unsigned long long> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#endif /* __STL_LONG_LONG */

__STL_TEMPLATE_NULL struct __type_traits<float> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

__STL_TEMPLATE_NULL struct __type_traits<double> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

# if !defined ( __STL_NO_LONG_DOUBLE )
__STL_TEMPLATE_NULL struct __type_traits<long double> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};
# endif

#if defined( __STL_CLASS_PARTIAL_SPECIALIZATION) && \
  ! defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template <class _Tp>
struct __type_traits<_Tp*> {
   typedef __true_type    has_trivial_default_constructor;
   typedef __true_type    has_trivial_copy_constructor;
   typedef __true_type    has_trivial_assignment_operator;
   typedef __true_type    has_trivial_destructor;
   typedef __true_type    is_POD_type;
};

#  define __STL_DEFINE_ARROW_OPERATOR  pointer operator->() const { return &(operator*()); }

#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// Important pointers specializations

# define __STL_TYPE_TRAITS_POD_SPECIALIZE(_Type) \
__STL_TEMPLATE_NULL \
struct __type_traits<_Type> { \
   typedef __true_type    has_trivial_default_constructor; \
   typedef __true_type    has_trivial_copy_constructor;    \
   typedef __true_type    has_trivial_assignment_operator; \
   typedef __true_type    has_trivial_destructor; \
   typedef __true_type    is_POD_type; \
};

// the following is a workaround for arrow operator problems
#  if defined  ( __SGI_STL_NO_ARROW_OPERATOR ) 

// Proxy -> operator workaround for compilers that produce
// type checking errors on unused ->() operators

template <class _Ref, class _Ptr>
struct __arrow_op_dispatch {
  _Ptr _M_ptr;
  __arrow_op_dispatch(_Ref __r) : _M_ptr(&__r) {}
  _Ptr operator ->() const { return _M_ptr; }
};

#   if defined (__STL_NO_PROXY_ARROW_OPERATOR)

// User wants to disable proxy -> operators
#    define __STL_DEFINE_ARROW_OPERATOR
#    define __STL_ARROW_SPECIALIZE_WITH_PTRS(_Tp)

#   else

struct __arrow_op_dummy { int _M_data ; };

#  define __STL_ARROW_SPECIALIZE(_Tp)  \
__STL_TEMPLATE_NULL struct __arrow_op_dispatch<_Tp&, _Tp*> { \
  __arrow_op_dispatch(_Tp&) {} \
  __arrow_op_dummy operator ->() const { return __arrow_op_dummy(); } \
};

# ifdef __SUNPRO_CC
#  define __STL_ARROW_SPECIALIZE_WITH_PTRS(_Tp) \
__STL_ARROW_SPECIALIZE(_Tp) \
__STL_ARROW_SPECIALIZE(const _Tp) \
__STL_ARROW_SPECIALIZE(_Tp*) \
__STL_ARROW_SPECIALIZE(_Tp* const) \
__STL_ARROW_SPECIALIZE(const _Tp*) \
__STL_ARROW_SPECIALIZE(_Tp**) \
__STL_ARROW_SPECIALIZE(const _Tp**) \
__STL_ARROW_SPECIALIZE(_Tp* const *) \
__STL_ARROW_SPECIALIZE(_Tp***) \
__STL_ARROW_SPECIALIZE(const _Tp***)
# else
#  define __STL_ARROW_SPECIALIZE_WITH_PTRS(_Tp) \
__STL_ARROW_SPECIALIZE(_Tp) \
__STL_ARROW_SPECIALIZE(const _Tp) \
__STL_ARROW_SPECIALIZE(_Tp*) \
__STL_ARROW_SPECIALIZE(const _Tp*) \
__STL_ARROW_SPECIALIZE(_Tp**) \
__STL_ARROW_SPECIALIZE(const _Tp**) \
__STL_ARROW_SPECIALIZE(_Tp* const *) \
__STL_ARROW_SPECIALIZE(_Tp***) \
__STL_ARROW_SPECIALIZE(const _Tp***)
# endif

#  define __STL_DEFINE_ARROW_OPERATOR __arrow_op_dispatch<reference, pointer> operator->() const \
 { return __arrow_op_dispatch<reference, pointer>(this->operator*()); }

#  endif /* __STL_NO_PROXY_ARROW_OPERATOR */
# else
// Compiler can handle generic -> operator.
#  define __STL_ARROW_SPECIALIZE_WITH_PTRS(_Tp)
#  ifdef __BORLANDC__
#   define __STL_DEFINE_ARROW_OPERATOR  pointer operator->() const { return &(*(*this)); }
#  else
#   define __STL_DEFINE_ARROW_OPERATOR  pointer operator->() const { return &(operator*()); }
#  endif

# endif /* __SGI_STL_NO_ARROW_OPERATOR */

# define __STL_TYPE_TRAITS_POD_SPECIALIZE_V(_Type) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(_Type*) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(const _Type*) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(_Type**) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(_Type* const *) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(const _Type**) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(_Type***) \
__STL_TYPE_TRAITS_POD_SPECIALIZE(const _Type***)

# define __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE(_Type) \
__STL_TYPE_TRAITS_POD_SPECIALIZE_V(_Type) \
__STL_ARROW_SPECIALIZE_WITH_PTRS(_Type)

#  if !defined ( __STL_NO_BOOL )
__STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( bool )
#  endif
__STL_TYPE_TRAITS_POD_SPECIALIZE_V(void)
  __STL_ARROW_SPECIALIZE_WITH_PTRS(void*)
  __STL_ARROW_SPECIALIZE_WITH_PTRS(void* const)
# ifndef __STL_NO_SIGNED_BUILTINS
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( signed char )
# endif
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( char )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( unsigned char )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( short )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( unsigned short )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( int )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( unsigned int )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( long )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( unsigned long )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( float )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( double )
#  if !defined ( __STL_NO_LONG_DOUBLE )
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( long double )
#  endif
#  if defined ( __STL_LONG_LONG)
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( long long )
#  endif
#if defined ( __STL_HAS_WCHAR_T ) && ! defined (__STL_WCHAR_T_IS_USHORT)
  __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE( wchar_t )
#  endif

# undef __STL_ARROW_SPECIALIZE
# undef __STL_ARROW_SPECIALIZE_WITH_PTRS
# undef __STL_TYPE_TRAITS_POD_PTRS_SPECIALIZE
  // # undef __STL_TYPE_TRAITS_POD_SPECIALIZE
# undef __STL_TYPE_TRAITS_POD_SPECIALIZE_V

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// The following could be written in terms of numeric_limits.  
// We're doing it separately to reduce the number of dependencies.

template <class _Tp> struct _Is_integer {
  typedef __false_type _Integral;
};

#ifndef __STL_NO_BOOL

__STL_TEMPLATE_NULL struct _Is_integer<bool> {
  typedef __true_type _Integral;
};

#endif /* __STL_NO_BOOL */

__STL_TEMPLATE_NULL struct _Is_integer<char> {
  typedef __true_type _Integral;
};

#ifndef __STL_NO_SIGNED_BUILTINS

__STL_TEMPLATE_NULL struct _Is_integer<signed char> {
  typedef __true_type _Integral;
};
#endif

__STL_TEMPLATE_NULL struct _Is_integer<unsigned char> {
  typedef __true_type _Integral;
};

#if defined ( __STL_HAS_WCHAR_T ) && ! defined (__STL_WCHAR_T_IS_USHORT)

__STL_TEMPLATE_NULL struct _Is_integer<wchar_t> {
  typedef __true_type _Integral;
};

#endif /* __STL_HAS_WCHAR_T */

__STL_TEMPLATE_NULL struct _Is_integer<short> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<unsigned short> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<int> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<unsigned int> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<long> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<unsigned long> {
  typedef __true_type _Integral;
};

#ifdef __STL_LONG_LONG

__STL_TEMPLATE_NULL struct _Is_integer<long long> {
  typedef __true_type _Integral;
};

__STL_TEMPLATE_NULL struct _Is_integer<unsigned long long> {
  typedef __true_type _Integral;
};

#endif /* __STL_LONG_LONG */

#endif /* __TYPE_TRAITS_H */

// Local Variables:
// mode:C++
// End:
