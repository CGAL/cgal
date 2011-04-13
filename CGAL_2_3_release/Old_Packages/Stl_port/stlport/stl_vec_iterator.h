/*
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

#ifndef __STLPORT_VEC_ITERATOR_H
# define __STLPORT_VEC_ITERATOR_H

#ifdef __STL_DEBUG

// _Vec_iter is being used by both vector and string


# if defined ( __STL_USE_ABBREVS )
#  define  _Vec_iter       _V__It
# endif


__STL_BEGIN_NAMESPACE
//============================================================

template <class _Tp>
bool __Vec_dereferenceable(const __owned_link& __that, const  _Tp* __ptr) {
    __stl_verbose_return(__that._Valid(), _StlMsg_INVALID_ITERATOR);
    _Tp* const * __start = (_Tp* const *)(__that._Owner()->_Owner());
    __stl_verbose_return((__ptr < *(__start+1)) && 
                         (__ptr >= *__start),
                          _StlMsg_NOT_DEREFERENCEABLE);    
    return true;
}

template <class _Tp>
bool __Vec_nonsingular(const __owned_link& __that, const  _Tp* __ptr) {
  __stl_verbose_return(__that._Valid(), _StlMsg_INVALID_ITERATOR);
  _Tp* const * __start = (_Tp* const *)(__that._Owner()->_Owner());
  __stl_verbose_return((__ptr <= *(__start+1)) && 
		       (__ptr >= *__start),
		       _StlMsg_SINGULAR_ITERATOR);    
    return true;
}

template <class _Tp, class _Traits>
struct _Vec_iter : public __owned_link {
public:
  typedef _Tp value_type;
  typedef typename _Traits::reference  reference;
  typedef typename _Traits::pointer    pointer;
  typedef ptrdiff_t difference_type;
  typedef random_access_iterator_tag iterator_category;
  pointer _M_iterator;
private:
  typedef _Vec_iter<_Tp, _Traits> _Self;
  typedef _Vec_iter<_Tp, _Nonconst_traits<_Tp> > _Nonconst_self;
  typedef _Vec_iter<_Tp, _Const_traits<_Tp> > _Const_self;
public:

  _Vec_iter() : __owned_link(0)  {}
  _Vec_iter(const __owned_list* __c, pointer __it) :
    __owned_link(__c), _M_iterator(__it) {}
  _Vec_iter(const _Nonconst_self& __it) :
    __owned_link(__it), _M_iterator(__it._M_iterator) {}
  ~_Vec_iter() {}
  reference operator*() const {
    __stl_debug_check(__Vec_dereferenceable(*this,_M_iterator));
    return *_M_iterator;
  }

  __STL_DEFINE_ARROW_OPERATOR
  
  _Self& operator++() {
    ++_M_iterator;
    __stl_debug_check(__Vec_nonsingular(*this,_M_iterator));
    return *this;
  }
  _Self operator++(int) {
    _Self __tmp = *this;
    ++_M_iterator;
    return __tmp;
  }
  _Self& operator--() {
    --_M_iterator;
    __stl_debug_check(__Vec_nonsingular(*this,_M_iterator));
    return *this;
  }
  _Self operator--(int) {
    _Self __tmp = *this;
    --_M_iterator;
    return __tmp;
  }
  difference_type operator-(const _Self& __y ) const {
    __stl_debug_check(__check_same_owner(*this, __y));
    return _M_iterator-__y._M_iterator;
  }
  _Self& operator+=(difference_type __n) {
    _M_iterator+=__n;
    __stl_debug_check(__Vec_nonsingular(*this,_M_iterator));
    return *this;
  }
  _Self& operator-=(difference_type __n) {
    return *this+=-__n;
  }
  _Self operator+(difference_type __n) const {
    _Self __tmp(*this);
    return __tmp += __n;
  }
  _Self operator-(difference_type __n) const {
    _Self __tmp(*this);
    return __tmp -= __n;
  }
  reference operator[](difference_type __n) const { return *(*this + __n); }
};

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator==(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));
  return __x._M_iterator==__y._M_iterator;
}

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator<(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
  __stl_debug_check(__check_same_owner(__x, __y));
  return __x._M_iterator < __y._M_iterator;
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator!=(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));
  return __x._M_iterator!=__y._M_iterator;
}

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator>(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
  return __y < __x;
}

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator>=(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
    return !(__x < __y);
}

template <class _Tp, class _Traits, class _Traits1>
inline bool 
operator<=(const _Vec_iter<_Tp, _Traits>& __x, const _Vec_iter<_Tp, _Traits1>& __y) {
    return !(__y < __x);
}

#else

template <class _Tp>
inline bool 
operator>(const _Vec_iter<_Tp, _Const_traits<_Tp> >& __x,
	  const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __y) { 
    return __y < __x;
}
template <class _Tp>
inline bool operator>=(const _Vec_iter<_Tp, _Const_traits<_Tp> >& __x,
		       const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __y) { 
    return !(__x < __y);
}
template <class _Tp>
inline bool operator<=(const _Vec_iter<_Tp, _Const_traits<_Tp> >& __x,
		       const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __y) { 
    return !(__y < __x);
}


template <class _Tp>
inline bool 
operator>(const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __x,
	  const _Vec_iter<_Tp, _Const_traits<_Tp> >& __y) { 
    return __y < __x;
}
template <class _Tp>
inline bool operator>=(const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __x,
		       const _Vec_iter<_Tp, _Const_traits<_Tp> >& __y) { 
    return !(__x < __y);
}
template <class _Tp>
inline bool operator<=(const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __x,
		       const _Vec_iter<_Tp, _Const_traits<_Tp> >& __y) { 
    return !(__y < __x);
}

template <class _Tp>
inline bool 
operator!=(const _Vec_iter<_Tp, _Const_traits<_Tp> >& __x, 
	   const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));
  return __x._M_iterator!=__y._M_iterator;
}

template <class _Tp>
inline bool 
operator!=(const _Vec_iter<_Tp, _Nonconst_traits<_Tp> >& __x, 
	   const _Vec_iter<_Tp, _Const_traits<_Tp> >& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));
  return __x._M_iterator!=__y._M_iterator;
}
#endif

template <class _Tp, class _Traits>
inline _Vec_iter<_Tp, _Traits> 
operator+(ptrdiff_t __n, const _Vec_iter<_Tp, _Traits>& __it) {
    _Vec_iter<_Tp, _Traits> __tmp(__it);
    return __tmp += __n;
}

#  if !defined (__STL_CLASS_PARTIAL_SPECIALIZATION)

template <class _Tp, class _Traits>
inline _Tp *
value_type(const _Vec_iter<_Tp, _Traits>&) { 
  return (_Tp*)0; 
}

template <class _Tp, class _Traits>
inline ptrdiff_t* distance_type(const  _Vec_iter<_Tp, _Traits>&) { return (ptrdiff_t*) 0; }

template <class _Tp, class _Traits>
inline random_access_iterator_tag iterator_category(const _Vec_iter<_Tp, _Traits>&) { 
    return random_access_iterator_tag();
}
#  endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// these exports are for basic_string
# if 0 // defined (__STL_USE_DECLSPEC)
__STL_EXPORT template class __STL_CLASS_DECLSPEC _Vec_iter<char, _Nonconst_traits<char> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC _Vec_iter<char, _Const_traits<char> >;
#  if defined (__STL_HAS_WCHAR_T)
__STL_EXPORT template class __STL_CLASS_DECLSPEC _Vec_iter<wchar_t, _Nonconst_traits<wchar_t> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC _Vec_iter<wchar_t, _Const_traits<wchar_t> >;
#  endif
# endif /* __STL_USE_DECLSPEC */


__STL_END_NAMESPACE

# endif /* __STL_DEBUG */

#endif /* INTERNAL_H */

// Local Variables:
// mode:C++
// End:

