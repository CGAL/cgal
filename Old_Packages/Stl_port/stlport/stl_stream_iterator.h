/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996-1998
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

/* NOTE: This is an internal header file, included by other STL headers.
 *   You should not attempt to use it directly.
 */

#ifndef __SGI_STL_INTERNAL_STREAM_ITERATOR_H
#define __SGI_STL_INTERNAL_STREAM_ITERATOR_H

#ifndef __SGI_STL_INTERNAL_ITERATOR_BASE_H
# include <stl_iterator_base.h>
#endif

// streambuf_iterators predeclarations must appear first
#ifndef __STLPORT_IOSFWD
# include <iosfwd>
#endif

// istream_iterator and ostream_iterator look very different if we're
// using new, templatized iostreams than if we're using the old cfront
// version.

#if defined (__STL_USE_NEW_IOSTREAMS) 

# if  (defined (__STL_USE_OWN_NAMESPACE) || defined (__STLPORT_NEW_IOSTREAMS) || !defined(__STL_USE_NAMESPACES))

# if defined (__STLPORT_NEW_IOSTREAMS)
#  include <istream>
# endif

__STL_BEGIN_NAMESPACE

#  ifndef __STL_LIMITED_DEFAULT_TEMPLATES
template <class _Tp, 
          class _CharT = __STL_DEFAULTCHAR, class _Traits = char_traits<_CharT>,
          class _Dist = ptrdiff_t> 
#   define __ISI_TMPL_HEADER_ARGUMENTS class _Tp, class _CharT, class _Traits, class _Dist
#   define __ISI_TMPL_ARGUMENTS _Tp, _CharT, _Traits, _Dist
#  else

#   if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS) && ! defined (__STL_DEFAULT_TYPE_PARAM)
#    define __ISI_TMPL_HEADER_ARGUMENTS class _Tp
#    define __ISI_TMPL_ARGUMENTS        _Tp
template <class _Tp>
#   else
#    define __ISI_TMPL_HEADER_ARGUMENTS class _Tp, class _Dist
#    define __ISI_TMPL_ARGUMENTS        _Tp, _Dist
template <class _Tp,__DFL_TYPE_PARAM(_Dist, ptrdiff_t)>
#   endif /* __STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS */

#  endif /* __STL_LIMITED_DEFAULT_TEMPLATES */

class istream_iterator {
# ifdef __STL_LIMITED_DEFAULT_TEMPLATES
  typedef char _CharT;
  typedef char_traits<char> _Traits;
#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS) && ! defined (__STL_DEFAULT_TYPE_PARAM)
  typedef ptrdiff_t _Dist;
#  endif
# endif
  typedef istream_iterator< __ISI_TMPL_ARGUMENTS > _Self;
public:
  typedef _CharT                         char_type;
  typedef _Traits                        traits_type;
  typedef basic_istream<_CharT, _Traits> istream_type;

  typedef input_iterator_tag             iterator_category;
  typedef _Tp                            value_type;
  typedef _Dist                          difference_type;
  typedef const _Tp*                     pointer;
  typedef const _Tp&                     reference;

  istream_iterator() : _M_stream(0), _M_ok(false) {}
  istream_iterator(istream_type& __s) : _M_stream(&__s) { _M_read(); }

  reference operator*() const { return _M_value; }
  pointer operator->() const { return &(operator*()); }

  _Self& operator++() { 
    _M_read(); 
    return *this;
  }
  _Self operator++(int)  {
    _Self __tmp = *this;
    _M_read();
    return __tmp;
  }

  bool _M_equal(const _Self& __x) const
    { return (_M_ok == __x._M_ok) && (!_M_ok || _M_stream == __x._M_stream); }

private:
  istream_type* _M_stream;
  _Tp _M_value;
  bool _M_ok;

  void _M_read() {
    _M_ok = (_M_stream && *_M_stream) ? true : false;
    if (_M_ok) {
      *_M_stream >> _M_value;
      _M_ok = *_M_stream ? true : false;
    }
  }
};

#ifndef __STL_LIMITED_DEFAULT_TEMPLATES
template <class _Tp,
          class _CharT = __STL_DEFAULTCHAR, class _Traits = char_traits<_CharT> >
#else
template <class _Tp>
#endif
class ostream_iterator {
# ifdef __STL_LIMITED_DEFAULT_TEMPLATES
  typedef char _CharT;
  typedef char_traits<char> _Traits;
  typedef ostream_iterator<_Tp> _Self;
# else
  typedef ostream_iterator<_Tp, _CharT, _Traits> _Self;
# endif
public:
  typedef _CharT                         char_type;
  typedef _Traits                        traits_type;
  typedef basic_ostream<_CharT, _Traits> ostream_type;

  typedef output_iterator_tag            iterator_category;
  typedef void                           value_type;
  typedef void                           difference_type;
  typedef void                           pointer;
  typedef void                           reference;

  ostream_iterator(ostream_type& __s) : _M_stream(&__s), _M_string(0) {}
  ostream_iterator(ostream_type& __s, const _CharT* __c) 
    : _M_stream(&__s), _M_string(__c)  {}
  _Self& operator=(const _Tp& __value) { 
    *_M_stream << __value;
    if (_M_string) *_M_stream << _M_string;
    return *this;
  }
  _Self& operator*() { return *this; }
  _Self& operator++() { return *this; } 
  _Self& operator++(int) { return *this; } 
private:
  ostream_type* _M_stream;
  const _CharT* _M_string;
};

# ifndef __STL_CLASS_PARTIAL_SPECIALIZATION
#  ifdef __STL_LIMITED_DEFAULT_TEMPLATES
template <class _Tp>
inline output_iterator_tag 
iterator_category(const ostream_iterator<_Tp>&) {
  return output_iterator_tag();
}
# else

template <class _Tp, class _CharT, class _Traits>
inline output_iterator_tag 
iterator_category(const ostream_iterator<_Tp, _CharT, _Traits>&) {
  return output_iterator_tag();
}
#  endif /* __STL_LIMITED_DEFAULT_TEMPLATES */
# endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */


# if defined (__STLPORT_NEW_IOSTREAMS)

// We do not read any characters until operator* is called.  The first
// time operator* is called, it calls getc.  Subsequent calls to getc 
// return a cached character, and calls to operator++ use snextc.  Before
// operator* or operator++ has been called, _M_is_initialized is false.
template<class _CharT, class _Traits>
class istreambuf_iterator
{
public:
  typedef _CharT                           char_type;
  typedef _Traits                          traits_type;
  typedef typename _Traits::int_type       int_type;
  typedef basic_streambuf<_CharT, _Traits> streambuf_type;
  typedef basic_istream<_CharT, _Traits>   istream_type;

  typedef input_iterator_tag               iterator_category;
  typedef _CharT                           value_type;
  typedef typename _Traits::off_type       difference_type;
  typedef const _CharT*                    pointer;
  typedef const _CharT&                    reference;

public:
  istreambuf_iterator(streambuf_type* __p = 0) { this->_M_init(__p); }
  istreambuf_iterator(istream_type& __is) { this->_M_init(__is.rdbuf()); }

  char_type operator*() const 
    { return _M_is_initialized ? _M_c : _M_dereference_aux(); }

  istreambuf_iterator<_CharT, _Traits>& operator++() { this->_M_nextc(); return *this; }

  istreambuf_iterator  operator++(int) {
    if (!_M_is_initialized)
      _M_postincr_aux();
    istreambuf_iterator<_CharT, _Traits> __tmp = *this;
    this->_M_nextc();
    return __tmp;
  }

  bool equal(const istreambuf_iterator<_CharT, _Traits>& __i) const {
    return this->_M_is_initialized && __i._M_is_initialized
      ? this->_M_eof == __i._M_eof
      : this->_M_equal_aux(__i);
  }

private:
  void _M_init(streambuf_type* __p) {
    _M_buf = __p;
    _M_eof = !__p;
    _M_is_initialized = _M_eof;
  }

  void _M_getc() const {
    int_type __c = _M_buf->sgetc();
    // SUNPRO has problems with mutable bitfields
# if !defined (__STL_NEED_MUTABLE) && ! defined (__SUNPRO_CC)
    _M_c = traits_type::to_char_type(__c);
    _M_eof = traits_type::eq_int_type(__c, traits_type::eof());
    _M_is_initialized = true;
# else
    typedef istreambuf_iterator<_CharT,_Traits> _Self;
    _Self* that = __CONST_CAST(_Self*, this);
    that->_M_c = traits_type::to_char_type(__c);
    that->_M_eof = traits_type::eq_int_type(__c, traits_type::eof());
    that->_M_is_initialized = true;
# endif
  }

  char_type _M_dereference_aux() const {
    this->_M_getc();
    return _M_c;
  }

  bool _M_equal_aux(const istreambuf_iterator<_CharT, _Traits>& __i) const {
    if (!this->_M_is_initialized)
      this->_M_getc();
    if (!__i._M_is_initialized)
      __i._M_getc();
    
    return this->_M_eof == __i._M_eof;
  }

  void _M_postincr_aux() {
    this->_M_getc();
  }

  void _M_nextc() {
    int_type __c = _M_buf->snextc();
    _M_c = traits_type::to_char_type(__c);    
    _M_eof = traits_type::eq_int_type(__c, traits_type::eof());
    _M_is_initialized = true;
  }

private:
  streambuf_type* _M_buf;
  mutable _CharT _M_c;
  mutable bool _M_eof : 1;
  mutable bool _M_is_initialized : 1;
};


template<class _CharT, class _Traits>
inline bool operator==(const istreambuf_iterator<_CharT, _Traits>& __x,
                       const istreambuf_iterator<_CharT, _Traits>& __y) {
  return __x.equal(__y);
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template<class _CharT, class _Traits>
inline bool operator!=(const istreambuf_iterator<_CharT, _Traits>& __x,
                       const istreambuf_iterator<_CharT, _Traits>& __y) {
  return !__x.equal(__y);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

// The default template argument is declared in iosfwd
template<class _CharT, class _Traits>
class ostreambuf_iterator
{
public:
  typedef _CharT                           char_type;
  typedef _Traits                          traits_type;
  typedef typename _Traits::int_type       int_type;
  typedef basic_streambuf<_CharT, _Traits> streambuf_type;
  typedef basic_ostream<_CharT, _Traits>   ostream_type;

  typedef output_iterator_tag              iterator_category;
  typedef void                             value_type;
  typedef void                             difference_type;
  typedef void                             pointer;
  typedef void                             reference;

public:
  ostreambuf_iterator(streambuf_type* __buf) __STL_NOTHROW : _M_buf(__buf), _M_ok(__buf!=0) {}
  ostreambuf_iterator(ostream_type& __o) __STL_NOTHROW
    : _M_buf(__o.rdbuf()), _M_ok(__o.rdbuf()!=0) {}

  ostreambuf_iterator<_CharT, _Traits>& operator=(char_type __c) {
    _M_ok = _M_ok && !traits_type::eq_int_type(_M_buf->sputc(__c),
                                               traits_type::eof());
    return *this;
  }    
  
  ostreambuf_iterator<_CharT, _Traits>& operator*()     { return *this; }
  ostreambuf_iterator<_CharT, _Traits>& operator++()    { return *this; }
  ostreambuf_iterator<_CharT, _Traits>& operator++(int) { return *this; }

  bool failed() const { return !_M_ok; }

private:
  streambuf_type* _M_buf;
  bool _M_ok;
};

# ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _CharT, class _Traits>
inline input_iterator_tag 
iterator_category(const istreambuf_iterator<_CharT, _Traits>&) {
  return input_iterator_tag();
}

template <class _CharT, class _Traits>
inline streamoff 
difference_type(const istreambuf_iterator<_CharT, _Traits>&) {
  typedef streamoff __off_type;
  return __off_type();
}

template <class _CharT, class _Traits>
inline _CharT* 
value_type(const istreambuf_iterator<_CharT, _Traits>&) {
  return (_CharT*)0;
}


template <class _CharT, class _Traits>
inline output_iterator_tag 
iterator_category(const ostreambuf_iterator<_CharT, _Traits>&) {
  return output_iterator_tag();
}

#   endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

#  endif /* __STLPORT_NEW_IOSTREAMS */

# endif /* namespaces */

# elif ! defined ( __STL_USE_NO_IOSTREAMS )

__STL_BEGIN_NAMESPACE

#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS) && ! defined (__STL_DEFAULT_TYPE_PARAM)
#  define __ISI_TMPL_HEADER_ARGUMENTS class _Tp
#  define __ISI_TMPL_ARGUMENTS        _Tp
template <class _Tp>
#  else
#  define __ISI_TMPL_HEADER_ARGUMENTS class _Tp, class _Dist
#  define __ISI_TMPL_ARGUMENTS        _Tp, _Dist
template <class _Tp, __DFL_TYPE_PARAM(_Dist, ptrdiff_t)>
#  endif
class istream_iterator {
protected:
  istream* _M_stream;
  _Tp _M_value;
  bool _M_end_marker;
  void _M_read() {
    _M_end_marker = (*_M_stream) ? true : false;
    // _M_end_marker = !_M_stream->fail();
    if (_M_end_marker) *_M_stream >> _M_value;
    _M_end_marker = (*_M_stream) ? true : false;
    //    _M_end_marker = !_M_stream->fail();
}
public:
  typedef input_iterator_tag  iterator_category;
  typedef _Tp                 value_type;
  typedef _Dist               difference_type;
  typedef const _Tp*          pointer;
  typedef const _Tp&          reference;

  istream_iterator() : _M_stream(&cin), _M_end_marker(false) {}
  istream_iterator(istream& __s) : _M_stream(&__s) { _M_read(); }
  reference operator*() const { return _M_value; }

  __STL_DEFINE_ARROW_OPERATOR

  istream_iterator< __ISI_TMPL_ARGUMENTS >& operator++() { 
    _M_read(); 
    return *this;
  }
  istream_iterator< __ISI_TMPL_ARGUMENTS > operator++(int)  {
    istream_iterator< __ISI_TMPL_ARGUMENTS > __tmp = *this;
    _M_read();
    return __tmp;
  }
  inline bool _M_equal(const istream_iterator< __ISI_TMPL_ARGUMENTS >& __y) const {
    return (_M_stream == __y._M_stream &&
	    _M_end_marker == __y._M_end_marker) ||
      _M_end_marker == false && __y._M_end_marker == false;
  }
};

template <class _Tp>
class ostream_iterator {
protected:
  ostream* _M_stream;
  const char* _M_string;
public:
  typedef output_iterator_tag iterator_category;
# ifdef __STL_CLASS_PARTIAL_SPECIALIZATION
  typedef void                value_type;
  typedef void                difference_type;
  typedef void                pointer;
  typedef void                reference;
# endif
  ostream_iterator(ostream& __s) : _M_stream(&__s), _M_string(0) {}
  ostream_iterator(ostream& __s, const char* __c) 
    : _M_stream(&__s), _M_string(__c)  {}
  ostream_iterator<_Tp>& operator=(const _Tp& __value) { 
    *_M_stream << __value;
    if (_M_string) *_M_stream << _M_string;
    return *this;
  }
  ostream_iterator<_Tp>& operator*() { return *this; }
  ostream_iterator<_Tp>& operator++() { return *this; } 
  ostream_iterator<_Tp>& operator++(int) { return *this; } 
};

#ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _Tp>
inline output_iterator_tag 
iterator_category(const ostream_iterator<_Tp>&) {
  return output_iterator_tag();
}

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

#endif /* __STL_USE_NEW_IOSTREAMS */

#ifndef __STL_USE_NO_IOSTREAMS

# if (defined (__STL_USE_OWN_NAMESPACE) || !defined(__STL_USE_NAMESPACES) || !defined(__STL_USE_NEW_IOSTREAMS) )

template < __ISI_TMPL_HEADER_ARGUMENTS >
inline bool 
operator==(const istream_iterator< __ISI_TMPL_ARGUMENTS >& __x,
           const istream_iterator< __ISI_TMPL_ARGUMENTS >& __y) {
  return __x._M_equal(__y);
}

#  ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template < __ISI_TMPL_HEADER_ARGUMENTS >
inline bool 
operator!=(const istream_iterator< __ISI_TMPL_ARGUMENTS >& __x,
           const istream_iterator< __ISI_TMPL_ARGUMENTS >& __y) {
  return !__x._M_equal(__y);
}

#  endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

# ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template < __ISI_TMPL_HEADER_ARGUMENTS >
inline input_iterator_tag 
iterator_category(const istream_iterator< __ISI_TMPL_ARGUMENTS >&)
{
  return input_iterator_tag();
}
template < __ISI_TMPL_HEADER_ARGUMENTS >
inline _Tp* 
value_type(const istream_iterator< __ISI_TMPL_ARGUMENTS >&) { return (_Tp*) 0; }

#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS) && ! defined (__STL_DEFAULT_TYPE_PARAM)
template < __ISI_TMPL_HEADER_ARGUMENTS >
inline ptrdiff_t* 
distance_type(const istream_iterator< __ISI_TMPL_ARGUMENTS >&) { return (ptrdiff_t*)0; }
#  else
template < __ISI_TMPL_HEADER_ARGUMENTS >
inline _Dist* 
distance_type(const istream_iterator< __ISI_TMPL_ARGUMENTS >&) { return (_Dist*)0; }
#  endif /* __STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS */

# endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

# endif /* (!defined(__STL_USE_NAMESPACES) || defined (__STL_USE_OWN_NAMESPACE)) */

__STL_END_NAMESPACE

#  undef __ISI_TMPL_HEADER_ARGUMENTS
#  undef __ISI_TMPL_ARGUMENTS

# endif /* NO_STREAMS */

#endif /* __SGI_STL_INTERNAL_STREAM_ITERATOR_H */

// Local Variables:
// mode:C++
// End:
