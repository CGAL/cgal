/*
 * Copyright (c) 1998
 * Silicon Graphics Computer Systems, Inc.
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

#ifndef __SGI_STL_BITSET_H
#define __SGI_STL_BITSET_H

// This implementation of bitset<> has a second template parameter,
// _WordT, which defaults to unsigned long.  *YOU SHOULD NOT USE
// THIS FEATURE*.  It is experimental, and it may be removed in
// future releases.

// A bitset of size N, using words of type _WordT, will have 
// N % (sizeof(_WordT) * CHAR_BIT) unused bits.  (They are the high-
// order bits in the highest word.)  It is a class invariant
// of class bitset<> that those unused bits are always zero.

// Most of the actual code isn't contained in bitset<> itself, but in the 
// base class _Base_bitset.  The base class works with whole words, not with
// individual bits.  This allows us to specialize _Base_bitset for the
// important special case where the bitset is only a single word.

// The C++ standard does not define the precise semantics of operator[].
// In this implementation the const version of operator[] is equivalent
// to test(), except that it does no range checking.  The non-const version
// returns a reference to a bit, again without doing any range checking.

#ifndef __SGI_STDEXCEPT
#  include <stdexcept>      
#endif

# ifndef __SGI_STL_INTERNAL_ALGOBASE_H
#  include <stl_algobase.h>
# endif

# ifndef __SGI_STL_INTERNAL_ALLOC_H
#  include <stl_alloc.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

# ifndef __SGI_STL_INTERNAL_UNINITIALIZED_H
#  include <stl_uninitialized.h>
# endif

# ifndef __STLPORT_STRING
#  include <string>
# endif

# ifndef __STLPORT_IOSTREAM
#  include <iostream>
# endif


#define __BITS_PER_WORDT(__wt) (CHAR_BIT*sizeof(__wt))
#define __BITSET_WORDS(__n,__wt) \
 ((__n) < 1 ? 1 : ((__n) + __BITS_PER_WORDT(__wt) - 1)/__BITS_PER_WORDT(__wt))

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1209
#endif

#if defined(__IBMCPP__)
// supress EDC3130: A constant is being used as a conditional expression
#pragma info(nocnd) 
#endif

// structure to aid in counting bits
template<bool __dummy> 
struct _Bit_count {
  static unsigned char _S_bit_count[256];
};

// Mapping from 8 bit unsigned integers to the index of the first one
// bit:
template<bool __dummy> 
struct _First_one {
  static unsigned char _S_first_one[256];
};

//
// Base class: general case.
//

template<size_t _Nw, class _WordT>
struct _Base_bitset {
  _WordT _M_w[_Nw];                // 0 is the least significant word.

  _Base_bitset( void ) { _M_do_reset(); }

  _Base_bitset(unsigned long __val);

  static size_t _S_whichword( size_t __pos ) {
    return __pos / __BITS_PER_WORDT(_WordT);
  }
  static size_t _S_whichbyte( size_t __pos ) {
    return (__pos % __BITS_PER_WORDT(_WordT)) / CHAR_BIT;
  }
  static size_t _S_whichbit( size_t __pos ) {
    return __pos % __BITS_PER_WORDT(_WordT);
  }
  static _WordT _S_maskbit( size_t __pos ) {
    return __STATIC_CAST(_WordT,1) << _S_whichbit(__pos);
  }

  _WordT& _M_getword(size_t __pos)       { return _M_w[_S_whichword(__pos)]; }
  _WordT  _M_getword(size_t __pos) const { return _M_w[_S_whichword(__pos)]; }

  _WordT& _M_hiword()       { return _M_w[_Nw - 1]; }
  _WordT  _M_hiword() const { return _M_w[_Nw - 1]; }

  void _M_do_and(const _Base_bitset<_Nw,_WordT>& __x) {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] &= __x._M_w[__i];
    }
  }

  void _M_do_or(const _Base_bitset<_Nw,_WordT>& __x) {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] |= __x._M_w[__i];
    }
  }

  void _M_do_xor(const _Base_bitset<_Nw,_WordT>& __x) {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] ^= __x._M_w[__i];
    }
  }

  void _M_do_left_shift(size_t __shift);

  void _M_do_right_shift(size_t __shift);

  void _M_do_flip() {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] = ~_M_w[__i];
    }
  }

  void _M_do_set() {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] = ~__STATIC_CAST(_WordT,0);
    }
  }

  void _M_do_reset() {
    for ( size_t __i = 0; __i < _Nw; __i++ ) {
      _M_w[__i] = 0;
    }
  }

  bool _M_is_equal(const _Base_bitset<_Nw,_WordT>& __x) const {
    for (size_t __i = 0; __i < _Nw; ++__i) {
      if (_M_w[__i] != __x._M_w[__i])
        return false;
    }
    return true;
  }

  bool _M_is_any() const {
    for ( size_t __i = 0; __i < __BITSET_WORDS(_Nw,_WordT); __i++ ) {
      if ( _M_w[__i] != __STATIC_CAST(_WordT,0) )
        return true;
    }
    return false;
  }

  size_t _M_do_count() const {
    size_t __result = 0;
    const unsigned char* __byte_ptr = (const unsigned char*)_M_w;
    const unsigned char* __end_ptr = (const unsigned char*)(_M_w+_Nw);

    while ( __byte_ptr < __end_ptr ) {
      __result += _Bit_count<true>::_S_bit_count[*__byte_ptr];
      __byte_ptr++;
    }
    return __result;
  }

  unsigned long _M_do_to_ulong() const; 

  // find first "on" bit
  size_t _M_do_find_first(size_t __not_found) const;

  // find the next "on" bit that follows "prev"
  size_t _M_do_find_next(size_t __prev, size_t __not_found) const;
};


//
// Base class: specialization for a single word.
//

# if defined (__STL_CLASS_PARTIAL_SPECIALIZATION) && \
! defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template<class _WordT>
struct _Base_bitset<1, _WordT> {
  _WordT _M_w;

  _Base_bitset( void ) { _M_do_reset(); }

  _Base_bitset(unsigned long __val) {
    _M_do_reset();
    const size_t __n = min(sizeof(unsigned long)*CHAR_BIT,
			   __BITS_PER_WORDT(_WordT));
    for(size_t __i = 0; __i < __n; ++__i, __val >>= 1)
      if ( __val & 0x1 )
	_M_w |= _S_maskbit(__i);
  }
  
  static size_t _S_whichword( size_t __pos ) {
    return __pos / __BITS_PER_WORDT(_WordT);
  }
  static size_t _S_whichbyte( size_t __pos ) {
    return (__pos % __BITS_PER_WORDT(_WordT)) / CHAR_BIT;
  }
  static size_t _S_whichbit( size_t __pos ) {
    return __pos % __BITS_PER_WORDT(_WordT);
  }
  static _WordT _S_maskbit( size_t __pos ) {
    return (__STATIC_CAST(_WordT,1)) << _S_whichbit(__pos);
  }

  _WordT& _M_getword(size_t)       { return _M_w; }
  _WordT  _M_getword(size_t) const { return _M_w; }

  _WordT& _M_hiword()       { return _M_w; }
  _WordT  _M_hiword() const { return _M_w; }

  void _M_do_and(const _Base_bitset<1,_WordT>& __x) { _M_w &= __x._M_w; }
  void _M_do_or(const _Base_bitset<1,_WordT>& __x)  { _M_w |= __x._M_w; }
  void _M_do_xor(const _Base_bitset<1,_WordT>& __x) { _M_w ^= __x._M_w; }
  void _M_do_left_shift(size_t __shift)     { _M_w <<= __shift; }
  void _M_do_right_shift(size_t __shift)    { _M_w >>= __shift; }
  void _M_do_flip()                       { _M_w = ~_M_w; }
  void _M_do_set()                        { _M_w = ~__STATIC_CAST(_WordT,0); }
  void _M_do_reset()                      { _M_w = 0; }

  bool _M_is_equal(const _Base_bitset<1,_WordT>& __x) const {
    return _M_w == __x._M_w;
  }
  bool _M_is_any() const {
    return _M_w != 0;
  }

  size_t _M_do_count() const {
    size_t __result = 0;
    const unsigned char* __byte_ptr = (const unsigned char*)&_M_w;
    const unsigned char* __end_ptr = ((const unsigned char*)&_M_w)+sizeof(_M_w);
    while ( __byte_ptr < __end_ptr ) {
      __result += _Bit_count<true>::_S_bit_count[*__byte_ptr];
      __byte_ptr++;
    }
    return __result;
  }

  unsigned long _M_do_to_ulong() const {
    if (sizeof(_WordT) <= sizeof(unsigned long))
        return __STATIC_CAST(unsigned long,_M_w);
    else {
      const _WordT __mask = __STATIC_CAST(_WordT,__STATIC_CAST(unsigned long,-1));
      if (_M_w & ~__mask) 
        __STL_THROW(overflow_error("bitset"));
      return __STATIC_CAST(unsigned long,_M_w);
    }
  }

  size_t _M_do_find_first(size_t __not_found) const;

  // find the next "on" bit that follows "prev"
  size_t _M_do_find_next(size_t __prev, size_t __not_found) const; 

};


#  if !defined( __STL_MEMBER_SPECIALIZATION_BUG )
//
// One last specialization: _M_do_to_ulong() and the constructor from
// unsigned long are very simple if the bitset consists of a single 
// word of type unsigned long.
//

__STL_TEMPLATE_NULL
inline unsigned long 
_Base_bitset<1, unsigned long>::_M_do_to_ulong() const { return _M_w; }

__STL_TEMPLATE_NULL
inline _Base_bitset<(size_t)1, unsigned long>::_Base_bitset(unsigned long __val) {
  _M_w = __val;
}
#  endif // __STL_MEMBER_SPECIALIZATION_BUG

# endif /* __STL_CLASS_PARTIAL_SPECIALIZATION  */

// ------------------------------------------------------------
// Helper class to zero out the unused high-order bits in the highest word.
#if defined( __STL_CLASS_PARTIAL_SPECIALIZATION) && \
! defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template <class _WordT, size_t _Extrabits> struct _Sanitize {
  static void _M_do_sanitize(_WordT& __val)
    { __val &= ~((~__STATIC_CAST(_WordT,0)) << _Extrabits); }
};

template <class _WordT> struct _Sanitize<_WordT, 0> {
  static void _M_do_sanitize(_WordT) {}
};

#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */

template <class _WordT, size_t _Extrabits> struct _Sanitize {
  static void _M_do_sanitize(_WordT& __val) {
    if (_Extrabits != 0)
      __val &= ~((~__STATIC_CAST(_WordT,0)) << _Extrabits);
  }
};

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// ------------------------------------------------------------
// Class bitset.
//   _Nb may be any nonzero number of type size_t.
//   Type _WordT may be any unsigned integral type.

typedef unsigned long __stl_ulong;

# if defined (__STL_DEFAULT_TYPE_PARAM)
template<size_t _Nb, class _WordT = __stl_ulong>
# else
#  define bitset __bitset
template<size_t _Nb, class _WordT>
# endif
class bitset : public _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT> 
{
private:
  typedef _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT> _Base;

  void _M_do_sanitize() {
    _Sanitize<_WordT,_Nb%__BITS_PER_WORDT(_WordT) >
      ::_M_do_sanitize(this->_M_hiword());
  }

public:
  // bit reference:
  struct reference {
  typedef _Base_bitset<_Nb,_WordT> _Bitset_base;
  typedef bitset<_Nb,_WordT> _Bitset;
    //    friend _Bitset;
    _WordT *_M_wp;
    size_t _M_bpos;

    // should be left undefined
    reference() {}

    reference( _Bitset& __b, size_t __pos ) {
      _M_wp = &__b._M_getword(__pos);
      _M_bpos = _Bitset_base::_S_whichbit(__pos);
    }

  public:
    ~reference() {}

    // for b[i] = __x;
    reference& operator=(bool __x) {
      if ( __x )
        *_M_wp |= _Bitset_base::_S_maskbit(_M_bpos);
      else
        *_M_wp &= ~_Bitset_base::_S_maskbit(_M_bpos);

      return *this;
    }

    // for b[i] = b[__j];
    reference& operator=(const reference& __j) {
      if ( (*(__j._M_wp) & _Bitset_base::_S_maskbit(__j._M_bpos)) )
        *_M_wp |= _Bitset_base::_S_maskbit(_M_bpos);
      else
        *_M_wp &= ~_Bitset_base::_S_maskbit(_M_bpos);

      return *this;
    }

    // flips the bit
    bool operator~() const { return (*(_M_wp) & _Bitset_base::_S_maskbit(_M_bpos)) == 0; }

    // for __x = b[i];
    operator bool() const { return (*(_M_wp) & _Bitset_base::_S_maskbit(_M_bpos)) != 0; }

    // for b[i].flip();
    reference& flip() {
      *_M_wp ^= _Bitset_base::_S_maskbit(_M_bpos);
      return *this;
    }
  };

  friend struct reference;

  // 23.3.5.1 constructors:
  bitset() {}
  bitset(unsigned long __val) : 
    _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT>(__val) {}
# ifdef __STL_MEMBER_TEMPLATES
  template<class _CharT, class _Traits, class _Alloc>
  explicit bitset(const basic_string<_CharT,_Traits,_Alloc>& __s,
                  size_t __pos = 0)
    : _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT>() 
  {
    if (__pos > __s.size()) 
      __STL_THROW(out_of_range("bitset"));
    _M_copy_from_string(__s, __pos,
                        basic_string<_CharT, _Traits, _Alloc>::npos);
  }
  template<class _CharT, class _Traits, class _Alloc>
  bitset(const basic_string<_CharT, _Traits, _Alloc>& __s,
          size_t __pos,
          size_t __n)
  : _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT>() 
  {
    if (__pos > __s.size()) 
      __STL_THROW(out_of_range("bitset"));
    _M_copy_from_string(__s, __pos, __n);
  }
#else /* __STL_MEMBER_TEMPLATES */
  explicit bitset(const string& __s,
                  size_t __pos = 0,
                  size_t __n = (size_t)-1) 
    : _Base_bitset<__BITSET_WORDS(_Nb,_WordT), _WordT>() 
  {
    if (__pos > __s.size()) 
      __STL_THROW(out_of_range("bitset"));
    _M_copy_from_string(__s, __pos, __n);
  }
#endif /* __STL_MEMBER_TEMPLATES */

  // 23.3.5.2 bitset operations:
  bitset<_Nb,_WordT>& operator&=(const bitset<_Nb,_WordT>& __rhs) {
    this->_M_do_and(__rhs);
    return *this;
  }

  bitset<_Nb,_WordT>& operator|=(const bitset<_Nb,_WordT>& __rhs) {
    this->_M_do_or(__rhs);
    return *this;
  }

  bitset<_Nb,_WordT>& operator^=(const bitset<_Nb,_WordT>& __rhs) {
    this->_M_do_xor(__rhs);
    return *this;
  }

  bitset<_Nb,_WordT>& operator<<=(size_t __pos) {
    this->_M_do_left_shift(__pos);
    this->_M_do_sanitize();
    return *this;
  }

  bitset<_Nb,_WordT>& operator>>=(size_t __pos) {
    this->_M_do_right_shift(__pos);
    this->_M_do_sanitize();
    return *this;
  }

  //
  // Extension:
  // Versions of single-bit set, reset, flip, test with no range checking.
  //

  bitset<_Nb,_WordT>& _Unchecked_set(size_t __pos) {
    this->_M_getword(__pos) |= _Base::_S_maskbit(__pos);
    return *this;
  }

  bitset<_Nb,_WordT>& _Unchecked_set(size_t __pos, int __val) {
    if (__val)
      this->_M_getword(__pos) |= _Base::_S_maskbit(__pos);
    else
      this->_M_getword(__pos) &= ~ _Base::_S_maskbit(__pos);

    return *this;
  }

  bitset<_Nb,_WordT>& _Unchecked_reset(size_t __pos) {
    this->_M_getword(__pos) &= ~ _Base::_S_maskbit(__pos);
    return *this;
  }

  bitset<_Nb,_WordT>& _Unchecked_flip(size_t __pos) {
    this->_M_getword(__pos) ^= _Base::_S_maskbit(__pos);
    return *this;
  }

  bool _Unchecked_test(size_t __pos) const {
    return (this->_M_getword(__pos) & _Base::_S_maskbit(__pos)) != __STATIC_CAST(_WordT,0);
  }

  // Set, reset, and flip.

  bitset<_Nb,_WordT>& set() {
    this->_M_do_set();
    this->_M_do_sanitize();
    return *this;
  }

  bitset<_Nb,_WordT>& set(size_t __pos) {
    if (__pos >= _Nb)
      __STL_THROW(out_of_range("bitset"));

    return _Unchecked_set(__pos);
  }

  bitset<_Nb,_WordT>& set(size_t __pos, int __val) {
    if (__pos >= _Nb)
      __STL_THROW(out_of_range("bitset"));

    return _Unchecked_set(__pos, __val);
  }

  bitset<_Nb,_WordT>& reset() {
    this->_M_do_reset();
    return *this;
  }

  bitset<_Nb,_WordT>& reset(size_t __pos) {
    if (__pos >= _Nb)
      __STL_THROW(out_of_range("bitset"));

    return _Unchecked_reset(__pos);
  }

  bitset<_Nb,_WordT>& flip() {
    this->_M_do_flip();
    this->_M_do_sanitize();
    return *this;
  }

  bitset<_Nb,_WordT>& flip(size_t __pos) {
    if (__pos >= _Nb)
      __STL_THROW(out_of_range("bitset"));

    return _Unchecked_flip(__pos);
  }

  bitset<_Nb,_WordT> operator~() const { 
    return bitset<_Nb,_WordT>(*this).flip();
  }

  // element access:
  //for b[i];
  reference operator[](size_t __pos) { return reference(*this,__pos); }
  bool operator[](size_t __pos) const { return _Unchecked_test(__pos); }

  unsigned long to_ulong() const { return this->_M_do_to_ulong(); }

#if defined (__STL_MEMBER_TEMPLATES) && \
    defined (__STL_EXPLICIT_FUNCTION_TMPL_ARGS)
  template <class _CharT, class _Traits, class _Alloc>
  basic_string<_CharT, _Traits, _Alloc> to_string() const {
    basic_string<_CharT, _Traits, _Alloc> __result;
    _M_copy_to_string(__result);
    return __result;
  }
#else
  string to_string() const {
    string __result;
    _M_copy_to_string(__result);
    return __result;
  }
#endif /* __STL_EXPLICIT_FUNCTION_TMPL_ARGS */

  size_t count() const { return this->_M_do_count(); }

  size_t size() const { return _Nb; }

  bool operator==(const bitset<_Nb,_WordT>& __rhs) const {
    return this->_M_is_equal(__rhs);
  }
  bool operator!=(const bitset<_Nb,_WordT>& __rhs) const {
    return !this->_M_is_equal(__rhs);
  }

  bool test(size_t __pos) const {
    if (__pos > _Nb)
      __STL_THROW(out_of_range("bitset"));

    return _Unchecked_test(__pos);
  }

  bool any() const { return this->_M_is_any(); }
  bool none() const { return !this->_M_is_any(); }

  bitset<_Nb,_WordT> operator<<(size_t __pos) const { 
    bitset<_Nb,_WordT> __result(*this);
    __result <<= __pos ;  return __result; 
  }
  bitset<_Nb,_WordT> operator>>(size_t __pos) const { 
    bitset<_Nb,_WordT> __result(*this);
    __result >>= __pos ;  return __result; 
  }

  //
  // EXTENSIONS: bit-find operations.  These operations are
  // experimental, and are subject to change or removal in future
  // versions.
  // 

  // find the index of the first "on" bit
  size_t _Find_first() const 
    { return this->_M_do_find_first(_Nb); }

  // find the index of the next "on" bit after prev
  size_t _Find_next( size_t __prev ) const 
    { return this->_M_do_find_next(__prev, _Nb); }

//
// Definitions of should-be non-inline member functions.
//
# if defined (__STL_MEMBER_TEMPLATES)
  template<class _CharT, class _Traits, class _Alloc>
    void _M_copy_from_string(const basic_string<_CharT,_Traits,_Alloc>& __s,
			     size_t __pos,
			     size_t __n)
#else
    void _M_copy_from_string(const string& __s,
			     size_t __pos,
			     size_t __n)
#endif
    {
      reset();
      size_t __tmp = _Nb;
      const size_t __Nbits = min(__tmp, min(__n, __s.size() - __pos));
      for ( size_t __i= 0; __i < __Nbits; ++__i) {
	switch(__s[__pos + __Nbits - __i - 1]) {
	case '0':
	  break;
	case '1':
	  set(__i);
	  break;
	default:
	  __STL_THROW(invalid_argument("bitset"));
	}
      }
    }
  
# if defined (__STL_MEMBER_TEMPLATES)
  template <class _CharT, class _Traits, class _Alloc>
    void _M_copy_to_string(basic_string<_CharT, _Traits, _Alloc>& __s) const
# else
    void _M_copy_to_string(string& __s) const
# endif
    {
      __s.assign(_Nb, '0');
      
      for (size_t __i = 0; __i < _Nb; ++__i) 
	if (_Unchecked_test(__i))
	  __s[_Nb - 1 - __i] = '1';
    }

# if defined (__STL_NON_TYPE_TMPL_PARAM_BUG)
  bitset<_Nb,_WordT> operator&(const bitset<_Nb,_WordT>& __y) {
    bitset<_Nb,_WordT> __result(*this);
    __result &= __y;
    return __result;
  }
  bitset<_Nb,_WordT> operator|(const bitset<_Nb,_WordT>& __y) {
    bitset<_Nb,_WordT> __result(*this);
    __result |= __y;
    return __result;
  }
  bitset<_Nb,_WordT> operator^(const bitset<_Nb,_WordT>& __y) {
    bitset<_Nb,_WordT> __result(*this);
    __result ^= __y;
    return __result;
  }
# endif 
};

// ------------------------------------------------------------

//
// 23.3.5.3 bitset operations:
//

# if ! defined (__STL_NON_TYPE_TMPL_PARAM_BUG)
template <size_t _Nb, class _WordT>
inline bitset<_Nb,_WordT> operator&(const bitset<_Nb,_WordT>& __x,
                                    const bitset<_Nb,_WordT>& __y) {
  bitset<_Nb,_WordT> __result(__x);
  __result &= __y;
  return __result;
}


template <size_t _Nb, class _WordT>
inline bitset<_Nb,_WordT> operator|(const bitset<_Nb,_WordT>& __x,
                                    const bitset<_Nb,_WordT>& __y) {
  bitset<_Nb,_WordT> __result(__x);
  __result |= __y;
  return __result;
}

template <size_t _Nb, class _WordT>
inline bitset<_Nb,_WordT> operator^(const bitset<_Nb,_WordT>& __x,
                                    const bitset<_Nb,_WordT>& __y) {
  bitset<_Nb,_WordT> __result(__x);
  __result ^= __y;
  return __result;
}

#if defined ( __STL_USE_NEW_IOSTREAMS )
template <class _CharT, class _Traits, size_t _Nb, class _WordT>
basic_istream<_CharT, _Traits>&
operator>>(basic_istream<_CharT, _Traits>& __is, bitset<_Nb,_WordT>& __x);


template <class _CharT, class _Traits, size_t _Nb, class _WordT>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& __os,
           const bitset<_Nb,_WordT>& __x);

#elif ! defined ( __STL_USE_NO_IOSTREAMS )

template <size_t _Nb, class _WordT>
istream&
operator>>(istream& __is, bitset<_Nb,_WordT>& __x);

template <size_t _Nb, class _WordT>
inline ostream& operator<<(ostream& __os, const bitset<_Nb,_WordT>& __x) {
  string __tmp;
  __x._M_copy_to_string(__tmp);
  return __os << __tmp;
}

#endif

# endif /* __STL_NON_TYPE_TMPL_PARAM_BUG */

#  undef  bitset

# if defined (__STL_DEFAULT_TYPE_PARAM)
#  define __bitset__ bitset
# else
#  define __bitset__ __bitset
// provide a "default" bitset adaptor
template <size_t _Nb>
class bitset : public __bitset__<_Nb, __stl_ulong>
{
public:
#   define _BS_SUPER __bitset__<_Nb, __stl_ulong >
  typedef _BS_SUPER  _Super;
  typedef typename _BS_SUPER::reference reference;                                 \
  __IMPORT_SUPER_COPY_ASSIGNMENT(bitset, bitset<_Nb>, _BS_SUPER)
  bitset() {}
  explicit bitset(const unsigned long __val) : _BS_SUPER(__val) { }
  explicit bitset(const string& __s,
		  size_t __pos = 0,
		  size_t __n = (size_t)-1) : _BS_SUPER(__s, __pos, __n) {}
  ~bitset() {}

#  if defined (__STL_BASE_MATCH_BUG)
  bool operator==(const bitset<_Nb>& __y) {
    return operator == ((const _Super&)*this,(const _Super&)__y);
  }
  bool operator<(const bitset<_Nb>& __y) {
    return operator < ((const _Super&)*this,(const _Super&)__y);
  }
#  endif /* __STL_BASE_MATCH_BUG */

};

#  undef _BS_SUPER

# endif /* __STL_DEFAULT_TYPE_PARAM */

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1209
#endif

__STL_END_NAMESPACE

#  undef __BITS_PER_WORDT
#  undef __BITSET_WORDS

#if defined(__IBMCPP__)
#pragma info(restore)
#endif

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_bitset.c>
# endif

#endif /* __SGI_STL_BITSET_H */


// Local Variables:
// mode:C++
// End:

