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

#ifndef __STLPORT_DEBUG_H
# define __STLPORT_DEBUG_H

# if defined (__STL_ASSERTIONS) || defined (__STL_DEBUG)

#ifndef __STL_CONFIG_H
# include <stl_config.h>
#endif

# if !defined(__STL_FILE__)
#   define __STL_FILE__ __FILE__
# endif
 
enum {
  _StlFormat_ERROR_RETURN,
  _StlFormat_ASSERTION_FAILURE,
  _StlFormat_VERBOSE_ASSERTION_FAILURE,
  _StlMsg_INVALID_ARGUMENT,
  _StlMsg_INVALID_CONTAINER,
  _StlMsg_EMPTY_CONTAINER,
  _StlMsg_ERASE_PAST_THE_END,
  _StlMsg_OUT_OF_BOUNDS,
  _StlMsg_NOT_OWNER,
  _StlMsg_INVALID_ITERATOR,
  _StlMsg_INVALID_LEFTHAND_ITERATOR,
  _StlMsg_INVALID_RIGHTHAND_ITERATOR,
  _StlMsg_DIFFERENT_OWNERS     ,
  _StlMsg_NOT_DEREFERENCEABLE  ,
  _StlMsg_INVALID_RANGE        ,
  _StlMsg_NOT_IN_RANGE_1       ,
  _StlMsg_NOT_IN_RANGE_2       ,
  _StlMsg_INVALID_ADVANCE      ,
  _StlMsg_SINGULAR_ITERATOR    ,
  // debug alloc messages
  _StlMsg_DBA_DELETED_TWICE    ,
  _StlMsg_DBA_NEVER_ALLOCATED  ,
  _StlMsg_DBA_TYPE_MISMATCH    ,
  _StlMsg_DBA_SIZE_MISMATCH    ,
  _StlMsg_DBA_UNDERRUN         ,
  _StlMsg_DBA_OVERRUN          ,
  _StlMsg_UNKNOWN
  /* _StlMsg_MAX */
};

/* have to hardcode that ;() */
# define _StlMsg_MAX 27

__STL_BEGIN_NAMESPACE

// This class is unique (not inherited from exception),
// to disallow catch in anything but (...)
struct __stl_debug_exception {
  // no members
};

class __owned_link;
class __owned_list;

template <class _Dummy>
struct __stl_debug_engine {

  // Basic routine to report any debug message
  // Use __STL_DEBUG_MESSAGE to override
  static void _Message(const char * format_str, ...);

  // Micsellanous function to report indexed error message
  static void _IndexedError(int __ind, const char* __f, int __l);

  // Basic assertion report mechanism.
  // Reports failed assertion via __stl_debug_message and calls _Terminate
  // if __STL_DEBUG_TERMINATE is specified, calls __stl_debug_terminate instead
  static void _Assert(const char* __expr, const char* __f, int __l);

  // The same, with additional diagnostics
  static void _VerboseAssert(const char* __expr, int __error_ind, const char* __f, int __l);

  // If exceptions are present, sends unique exception
  // If not, calls abort() to terminate
  // Use __STL_DEBUG_TERMINATE to override
  static void _Terminate();

  // owned_list/link delegate non-inline functions here

  static bool _Check_same_owner( const __owned_link& __i1, 
				 const __owned_link& __i2);
  static bool _Check_same_owner_or_null( const __owned_link& __i1, 
					 const __owned_link& __i2);
  static bool _Check_if_owner( const __owned_list*, const __owned_link&);

  static void _Verify(const __owned_list*);
  
  static void _Swap_owners(__owned_list&, __owned_list& , bool __swap_roots);
 
  static void _Invalidate_all(__owned_list*);
  
  static void _M_detach(__owned_list*, __owned_link*);

  static void _M_attach(__owned_list*, __owned_link*);

  // debug messages and formats
  static const char* _Message_table[_StlMsg_MAX];
};

// This is a hook to instantiate STLport exports in a designated DLL
# if defined (__STL_USE_DECLSPEC)
__STL_EXPORT template struct __STL_CLASS_DECLSPEC __stl_debug_engine<bool>;
# endif /* __STL_USE_DECLSPEC */

typedef __stl_debug_engine<bool> __stl_debugger;

__STL_END_NAMESPACE

#   define __stl_assert(expr) \
     do { if (!(expr)) { \
             __STLPORT_STD::__stl_debugger::_Assert( # expr, __STL_FILE__, __LINE__); \
             } \
        } while (0)

# endif /* __STL_ASSERTIONS || __STL_DEBUG */


// this section is for __STL_DEBUG only 
#if defined ( __STL_DEBUG )

#  define __stl_verbose_assert(expr,__diag_num) \
    do { if (!(expr)) { __STLPORT_STD::__stl_debug_engine<bool>::_VerboseAssert\
                                 ( # expr,  __diag_num, __STL_FILE__, __LINE__ ); \
         } \
       } while(0)
#  define __stl_debug_check           __stl_assert
#  define __stl_debug_do(expr)        expr

#if defined(__MRC__)		//*TY 01/15/1999 - MrCpp optimizer confuses on this debugging code
#pragma options opt none
#endif

#  define __stl_verbose_return(__expr,__diag_num) if (!(__expr)) { \
       __STLPORT_STD::__stl_debug_engine<bool>::_IndexedError(__diag_num, __FILE__ , __LINE__); \
       return false; }


#if defined (__STL_THREADS) && ! defined (__SGI_STL_INTERNAL_THREADS_H)
# include <stl_threads.h>
#endif

#ifndef __SGI_STL_INTERNAL_ITERATOR_BASE_H
# include <stl_iterator_base.h>
#endif

__STL_BEGIN_NAMESPACE

//=============================================================
template <class __Iterator>
inline bool __valid_range(const __Iterator& __i1 ,const __Iterator& __i2, 
                          random_access_iterator_tag) { 
    return __i1<=__i2;
}

template <class __Iterator>
inline bool __valid_range(const __Iterator& __i1 ,const __Iterator& __i2, 
                          bidirectional_iterator_tag) { 
    // check if comparable
    bool __dummy(__i1==__i2);
    return (__dummy==__dummy); 
}

template <class __Iterator>
inline bool __valid_range(const __Iterator& __i1 ,const __Iterator& __i2, forward_iterator_tag) { 
    // check if comparable
    bool __dummy(__i1==__i2);
    return (__dummy==__dummy); 
}

template <class __Iterator>
inline bool __valid_range(const __Iterator&,const __Iterator&, input_iterator_tag) { 
    return false; 
}

template <class __Iterator>
inline bool __valid_range(const __Iterator&,const __Iterator&, output_iterator_tag) { 
    return false; 
}

template <class __Iterator>
inline bool __valid_range(const __Iterator& __i1, const __Iterator& __i2) { 
    return __valid_range(__i1,__i2,__ITERATOR_CATEGORY(__i1));
}

// Note : that means in range [i1, i2].
template <class __Iterator>
inline bool __in_range(const __Iterator& __it, const __Iterator& __i1,
		       const __Iterator& __i2) { 
    return __valid_range(__i1,__it,__ITERATOR_CATEGORY(__i1)) && 
        __valid_range(__it,__i2,__ITERATOR_CATEGORY(__it));
}

template <class __Iterator>
inline bool __in_range(const __Iterator& __first, const __Iterator& __last, 
                       const __Iterator& __start, const __Iterator& __finish) { 
    return __valid_range(__first,__last,__ITERATOR_CATEGORY(__first)) &&
        __valid_range(__start,__first,__ITERATOR_CATEGORY(__first)) && 
        __valid_range(__last,__finish,__ITERATOR_CATEGORY(__last));
}

//==========================================================


class __STL_CLASS_DECLSPEC __owned_link {
private:
  typedef __owned_link  _Self;
  typedef _Self* _Self_ptr;
  typedef __owned_list  _Owner_type;
  typedef __owned_list* _Owner_ptr;
public:
  inline __owned_link(const _Owner_type* __c);
  __owned_link() : _M_owner(0) {}
  
  const _Owner_type* _Owner() const { 
    return _M_owner; 
  }
  void _Set_owner(const _Owner_type* __o)  { 
    _M_owner= __o; 
  }
  bool  _Valid() const { 
    return _M_owner !=0; 
  }
  inline __owned_link(const _Self& __rhs);
  inline _Self& operator=(const _Self& __rhs);
  inline ~__owned_link();
  void _Invalidate() { _M_owner=0; }
  void _Link_to_self() { _M_next= 0; }
  _Self* _Next() { return _M_next; }
  const _Self* _Next() const { return _M_next; }
public:
  _Self* _M_next;
  const _Owner_type* _M_owner;
};


class __STL_CLASS_DECLSPEC __owned_list {
private:
  typedef __owned_list _Self;
  typedef _Self* _Self_ptr;
  typedef __owned_link  _Node;
  typedef _Node* _Node_ptr;
  typedef void* _Owner_ptr;
  typedef const void* _Owner_const_ptr;
  // should never be called, should be left undefined,
  // but some compilers complain about it ;(
  __owned_list(const _Self&){}
  void operator=(const _Self&) {}
public:
  __owned_list() {
  }
  ~__owned_list() {
    _Invalidate_all();
    // that prevents detach
    _M_node._Invalidate();
  }
  _Owner_const_ptr _Owner() const { 
    return (_Owner_const_ptr)_M_node._M_owner; 
  }
  _Owner_ptr _Owner() { 
    return (_Owner_ptr)_M_node._M_owner; 
  }
  bool  _Valid() const { 
    return _M_node._M_owner!=0; 
  }
  void _Invalidate() { _M_node._M_owner=0; }
  void _Safe_init(const void* __c) { 
    _M_node._M_owner = (_Self_ptr)__c; 
    _M_node._M_next=0;
  }
  void _Orphan(_Node& __it) {
    if ( __it._M_owner && __it._M_next != &__it ) {
      __stl_debugger::_M_detach(this, &__it);
    }
    __it._M_next = &__it;
    __it._M_owner = this;
  }
  
  _Node* _First() { return _M_node._Next(); }
  _Node* _Last() { return 0 ; }
  const _Node* _First() const { return (_Node*)_M_node._M_next; }
  const _Node* _Last() const { return 0 ;}
  
  void _Verify() const {
    __stl_debugger::_Verify(this); 
  }
  
  void _Swap_owners(_Self& __y, bool __swap_roots =false) { 
    __stl_debugger::_Swap_owners(*this, __y, __swap_roots); 
  }
 
  void _Invalidate_all() { 
    __stl_debugger::_Invalidate_all(this);
  }
  
  void _M_detach(const _Node* __node) const { 
    __stl_debugger::_M_detach(__CONST_CAST(__owned_list*, this), __CONST_CAST(__owned_link*, __node));
  }
  
  mutable _Node                     _M_node; 
# ifdef __STL_THREADS
  mutable _STL_mutex_lock           _M_lock;
# endif
  
private:
  void _M_attach(const _Node* __n) const { 
    __stl_debugger::_M_attach(__CONST_CAST(__owned_list*, this), __CONST_CAST(__owned_link*,__n)); 
  }
  friend class __owned_link;
  friend struct __stl_debug_engine<bool>;
};


//==========================================================

inline __owned_link::__owned_link(const __owned_list* __c): _M_owner(0) {
  __c->_M_attach(this);
}

inline __owned_link::__owned_link(const __owned_link& __rhs) : _M_owner(0) {
  __rhs._M_owner->_M_attach(this);
}

inline __owned_link& __owned_link::operator=(const __owned_link& __rhs) {
  const _Owner_type* __new_owner = __rhs._Owner();
  const _Owner_type* __old_owner = _Owner();
  if ( __old_owner != __new_owner ) {
    __old_owner->_M_detach(this);
    __new_owner->_M_attach(this);
  }
  return *this;
}

inline __owned_link::~__owned_link() {
  _M_owner->_M_detach(this);
}


//==========================================================

// forward declaratioins

template <class _Iterator>
bool __check_range(const _Iterator&, const _Iterator&);
template <class _Iterator>
bool __check_range(const _Iterator&, 
                   const _Iterator&, const _Iterator&);
template <class _Iterator>
bool __check_range(const _Iterator&, const _Iterator& , 
                   const _Iterator&, const _Iterator& );

template <class _Iterator>
void __invalidate_range(const __owned_list* __base, 
                        const _Iterator& __first,
                        const _Iterator& __last);

template <class _Iterator>
void __invalidate_iterator(const __owned_list* __base, 
                           const _Iterator& __it);

//============================================================

inline bool 
__check_same_owner( const __owned_link& __i1, const __owned_link& __i2) {
  return __stl_debugger::_Check_same_owner(__i1,__i2);
}
inline bool 
__check_same_owner_or_null( const __owned_link& __i1, const __owned_link& __i2) {
  return __stl_debugger::_Check_same_owner_or_null(__i1,__i2);
}

template <class Iterator>
inline bool __check_if_owner( const __owned_list* __owner, 
		       const Iterator& __it) {
  return __stl_debugger::_Check_if_owner(__owner, (const __owned_link&)__it);
}

__STL_END_NAMESPACE

# endif /* __STL_DEBUG */

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stldebug.c>
# endif

#endif /* DEBUG_H */

// Local Variables:
// mode:C++
// End:
