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

# ifndef __STLPORT_DEBUG_C
#  define  __STLPORT_DEBUG_C

//==========================================================
// .c section
//  owned_list non-inline methods and global functions 
//==========================================================

#if defined ( __STL_ASSERTIONS )

# ifndef __STL_ASSERT_MSG_TRAILER
#  define __STL_ASSERT_MSG_TRAILER
# endif

// dwa 12/30/98 - if __STL_DEBUG_MESSAGE is defined, the user can supply own definition.
# if !defined( __STL_DEBUG_MESSAGE )
#   define __stl_debug_message __STLPORT_STD::__stl_debug_engine<bool>::_Message
# else
    extern  void __stl_debug_message(const char * format_str, ...);
# endif

// fbp: if __STL_DEBUG_TERMINATE is defined, the user can supply own definition.
# if !defined( __STL_DEBUG_TERMINATE )
#   define __stl_debug_terminate __STLPORT_STD::__stl_debug_engine<bool>::_Terminate
# else
    extern  void __stl_debug_terminate(void);
# endif

__STL_BEGIN_NAMESPACE
# define __STL_MESSAGE_TABLE_BODY = { \
"\n%s:%d STL error: %s\n" , \
"%s:%d STL assertion failure : %s\n" __STL_ASSERT_MSG_TRAILER, \
"\n%s:%d STL error : %s\n%s:%d STL assertion failure:     %s \n" __STL_ASSERT_MSG_TRAILER, \
"Invalid argument to operation (see operation documentation)",                  \
"Taking an iterator out of destroyed (or otherwise corrupted) container",       \
"Trying to extract an object out from empty container",\
"Past-the-end iterator could not be erased",  \
"Index out of bounds",  \
"Container doesn't own the iterator",  \
"Uninitialized or invalidated (by mutating operation) iterator used",  \
"Uninitialized or invalidated (by mutating operation) lefthand iterator in expression",  \
"Uninitialized or invalidated (by mutating operation) righthand iterator in expression",  \
"Iterators used in expression are from different owners",  \
"Iterator could not be dereferenced (past-the-end ?)",  \
"Range [first,last) is invalid",  \
"Iterator is not in range [first,last)",  \
"Range [first,last) is not in range [start,finish)",  \
"The advance would produce invalid iterator",  \
"Iterator is singular (advanced beyond the bounds ?)",  \
"Memory block deallocated twice",  \
"Deallocating a block that was never allocated",  \
"Deallocating a memory block allocated for another type",  \
"Size of block passed to deallocate() doesn't match block size",  \
"Pointer underrun - safety margin at front of memory block overwritten",  \
"Pointer overrrun - safety margin at back of memory block overwritten",   \
"Unknown problem" \
  }

# if ( __STL_STATIC_TEMPLATE_DATA > 0 )
template <class _Dummy>
const char* __stl_debug_engine<_Dummy>::_Message_table[_StlMsg_MAX]  __STL_MESSAGE_TABLE_BODY;

# else
__DECLARE_INSTANCE(const char*, __stl_debug_engine<bool>::_Message_table[_StlMsg_MAX],
		   __STL_MESSAGE_TABLE_BODY);

# endif
__STL_END_NAMESPACE

// abort()
#    include <cstdlib>

#  if !defined( __STL_DEBUG_MESSAGE )

#    include <cstdarg>
#    include <cstdio>

__STL_BEGIN_NAMESPACE


template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_Message(const char * __format_str, ...)
{

#  if defined(__MSL__) && __MSL__ <= 0x4011
	using __STL_VENDOR_STD::__files;
#  endif
	va_list __args;
	va_start( __args, __format_str );

# if defined (__STL_WINCE)

	TCHAR __buffer[512];
	wvsprintf(__buffer, _T(__format_str), __args);
	__STL_WINCE_TRACE(__buffer);
# else
	__STL_VENDOR_CSTD::vfprintf(stderr, __format_str, __args);

# endif /* WINCE */

# ifdef __STL_DEBUG_MESSAGE_POST
	__STL_DEBUG_MESSAGE_POST
# endif

       va_end(__args);

}

__STL_END_NAMESPACE

#  endif /* __STL_DEBUG_MESSAGE */

__STL_BEGIN_NAMESPACE


template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_IndexedError(int __error_ind, const char* __f, int __l)
{
  __stl_debug_message(_Message_table[_StlFormat_ERROR_RETURN], 
		      __f, __l, _Message_table[__error_ind]);
}

template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_VerboseAssert(const char* __expr, int __error_ind, const char* __f, int __l)
{
# if defined (_MFC_VER)
TRACE(_T(_Message_table[_StlFormat_VERBOSE_ASSERTION_FAILURE]), \
      __FILE__, __LINE__, _Message_table[__error_ind], __f, __l, # __expr );
    ASSERT(0);
# else
  __stl_debug_message(_Message_table[_StlFormat_VERBOSE_ASSERTION_FAILURE],
		      __f, __l, _Message_table[__error_ind], __f, __l, __expr);
# endif
  __stl_debug_terminate();
}

template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_Assert(const char* __expr, const char* __f, int __l)
{
# if defined (_MFC_VER)
  TRACE(_T(_Message_table[_StlFormat_ASSERTION_FAILURE]), __f, __l, __expr );
  ASSERT(0);
# else
  __stl_debug_message(_Message_table[_StlFormat_ASSERTION_FAILURE],__f, __l, __expr);
  __stl_debug_terminate();
# endif
}

// if exceptions are present, sends unique exception
// if not, calls abort() to terminate
template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_Terminate()
{
# if defined (__STL_USE_EXCEPTIONS) && ! defined (__STL_NO_DEBUG_EXCEPTIONS)
  throw __stl_debug_exception();
# elif defined (__STL_WINCE)
  TerminateProcess(GetCurrentProcess(), 0);
# else
  abort();
# endif
}

__STL_END_NAMESPACE

# endif /* __STL_ASSERTIONS */

#ifdef __STL_DEBUG

__STL_BEGIN_NAMESPACE

# ifdef __STL_THREADS
#  ifndef __STL_NEED_MUTABLE 
#   define __STL_ACQUIRE_LOCK(_Lock) _Lock._M_acquire_lock();
#   define __STL_RELEASE_LOCK(_Lock) _Lock._M_release_lock();
#  else
#   define __STL_ACQUIRE_LOCK(_Lock) ((_STL_mutex_lock&)_Lock)._M_acquire_lock();
#   define __STL_RELEASE_LOCK(_Lock) ((_STL_mutex_lock&)_Lock)._M_release_lock();
#  endif /* __STL_NEED_MUTABLE */
# else
#  define __STL_ACQUIRE_LOCK(_Lock)
#  define __STL_RELEASE_LOCK(_Lock)
# endif /* __STL_THREADS */



// [ i1, i2)
template <class _Iterator>
inline bool __in_range_aux(const _Iterator& __it, const _Iterator& __first,
                           const _Iterator& __last, random_access_iterator_tag) {
    return ( __it >= __first && 
             __it < __last);
}

template <class _Iterator1, class _Iterator>
# if defined (__STL_MSVC) && (__STL_MSVC >= 1100)
inline bool __in_range_aux(_Iterator1 __it, const _Iterator& __first,
# else
inline bool __in_range_aux(const _Iterator1& __it, const _Iterator& __first,
# endif
                           const _Iterator& __last, forward_iterator_tag) {
    _Iterator1 __i(__first);
    for (;  __i != __last && __i != __it; ++__i);
    return (__i!=__last);
}

# if defined (__STL_NONTEMPL_BASE_MATCH_BUG) /* OBSOLETE by inheritance */
template <class _Iterator1, class _Iterator>
inline bool __in_range_aux(const _Iterator1& __it, const _Iterator& __first,
                           const _Iterator& __last, bidirectional_iterator_tag) {
    _Iterator1 __i(__first);
    for (;  __i != __last && __i != __it; ++__i);
    return (__i !=__last);
}
# endif

//==========================================================
//  owned_list non-inline methods 
//==========================================================

template <class _Dummy>
void __stl_debug_engine<_Dummy>::_Invalidate_all(__owned_list* __l) {
  // crucial
  if (__l->_M_node._M_owner) {
    for (__owned_link*  __position = (__owned_link*)__l->_M_node._M_next; 
	 __position != 0; __position= (__owned_link*)__position->_M_next) {
      __position->_M_owner=0;
    }
    __l->_M_node._M_next =0;
  }
}

template <class _Dummy>
void __stl_debug_engine<_Dummy>::_Verify(const __owned_list* __l) {
  __stl_assert(__l->_M_node._Owner() != 0);
  for (__owned_link* __position = (__owned_link*)__l->_M_node._M_next; 
       __position != 0; __position= (__owned_link*)__position->_M_next) {
    __stl_assert(__position->_Owner()== __l);
  }
}

template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_Swap_owners(__owned_list& __x, __owned_list& __y, bool __swap_roots) {
  __x._Invalidate_all();
  __y._Invalidate_all();
  if (__swap_roots) {
    const __owned_list* __tmp = __x._M_node._M_owner;
    __x._M_node._M_owner=__y._M_node._M_owner;
    __y._M_node._M_owner=__tmp;
  }
}

template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_M_detach(__owned_list* __l, __owned_link* __c_node) {
  if (__l  != 0) {

    __stl_verbose_assert(__l->_Owner()!=0, _StlMsg_INVALID_CONTAINER);

    __STL_ACQUIRE_LOCK(__l->_M_lock)

    __owned_link* __prev, *__next;
   
    for (__prev = &__l->_M_node; (__next = __prev->_M_next) != __c_node; 
	 __prev = __next) {}
      
    __prev->_M_next = __c_node->_M_next;
    __c_node->_M_owner=0;

    __STL_RELEASE_LOCK(__l->_M_lock)
  }
}

template <class _Dummy>
void 
__stl_debug_engine<_Dummy>::_M_attach(__owned_list* __l, __owned_link* __c_node) {
  if (__l ==0) {
    (__c_node)->_M_owner = 0;    
  } else {
    __stl_verbose_assert(__l->_Owner()!=0, _StlMsg_INVALID_CONTAINER);
    __STL_ACQUIRE_LOCK(__l->_M_lock)
    __c_node->_M_owner = __l;
    __c_node->_M_next = __l->_M_node._M_next;
    __l->_M_node._M_next = __c_node;
    __STL_RELEASE_LOCK(__l->_M_lock)
  }
}

template <class _Dummy>
bool __stl_debug_engine<_Dummy>::_Check_same_owner( const __owned_link& __i1, 
						    const __owned_link& __i2)
{
  __stl_verbose_return(__i1._Valid(), _StlMsg_INVALID_LEFTHAND_ITERATOR);
  __stl_verbose_return(__i2._Valid(), _StlMsg_INVALID_RIGHTHAND_ITERATOR);
  __stl_verbose_return((__i1._Owner()==__i2._Owner()), _StlMsg_DIFFERENT_OWNERS);
  return true;
}

template <class _Dummy>
bool __stl_debug_engine<_Dummy>::_Check_same_owner_or_null( const __owned_link& __i1, 
							    const __owned_link& __i2)
{
  __stl_verbose_return(__i1._Owner()==__i2._Owner(), _StlMsg_DIFFERENT_OWNERS);
  return true;
}

template <class _Dummy>
bool __stl_debug_engine<_Dummy>::_Check_if_owner( const __owned_list * __l, const __owned_link& __it)
{
  const __owned_list* __owner_ptr = __it._Owner();
  __stl_verbose_return(__owner_ptr!=0, _StlMsg_INVALID_ITERATOR);
  __stl_verbose_return(__l==__owner_ptr, _StlMsg_NOT_OWNER);
  return true;
}

//==========================================================
//  global non-inline functions
//==========================================================

template <class _Iterator>
bool __check_range(const _Iterator& __first, const _Iterator& __last) {
    __stl_verbose_return(__valid_range(__first,__last), _StlMsg_INVALID_RANGE );
    return true;
}

template <class _Iterator>
bool __check_range(const _Iterator& __it, 
                   const _Iterator& __start, const _Iterator& __finish) {
    __stl_verbose_return(__in_range(__it,__start, __finish), 
                         _StlMsg_NOT_IN_RANGE_1);
    return true;
}

template <class _Iterator>
bool __check_range(const _Iterator& __first, const _Iterator& __last, 
                   const _Iterator& __start, const _Iterator& __finish) {
    __stl_verbose_return(__in_range(__first, __last, __start, __finish), 
                         _StlMsg_NOT_IN_RANGE_2);
    return true;
}

//===============================================================

template <class _Iterator>
void __invalidate_range(const __owned_list* __base, 
                        const _Iterator& __first,
                        const _Iterator& __last)
{
    typedef _Iterator* _Safe_iterator_ptr;
    typedef __owned_link _L_type;
    __STL_ACQUIRE_LOCK(__base->_M_lock)
    _L_type* __pos;
    _L_type* __prev;

    for (__prev = (_L_type*)&__base->_M_node, __pos= (_L_type*)__prev->_M_next; 
         __pos!=0;) {	    
        if ((!(&__first == (_Iterator*)__pos || &__last == (_Iterator*)__pos))
            &&  __in_range_aux(
			       *(_Iterator*)__pos,
			       __first,
			       __last,
			       __ITERATOR_CATEGORY(__first))) {
	  __pos->_M_owner = 0;
	  __pos = (_L_type*) (__prev->_M_next = __pos->_M_next);
	}
	else {
	  __prev = __pos;
	  __pos=(_L_type*)__pos->_M_next;
	}
    }
    __STL_RELEASE_LOCK(__base->_M_lock)    
}

template <class _Iterator>
void __invalidate_iterator(const __owned_list* __base, 
			   const _Iterator& __it)
{
    typedef __owned_link   _L_type;
    _L_type*  __position, *__prev;
    __STL_ACQUIRE_LOCK(__base->_M_lock)
    for (__prev = (_L_type*)&__base->_M_node, __position = (_L_type*)__prev->_M_next; 
         __position!= 0;) {
      // this requires safe iterators to be derived from __owned_link
       if ((__position != (_L_type*)&__it) && *((_Iterator*)__position)==__it) {
	    __position->_M_owner = 0;
	    __position = (_L_type*) (__prev->_M_next = __position->_M_next);
        }
       else {
	 __prev = __position;
	 __position=(_L_type*)__position->_M_next;
       }
    }
    __STL_RELEASE_LOCK(__base->_M_lock)
}

__STL_END_NAMESPACE

#endif /* __STL_DEBUG */

#endif

// Local Variables:
// mode:C++
// End:

