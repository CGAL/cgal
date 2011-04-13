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
#ifndef __STL_ALLOC_C
#define __STL_ALLOC_C

// Specialised debug form of malloc which does not provide "false"
// memory leaks when run with debug CRT libraries.
#if defined(__STL_MSVC) && __STL_MSVC>=1020 && defined(_DEBUG)
#  include <crtdbg.h>
#  define   __STL_CHUNK_MALLOC(s)         _malloc_dbg(s, _CRT_BLOCK, __FILE__, __LINE__)
#else	// !_DEBUG
# ifdef __STL_NODE_ALLOC_USE_MALLOC
#  define   __STL_CHUNK_MALLOC(s)         ::malloc(s)
# else
#  define   __STL_CHUNK_MALLOC(s)         __stl_new(s)
# endif
#endif	// !_DEBUG

__STL_BEGIN_NAMESPACE

template <int __inst>
void * __malloc_alloc<__inst>::_S_oom_malloc(size_t __n)
{
  __oom_handler_type __my_malloc_handler;
  void * __result;

  for (;;) {
    __my_malloc_handler = __oom_handler;
    if (0 == __my_malloc_handler) { __THROW_BAD_ALLOC; }
    (*__my_malloc_handler)();
    __result = malloc(__n);
    if (__result) return(__result);
  }
#if defined(__STL_NEED_UNREACHABLE_RETURN)
  return 0;
#endif

}

template <int __inst>
void* __malloc_alloc<__inst>::_S_oom_realloc(void* __p, size_t __n)
{
  __oom_handler_type __my_malloc_handler;
  void* __result;

  for (;;) {
    __my_malloc_handler = __oom_handler;
    if (0 == __my_malloc_handler) { __THROW_BAD_ALLOC; }
    (*__my_malloc_handler)();
    __result = realloc(__p, __n);
    if (__result) return(__result);
  }
#if defined(__STL_NEED_UNREACHABLE_RETURN)
  return 0;
#endif
}


# ifdef __STL_DEBUG_ALLOC

template <class _Alloc>
void * __debug_alloc<_Alloc>::allocate(size_t __n) {
  size_t __real_n = __n + __extra_before_chunk() + __extra_after_chunk();
  __alloc_header *__result = (__alloc_header *)__allocator_type::allocate(__real_n);
  memset((char*)__result, __shred_byte, __real_n*sizeof(value_type));
  __result->__magic = __magic;
  __result->__type_size = sizeof(value_type);
  __result->_M_size = __n;
  return ((char*)__result) + (long)__extra_before;
}

template <class _Alloc>
void __debug_alloc<_Alloc>::deallocate(void *__p, size_t __n) {
  __alloc_header * __real_p = (__alloc_header*)((char *)__p -(long)__extra_before);
  // check integrity
  __stl_verbose_assert(__real_p->__magic != __deleted_magic, _StlMsg_DBA_DELETED_TWICE);
  __stl_verbose_assert(__real_p->__magic == __magic, _StlMsg_DBA_NEVER_ALLOCATED);
  __stl_verbose_assert(__real_p->__type_size == 1,
		       _StlMsg_DBA_TYPE_MISMATCH);
  __stl_verbose_assert(__real_p->_M_size == __n, _StlMsg_DBA_SIZE_MISMATCH);
  // check pads on both sides
  unsigned char* __tmp;
  for (__tmp= (unsigned char*)(__real_p+1); __tmp < (unsigned char*)__p; __tmp++)
    __stl_verbose_assert(*__tmp==__shred_byte, _StlMsg_DBA_UNDERRUN);
  size_t __real_n= __n + __extra_before_chunk() + __extra_after_chunk();

  for (__tmp= ((unsigned char*)__p)+__n*sizeof(value_type); 
       __tmp < ((unsigned char*)__real_p)+__real_n ; __tmp++)
    __stl_verbose_assert(*__tmp==__shred_byte, _StlMsg_DBA_OVERRUN);
  // that may be unfortunate, just in case
  __real_p->__magic=__deleted_magic;
  memset((char*)__p, __shred_byte, __n*sizeof(value_type));
  __allocator_type::deallocate(__real_p, __real_n);
}

# if 0
template <class _Alloc>
void * 
__debug_alloc<_Alloc>::reallocate(void *__p, size_t __old_sz, size_t __new_sz) {
  __alloc_header * __real_p = (__alloc_header*)((char *)__p - (long)__extra_before);
  size_t __extra = __extra_before_chunk() + __extra_after_chunk();
  __stl_verbose_assert(__real_p->__magic != __deleted_magic, _StlMsg_DBA_DELETED_TWICE);
  __stl_verbose_assert(__real_p->__magic == __magic, _StlMsg_DBA_NEVER_ALLOCATED);
  __stl_verbose_assert(__real_p->__type_size == sizeof(value_type), 
		       _StlMsg_DBA_TYPE_MISMATCH);
  __stl_verbose_assert(__real_p->_M_size == __old_sz, _StlMsg_DBA_SIZE_MISMATCH);
  __real_p = (__alloc_header*)__allocator_type::reallocate(__real_p, __old_sz + __extra, 
							 __new_sz + __extra);
  __real_p->_M_size = __new_sz;
  return ((char*)__real_p) + (long)__extra_before;
}
# endif
#endif


/* We allocate memory in large chunks in order to avoid fragmenting     */
/* the malloc heap too much.                                            */
/* We assume that size is properly aligned.                             */
/* We hold the allocation lock.                                         */
template <bool __threads, int __inst>
char*
__node_alloc<__threads, __inst>::_S_chunk_alloc(size_t _p_size, 
						int& __nobjs)
{
  char* __result;
  size_t __total_bytes = _p_size * __nobjs;
  size_t __bytes_left = _S_end_free - _S_start_free;

  if (__bytes_left >= __total_bytes) {
    __result = _S_start_free;
    _S_start_free += __total_bytes;
    return(__result);
  } else if (__bytes_left >= _p_size) {
    __nobjs = (int)(__bytes_left/_p_size);
    __total_bytes = _p_size * __nobjs;
    __result = _S_start_free;
    _S_start_free += __total_bytes;
    return(__result);
  } else {
    size_t __bytes_to_get = 
      2 * __total_bytes + _S_round_up(_S_heap_size >> 4);
    // Try to make use of the left-over piece.
    if (__bytes_left > 0) {
      _Obj* __STL_VOLATILE* __my_free_list =
	_S_free_list + _S_freelist_index(__bytes_left);

      ((_Obj*)_S_start_free) -> _M_free_list_link = *__my_free_list;
      *__my_free_list = (_Obj*)_S_start_free;
    }
    _S_start_free = (char*)__STL_CHUNK_MALLOC(__bytes_to_get);
    if (0 == _S_start_free) {
      size_t __i;
      _Obj* __STL_VOLATILE* __my_free_list;
      _Obj* __p;
      // Try to make do with what we have.  That can't
      // hurt.  We do not try smaller requests, since that tends
      // to result in disaster on multi-process machines.
      for (__i = _p_size; __i <= (size_t)_MAX_BYTES; __i += (size_t)_ALIGN) {
	__my_free_list = _S_free_list + _S_freelist_index(__i);
	__p = *__my_free_list;
	if (0 != __p) {
	  *__my_free_list = __p -> _M_free_list_link;
	  _S_start_free = (char*)__p;
	  _S_end_free = _S_start_free + __i;
	  return(_S_chunk_alloc(_p_size, __nobjs));
	  // Any leftover piece will eventually make it to the
	  // right free list.
	}
      }
      _S_end_free = 0;	// In case of exception.
      _S_start_free = (char*)__STL_CHUNK_MALLOC(__bytes_to_get);
    /*
      (char*)malloc_alloc::allocate(__bytes_to_get);
      */

      // This should either throw an
      // exception or remedy the situation.  Thus we assume it
      // succeeded.
    }
    _S_heap_size += __bytes_to_get;
    _S_end_free = _S_start_free + __bytes_to_get;
    return(_S_chunk_alloc(_p_size, __nobjs));
  }
}


/* Returns an object of size __n, and optionally adds to size __n free list.*/
/* We assume that __n is properly aligned.                                */
/* We hold the allocation lock.                                         */
template <bool __threads, int __inst>
void*
__node_alloc<__threads, __inst>::_S_refill(size_t __n)
{
  int __nobjs = 20;
  __n = _S_round_up(__n);
  char* __chunk = _S_chunk_alloc(__n, __nobjs);
  _Obj* __STL_VOLATILE* __my_free_list;
  _Obj* __result;
  _Obj* __current_obj;
  _Obj* __next_obj;
  int __i;

  if (1 == __nobjs) return(__chunk);
  __my_free_list = _S_free_list + _S_freelist_index(__n);

  /* Build free list in chunk */
  __result = (_Obj*)__chunk;
  *__my_free_list = __next_obj = (_Obj*)(__chunk + __n);
  for (__i = 1; ; __i++) {
    __current_obj = __next_obj;
    __next_obj = (_Obj*)((char*)__next_obj + __n);
    if (__nobjs - 1 == __i) {
      __current_obj -> _M_free_list_link = 0;
      break;
    } else {
      __current_obj -> _M_free_list_link = __next_obj;
    }
  }
  return(__result);
}

# if 0
template <bool threads, int inst>
void*
__node_alloc<threads, inst>::reallocate(void* __p,
                                                    size_t __old_sz,
                                                    size_t __new_sz)
{
  void* __result;
  size_t __copy_sz;

  if (__old_sz > (size_t) _MAX_BYTES && __new_sz > (size_t) _MAX_BYTES) {
    return(realloc(__p, __new_sz));
  }
  if (_S_round_up(__old_sz) == _S_round_up(__new_sz)) return(__p);
  __result = allocate(__new_sz);
  __copy_sz = __new_sz > __old_sz? __old_sz : __new_sz;
  memcpy(__result, __p, __copy_sz);
  deallocate(__p, __old_sz);
  return(__result);
}
# endif


# if ( __STL_STATIC_TEMPLATE_DATA > 0 )

// malloc_alloc out-of-memory handling
template <int __inst>
__oom_handler_type __malloc_alloc<__inst>::__oom_handler=(__oom_handler_type)0 ;

#ifdef __STL_THREADS
    template <bool __threads, int __inst>
    _STL_mutex_base
    _Node_Alloc_Lock<__threads, __inst>::_S_lock __STL_MUTEX_INITIALIZER;
#endif

template <bool __threads, int __inst>
_Node_alloc_obj * __STL_VOLATILE
__node_alloc<__threads, __inst>::_S_free_list[_NFREELISTS]
= {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// The 16 zeros are necessary to make version 4.1 of the SunPro
// compiler happy.  Otherwise it appears to allocate too little
// space for the array.

template <bool __threads, int __inst>
char *__node_alloc<__threads, __inst>::_S_start_free = 0;

template <bool __threads, int __inst>
char *__node_alloc<__threads, __inst>::_S_end_free = 0;

template <bool __threads, int __inst>
size_t __node_alloc<__threads, __inst>::_S_heap_size = 0;


# else /* ( __STL_STATIC_TEMPLATE_DATA > 0 ) */

__DECLARE_INSTANCE(__oom_handler_type, __malloc_alloc<0>::__oom_handler, =0);

__DECLARE_INSTANCE(char *, __single_client_alloc::_S_start_free,=0);
__DECLARE_INSTANCE(char *, __single_client_alloc::_S_end_free,=0);
__DECLARE_INSTANCE(size_t, __single_client_alloc::_S_heap_size,=0);
__DECLARE_INSTANCE(_Node_alloc_obj * __STL_VOLATILE,
                   __single_client_alloc::_S_free_list[_NFREELISTS],
                   ={0});
__DECLARE_INSTANCE(char *, __multithreaded_alloc::_S_start_free,=0);
__DECLARE_INSTANCE(char *, __multithreaded_alloc::_S_end_free,=0);
__DECLARE_INSTANCE(size_t, __multithreaded_alloc::_S_heap_size,=0);
__DECLARE_INSTANCE(_Node_alloc_obj * __STL_VOLATILE,
                   __multithreaded_alloc::_S_free_list[_NFREELISTS],
                   ={0});
#   ifdef __STL_THREADS
__DECLARE_INSTANCE(_STL_mutex_base,
                   __single_client_alloc::_S_node_allocator_lock,
                   __STL_MUTEX_INITIALIZER);
__DECLARE_INSTANCE(_STL_mutex_base,
                   __multithreaded_alloc::_S_node_allocator_lock,
                   __STL_MUTEX_INITIALIZER);
#   endif
#  endif /* __STL_STATIC_TEMPLATE_DATA */

__STL_END_NAMESPACE

#endif /*  __STL_ALLOC_C */

// Local Variables:
// mode:C++
// End:
