// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_MEMORY_H
#define CGAL_RS_MEMORY_H

#include <CGAL/RS/basic.h>
#include <gmp.h>
#include <cstdlib>
#include <rs_exports.h>
#ifdef CGAL_USE_OLD_RS3
extern "C"{
#include <rs_gc.h>
}
#endif

namespace CGAL{

extern "C"{
extern void* RS_gmpalloc(size_t);
extern void* RS_gmprealloc(void*,size_t,size_t);
extern void RS_gmpfree(void*,size_t);
}

#ifdef CGAL_USE_OLD_RS3
inline void * rs3_rs_alloc(size_t s) {
  return(rs_alloc(s,&(rs_default_gc_session.default_heap[0])));
}

inline void * rs3_rs_realloc(void * p,size_t s) {
  assert(1==0);
  return(NULL);
}

inline void rs3_rs_free(void * p) {
  assert(1==0);
}

inline void * rs3_rs_gmp_realloc(void * old_p,size_t old_size,size_t new_size){
  return(rs_realloc(old_p,((RS_UI)old_size),((RS_UI)new_size),&(rs_default_gc_session.default_heap[0])));
}

inline void rs3_rs_gmp_free(void * p,size_t s) {}
#endif  // CGAL_USE_OLD_RS3

//--------------------------------------------------
// extern void * (*__cgalrs_allocate_func) (size_t);
// extern void * (*__cgalrs_reallocate_func) (void *, size_t, size_t);
// extern void   (*__cgalrs_free_func) (void *, size_t);
//--------------------------------------------------
inline void* __cgalrs_default_allocate(size_t s){
        return malloc(s);
}

inline void* __cgalrs_default_reallocate(void *a,size_t o,size_t n){
        return realloc(a,n);
}

inline void __cgalrs_default_free(void *a,size_t s){
        return free(a);
}

CGALRS_THREAD_ATTR void * (*__cgalrs_allocate_func) (size_t) =
                __cgalrs_default_allocate;

CGALRS_THREAD_ATTR void * (*__cgalrs_reallocate_func) (void *, size_t, size_t) =
                __cgalrs_default_reallocate;

CGALRS_THREAD_ATTR void   (*__cgalrs_free_func) (void *, size_t) =
                __cgalrs_default_free;

inline void __cgalrs_dummy_free(void *p,size_t s){}

inline void cgalrs_set_memory_functions(
                        void *(*alloc_func) (size_t),
                        void *(*realloc_func) (void *, size_t, size_t),
                        void (*free_func) (void *, size_t)){
        __cgalrs_allocate_func=
                alloc_func?alloc_func:__cgalrs_default_allocate;
        __cgalrs_reallocate_func=
                realloc_func?realloc_func:__cgalrs_default_reallocate;
        __cgalrs_free_func=
                free_func?free_func:__cgalrs_default_free;
}

inline void cgalrs_get_memory_functions(
                        void *(**alloc_func) (size_t),
                        void *(**realloc_func) (void *, size_t, size_t),
                        void (**free_func) (void *, size_t)){
        if(alloc_func!=NULL)
                *alloc_func=__cgalrs_allocate_func;
        if(realloc_func!=NULL)
                *realloc_func=__cgalrs_reallocate_func;
        if(free_func!=NULL)
                *free_func=__cgalrs_free_func;
}

} // namespace CGAL

#endif  // CGAL_RS_MEMORY_H
