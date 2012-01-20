// Copyright (c) 2007 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_RS__PAGEALLOC_H
#define CGAL_RS__PAGEALLOC_H

#include <cstdlib>
#include <cstring>

namespace CGAL{
namespace RS_MGCD{

#define CGALRS_PAGESIZE        4194304
#define CGALRS_TABLESIZE       2048
#define CGALRS_PAGES           8
#define CGALRS_VOIDSCAST       unsigned long

struct pinfo{
    void *start;
    size_t size;
};

CGALRS_THREAD_ATTR void**  pages_startptr;
CGALRS_THREAD_ATTR size_t  pages_max;//=CGALRS_PAGES;
CGALRS_THREAD_ATTR size_t  pages_allocated;
CGALRS_THREAD_ATTR size_t  pages_current;

CGALRS_THREAD_ATTR size_t  page_remainingbytes;
CGALRS_THREAD_ATTR void*   page_currentptr;

CGALRS_THREAD_ATTR struct pinfo    *nodes_allocated;
CGALRS_THREAD_ATTR size_t          nodes_total;
CGALRS_THREAD_ATTR size_t          nodes_assigned;

class Page_alloc{

    protected:
        static
        void* meminit(){
            pages_startptr=(void**)malloc(CGALRS_PAGES*sizeof(void*));
            pages_startptr[0]=malloc(CGALRS_PAGESIZE);
            pages_allocated=1;
            pages_current=0;
            page_remainingbytes=CGALRS_PAGESIZE;
            page_currentptr=pages_startptr[0];
            nodes_total=CGALRS_TABLESIZE;
            nodes_allocated=
                (struct pinfo*)malloc(CGALRS_TABLESIZE*sizeof(struct pinfo));
            nodes_assigned=0;
            return page_currentptr;
        };

        static
        void* newpage(){
            void *r;
            if(pages_allocated>pages_current+1){
                ++pages_current;
                r=pages_startptr[pages_current];
                page_currentptr=r;
                page_remainingbytes=CGALRS_PAGESIZE;
                return r;
            }
            // iso c++ forbids to initialize a static member (pages_max),
            // so we have to start using pages_max when the amount of
            // allocated pages reaches the value CGALRS_PAGES (this is not of
            // course the cleanest way to do it)
            if(pages_allocated==CGALRS_PAGES)
                pages_max=2*CGALRS_PAGES;
            else
                pages_max=0;
            if(pages_allocated==pages_max){
                pages_max*=2;
                pages_startptr=
                    (void**)realloc(pages_startptr,pages_max*sizeof(void*));
            }
            r=malloc(CGALRS_PAGESIZE);
            pages_startptr[pages_allocated]=r;
            page_currentptr=r;
            ++pages_allocated;
            page_remainingbytes=CGALRS_PAGESIZE;
            return r;
        };

        // if they ask us to allocate a memory size bigger than the page,
        // we are lost (we could in that case make a bigger page)
        static
        void* palloc(size_t size){
            void* r;
            if(size>page_remainingbytes){
                newpage();
                return palloc(size);
            }
            if(nodes_assigned==nodes_total){
                nodes_total*=2;
                nodes_allocated=(struct pinfo*)realloc
                    (nodes_allocated,nodes_total*sizeof(struct pinfo));
            }
            page_remainingbytes-=size;
            r=page_currentptr;
            page_currentptr=(void*)((CGALRS_VOIDSCAST)page_currentptr+size);
            // c++ does not support nodes_allocated[nodes_assigned]={r,s}
            nodes_allocated[nodes_assigned].start=r;
            nodes_allocated[nodes_assigned].size=size;
            ++nodes_assigned;
            return r;
        };

        static
        void* prealloc(void *ptr,size_t size){
            void *dest;
            size_t i=0;
            while(nodes_allocated[i].start!=ptr)
                ++i;
            if(nodes_allocated[i].size<size){
                dest=palloc(size);
                nodes_allocated[i].start=dest;
                nodes_allocated[i].size=size;
                return memcpy(dest,ptr,nodes_allocated[i].size);
            }
            return ptr;
        };

        #define CGALRS_PFREE(X)        {}
        //void CGALRS_PFREE(void* ptr){
        //  size_t i=0;
        //  while(nodes_allocated[i].start!=ptr)
        //      ++i;
        //  nodes_allocated[i].start=0;
        //  return;
        //};

        static
        void* memclear(){
            pages_current=0;
            page_remainingbytes=CGALRS_PAGESIZE;
            page_currentptr=pages_startptr[0];
            nodes_assigned=0;
            return page_currentptr;
        }

        static
        void memrelease(){
            size_t i;
            for(i=0;i<pages_allocated;++i)
                free(pages_startptr[i]);
            return;
        };

}; // class Page_alloc

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__PAGEALLOC_H
