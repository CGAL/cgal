// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/kernel_event_support.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================


#if !defined(CGAL_KERNEL_EVENTS)
#define CGAL_KERNEL_EVENTS

#include <CGAL/basic.h>
#include <CGAL/event.h>
#include <map>
#include <cstring>


CGAL_BEGIN_NAMESPACE

template<class T>
struct kernel_event {

static CGAL::event   EVENT;

template<class O1,
         class A1>
void operator()(const O1& func1, const A1& arg1)
{
#if defined(DEBUG)
   std::cout << "occur (1)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&>(EVENT, func1, arg1);
}

template<class O1,
         class A1,class A2>
void operator()(const O1& func1,const A1& arg1,const A2& arg2) const
{
#if defined(DEBUG)
   std::cout << "occur (2)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&>(EVENT, func1, arg1, arg2);
}

template<class O1,
         class A1,class A2,class A3>
void operator()(const O1& func1, 
                const A1& arg1,const A2& arg2,const A3& arg3) const
{
#if defined(DEBUG)
   std::cout << "occur (3)- " << typeid(func1).name() << "\n";
   std::cout << arg1 << " " << arg2 << " " << arg3 <<  "\n";
#endif
   CGAL::occur<const O1&, const A1&,const A2& ,const A3&>(kernel_event<T>::EVENT, func1, arg1, arg2, arg3);
}

template<class O1,
         class A1,class A2,class A3,class A4>
void operator()(const O1& func1,
                const A1& arg1, const A2& arg2,const A3& arg3, 
		const A4& arg4) const
{
#if defined(DEBUG)
   std::cout << "occur (4)- " << typeid(func1).name() << "\n";
#endif
   // add O1
   CGAL::occur<const O1&, const A1&,const A2&,const A3&,const A4&>(EVENT, func1, arg1, arg2, arg3, arg4);
}

template<class O1,
         class A1,class A2,class A3,class A4,class A5>
void operator()(const O1& func1, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5) const
{
#if defined(DEBUG)
   std::cout << "occur (5)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&,const A3&,const A4&,const A5&>(EVENT, func1, arg1, arg2, arg3, arg4, arg5);
}

template<class O1,
         class A1,class A2,class A3,class A4,class A5,class A6>
void operator()(const O1& func1, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6) const
{
#if defined(DEBUG)
   std::cout << "occur (6)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&,const A3&,const A4&,const A5&,const A6&>(EVENT, func1, arg1, arg2, arg3, arg4, arg5, arg6);
}

template<class O1,
         class A1,class A2,class A3,class A4,class A5,class A6,class A7>
void operator()(const O1& func1,
                const A1& arg1,const A2& arg2,const A3& arg3, 
		const A4& arg4,const A5& arg5,const A6& arg6,
		const A7& arg7) const
{
#if defined(DEBUG)
   std::cout << "occur (7)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&,const A3&,const A4&,const A5&,const A6&,const A7&>(EVENT, func1, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
}

template<class O1,
         class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8>
void operator()(const O1& func1, 
                const A1& arg1,const A2& arg2,const A3& arg3, 
		const A4& arg4,const A5& arg5,const A6& arg6,
		const A7& arg7,const A8& arg8) const
{
#if defined(DEBUG)
   std::cout << "occur (8)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&,const A3&,const A4&,const A5&,const A6&,const A7&,const A8&>(EVENT, func1, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
}

template<class O1,
         class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8,class A9>
void operator()(const O1& func1, 
                const A1& arg1,const A2& arg2,const A3& arg3, 
		const A4& arg4,const A5& arg5,const A6& arg6,
		const A7& arg7,const A8& arg8,const A9& arg9) const
{
#if defined(DEBUG)
   std::cout << "occur (9)- " << typeid(func1).name() << "\n";
#endif
   CGAL::occur<const O1&,const A1&,const A2&,const A3&,const A4&,const A5&,const A6&,const A7&,const A8&,const A9&>(EVENT, func1, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
}

};

template<class T>
CGAL::event kernel_event<T>::EVENT;

CGAL_END_NAMESPACE

#endif
