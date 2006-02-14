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
// file          : include/CGAL/Do_nothing.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================

#ifndef CGAL_DO_NOTHING_H
#define CGAL_DO_NOTHING_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

struct Do_nothing {

template<class O1,class O2,class A1>
void operator()(const O1& func1, const O2& func2, const A1& arg1)
{ }

template<class O1,class O2,class A1,class A2>
void operator()(const O1& func1, const O2& func2, const A1& arg1, const A2& arg2) const
{ }

template<class O1,class O2,class A1,class A2,class A3>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6,class A7>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6,
		const A7& arg7) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6,
		const A7& arg7, const A8& arg8) const
{ }

template<class O1,class O2,class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8,class A9>
void operator()(const O1& func1, const O2& func2, 
                const A1& arg1, const A2& arg2, const A3& arg3, 
		const A4& arg4, const A5& arg5, const A6& arg6,
		const A7& arg7, const A8& arg8, const A9& arg9) const
{ }

};

CGAL_END_NAMESPACE

#endif 

