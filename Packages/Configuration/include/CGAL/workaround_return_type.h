// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-05 $
// release_date  : $CGAL_Date: 1997/12/17 $
//
// file          : include/CGAL/workaround_return_type.h
// source        :
// revision      :
// revision_date :
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Utrecht University
//
// ============================================================================


#ifndef CGAL_WORKAROUND_RETURN_TYPE_H
#define CGAL_WORKAROUND_RETURN_TYPE_H

/*
This is a workaround for compilers that do not allow to return a type derived
in a special way from a template parameter:
template <class A>
A::X foo(A a) {}

gcc 2.7.2 and a Sun compiler have this feature.
*/

#ifdef CGAL_CFG_RETURN_TYPE_BUG_1
template <class R>
struct R_FT_Return_type
{
    R_FT_Return_type(typename R::FT xt) {_xt = xt; }
    operator typename R::FT() const { return _xt; }
private:
   typename R::FT _xt;
};

#define R_FT_return(template_parameter) R_FT_Return_type<template_parameter>

template <class R>
struct R_RT_Return_type
{
    R_RT_Return_type(typename R::RT xt) {_xt = xt; }
    operator typename R::RT() const { return _xt; }
private:
    typename R::RT _xt;
};

#define R_RT_return(template_parameter) R_RT_Return_type<template_parameter>

#else
#define R_FT_return(template_parameter) typename template_parameter::FT

#define R_RT_return(template_parameter) typename template_parameter::RT
#endif

#endif // CGAL_WORKAROUND_RETURN_TYPE_H

