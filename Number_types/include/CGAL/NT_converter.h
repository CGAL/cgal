// Copyright (c) 2001-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_NT_CONVERTER_H
#define CGAL_NT_CONVERTER_H

#include <functional>
#include <CGAL/number_type_config.h>
#include <CGAL/number_utils.h>

template <bool> class Interval_nt;

namespace CGAL {

// A number type converter usable as default, using the conversion operator.

template < class NT1, class NT2 >
struct NT_converter
  : public CGAL::cpp98::unary_function< NT1, NT2 >
{
    NT2
    operator()(const NT1 &a) const
    {
        return NT2(a);
    }
};

// Some specializations for :
// - double to call to_double().
// - Interval_nt<> to call to_interval().
// - NT1 == NT2 to return a reference instead of copying.

template < class NT1 >
struct NT_converter < NT1, NT1 >
  : public CGAL::cpp98::unary_function< NT1, NT1 >
{
    const NT1 &
    operator()(const NT1 &a) const
    {
        return a;
    }
};

template < class NT1 >
struct NT_converter < NT1, double >
  : public CGAL::cpp98::unary_function< NT1, double >
{
    double
    operator()(const NT1 &a) const
    {
        return CGAL_NTS to_double(a);
    }
};

template <>
struct NT_converter < double, double >
  : public CGAL::cpp98::unary_function< double, double >
{
    const double &
    operator()(const double &a) const
    {
        return a;
    }
};

template < class NT1, bool b >
struct NT_converter < NT1, Interval_nt<b> >
  : public CGAL::cpp98::unary_function< NT1, Interval_nt<b> >
{
    Interval_nt<b>
    operator()(const NT1 &a) const
    {
        return CGAL_NTS to_interval(a);
    }
};

template < bool b >
struct NT_converter < Interval_nt<b>, Interval_nt<b> >
  : public CGAL::cpp98::unary_function< Interval_nt<b>, Interval_nt<b> >
{
    const Interval_nt<b> &
    operator()(const Interval_nt<b> &a) const
    {
        return a;
    }
};

} //namespace CGAL

#endif // CGAL_NT_CONVERTER_H
