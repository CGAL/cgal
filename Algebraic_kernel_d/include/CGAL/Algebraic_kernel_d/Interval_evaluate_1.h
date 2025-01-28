// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_INTERVAL_EVALUATE_1
#define CGAL_INTERVAL_EVALUATE_1 1

#include <iterator>

#include <CGAL/basic.h>
#include <boost/numeric/interval.hpp>
#include <CGAL/algorithm.h>
#include <CGAL/array.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Polynomial_traits_d.h>

namespace CGAL {

namespace internal {

template<typename Polynomial_1, typename Bound>
struct Interval_evaluate_1 : public CGAL::cpp98::binary_function
<Polynomial_1,std::pair<Bound,Bound>,
  std::pair<typename CGAL::Coercion_traits<typename
     CGAL::Polynomial_traits_d<Polynomial_1>::Coefficient_type,Bound>::Type,
  typename CGAL::Coercion_traits<typename
     CGAL::Polynomial_traits_d<Polynomial_1>::Coefficient_type,Bound>::Type> > {

  typedef CGAL::Polynomial_traits_d< Polynomial_1 > PT_1;

  typedef typename PT_1::Innermost_coefficient_type Innermost_coefficient_type;

  typedef CGAL::Coercion_traits< Innermost_coefficient_type, Bound > CT;

  typedef typename CT::Type Coercion_type;

  typedef std::pair< Coercion_type, Coercion_type > result_type;

  result_type operator()(const Polynomial_1& p,
                         const std::pair< Bound, Bound >& b) const {
    return this->operator()(p, CGAL::make_array(b.first, b.second));
  }

  result_type operator()(const Polynomial_1& p,
                         const std::array< Bound, 2 >& b) const {

    typename CT::Cast cast;

    typedef ::boost::numeric::interval< Coercion_type > Coercion_interval;

    typedef typename PT_1::Coefficient_const_iterator
      Coefficient_const_iterator;

    Coercion_interval ix(cast(b[0]), cast(b[1]));

    typedef typename PT_1::Coefficient_const_iterator_range
      Coefficient_const_iterator_range;

    Coefficient_const_iterator_range range =
      typename PT_1::Construct_coefficient_const_iterator_range()(p);

    Coefficient_const_iterator it = std::prev(range.second);

    Coercion_interval res(cast(*it));

    Coefficient_const_iterator p_begin = range.first;
    while(it != p_begin) {
      it--;
      res = res * ix + Coercion_interval(cast(*it));
    }
    return std::make_pair(res.lower(),res.upper());
  }

};

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERVAL_EVALUATE_1
