// Copyright (c) 2006-2010 Max-Planck-Institute Saarbruecken (Germany).
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


#ifndef CGAL_INTERVAL_EVALUATE_2
#define CGAL_INTERVAL_EVALUATE_2 1

#include <iterator>

#include <CGAL/basic.h>
#include <boost/numeric/interval.hpp>
#include <CGAL/algorithm.h>
#include <CGAL/array.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_1.h>

namespace CGAL {

namespace internal {

template<typename Polynomial_2, typename Bound>
struct Interval_evaluate_2 : public CGAL::cpp98::binary_function
<Polynomial_2,std::array<Bound,4>,
      std::pair<typename CGAL::Coercion_traits<typename CGAL::Polynomial_traits_d<Polynomial_2>::Innermost_coefficient_type,Bound>::Type,
                typename CGAL::Coercion_traits<typename CGAL::Polynomial_traits_d<Polynomial_2>::Innermost_coefficient_type,Bound>::Type> > {

public:

  typedef CGAL::Polynomial_traits_d< Polynomial_2 > PT_2;

  typedef typename PT_2::Innermost_coefficient_type Innermost_coefficient_type;

  typedef CGAL::Coercion_traits< Innermost_coefficient_type, Bound > CT;

  typedef typename CT::Type Coercion_type;

  typedef std::pair< Coercion_type, Coercion_type > result_type;

  result_type operator()(const Polynomial_2& p,
                         const std::array< Bound, 4 >& b) const {

    typename CT::Cast cast;

    typedef ::boost::numeric::interval< Coercion_type > Coercion_interval;

    typedef typename PT_2::Coefficient_const_iterator
      Coefficient_const_iterator;

    typedef typename PT_2::Coefficient_const_iterator_range
      Coefficient_const_iterator_range;

    typedef typename PT_2::Coefficient_type Polynomial_1;

    CGAL::internal::Interval_evaluate_1< Polynomial_1,Bound >
      interval_evaluate_1;

    typedef typename CGAL::internal::Interval_evaluate_1< Polynomial_1,Bound >::
      result_type Interval_result_type;

    std::pair< Bound, Bound > x_pair = std::make_pair(b[0],b[1]);

    Coercion_interval iy(cast(b[2]),cast(b[3]));

    // CGAL::Polynomial does not provide Coercion_traits for number
    // types => therefore evaluate manually
    Coefficient_const_iterator_range range =
      typename PT_2::Construct_coefficient_const_iterator_range()(p);

    Coefficient_const_iterator it = std::prev(range.second);

    Interval_result_type initial_pair = interval_evaluate_1(*it,x_pair);
    Coercion_interval res(initial_pair.first,initial_pair.second);

    Coefficient_const_iterator p_begin = range.first;

    while((it) != p_begin) {
      it--;
      Interval_result_type curr_iv = interval_evaluate_1(*it,x_pair);
      res = res * iy + Coercion_interval(curr_iv.first,curr_iv.second);
    }
    return std::make_pair(res.lower(),res.upper());
  }

};

} // namespace internal


} // namespace CGAL

#endif // CGAL_INTERVAL_EVALUATE_2
