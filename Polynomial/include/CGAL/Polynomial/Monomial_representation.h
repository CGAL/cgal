// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_POLYNOMIAL_MONOMIAL_REPRESENTAION_H
#define CGAL_POLYNOMIAL_MONOMIAL_REPRESENTAION_H

#include <CGAL/Exponent_vector.h>
#include <CGAL/Polynomial/misc.h>


namespace CGAL {

namespace internal{
template <typename Polynomial> struct Monomial_representation;

// Polynomial muss be at least univariate !
template <typename Coeff_ >
struct Monomial_representation<Polynomial<Coeff_> >{
private:
  typedef typename Innermost_coefficient_type<Polynomial<Coeff_> >::Type
  Innermost_coefficient;

  // Polynomial is univariate
  // final creation of pair<Exponent_vector,Innermost_coefficient>
  template <typename Polynomial, typename OutputIterator>
  OutputIterator
  create_mrep(const Polynomial& p, OutputIterator oit , Exponent_vector& ivec, Tag_true) const {
    int degree = 0;
    for(typename Polynomial::const_iterator it = p.begin(); it != p.end(); it++){
      ivec[0] = degree;
      if(!CGAL::is_zero(*it))
        *oit++ = std::make_pair(ivec,*it);
      degree++;
    }
    ivec[0]=0;
    return oit;
  }

  // polynomial is multivariate
  // define correct exponent for dimension and recurse
  template <typename Polynomial, typename OutputIterator>
  OutputIterator
  create_mrep(const Polynomial& p, OutputIterator oit , Exponent_vector& ivec, Tag_false) const {
    if(CGAL::is_zero(p)) return oit;
    static const int dim = Dimension<Polynomial>::value ;
    int degree = 0;
    for(typename Polynomial::const_iterator it = p.begin(); it != p.end(); it++){
      ivec[dim-1] = degree;
      oit = create_mrep(*it,oit,ivec,CGAL::Boolean_tag<1 == dim-1>());
      degree++;
    }
    ivec[dim-1] = 0;
    return oit;
  }

public:
  template <typename OutputIterator>
  OutputIterator operator()(const Polynomial<Coeff_>& p, OutputIterator oit) const {
    typedef Polynomial<Coeff_> Polynomial;
    typedef CGAL::Boolean_tag<1 == Dimension<Polynomial>::value> Is_univariate;
    CGAL::Exponent_vector ivec((std::vector<int>)(Dimension<Polynomial>::value));
    if(p.is_zero()){
      *oit++ = std::make_pair(ivec,Innermost_coefficient(0));
      return oit;
    }
    return create_mrep(p, oit, ivec, Is_univariate());
  }
};

} // namespace internal
} //namespace CGAL

#endif //CGAL_POLYNOMIAL_MONOMIAL_REPRESENTAION_H
