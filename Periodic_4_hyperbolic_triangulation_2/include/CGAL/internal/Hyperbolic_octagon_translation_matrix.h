// Copyright (c) 2010-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov
//

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/internal/Exact_complex.h>
#include <CGAL/number_utils.h>

#include <iostream>
#include <vector>

namespace CGAL {

template <class ECplx>
class Hyperbolic_octagon_translation_matrix
{
public:
  typedef ECplx                              Matrix_element;
  typedef typename ECplx::NT                 NT;

private:
  typedef Hyperbolic_octagon_translation_matrix<ECplx>   Self;

protected:
  Matrix_element  _alpha;
  Matrix_element  _beta;

public:
  static void generators(std::vector<Self>& gens)
  {
    // Create matrix coefficients
    NT sq2 = CGAL::sqrt(NT(2));
    NT xi  = NT(1) + sq2;
    NT rxi = CGAL::sqrt(xi);

    // All matrices have the same _alpha
    ECplx alpha(xi,NT(0));

    // This vector holds the different _betas
    std::vector< ECplx > beta;

    beta.push_back(ECplx(sq2 * rxi, 0));
    beta.push_back(ECplx(rxi, rxi));
    beta.push_back(ECplx(0, sq2 * rxi ));
    beta.push_back(ECplx(-rxi, rxi));
    beta.push_back(ECplx(-sq2 * rxi, 0));
    beta.push_back(ECplx(-rxi, -rxi));
    beta.push_back(ECplx(0, -sq2 * rxi));
    beta.push_back(ECplx(rxi, -rxi));

    for(int i=0; i<8; ++i)
      gens.push_back(Self(alpha, beta[i]));
  }

  Hyperbolic_octagon_translation_matrix() : _alpha(1, 0), _beta(0, 0) {}
  Hyperbolic_octagon_translation_matrix(const Matrix_element& _a, const Matrix_element& _b)
    : _alpha(_a), _beta(_b)
  {}

  Matrix_element alpha() const { return _alpha; }
  Matrix_element beta() const { return _beta; }

  Self operator*(const Self& rh) const
  {
    return Self(_alpha*rh.alpha() + _beta*rh.beta().conj(),
                _alpha*rh.beta() + _beta*rh.alpha().conj());
  }

  Self inverse() const
  {
    return Self(_alpha.conj(), -_beta);
  }

  NT trace() const {
    return NT(2)*_alpha.real();
  }

  NT length() const
  {
    NT tr = _alpha.real();
    if(tr < 0)
      tr = -tr;

    return NT(2)*acosh(tr);
  }

  NT det() const
  {
    return _alpha.square_modulus() - _beta.square_modulus();
  }

  Matrix_element operator()(Matrix_element z) const
  {
    return (_alpha*z + _beta)/(_beta.conj()*z + _alpha.conj());
  }
};

template<class ECplx >
bool operator==(const Hyperbolic_octagon_translation_matrix<ECplx>& lh,
                const Hyperbolic_octagon_translation_matrix<ECplx>& rh) {
  return (lh.alpha() == rh.alpha() && lh.beta() == rh.beta());
}

template<class ECplx>
std::ostream& operator<<(std::ostream& os, const Hyperbolic_octagon_translation_matrix<ECplx>& m) {
  os << m.alpha() << ", " << m.beta();
  return os;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
