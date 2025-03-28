// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallée (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

#ifndef CGAL_COMPLEX_NUMBER_H
#define CGAL_COMPLEX_NUMBER_H

#include <fstream>
#include <sstream>
#include <utility>

namespace CGAL {

/*
Templated by a field FT. Represents a complex number over FT.
*/
template <class FT>
class Complex_number
{
  typedef Complex_number<FT> Self;

  FT real_, imag_;

public:
  Complex_number(const FT& real_part)
    : real_(real_part),
      imag_(0)
  {}

  Complex_number(const FT& real_part, const FT& imaginary_part)
    : real_(real_part),
      imag_(imaginary_part)
  {}

  Complex_number()
    : Complex_number(FT(0), FT(0))
  {}

  template<class U,class V>
  Complex_number(U&& real_part, V&& imaginary_part)
    : real_(std::forward<U>(real_part)),
      imag_(std::forward<V>(imaginary_part))
  {}

  void real(const FT& real_part) {
    real_ = real_part;
  }

  void imag(const FT& imaginary_part) {
    imag_ = imaginary_part;
  }

  FT real() const {
    return real_;
  }

  FT imag() const {
    return imag_;
  }

  Self& operator+=(const Self& other);
  Self& operator-=(const Self& other);
  Self& operator*=(const Self& other);
  Self& operator/=(const Self& other);

  // These member versions are not working ?
  /* Self operator+(const Self& z) const; */
  /* Self operator-(const Self& z) const; */

  // Hidden friends
  friend Self operator+(const Self& z) {
    return z;
  }

  friend Self operator-(const Self& z) {
    return Self(-z.real_,-z.imag_);
  }

  friend bool operator==(const Self& z1, const Self& z2) {
    return (z1.real_==z2.real_ && z1.imag_==z2.imag_);
  }

  friend bool operator!=(const Self& z1, const Self& z2) {
    return !operator==(z1, z2);
  }

  friend Self operator+(const Self& z1, const Self& z2) {
    return Self(z1.real_+z2.real_, z1.imag_+z2.imag_);
  }

  friend Self operator-(const Self& z1, const Self& z2) {
    return Self(z1.real_-z2.real_, z1.imag_-z2.imag_);
  }

  friend Self operator*(const Self& z1, const Self& z2) {
    return Self(z1.real_*z2.real_-z1.imag_*z2.imag_, z1.real_*z2.imag_+z1.imag_*z2.real_);
  }

  friend Self operator/(const Self& z1, const Self& z2) {
    FT m2 = norm(z2);
    return Self(z1.real_/m2, z1.imag_/m2)*conj(z2);
  }

  friend std::ostream& operator<<(std::ostream& s, const Self& z) {
    s << z.real_ << std::endl << z.imag_ << std::endl;
    return s;
  }

  friend void operator>>(std::istream& s, Self& z) {
    FT ft;
    s >> ft;
    z.real(ft);
    s >> ft;
    z.imag(ft);
  }
};

////////////////////////////////////////////////////////////////////////////////
template<class FT>
Complex_number<FT>& Complex_number<FT>::operator+=(const Complex_number<FT>& other)
{
  real_ += other.real();
  imag_ += other.imag();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator-=(const Complex_number<FT>& other)
{
  real_ -= other.real();
  imag_ -= other.imag();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator*=(const Complex_number<FT>& other)
{
  real_ = real_*other.real() - imag_*other.imag();
  imag_ = real_*other.imag() + imag_*other.real();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator/=(const Complex_number<FT>& other)
{
  FT m2 = norm(other);
  real_ /= m2;
  imag_ /= m2;
  this *= conj(other);
  return *this;
}

template<class FT>
FT norm(const Complex_number<FT>& z)
{
  return z.real()*z.real() + z.imag()*z.imag();
}

template<class FT>
Complex_number<FT> conj(const Complex_number<FT>& z)
{
  return Complex_number<FT>(z.real(), -z.imag());
}

} // namespace CGAL

#endif // CGAL_COMPLEX_NUMBER_H
