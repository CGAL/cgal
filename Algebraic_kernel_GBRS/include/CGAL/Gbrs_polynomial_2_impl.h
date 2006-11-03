#ifndef CGAL_GBRS_POLYNOMIAL_2_IMPL_H
#define CGAL_GBRS_POLYNOMIAL_2_IMPL_H

CGAL_BEGIN_NAMESPACE

template <class T>
Rational_polynomial_2 Rational_polynomial_2::operator*(const T &n)const{
	Rational_polynomial_2 r(*this);
	return (r*=n);
};

template <class T> inline Rational_polynomial_2 operator*(const T &n,
		const Rational_polynomial_2 &p){
	return (p*n);
};

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_POLYNOMIAL_2_IMPL_H
