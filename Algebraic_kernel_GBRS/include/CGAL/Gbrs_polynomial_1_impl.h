#ifndef CGAL_GBRS_POLYNOMIAL_1_IMPL_H
#define CGAL_GBRS_POLYNOMIAL_1_IMPL_H

CGAL_BEGIN_NAMESPACE

template <class T>
Rational_polynomial_1 Rational_polynomial_1::operator* (const T &n) const {
	Rational_polynomial_1 r (*this);
	r.scale (n);
	return r;
};

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_POLYNOMIAL_1_IMPL_H
