#ifndef QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H
#define QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H

CGAL_BEGIN_NAMESPACE

template <typename Homogeneous>
typename Homogeneous::Point_3
quotient_coordinates_to_homogeneous_point(
				  typename Homogeneous::FT x,
				  typename Homogeneous::FT y,
				  typename Homogeneous::FT z) {
  typedef typename Homogeneous::Point_3 Point_3;
  if( (x.denominator() == y.denominator()) && 
      (x.denominator() == z.denominator())) {
    Point_3 p( x.numerator(),
	       y.numerator(),
	       z.numerator(),
	       x.denominator());
    return normalized(p);
  }
  else {
    Point_3 p( x.numerator()   * y.denominator() * z.denominator(),
	       x.denominator() * y.numerator()   * z.denominator(),
	       x.denominator() * y.denominator() * z.numerator(),
	       x.denominator() * y.denominator() * z.denominator());
    return normalized(p);
  }
}

CGAL_END_NAMESPACE

#endif // QUOTIENT_COORDINATES_TO_HOMOGENEOUS_POINT_H
