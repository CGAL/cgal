#ifndef CGAL_SVD_DO_INTERSECT_C2_H
#define CGAL_SVD_DO_INTERSECT_C2_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

//---------------------------------------------------------------------

template<class RT>
std::pair<int,int>
svd_do_intersect_C2(const RT& x1, const RT& y1,
		    const RT& x2, const RT& y2,
		    const RT& x3, const RT& y3,
		    const RT& x4, const RT& y4)
{
  RT delta = -det2x2_by_formula(x2 - x1, x4 - x3, y2 - y1, y4 - y3);

  Sign s = CGAL::sign( delta );
  if ( s != CGAL::ZERO ) {
    return svd_do_intersect_non_parallel_C2(x1, y1, x2, y2,
					    x3, y3, x4, y4, delta);
  } else {
    return svd_do_intersect_parallel_C2(x1, y1, x2, y2,
					x3, y3, x4, y4, delta);
  }
}


//---------------------------------------------------------------------


template<class RT>
std::pair<int,int>
svd_do_intersect_non_parallel_C2(const RT& x1, const RT& y1,
				 const RT& x2, const RT& y2,
				 const RT& x3, const RT& y3,
				 const RT& x4, const RT& y4,
				 const RT& D)
{
  RT Dt = -det2x2_by_formula(x3 - x1, x4 - x3,
			     y3 - y1, y4 - y3);

  RT Ds = det2x2_by_formula(x2 - x1, x3 - x1,
			    y2 - y1, y3 - y1);

  Sign s_D = CGAL::sign( D );
  Sign s_Dt = CGAL::sign( Dt );
  Sign s_Ds = CGAL::sign( Ds );

  Sign s_tdiff = CGAL::sign(Dt - D);
  Sign s_sdiff = CGAL::sign(Ds - D);

  Sign s_t = Sign(s_Dt * s_D);
  Sign s_s = Sign(s_Ds * s_D);

  Sign s_t_minus_1 = Sign(s_tdiff * s_D);
  Sign s_s_minus_1 = Sign(s_sdiff * s_D);

  if ( s_t == CGAL::NEGATIVE || s_t_minus_1 == CGAL::POSITIVE ||
       s_s == CGAL::NEGATIVE || s_s_minus_1 == CGAL::POSITIVE ) {
    //  t < 0 or t > 1 or s < 0 or s > 1
    return std::pair<int,int>(3,3);
  }

  int it(0), is(0);
  if ( s_t == CGAL::ZERO ) {
    it = 0;
  } else if ( s_t_minus_1 == CGAL::ZERO ) {
    it = 1;
  } else {
    it = 2;
  }
  if ( s_s == CGAL::ZERO ) {
    is = 0;
  } else if ( s_s_minus_1 == CGAL::ZERO ) {
    is = 1;
  } else {
    is = 2;
  }
  return std::pair<int,int>(it, is);
}




//---------------------------------------------------------------------

template<class RT>
std::pair<int,int>
svd_do_intersect_parallel_C2(const RT& x1, const RT& y1,
			     const RT& x2, const RT& y2,
			     const RT& x3, const RT& y3,
			     const RT& x4, const RT& y4,
			     const RT& D)
{
  RT D1 = det2x2_by_formula(x2 - x1, x3 - x1,
			    y2 - y1, y3 - y1);

  if ( CGAL::sign( D1 ) != CGAL::ZERO ) {
    return std::pair<int,int>(3,3);
  }

  RT Dt3, Dt4, Dt;
  if ( CGAL::compare(x2, x1) != CGAL::EQUAL ) {
    Dt  = x2 - x1;
    Dt3 = x3 - x1;
    Dt4 = x4 - x1;
  } else {
    Dt  = y2 - y1;
    Dt3 = y3 - y1;
    Dt4 = y4 - y1;
  }

  Sign s_Dt = CGAL::sign( Dt );
  Sign s_Dt3 = CGAL::sign( Dt3 );
  Sign s_Dt4 = CGAL::sign( Dt4 );

  Sign s_t3 = Sign(s_Dt3 * s_Dt);
  Sign s_t4 = Sign(s_Dt4 * s_Dt);

  Sign s_t3diff = CGAL::sign( Dt3 - Dt );
  Sign s_t4diff = CGAL::sign( Dt4 - Dt );

  Sign s_t3_minus_1 = Sign(s_t3diff * s_Dt);
  Sign s_t4_minus_1 = Sign(s_t4diff * s_Dt);

  if ( (s_t3 == CGAL::NEGATIVE || s_t3_minus_1 == CGAL::POSITIVE) &&
       (s_t4 == CGAL::NEGATIVE || s_t4_minus_1 == CGAL::POSITIVE) ) {
    //  (t3 < 0 or t3 > 1) and (t4 < 0 or t4 > 1)
    return std::pair<int,int>(3,3); // no intersection
  }

  int it3(0), it4(0);
  if ( s_t3 == CGAL::ZERO ) { // t3 == 0
    it3 = 0;
  } else if ( s_t3_minus_1 == CGAL::ZERO ) { // t3 == 1
    it3 = 1;
  } else if ( s_t3 == CGAL::POSITIVE &&
	      s_t3_minus_1 == CGAL::NEGATIVE ) { // 0 < t3 < 1
    it3 = 2;
  } else { // t3 < 0 or t3 > 1
    it3 = 3;
  }

  if ( s_t4 == CGAL::ZERO ) { // t4 == 0
    it4 = 0;
  } else if ( s_t4_minus_1 == CGAL::ZERO ) { // t4 == 1
    it4 = 1;
  } else if ( s_t4 == CGAL::POSITIVE &&
	      s_t4_minus_1 == CGAL::NEGATIVE ) { // 0 < t4 < 1
    it4 = 2;
  } else { // t4 < 0 or t4 > 1
    it4 = 3;
  }

  if ( it3 < 2 && it4 < 2 ) { // segments are identical
    return std::pair<int,int>(4,4);
  } else if ( it3 < 2 && it4 == 3 ) { // segments intersect at p1 or p2
    return std::pair<int,int>(it3,it3);
  } else if ( it3 == 3 && it4 < 2 ) { // segments intersect at p1 or p2
    return std::pair<int,int>(it4,it4);
  } else {
    // MK: this case has to be further investigating to produce finer
    //     answers wrt the exact configuration
    return std::pair<int,int>(5,5);
  }
  
}


//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class K>
class Svd_do_intersect_C2
{
public:
  //  typedef typename K::Point_2      Point_2;
  typedef typename K::Segment_2    Segment_2;

  //  typedef std::pair<Intersection_type,Intersection_type>  result_type;
  typedef std::pair<int,int>       result_type;

  enum
    { FIRST_ENDPOINT = 0, SECOND_ENDPOINT, INTERIOR_POINT,
      NO_INTERSECTION, IDENTICAL, SUBSEGMENT
    } Intersection_type;

public:
  result_type
  operator()(const Segment_2& s1, const Segment_2& s2) const
  {
    std::pair<int,int> res =
      svd_do_intersect_C2( s1.source().x(), s1.source().y(),
			   s1.target().x(), s1.target().y(),
			   s2.source().x(), s2.source().y(),
			   s2.target().x(), s2.target().y() );

    return res;

    //    Intersection_type it1 = res.first;
    //    Intersection_type it2 = res.second;

    //    return result_type(it1, it2);
  }
};

//---------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_DO_INTERSECT_C2_H
