#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAIS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAIS_2_H

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Weighted_point_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>


template < class R , class W >
class CGAL_Regular_triangulation_euclidean_traits_2
  : public CGAL_Triangulation_euclidean_traits_2<R>
{
public:
  typedef R Rep;
  typedef W Weight;
  typedef CGAL_Triangulation_euclidean_traits_2 <R>  Traits;
  typedef Traits::Point  Bare_point;
  typedef Traits::Segment Segment;
  typedef Traits::Triangle  Triangle;
  typedef Traits::Line Line;
  typedef Traits::Direction  Direction;
  typedef Traits::Ray  Ray;
  typedef CGAL_Weighted_point_2<Bare_point, W> Weighted_point;
  typedef Weighted_point Point;

 
  CGAL_Oriented_side power_test(const Weighted_point &p,
				const Weighted_point &q,
				const Weighted_point &r,
				const Weighted_point &s) const
  {
    if (compare(p,s)) return CGAL_ON_ORIENTED_BOUNDARY;
    if (compare(q,s)) return CGAL_ON_ORIENTED_BOUNDARY;
    if (compare(r,s)) return CGAL_ON_ORIENTED_BOUNDARY;

    return CGAL_power_test(p, q, r, s);
  }
	
  CGAL_Oriented_side power_test(const Weighted_point &p,
				const Weighted_point &q,
				const Weighted_point &s) const
  {
    if (compare(p,s)) return CGAL_ON_ORIENTED_BOUNDARY;
    if (compare(q,s)) return CGAL_ON_ORIENTED_BOUNDARY;
    return CGAL_power_test(p, q, s);
  }

};

// power test for 2 dimensions triangulation

template <class FT,class Weight >
CGAL_Oriented_side
CGAL_power_test(const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &p,
                const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &q,
                const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &r,
                const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &test)
{
	FT FT0(0);
	FT FT1(1);
	const FT &px = p.x();
	const FT &py = p.y();
	const FT &pw = p.weight();
	const FT &qx = q.x();
	const FT &qy = q.y();
	const FT &qw = q.weight();
	const FT &rx = r.x();
	const FT &ry = r.y();
	const FT &rw = r.weight();
	const FT &tx = test.x();
	const FT &ty = test.y();
	const FT &tw = test.weight();

	CGAL_kernel_precondition( ! CGAL_collinear(p,q,r) );

	FT det = CGAL_det4x4_by_formula(px, py, px*px + py*py -pw*pw, FT1,
	                                qx, qy, qx*qx + qy*qy -qw*qw, FT1,
	                                rx, ry, rx*rx + ry*ry -rw*rw, FT1,
	                                tx, ty, tx*tx + ty*ty -tw*tw, FT1);

	return (det<FT0) ? CGAL_ON_NEGATIVE_SIDE  : 
	  ((det==FT0) ? CGAL_ON_ORIENTED_BOUNDARY : CGAL_ON_POSITIVE_SIDE);
}


template < class RT,class Weight >
CGAL_Oriented_side
CGAL_power_test( 
const CGAL_Weighted_point_2<
CGAL_Point_2<
CGAL_Homogeneous<RT> >, 
Weight>& q,
const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& r,
const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& s,
const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& t)
{
	const RT& qhx = q.hx();
	const RT& qhy = q.hy();
	const RT& qhw = q.hw();
	const RT& qhwt = q.weight();
	const RT& rhx = r.hx();
	const RT& rhy = r.hy();
	const RT& rhw = r.hw();
	const RT& rhwt = r.weight();
	const RT& shx = s.hx();
	const RT& shy = s.hy();
	const RT& shw = s.hw();
	const RT& shwt = s.weight();
	const RT& thx = t.hx();
	const RT& thy = t.hy();
	const RT& thw = t.hw();
	const RT& thwt = t.weight();
	const RT  RT0 = RT(0);

	CGAL_kernel_precondition( ! CGAL_collinear(q,r,s) );

	// compute sign of      |qx  qy  qx^2+qy^2-qwt^2  1 |   | a b c d |
	//                      |        --   r   --        | = | e f g h |
	//     determinant      |        --   s   --        | = | i j k l |
	//                      |        --   t   --        |   | m n o p |
	//           where

	RT a = qhx*qhw;
	RT b = qhy*qhw;
	RT c = qhx*qhx + qhy*qhy - qhwt*qhwt;
	RT d = qhw*qhw;

	RT e = rhx*rhw;
	RT f = rhy*rhw;
	RT g = rhx*rhx + rhy*rhy - rhwt*rhwt;
	RT h = rhw*rhw;

	RT i = shx*shw;
	RT j = shy*shw;
	RT k = shx*shx + shy*shy - shwt*shwt;
	RT l = shw*shw;

	RT m = thx*thw;
	RT n = thy*thw;
	RT o = thx*thx + thy*thy - thwt*thwt;
	RT p = thw*thw;

	RT det = CGAL_det4x4_by_formula(a,b,c,d,
	                                e,f,g,h,
	                                i,j,k,l,
	                                m,n,o,p);


	return (det<RT0) ? CGAL_ON_NEGATIVE_SIDE
	                 : ((det==RT0) ? CGAL_ON_ORIENTED_BOUNDARY : CGAL_ON_POSITIVE_SIDE);
}



// power test for 1 dimension triangulation

template <class FT,class Weight >
CGAL_Oriented_side
CGAL_power_test(const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &p,
                const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &q,
                const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Cartesian<FT> >, Weight> &test)
{
	FT FT0(0);
	FT FT1(1);
	FT det;
	const FT &px = p.x();
	const FT &py = p.y();
	const FT &pw = p.weight();
	const FT &qx = q.x();
	const FT &qy = q.y();
	const FT &qw = q.weight();
	const FT &tx = test.x();
	const FT &ty = test.y();
	const FT &tw = test.weight();

	CGAL_kernel_precondition( CGAL_collinear(p.point(),q.point(),test.point()) );
	CGAL_kernel_precondition( p.point() != q.point() );
	if (px!=qx)
	{	
		det = CGAL_det3x3_by_formula(px, px*px + py*py -pw*pw, FT1,
		                             qx, qx*qx + qy*qy -qw*qw, FT1,
		                             tx, tx*tx + ty*ty -tw*tw, FT1);
		if (px<qx)
			det=-det;
	}
	else
	{
		det = CGAL_det3x3_by_formula(py, px*px + py*py -pw*pw, FT1,
		                             qy, qx*qx + qy*qy -qw*qw, FT1,
		                             ty, tx*tx + ty*ty -tw*tw, FT1);
		if (py<qy)
			det=-det;
	}
	return (det<FT0) ? CGAL_ON_NEGATIVE_SIDE
	                 : ((det==FT0) ? CGAL_ON_ORIENTED_BOUNDARY : CGAL_ON_POSITIVE_SIDE);
};



template < class RT,class Weight >
CGAL_Oriented_side
CGAL_power_test( const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& q,
                 const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& r,
                 const CGAL_Weighted_point_2<CGAL_Point_2<CGAL_Homogeneous<RT> >, Weight>& t)
{
	const RT& qhx = q.hx();
	const RT& qhy = q.hy();
	const RT& qhw = q.hw();
	const RT& qhwt = q.weight();
	const RT& rhx = r.hx();
	const RT& rhy = r.hy();
	const RT& rhw = r.hw();
	const RT& rhwt = r.weight();
	const RT& shx = s.hx();
	const RT& shy = s.hy();
	const RT& shw = s.hw();
	const RT& shwt = s.weight();
	const RT& thx = t.hx();
	const RT& thy = t.hy();
	const RT& thw = t.hw();
	const RT& thwt = t.weight();
	const RT  RT0 = RT(0);

	CGAL_kernel_precondition( CGAL_collinear(q,r,t) );
	CGAL_kernel_precondition( q.point() != r.point() );
	if (qhx != rhx )
	{	RT a = qhx*qhw;
		RT e = rhx*rhw;
		RT i = thx*thw;
	}
	else
	{	RT a = qhy*qhw;
		RT e = rhy*rhw;
		RT i = thy*thw;
	}

	RT c = qhx*qhx + qhy*qhy - qhwt*qhwt;
	RT d = qhw*qhw;

	RT g = rhx*rhx + rhy*rhy - rhwt*rhwt;
	RT h = rhw*rhw;

	RT k = thx*thx + thy*thy - thwt*thwt;
	RT l = thw*thw;

	RT det = CGAL_det4x4_by_formula(a,c,d,
	                                e,g,h,
	                                i,k,l);

	if (a<e)
		det=-det;

	return (det<RT0) ? CGAL_ON_NEGATIVE_SIDE
	                 : ((det==RT0) ? CGAL_ON_ORIENTED_BOUNDARY : CGAL_ON_POSITIVE_SIDE);
}

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAIS_2_H
