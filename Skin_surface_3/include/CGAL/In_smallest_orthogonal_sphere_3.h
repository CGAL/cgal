
#ifndef CGAL_IN_SMALLEST_ORTHOGONAL_SPHEREC3
#define CGAL_IN_SMALLEST_ORTHOGONAL_SPHEREC3

namespace CGAL {
//return the sign of the power test of weighted point (sx,sy,sz,sw)
//with respect to the smallest sphere orthogonal to
//p,q,r
template< class FT >
CGAL_MEDIUM_INLINE
Sign
in_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{
 // Translate p to origin and compute determinants
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;

  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;

  FT qq = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
  FT rr = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + CGAL_NTS square(rpz);
  FT qr = qpx*rpx + qpy*rpy + qpz*rpz;

  FT qpw = qq - qw + pw ;
  FT rpw = rr - rw + pw ;

  FT den = determinant(qq,qr,
			     qr,rr);
  FT detq = determinant(qpw,qr,
			      rpw,rr);
  FT detr = determinant(qq,qpw,
			      qr,rpw);

  // Smallest  smallest orthogonal sphere center
  // c =  detq/2*den  q + detr/2*den  r (origin at p)
  // square radius  c^2 - pw

  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;
  FT ss = CGAL_NTS square(spx) + CGAL_NTS square(spy) + CGAL_NTS  square(spz);
  FT sq = spx*qpx + spy*qpy + spz*qpz;
  FT sr = spx*rpx + spy*rpy + spz*rpz;

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );
  // return  sign of (c- s)^2 - (c^2 - pw) - sw    note that den >= 0 -
  return CGAL_NTS sign( den*(ss - sw + pw)- detq*sq - detr*sr);
}



// return the sign of the power test of weighted point (rx,ry,rz,rw)
 // with respect to the smallest sphere orthogoanal to
// p,q
template< class FT >
Sign
in_smallest_orthogonal_sphereC3(
 const FT &px, const FT &py, const FT &pz, const FT  &pw,
 const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
 const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  FT FT2(2);
  FT FT4(4);
  FT dpx = px - qx;
  FT dpy = py - qy;
  FT dpz = pz - qz;
  FT dpw = pw - qw;
  FT dp2 = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) + CGAL_NTS square(dpz);
  FT drx = rx - (px + qx)/FT2;
  FT dry = ry - (py + qy)/FT2;
  FT drz = rz - (pz + qz)/FT2;
  FT drw = rw - (pw + qw)/FT2;
  FT dr2 = CGAL_NTS square(drx) + CGAL_NTS square(dry) + CGAL_NTS square(drz);
  FT dpr = dpx*drx + dpy*dry +dpz*drz;
  return CGAL_NTS sign (dr2 - dp2/FT4 + dpr*dpw/dp2 - drw );
}




template< typename K >
class In_smallest_orthogonal_sphere_3
{
public:
  typedef typename K::Weighted_point_3               Weighted_point_3;
  typedef typename K::Sign                           Sign;

  typedef Sign             result_type;

  Sign operator() ( const Weighted_point_3 & p,
		    const Weighted_point_3 & q,
		    const Weighted_point_3 & r,
		    const Weighted_point_3 & s,
		    const Weighted_point_3 & t) const
  {
    K traits;
    typename K::Power_side_of_oriented_power_sphere_3 power_test_3 = traits.power_side_of_oriented_power_sphere_3_object();
    typename K::Orientation_3  orientation = traits.orientation_3_object();
    typename K::Orientation o = orientation(p,q,r,s);
    typename K::Oriented_side os = power_test_3(p,q,r,s,t);
    CGAL_triangulation_assertion( o != COPLANAR);
    // the minus sign below is due to the fact that power_test_3
    // return in fact minus the 5x5 determinant of lifted (p,q,r,s.t)
    return - o * os;
  }

  Sign operator() ( const Weighted_point_3 & p,
		    const Weighted_point_3 & q,
		    const Weighted_point_3 & r,
		    const Weighted_point_3 & s) const
  {
    return in_smallest_orthogonal_sphereC3(
                              p.x(), p.y(), p.z(), p.weight(),
			      q.x(), q.y(), q.z(), q.weight(),
			      r.x(), r.y(), r.z(), r.weight(),
			      s.x(), s.y(), s.z(), s.weight());
  }

  Sign operator() ( const Weighted_point_3 & p,
		    const Weighted_point_3 & q,
		    const Weighted_point_3 & r) const
  {
    return in_smallest_orthogonal_sphereC3(
                              p.x(), p.y(), p.z(), p.weight(),
			      q.x(), q.y(), q.z(), q.weight(),
			      r.x(), r.y(), r.z(), r.weight());
  }

  Sign operator() ( const Weighted_point_3 & p,
		    const Weighted_point_3 & q) const
  {
    return CGAL_NTS sign( CGAL_NTS square(p.x()-q.x()) +
			  CGAL_NTS square(p.y()-q.y()) +
			  CGAL_NTS square(p.z()-q.z()) +
			  p.weight() - q.weight());
  }

};

}

#endif
