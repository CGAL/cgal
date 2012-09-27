
#ifndef CGAL_PREDICATES_ALGEBRAIC_KERNEL_FOR_SPHERES_NTC3_H
#define CGAL_PREDICATES_ALGEBRAIC_KERNEL_FOR_SPHERES_NTC3_H

#include <CGAL/Root_of_2.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Interval_nt.h>

namespace CGAL {

// Because of the way the filtered predicates generator script works,
// cmp_dist_to_pointC3() must be defined _before_ ths following one.
// on a un cone avec apex sur (0,0,0)
// l'autre version avec apex sur x, y, z est straight-forward
template <class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Same_uncertainty_nt<Orientation, FT>::type
in_coneC3(const FT &px, const FT &py, const FT &pz, 
          const FT &qx, const FT &qy, const FT &qz, 
          const FT &sx, const FT &sy, const FT &sz, 
          const FT &tx, const FT &ty, const FT &tz)
{
	const FT A1 = determinant(qx, qy, qz, 
	                    sx, sy, sz,
	                    tx, ty, tz);

	const FT A3 = determinant(px, py, pz,
	                    qx, qy, qz, 
	                    tx, ty, tz);

	const FT A2 = -determinant(px, py, pz,
	                    sx, sy, sz, 
                      tx, ty, tz);

	const FT A4 = -determinant(px, py, pz,
	                    qx, qy, qz, 
	                    sx, sy, sz);
	
	const FT G1 = CGAL_NTS square(px) + CGAL_NTS square(py) + CGAL_NTS square(pz);
	const FT G2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz);
	const FT G3 = CGAL_NTS square(sx) + CGAL_NTS square(sy) + CGAL_NTS square(sz);
	const FT G4 = CGAL_NTS square(tx) + CGAL_NTS square(ty) + CGAL_NTS square(tz);
	
	// now compare A1 sqrt G1 + A3 sqrt G3 + A2 sqrt G2 + A4 sqrt G4 with 0 (evaluate the sign)

  Sign s1, s2, t, r1, r2;	
	bool unknown1 = true;
	bool unknown2 = true;
	
  // first evaluate the sign of A1 + A3 sqrt G3 
	s1 = sign(A1);
	t = sign(A3);
	if(s1 == ZERO && t == ZERO) { r1 = ZERO; unknown1 = false; }
	else if(s1 != NEGATIVE && t != NEGATIVE) { r1 = POSITIVE; unknown1 = false; }
	else if(s1 != POSITIVE && t != POSITIVE) { r1 = NEGATIVE; unknown1 = false; }

  // second evaluate the sign of A2 sqrt G2 + A4 sqrt G4 
	s2 = sign(A2);
	t = sign(A4);
	if(s2 == ZERO && t == ZERO) { r2 = ZERO; unknown2 = false; }
	else if(s2 != NEGATIVE && t != NEGATIVE) { r2 = POSITIVE; unknown2 = false; }
	else if(s2 != POSITIVE && t != POSITIVE) { r2 = NEGATIVE; unknown2 = false; }
	
	FT sqrA1G1, sqrA2G2, sqrA3G3, sqrA4G4;
	
	if(unknown1) {
		sqrA1G1 = (CGAL_NTS square(A1)) * G1;
		sqrA3G3 = ((CGAL_NTS square(A3)) * G3);
		r1 = s1 * CGAL_NTS compare(sqrA1G1,  sqrA3G3);
	}
	if(unknown2) {
		sqrA2G2 = ((CGAL_NTS square(A2)) * G2);
		sqrA4G4 = ((CGAL_NTS square(A4)) * G4);	
	  r2 = s2 * CGAL_NTS compare(sqrA2G2,  sqrA4G4);
  }

  // third solve the cases where we already know the answers
	if(r1 == ZERO && r2 == ZERO) return DEGENERATE;
	else if(r1 != NEGATIVE && r2 != NEGATIVE) return COUNTERCLOCKWISE;
	else if(r1 != POSITIVE && r2 != POSITIVE) return CLOCKWISE;
	
	if(!unknown1) {
		sqrA1G1 = (CGAL_NTS square(A1)) * G1;
		sqrA3G3 = ((CGAL_NTS square(A3)) * G3);
		unknown1 = true;		
	}
	
	if(!unknown2) {
		sqrA2G2 = ((CGAL_NTS square(A2)) * G2);
		sqrA4G4 = ((CGAL_NTS square(A4)) * G4);
	}
	
	// we need to compare |A1 sqrt G1 + A3 sqrt G3| with |A2 sqrt G2 + A4 sqrt G4|
	// which is the same to compare (A1^2 G1 + A3^2 G3 + 2 A1 A3 sqrt G1 G3) with (A2^2 G2 + A4^2 G4 + 2 A2 A4 sqrt G2 G4)
	const FT B1 = sqrA1G1 + sqrA3G3;
	const FT B3 = 2 * A1 * A3;
	const FT H3 = (G1 * G3);
	const FT B2 = (sqrA2G2 + sqrA4G4); 
	const FT B4 = 2 * A2 * A4;
	const FT H4 = (G2 * G4);	

	typename Root_of_traits< FT >::RootOf_2 xa = make_root_of_2(B1, B3, H3);
	typename Root_of_traits< FT >::RootOf_2 xb = make_root_of_2(B2, B4, H4);

  return enum_cast<Orientation>( r1 * CGAL_NTS compare(xa, xb) );
  
}

// Optimizing for interval computation
template < >
CGAL_KERNEL_MEDIUM_INLINE
Same_uncertainty_nt<Orientation, Interval_nt_advanced>::type
in_coneC3(const Interval_nt_advanced &px, const Interval_nt_advanced &py, const Interval_nt_advanced &pz, 
          const Interval_nt_advanced &qx, const Interval_nt_advanced &qy, const Interval_nt_advanced &qz, 
          const Interval_nt_advanced &sx, const Interval_nt_advanced &sy, const Interval_nt_advanced &sz, 
          const Interval_nt_advanced &tx, const Interval_nt_advanced &ty, const Interval_nt_advanced &tz)
{
	const Interval_nt_advanced A1 = determinant(qx, qy, qz, sx, sy, sz, tx, ty, tz);
  const Interval_nt_advanced A3 = determinant(px, py, pz, qx, qy, qz, tx, ty, tz);
  const Interval_nt_advanced A2 = -determinant(px, py, pz, sx, sy, sz, tx, ty, tz);
  const Interval_nt_advanced A4 = -determinant(px, py, pz, qx, qy, qz, sx, sy, sz);

  const Interval_nt_advanced G1 = CGAL_NTS square(px) + CGAL_NTS square(py) + CGAL_NTS square(pz);
  const Interval_nt_advanced G2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz);
  const Interval_nt_advanced G3 = CGAL_NTS square(sx) + CGAL_NTS square(sy) + CGAL_NTS square(sz);
  const Interval_nt_advanced G4 = CGAL_NTS square(tx) + CGAL_NTS square(ty) + CGAL_NTS square(tz);

  const Interval_nt_advanced res = (CGAL_NTS sqrt(G1)) * A1 +
               (CGAL_NTS sqrt(G2)) * A2 +
               (CGAL_NTS sqrt(G3)) * A3 +
							 (CGAL_NTS sqrt(G4)) * A4;
	
	return sign(res);
}

//supposing center at (0,0,0)
template < class FT >
inline
typename Compare<FT>::result_type
cmp_dist_to_point_projected_on_sphereC3(const FT &px, const FT &py, const FT &pz,
                                        const FT &qx, const FT &qy, const FT &qz,
                                        const FT &rx, const FT &ry, const FT &rz)
{
	const FT vpx = px;
	const FT vpy = py;
	const FT vpz = pz;
	const FT vqx = qx;
	const FT vqy = qy;
	const FT vqz = qz;
	const FT vrx = rx;
	const FT vry = ry;
	const FT vrz = rz;
	const FT dot_pq = vpx * vqx + vpy * vqy + vpz * vqz;
	const FT dot_pr = vpx * vrx + vpy * vry + vpz * vrz;
	return CGAL_NTS compare(dot_pr, dot_pq); 
	// the bigger the dot product, the shorter the length of the arc
}


} //namespace CGAL

#endif // CGAL_PREDICATES_ALGEBRAIC_KERNEL_FOR_SPHERES_NTC3_H
