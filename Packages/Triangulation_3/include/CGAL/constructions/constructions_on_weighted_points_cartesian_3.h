#ifndef CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
#define CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H

CGAL_BEGIN_NAMESPACE 

template <class FT>
void
determinants_for_weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &num_x,  FT &num_y, FT &num_z, FT& den)
{
  // translate origin to p
  // and compute determinants for weighted_circumcenter and
  // circumradius
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + 
           CGAL_NTS square(rpz) - rw + pw;
  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;
  FT sp2 = CGAL_NTS square(spx) + CGAL_NTS square(spy) + 
           CGAL_NTS square(spz) - sw + pw;

  num_x = det3x3_by_formula(qpy,qpz,qp2,
			    rpy,rpz,rp2,
			    spy,spz,sp2);
  num_y = det3x3_by_formula(qpx,qpz,qp2,
			    rpx,rpz,rp2,
			    spx,spz,sp2);
  num_z = det3x3_by_formula(qpx,qpy,qp2,
			    rpx,rpy,rp2,
			    spx,spy,sp2);
  den   = det3x3_by_formula(qpx,qpy,qpz,
			    rpx,rpy,rpz,
			    spx,spy,spz);
}


template < class FT>
void
weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &x, FT &y, FT &z)
{
  // this function  compute the weighted circumcenter point only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x,  num_y, num_z,den);

  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}

template < class FT>
void
weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &x, FT &y, FT &z, FT &w)
{
  // this function  compute the weighted circumcenter point 
  // and the squared weighted circumradius
  
  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x,  num_y, num_z, den);

  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
  
  w = (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
      *CGAL_NTS square(inv) - pw;
}


template< class FT >
FT
squared_radius_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{

  // this function  compute the squared weighted circumradius only
  
  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x, num_y, num_z,den);

  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

   return
    (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
    *CGAL_NTS square(inv) - pw;
}


template <class FT>
void
determinants_for_weighted_circumcenterC3(
	        const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		FT &num_x,  FT &num_y, FT &num_z, FT &den)
{
  // translate origin to p
  // and compute determinants for weighted_circumcenter and
  // circumradius

  // Translate s to origin to simplify the expression.
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + 
           CGAL_NTS square(rpz) - rw + pw;

  FT sx = qpy*rpz-qpz*rpy;
  FT sy = qpz*rpx-qpx*rpz;
  FT sz = qpx*rpy-qpy*rpx;

  // The following determinants can be developped and simplified.
  //
  // FT num_x = det3x3_by_formula(qpy,qpz,qp2,
  //                              rpy,rpz,rp2,
  //                              sy,sz,FT(0));
  // FT num_y = det3x3_by_formula(qpx,qpz,qp2,
  //                              rpx,rpz,rp2,
  //                              sx,sz,FT(0));
  // FT num_z = det3x3_by_formula(qpx,qpy,qp2,
  //                              rpx,rpy,rp2,
  //                              sx,sy,FT(0));

  num_x = qp2 * det2x2_by_formula(rpy,rpz,sy,sz)
        - rp2 * det2x2_by_formula(qpy,qpz,sy,sz);

  num_y = qp2 * det2x2_by_formula(rpx,rpz,sx,sz)
	- rp2 * det2x2_by_formula(qpx,qpz,sx,sz);

  num_z = qp2 * det2x2_by_formula(rpx,rpy,sx,sy)
	- rp2 * det2x2_by_formula(qpx,qpy,sx,sy);

  den   = det3x3_by_formula(qpx,qpy,qpz,
			    rpx,rpy,rpz,
			    sx,sy,sz);
}

template < class FT >
void
weighted_circumcenterC3( 
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		  FT &x, FT &y, FT &z)
{
  // this function  compute the weighted circumcenter point only

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_kernel_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}


template < class FT >
void
weighted_circumcenterC3( 
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		  FT &x, FT &y, FT &z, FT &w)
{
  // this function  compute the weighted circumcenter and
  // the weighted squared circumradius

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_kernel_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;

  w = (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
      *CGAL_NTS square(inv)  - pw;
}


template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  // this function  compute the weighted squared circumradius only

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_kernel_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

 return
    (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
     *CGAL_NTS square(inv)  - pw;
}



template < class FT >
void
weighted_circumcenterC3( 
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  FT &x, FT &y, FT &z)
{
// this function  compute the weighted circumcenter point only
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;
  
  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;
} 

 
template < class FT >
void
weighted_circumcenterC3( 
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  FT &x, FT &y, FT &z, FT &w)
{
 // this function  compute the weighted circumcenter point and
  // the weighted circumradius
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;
  
  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;

  w = CGAL_NTS square(alpha)*qp2 - pw;
} 


template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{ 
  // this function  computes
  // the weighted circumradius only
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;
  
  return  CGAL_NTS square(alpha)*qp2 - pw;
}


template< class FT >
FT
power_productC3( 
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{ 
  // computes the power product of two weighted points
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + 
           CGAL_NTS square(qpz);
  return qp2 - pw - qw ;
}

template < class RT , class We>
void
radical_axisC3(const RT &px, const RT &py, const RT &pz, const We &pw,
	       const RT &qx, const RT &qy, const RT &qz, const We &qw,
	       const RT &rx, const RT &ry, const RT &rz, const We &rw,
	       RT &a, RT &b, RT& c )
{
  RT dqx=qx-px, dqy=qy-py, dqz=qz-pz, drx=rx-px, dry=ry-py, drz=rz-pz;

  //il manque des tests...

  a= RT(1)*det2x2_by_formula(dqy, dqz, dry, drz);
  b= - RT(1)*det2x2_by_formula(dqx, dqz, drx, drz);
  c= RT(1)*det2x2_by_formula(dqx, dqy, drx, dry);
}





 // I will use this to test if the radial axis of three spheres
  // intersect the triangle formed by the centers.
//   // resolution of the system (where we note c the center)
//   // |       dc^2 = cw + rw
//   // |  (dp-dc)^2 = pw + cw
//   // |  (dq-dc)^2 = qw + cw
//   // |         dc = Lamdba*dp + Mu*dq

//   FT FT2(2);
//   FT dpx = px-rx;
//   FT dpy = py-ry;
//   FT dpz = pz-rz;
//   FT dp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS  square(dpz);
//   FT dpp = dp-pw+rw;
//   FT dqx = qx-rx;
//   FT dqy = qy-ry;
//   FT dqz = qz-rz;
//   FT dq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz);
//   FT dqq = dq-qw+rw;
//   FT dpdq = dpx*dqx+dpy*dqy+dpz*dqz;
//   FT denom = FT2*(dp*dq-CGAL_NTS square(dpdq));
//   FT Lambda = (dpp*dq-dqq*dpdq)/denom;
//   FT Mu = (dqq*dp-dpp*dpdq)/denom;

//   return (CGAL_NTS square(Lambda)*dp+CGAL_NTS square(Mu)*dq
// 	  +FT2*Lambda*Mu*dpdq - rw);




CGAL_END_NAMESPACE
#endif //CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
