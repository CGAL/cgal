#ifndef CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
#define CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H

CGAL_BEGIN_NAMESPACE 

template < class FT>
void
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
		       const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		       const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		       const FT &sx, const FT &sy, const FT &sz, const FT &sw,
		       FT &x, FT &y, FT &z)
{
  FT psx = px - sx;
  FT psy = py - sy;
  FT psz = pz - sz;
  FT qsx = qx - sx;
  FT qsy = qy - sy;
  FT qsz = qz - sz;
  FT rsx = rx - sx;
  FT rsy = ry - sy;
  FT rsz = rz - sz;
  FT p2w = CGAL_NTS square(px)+CGAL_NTS square(py)+CGAL_NTS  square(pz) - pw;
  FT q2w = CGAL_NTS square(qx)+CGAL_NTS square(qy)+CGAL_NTS  square(qz) - qw;
  FT r2w = CGAL_NTS square(rx)+CGAL_NTS square(ry)+CGAL_NTS  square(rz) - rw;
  FT s2w = CGAL_NTS square(sx)+CGAL_NTS square(sy)+CGAL_NTS  square(sz) - sw;
  FT psw = p2w - s2w;
  FT qsw = q2w - s2w;
  FT rsw = r2w - s2w;

  FT Dx = det3x3_by_formula(psy,qsy,rsy,
			    psz,qsz,rsz,
			    psw,qsw,rsw);
  FT Dy = det3x3_by_formula(psx,qsx,rsx,
			    psz,qsz,rsz,
			    psw,qsw,rsw); 
  FT Dz = det3x3_by_formula(psx,qsx,rsx,
			    psy,qsy,rsy,
			    psw,qsw,rsw);
  FT den = det3x3_by_formula(psx,qsx,rsx,
			     psy,qsy,rsy,
			     psz,qsz,rsz);
  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
  FT inv = FT(1)/(FT(2) * den);
  x = Dx * inv;
  y = - Dy * inv;
  z = Dz *inv;
  return;
}


template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{
  FT FT4(4);
  FT dpx = px-sx;
  FT dpy = py-sy;
  FT dpz = pz-sz;
  FT dpp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS square(dpz)
           -pw+sw;
  FT dqx = qx-sx;
  FT dqy = qy-sy;
  FT dqz = qz-sz;
  FT dqq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz)
           -qw+sw;
  FT drx = rx-sx;
  FT dry = ry-sy;
  FT drz = rz-sz;
  FT drr = CGAL_NTS square(drx)+CGAL_NTS square(dry)+CGAL_NTS square(drz)
           -rw+sw;

  FT det0 = det3x3_by_formula(dpx,dpy,dpz,dqx,dqy,dqz,drx,dry,drz);
  
  FT det1 = det3x3_by_formula(dpp,dpy,dpz,dqq,dqy,dqz,drr,dry,drz);
  FT det2 = det3x3_by_formula(dpx,dpp,dpz,dqx,dqq,dqz,drx,drr,drz);
  FT det3 = det3x3_by_formula(dpx,dpy,dpp,dqx,dqy,dqq,drx,dry,drr);

  return
    (CGAL_NTS square(det1)+CGAL_NTS square(det2)+CGAL_NTS square(det3))/
    (FT4*CGAL_NTS square(det0)) - sw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  // resolution of the system (where we note c the center)
  // |       dc^2 = cw + rw
  // |  (dp-dc)^2 = pw + cw
  // |  (dq-dc)^2 = qw + cw
  // |         dc = Lamdba*dp + Mu*dq

  FT FT2(2);
  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dpz = pz-rz;
  FT dp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS square(dpz);
  FT dpp = dp-pw+rw;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dqz = qz-rz;
  FT dq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz);
  FT dqq = dq-qw+rw;
  FT dpdq = dpx*dqx+dpy*dqy+dpz*dqz;
  FT denom = FT2*(dp*dq-CGAL_NTS square(dpdq));
  FT Lambda = (dpp*dq-dqq*dpdq)/denom;
  FT Mu = (dqq*dp-dpp*dpdq)/denom;

  return (CGAL_NTS square(Lambda)*dp+CGAL_NTS square(Mu)*dq
	  +FT2*Lambda*Mu*dpdq - rw);
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{ 
  FT FT4(4);
  FT dp = CGAL_NTS square(px-qx)+CGAL_NTS square(py-qy)+CGAL_NTS square(pz-qz);

  return (CGAL_NTS square(dp-pw+qw)/(FT4*dp)-qw);
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





CGAL_END_NAMESPACE
#endif //CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
