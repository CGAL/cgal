inline
bool
collinearC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2)
{
  typedef Static_filter_error FT;

  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dpz = pz-rz;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dqz = qz-rz;
  return (sign_of_determinant2x2_SAF(dpx,dqx,dpy,dqy,
		epsilon_0) == ZERO)
      && (sign_of_determinant2x2_SAF(dpx,dqx,dpz,dqz,
		epsilon_1) == ZERO)
      && (sign_of_determinant2x2_SAF(dpy,dqy,dpz,dqz,
		epsilon_2) == ZERO);
}

inline
bool
collinearC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2)
{
  typedef Restricted_double FT;

  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dpz = pz-rz;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dqz = qz-rz;
  return (sign_of_determinant2x2_SAF(dpx,dqx,dpy,dqy,
		epsilon_0) == ZERO)
      && (sign_of_determinant2x2_SAF(dpx,dqx,dpz,dqz,
		epsilon_1) == ZERO)
      && (sign_of_determinant2x2_SAF(dpy,dqy,dpz,dqz,
		epsilon_2) == ZERO);
}

inline
Orientation
orientationC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    const Static_filter_error &sx,
    const Static_filter_error &sy,
    const Static_filter_error &sz,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return Orientation(sign_of_determinant3x3_SAF(qx-px,rx-px,sx-px,
                                            qy-py,ry-py,sy-py,
                                            qz-pz,rz-pz,sz-pz,
		epsilon_0));
}

inline
Orientation
orientationC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const Restricted_double &sx,
    const Restricted_double &sy,
    const Restricted_double &sz,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return Orientation(sign_of_determinant3x3_SAF(qx-px,rx-px,sx-px,
                                            qy-py,ry-py,sy-py,
                                            qz-pz,rz-pz,sz-pz,
		epsilon_0));
}

inline
Oriented_side
side_of_oriented_sphereC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    const Static_filter_error &sx,
    const Static_filter_error &sy,
    const Static_filter_error &sz,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tz,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  FT ptx = px - tx;
  FT pty = py - ty;
  FT ptz = pz - tz;
  FT pt2 = square(ptx) + square(pty) + square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = square(qtx) + square(qty) + square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = square(rtx) + square(rty) + square(rtz);
  FT stx = sx - tx;
  FT sty = sy - ty;
  FT stz = sz - tz;
  FT st2 = square(stx) + square(sty) + square(stz);
  return Oriented_side(sign_of_determinant4x4_SAF(ptx,pty,ptz,pt2,
                                              rtx,rty,rtz,rt2,
                                              qtx,qty,qtz,qt2,
                                              stx,sty,stz,st2,
		epsilon_0));
}

inline
Oriented_side
side_of_oriented_sphereC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const Restricted_double &sx,
    const Restricted_double &sy,
    const Restricted_double &sz,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tz,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  FT ptx = px - tx;
  FT pty = py - ty;
  FT ptz = pz - tz;
  FT pt2 = square(ptx) + square(pty) + square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = square(qtx) + square(qty) + square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = square(rtx) + square(rty) + square(rtz);
  FT stx = sx - tx;
  FT sty = sy - ty;
  FT stz = sz - tz;
  FT st2 = square(stx) + square(sty) + square(stz);
  return Oriented_side(sign_of_determinant4x4_SAF(ptx,pty,ptz,pt2,
                                              rtx,rty,rtz,rt2,
                                              qtx,qty,qtz,qt2,
                                              stx,sty,stz,st2,
		epsilon_0));
}

inline
Bounded_side
side_of_bounded_sphereC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    const Static_filter_error &sx,
    const Static_filter_error &sy,
    const Static_filter_error &sz,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tz,
    double & epsilon_0,
    double & epsilon_1)
{
  typedef Static_filter_error FT;

  Oriented_side s = side_of_oriented_sphereC3_SAF(px, py, pz,
                                              qx, qy, qz,
                                              rx, ry, rz,
                                              sx, sy, sz,
                                              tx, ty, tz,
		epsilon_0);
  Orientation o = orientationC3_SAF(px, py, pz,
                                qx, qy, qz,
                                rx, ry, rz,
                                sx, sy, sz,
		epsilon_1);
  return Bounded_side(s * o);
}

inline
Bounded_side
side_of_bounded_sphereC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const Restricted_double &sx,
    const Restricted_double &sy,
    const Restricted_double &sz,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tz,
    const double & epsilon_0,
    const double & epsilon_1)
{
  typedef Restricted_double FT;

  Oriented_side s = side_of_oriented_sphereC3_SAF(px, py, pz,
                                              qx, qy, qz,
                                              rx, ry, rz,
                                              sx, sy, sz,
                                              tx, ty, tz,
		epsilon_0);
  Orientation o = orientationC3_SAF(px, py, pz,
                                qx, qy, qz,
                                rx, ry, rz,
                                sx, sy, sz,
		epsilon_1);
  return Bounded_side(s * o);
}

inline
Comparison_result
cmp_dist_to_pointC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(squared_distanceC3(px,py,pz,qx,qy,qz),
                       squared_distanceC3(px,py,pz,rx,ry,rz),
		epsilon_0);
}

inline
Comparison_result
cmp_dist_to_pointC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(squared_distanceC3(px,py,pz,qx,qy,qz),
                       squared_distanceC3(px,py,pz,rx,ry,rz),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_planeC3_SAF(
    const Static_filter_error &pa,
    const Static_filter_error &pb,
    const Static_filter_error &pc,
    const Static_filter_error &pd,
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(scaled_distance_to_planeC3(pa,pb,pc,pd,px,py,pz),
                       scaled_distance_to_planeC3(pa,pb,pc,pd,qx,qy,qz),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_planeC3_SAF(
    const Restricted_double &pa,
    const Restricted_double &pb,
    const Restricted_double &pc,
    const Restricted_double &pd,
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(scaled_distance_to_planeC3(pa,pb,pc,pd,px,py,pz),
                       scaled_distance_to_planeC3(pa,pb,pc,pd,qx,qy,qz),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_planeC3_SAF(
    const Static_filter_error &ppx,
    const Static_filter_error &ppy,
    const Static_filter_error &ppz,
    const Static_filter_error &pqx,
    const Static_filter_error &pqy,
    const Static_filter_error &pqz,
    const Static_filter_error &prx,
    const Static_filter_error &pry,
    const Static_filter_error &prz,
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      px,py,pz),
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      qx,qy,qz) ,
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_planeC3_SAF(
    const Restricted_double &ppx,
    const Restricted_double &ppy,
    const Restricted_double &ppz,
    const Restricted_double &pqx,
    const Restricted_double &pqy,
    const Restricted_double &pqz,
    const Restricted_double &prx,
    const Restricted_double &pry,
    const Restricted_double &prz,
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      px,py,pz),
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      qx,qy,qz) ,
		epsilon_0);
}

