inline
Comparison_result
compare_xC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &l1a,
    const Static_filter_error &l1b,
    const Static_filter_error &l1c,
    const Static_filter_error &l2a,
    const Static_filter_error &l2b,
    const Static_filter_error &l2c,
    double & epsilon_0,
    double & epsilon_1)
{
  typedef Static_filter_error FT;

  Sign sign1 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_0);
  Sign sign2 = sign_of_determinant3x3_SAF(l1a, l1b, l1c,
                                      l2a, l2b, l2c,
                                      -FT(1), FT(0), px,
		epsilon_1);
  CGAL_kernel_assertion( sign1 != 0 );
  return Comparison_result (sign1 * sign2);
}

inline
Comparison_result
compare_xC2_SAF(
    const Restricted_double &px,
    const Restricted_double &l1a,
    const Restricted_double &l1b,
    const Restricted_double &l1c,
    const Restricted_double &l2a,
    const Restricted_double &l2b,
    const Restricted_double &l2c,
    const double & epsilon_0,
    const double & epsilon_1)
{
  typedef Restricted_double FT;

  Sign sign1 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_0);
  Sign sign2 = sign_of_determinant3x3_SAF(l1a, l1b, l1c,
                                      l2a, l2b, l2c,
                                      -FT(1), FT(0), px,
		epsilon_1);
  CGAL_kernel_assertion( sign1 != 0 );
  return Comparison_result (sign1 * sign2);
}

inline
Comparison_result
compare_xC2_SAF(
    const Static_filter_error &l1a,
    const Static_filter_error &l1b,
    const Static_filter_error &l1c,
    const Static_filter_error &l2a,
    const Static_filter_error &l2b,
    const Static_filter_error &l2c,
    const Static_filter_error &h1a,
    const Static_filter_error &h1b,
    const Static_filter_error &h1c,
    const Static_filter_error &h2a,
    const Static_filter_error &h2b,
    const Static_filter_error &h2c,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2)
{
  typedef Static_filter_error FT;

  Sign sign1 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_0);
  Sign sign2 = sign_of_determinant2x2_SAF(h1a, h1b, h2a, h2b,
		epsilon_1);
  
  
  FT FT0(0);
  Sign sign3 = sign_of_determinant4x4_SAF(l1a, l1b, FT0, l1c,
                                      l2a, l2b, FT0, l2c,
                                      h1a, FT0, h1b, h1c,
                                      h2a, FT0, h2b, h2c,
		epsilon_2);
  CGAL_kernel_assertion( (sign1 != 0) && (sign2 != 0) );
  return Comparison_result (- sign1 * sign2 * sign3);
}

inline
Comparison_result
compare_xC2_SAF(
    const Restricted_double &l1a,
    const Restricted_double &l1b,
    const Restricted_double &l1c,
    const Restricted_double &l2a,
    const Restricted_double &l2b,
    const Restricted_double &l2c,
    const Restricted_double &h1a,
    const Restricted_double &h1b,
    const Restricted_double &h1c,
    const Restricted_double &h2a,
    const Restricted_double &h2b,
    const Restricted_double &h2c,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2)
{
  typedef Restricted_double FT;

  Sign sign1 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_0);
  Sign sign2 = sign_of_determinant2x2_SAF(h1a, h1b, h2a, h2b,
		epsilon_1);
  
  
  FT FT0(0);
  Sign sign3 = sign_of_determinant4x4_SAF(l1a, l1b, FT0, l1c,
                                      l2a, l2b, FT0, l2c,
                                      h1a, FT0, h1b, h1c,
                                      h2a, FT0, h2b, h2c,
		epsilon_2);
  CGAL_kernel_assertion( (sign1 != 0) && (sign2 != 0) );
  return Comparison_result (- sign1 * sign2 * sign3);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &la,
    const Static_filter_error &lb,
    const Static_filter_error &lc,
    double & epsilon_0,
    double & epsilon_1)
{
  typedef Static_filter_error FT;

  Sign sign1 = CGAL::sign_SAF(lb,
		epsilon_0);
  Sign sign2 = CGAL::sign_SAF(la*px + lb*py + lc,
		epsilon_1);
  CGAL_kernel_assertion( sign1 != 0 );
  return Comparison_result (sign1 * sign2);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &la,
    const Restricted_double &lb,
    const Restricted_double &lc,
    const double & epsilon_0,
    const double & epsilon_1)
{
  typedef Restricted_double FT;

  Sign sign1 = CGAL::sign_SAF(lb,
		epsilon_0);
  Sign sign2 = CGAL::sign_SAF(la*px + lb*py + lc,
		epsilon_1);
  CGAL_kernel_assertion( sign1 != 0 );
  return Comparison_result (sign1 * sign2);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &l1a,
    const Static_filter_error &l1b,
    const Static_filter_error &l1c,
    const Static_filter_error &l2a,
    const Static_filter_error &l2b,
    const Static_filter_error &l2c,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2)
{
  typedef Static_filter_error FT;

  Sign sign1 = CGAL::sign_SAF(l1b,
		epsilon_0);
  Sign sign2 = CGAL::sign_SAF(l2b,
		epsilon_1);
  Sign sign3 = sign_of_determinant2x2_SAF(l1a*px+l1c,l2a*px+l2c,l1b,l2b,
		epsilon_2);
  CGAL_kernel_assertion( (sign1 != 0) && (sign2 != 0) );
  return Comparison_result (- sign1 * sign2 * sign3);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Restricted_double &px,
    const Restricted_double &l1a,
    const Restricted_double &l1b,
    const Restricted_double &l1c,
    const Restricted_double &l2a,
    const Restricted_double &l2b,
    const Restricted_double &l2c,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2)
{
  typedef Restricted_double FT;

  Sign sign1 = CGAL::sign_SAF(l1b,
		epsilon_0);
  Sign sign2 = CGAL::sign_SAF(l2b,
		epsilon_1);
  Sign sign3 = sign_of_determinant2x2_SAF(l1a*px+l1c,l2a*px+l2c,l1b,l2b,
		epsilon_2);
  CGAL_kernel_assertion( (sign1 != 0) && (sign2 != 0) );
  return Comparison_result (- sign1 * sign2 * sign3);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Static_filter_error &l1a,
    const Static_filter_error &l1b,
    const Static_filter_error &l1c,
    const Static_filter_error &l2a,
    const Static_filter_error &l2b,
    const Static_filter_error &l2c,
    const Static_filter_error &ha,
    const Static_filter_error &hb,
    const Static_filter_error &hc,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2,
    double & epsilon_3)
{
  typedef Static_filter_error FT;

  Sign sign0 = sign_of_determinant2x2_SAF(l1a,l1b,l2a,l2b,
		epsilon_0);
  Sign sign1 = sign_of_determinant3x3_SAF(ha,hb,hc,l1a,l1b,l1c,l2a,l2b,l2c,
		epsilon_1);
  CGAL_kernel_assertion( (sign0 != ZERO) && (sign_SAF(hb,
		epsilon_2) != ZERO) );
  return Comparison_result (sign0 * CGAL::sign_SAF(hb,
		epsilon_3) * sign1);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Restricted_double &l1a,
    const Restricted_double &l1b,
    const Restricted_double &l1c,
    const Restricted_double &l2a,
    const Restricted_double &l2b,
    const Restricted_double &l2c,
    const Restricted_double &ha,
    const Restricted_double &hb,
    const Restricted_double &hc,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2,
    const double & epsilon_3)
{
  typedef Restricted_double FT;

  Sign sign0 = sign_of_determinant2x2_SAF(l1a,l1b,l2a,l2b,
		epsilon_0);
  Sign sign1 = sign_of_determinant3x3_SAF(ha,hb,hc,l1a,l1b,l1c,l2a,l2b,l2c,
		epsilon_1);
  CGAL_kernel_assertion( (sign0 != ZERO) && (sign_SAF(hb,
		epsilon_2) != ZERO) );
  return Comparison_result (sign0 * CGAL::sign_SAF(hb,
		epsilon_3) * sign1);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Static_filter_error &l1a,
    const Static_filter_error &l1b,
    const Static_filter_error &l1c,
    const Static_filter_error &l2a,
    const Static_filter_error &l2b,
    const Static_filter_error &l2c,
    const Static_filter_error &h1a,
    const Static_filter_error &h1b,
    const Static_filter_error &h1c,
    const Static_filter_error &h2a,
    const Static_filter_error &h2b,
    const Static_filter_error &h2c,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2,
    double & epsilon_3)
{
  typedef Static_filter_error FT;

  FT FT0(0);
  Sign s1 = lexicographical_sign_SAF(h1b, -h1a,
		epsilon_0);
  Sign s2 = lexicographical_sign_SAF(h2b, -h2a,
		epsilon_1);
  Sign s3 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_2);
  Sign s4 = sign_of_determinant4x4_SAF(h2a, h2b, FT0, h2c,
                                   l1a, FT0, l1b, l1c,
                                   l2a, FT0, l2b, l2c,
                                   h1a, h1b, FT0, h1c,
		epsilon_3);
  return Comparison_result (s1 * s2 * s3 * s4);
}

inline
Comparison_result
compare_y_at_xC2_SAF(
    const Restricted_double &l1a,
    const Restricted_double &l1b,
    const Restricted_double &l1c,
    const Restricted_double &l2a,
    const Restricted_double &l2b,
    const Restricted_double &l2c,
    const Restricted_double &h1a,
    const Restricted_double &h1b,
    const Restricted_double &h1c,
    const Restricted_double &h2a,
    const Restricted_double &h2b,
    const Restricted_double &h2c,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2,
    const double & epsilon_3)
{
  typedef Restricted_double FT;

  FT FT0(0);
  Sign s1 = lexicographical_sign_SAF(h1b, -h1a,
		epsilon_0);
  Sign s2 = lexicographical_sign_SAF(h2b, -h2a,
		epsilon_1);
  Sign s3 = sign_of_determinant2x2_SAF(l1a, l1b, l2a, l2b,
		epsilon_2);
  Sign s4 = sign_of_determinant4x4_SAF(h2a, h2b, FT0, h2c,
                                   l1a, FT0, l1b, l1c,
                                   l2a, FT0, l2b, l2c,
                                   h1a, h1b, FT0, h1c,
		epsilon_3);
  return Comparison_result (s1 * s2 * s3 * s4);
}

inline
Comparison_result
compare_deltax_deltayC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &qx,
    const Static_filter_error &ry,
    const Static_filter_error &sy,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

    return CGAL::compare_SAF(abs(px-qx), abs(ry-sy),
		epsilon_0);
}

inline
Comparison_result
compare_deltax_deltayC2_SAF(
    const Restricted_double &px,
    const Restricted_double &qx,
    const Restricted_double &ry,
    const Restricted_double &sy,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

    return CGAL::compare_SAF(abs(px-qx), abs(ry-sy),
		epsilon_0);
}

inline
Orientation
orientationC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return sign_of_determinant2x2_SAF(px-rx,py-ry,qx-rx,qy-ry,
		epsilon_0);
}

inline
Orientation
orientationC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return sign_of_determinant2x2_SAF(px-rx,py-ry,qx-rx,qy-ry,
		epsilon_0);
}

inline
Oriented_side
side_of_oriented_circleC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  
  
  
  
  
  
  FT ptx = px-tx;
  FT pty = py-ty;
  FT qtx = qx-tx;
  FT qty = qy-ty;
  FT rtx = rx-tx;
  FT rty = ry-ty;
  return Oriented_side(
           sign_of_determinant3x3_SAF(ptx, pty, square(ptx) + square(pty),
                                  qtx, qty, square(qtx) + square(qty),
                                  rtx, rty, square(rtx) + square(rty),
		epsilon_0));
}

inline
Oriented_side
side_of_oriented_circleC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  
  
  
  
  
  
  FT ptx = px-tx;
  FT pty = py-ty;
  FT qtx = qx-tx;
  FT qty = qy-ty;
  FT rtx = rx-tx;
  FT rty = ry-ty;
  return Oriented_side(
           sign_of_determinant3x3_SAF(ptx, pty, square(ptx) + square(pty),
                                  qtx, qty, square(qtx) + square(qty),
                                  rtx, rty, square(rtx) + square(rty),
		epsilon_0));
}

inline
Bounded_side
side_of_bounded_circleC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    double & epsilon_0,
    double & epsilon_1)
{
  typedef Static_filter_error FT;

  Oriented_side s = side_of_oriented_circleC2_SAF(px,py,qx,qy,rx,ry,tx,ty,
		epsilon_0);
  Orientation o = orientationC2_SAF(px,py,qx,qy,rx,ry,
		epsilon_1);

  return Bounded_side (s * o);
}

inline
Bounded_side
side_of_bounded_circleC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const double & epsilon_0,
    const double & epsilon_1)
{
  typedef Restricted_double FT;

  Oriented_side s = side_of_oriented_circleC2_SAF(px,py,qx,qy,rx,ry,tx,ty,
		epsilon_0);
  Orientation o = orientationC2_SAF(px,py,qx,qy,rx,ry,
		epsilon_1);

  return Bounded_side (s * o);
}

inline
Comparison_result
cmp_dist_to_pointC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(squared_distanceC2(px,py,qx,qy),
                       squared_distanceC2(px,py,rx,ry),
		epsilon_0);
}

inline
Comparison_result
cmp_dist_to_pointC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(squared_distanceC2(px,py,qx,qy),
                       squared_distanceC2(px,py,rx,ry),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_lineC2_SAF(
    const Static_filter_error &la,
    const Static_filter_error &lb,
    const Static_filter_error &lc,
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(scaled_distance_to_directionC2(la,lb,px,py),
                       scaled_distance_to_directionC2(la,lb,qx,qy),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_lineC2_SAF(
    const Restricted_double &la,
    const Restricted_double &lb,
    const Restricted_double &lc,
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(scaled_distance_to_directionC2(la,lb,px,py),
                       scaled_distance_to_directionC2(la,lb,qx,qy),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_lineC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &sx,
    const Static_filter_error &sy,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::compare_SAF(scaled_distance_to_lineC2(px,py,qx,qy,rx,ry),
                       scaled_distance_to_lineC2(px,py,qx,qy,sx,sy),
		epsilon_0);
}

inline
Comparison_result
cmp_signed_dist_to_lineC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &sx,
    const Restricted_double &sy,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::compare_SAF(scaled_distance_to_lineC2(px,py,qx,qy,rx,ry),
                       scaled_distance_to_lineC2(px,py,qx,qy,sx,sy),
		epsilon_0);
}

