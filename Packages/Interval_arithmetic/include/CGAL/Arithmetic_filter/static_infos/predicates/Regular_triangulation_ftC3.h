inline
Oriented_side
power_testC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &pwt,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &qwt,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    const Static_filter_error &rwt,
    const Static_filter_error &sx,
    const Static_filter_error &sy,
    const Static_filter_error &sz,
    const Static_filter_error &swt,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tz,
    const Static_filter_error &twt,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
    FT dsx = sx - tx;
    FT dsy = sy - ty;
    FT dsz = sz - tz;
    FT dst = square(dsx) + square(dsy) + square(dsz) - swt + twt;

    return Oriented_side( - sign_of_determinant4x4_SAF(dpx, dpy, dpz, dpt,
						   dqx, dqy, dqz, dqt,
						   drx, dry, drz, drt,
						   dsx, dsy, dsz, dst,
		epsilon_0));
}

inline
Oriented_side
power_testC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &pwt,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &qwt,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const Restricted_double &rwt,
    const Restricted_double &sx,
    const Restricted_double &sy,
    const Restricted_double &sz,
    const Restricted_double &swt,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tz,
    const Restricted_double &twt,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
    FT dsx = sx - tx;
    FT dsy = sy - ty;
    FT dsz = sz - tz;
    FT dst = square(dsx) + square(dsy) + square(dsz) - swt + twt;

    return Oriented_side( - sign_of_determinant4x4_SAF(dpx, dpy, dpz, dpt,
						   dqx, dqy, dqz, dqt,
						   drx, dry, drz, drt,
						   dsx, dsy, dsz, dst,
		epsilon_0));
}

inline
Oriented_side
power_testC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &pwt,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &qwt,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rz,
    const Static_filter_error &rwt,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tz,
    const Static_filter_error &twt,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2,
    double & epsilon_3,
    double & epsilon_4,
    double & epsilon_5)
{
  typedef Static_filter_error FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
    Sign cmp;

    
    cmp = sign_of_determinant3x3_SAF(dpx, dpy, dpt,
		                 dqx, dqy, dqt,
				 drx, dry, drt,
		epsilon_0);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2_SAF(px-rx, py-ry,
		                                          qx-rx, qy-ry,
		epsilon_1));

    
    cmp = sign_of_determinant3x3_SAF(dpx, dpz, dpt,
		                 dqx, dqz, dqt,
				 drx, drz, drt,
		epsilon_2);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2_SAF(px-rx, pz-rz,
		                                          qx-rx, qz-rz,
		epsilon_3));

    
    cmp = sign_of_determinant3x3_SAF(dpy, dpz, dpt,
		                 dqy, dqz, dqt,
				 dry, drz, drt,
		epsilon_4);
    return Oriented_side(cmp * sign_of_determinant2x2_SAF(py-ry, pz-rz,
		                                      qy-ry, qz-rz,
		epsilon_5));
}

inline
Oriented_side
power_testC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &pwt,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &qwt,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rz,
    const Restricted_double &rwt,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tz,
    const Restricted_double &twt,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2,
    const double & epsilon_3,
    const double & epsilon_4,
    const double & epsilon_5)
{
  typedef Restricted_double FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
    Sign cmp;

    
    cmp = sign_of_determinant3x3_SAF(dpx, dpy, dpt,
		                 dqx, dqy, dqt,
				 drx, dry, drt,
		epsilon_0);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2_SAF(px-rx, py-ry,
		                                          qx-rx, qy-ry,
		epsilon_1));

    
    cmp = sign_of_determinant3x3_SAF(dpx, dpz, dpt,
		                 dqx, dqz, dqt,
				 drx, drz, drt,
		epsilon_2);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2_SAF(px-rx, pz-rz,
		                                          qx-rx, qz-rz,
		epsilon_3));

    
    cmp = sign_of_determinant3x3_SAF(dpy, dpz, dpt,
		                 dqy, dqz, dqt,
				 dry, drz, drt,
		epsilon_4);
    return Oriented_side(cmp * sign_of_determinant2x2_SAF(py-ry, pz-rz,
		                                      qy-ry, qz-rz,
		epsilon_5));
}

inline
Oriented_side
power_testC3_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pz,
    const Static_filter_error &pwt,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qz,
    const Static_filter_error &qwt,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tz,
    const Static_filter_error &twt,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2,
    double & epsilon_3,
    double & epsilon_4,
    double & epsilon_5)
{
  typedef Static_filter_error FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    Comparison_result cmp;

    
    cmp = CGAL::compare_SAF(px, qx,
		epsilon_0);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpx, dpt, dqx, dqt,
		epsilon_1));

    
    cmp = CGAL::compare_SAF(py, qy,
		epsilon_2);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpy, dpt, dqy, dqt,
		epsilon_3));

    
    cmp = CGAL::compare_SAF(pz, qz,
		epsilon_4);
    return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpz, dpt, dqz, dqt,
		epsilon_5));
}

inline
Oriented_side
power_testC3_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pz,
    const Restricted_double &pwt,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qz,
    const Restricted_double &qwt,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tz,
    const Restricted_double &twt,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2,
    const double & epsilon_3,
    const double & epsilon_4,
    const double & epsilon_5)
{
  typedef Restricted_double FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
    Comparison_result cmp;

    
    cmp = CGAL::compare_SAF(px, qx,
		epsilon_0);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpx, dpt, dqx, dqt,
		epsilon_1));

    
    cmp = CGAL::compare_SAF(py, qy,
		epsilon_2);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpy, dpt, dqy, dqt,
		epsilon_3));

    
    cmp = CGAL::compare_SAF(pz, qz,
		epsilon_4);
    return Oriented_side(cmp * sign_of_determinant2x2_SAF(dpz, dpt, dqz, dqt,
		epsilon_5));
}

