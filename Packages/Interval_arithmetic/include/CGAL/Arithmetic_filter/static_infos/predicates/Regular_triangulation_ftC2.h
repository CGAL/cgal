inline
Oriented_side
power_testC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pwt,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qwt,
    const Static_filter_error &rx,
    const Static_filter_error &ry,
    const Static_filter_error &rwt,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &twt,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

    

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = square(dpx) + square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = square(dqx) + square(dqy) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = square(drx) + square(dry) - rwt + twt;

    return Oriented_side(sign_of_determinant3x3_SAF(dpx, dpy, dpz,
                                                dqx, dqy, dqz,
                                                drx, dry, drz,
		epsilon_0));
}

inline
Oriented_side
power_testC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pwt,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qwt,
    const Restricted_double &rx,
    const Restricted_double &ry,
    const Restricted_double &rwt,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &twt,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

    

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = square(dpx) + square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = square(dqx) + square(dqy) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = square(drx) + square(dry) - rwt + twt;

    return Oriented_side(sign_of_determinant3x3_SAF(dpx, dpy, dpz,
                                                dqx, dqy, dqz,
                                                drx, dry, drz,
		epsilon_0));
}

inline
Oriented_side
power_testC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pwt,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qwt,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &twt,
    double & epsilon_0,
    double & epsilon_1,
    double & epsilon_2,
    double & epsilon_3)
{
  typedef Static_filter_error FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = square(dpx) + square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = square(dqx) + square(dqy) - qwt + twt;

    
    Comparison_result cmpx = CGAL::compare_SAF(px, qx,
		epsilon_0);
    if (cmpx != EQUAL)
	return Oriented_side(cmpx * sign_of_determinant2x2_SAF(dpx, dpz, dqx, dqz,
		epsilon_1));

    
    Comparison_result cmpy = CGAL::compare_SAF(py, qy,
		epsilon_2);
    return Oriented_side(cmpy * sign_of_determinant2x2_SAF(dpy, dpz, dqy, dqz,
		epsilon_3));
}

inline
Oriented_side
power_testC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pwt,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qwt,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &twt,
    const double & epsilon_0,
    const double & epsilon_1,
    const double & epsilon_2,
    const double & epsilon_3)
{
  typedef Restricted_double FT;

    
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = square(dpx) + square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = square(dqx) + square(dqy) - qwt + twt;

    
    Comparison_result cmpx = CGAL::compare_SAF(px, qx,
		epsilon_0);
    if (cmpx != EQUAL)
	return Oriented_side(cmpx * sign_of_determinant2x2_SAF(dpx, dpz, dqx, dqz,
		epsilon_1));

    
    Comparison_result cmpy = CGAL::compare_SAF(py, qy,
		epsilon_2);
    return Oriented_side(cmpy * sign_of_determinant2x2_SAF(dpy, dpz, dqy, dqz,
		epsilon_3));
}

