inline
Oriented_side
in_smallest_orthogonalcircleC2_SAF(
    const Static_filter_error &px,
    const Static_filter_error &py,
    const Static_filter_error &pw,
    const Static_filter_error &qx,
    const Static_filter_error &qy,
    const Static_filter_error &qw,
    const Static_filter_error &tx,
    const Static_filter_error &ty,
    const Static_filter_error &tw,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  FT dpx = px-qx;
  FT dpy = py-qy;
  FT dtx = tx-qx;
  FT dty = ty-qy;
  FT dpz = square(dpx)+square(dpy);
 
  return Oriented_side (sign_SAF((square(dtx)+square(dty)-tw+qw)*dpz
			     -(dpz-pw+qw)*(dpx*dtx+dpy*dty),
		epsilon_0));
}

inline
Oriented_side
in_smallest_orthogonalcircleC2_SAF(
    const Restricted_double &px,
    const Restricted_double &py,
    const Restricted_double &pw,
    const Restricted_double &qx,
    const Restricted_double &qy,
    const Restricted_double &qw,
    const Restricted_double &tx,
    const Restricted_double &ty,
    const Restricted_double &tw,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  FT dpx = px-qx;
  FT dpy = py-qy;
  FT dtx = tx-qx;
  FT dty = ty-qy;
  FT dpz = square(dpx)+square(dpy);
 
  return Oriented_side (sign_SAF((square(dtx)+square(dty)-tw+qw)*dpz
			     -(dpz-pw+qw)*(dpx*dtx+dpy*dty),
		epsilon_0));
}

