#ifndef CHR_MOEBIUS_DIAGRAM_REGULAR_FTC2_H
#define CHR_MOEBIUS_DIAGRAM_REGULAR_FTC2_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
void
inline
moebius_to_regularC2 (const FT &x, const FT &y, const FT &l, const FT &m,
		     FT &x_, FT &y_, FT &z_, FT &w_)
{ 
  x_ = l * x;
  y_ = l * y;
  z_ = l / FT (-2);
  FT l2 = l * l;
  FT p2 = x * x + y * y;
  w_ = l2 * p2 + l2 / FT (4) - l * p2 + m;
  //  TRACE ("      ("<<x<<","<<y<<", "<<l<<","<<m<<") -> ("<<x_<<","<<y_<<","<<z_<<", "<<w_<<")\n");
}

template <class RT>
Orientation
moebius_orientationC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		      const RT &qx, const RT &qy, const RT &ql, const RT &qm,
		      const RT &rx, const RT &ry, const RT &rl, const RT &rm,
		      const RT &sx, const RT &sy, const RT &sl, const RT &sm)
{
  const RT plx = pl * px, ply = pl * py;

  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const RT sx_ = sl * sx - plx, sy_ = sl * sy - ply, sl_ = pl - sl;

  return Orientation (sign_of_determinant3x3 (qx_,qy_,ql_,
					      rx_,ry_,rl_,
					      sx_,sy_,sl_));
}

template <class RT>
Orientation
moebius_coplanar_orientationC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
			       const RT &qx, const RT &qy, const RT &ql, const RT &qm,
			       const RT &rx, const RT &ry, const RT &rl, const RT &rm)
{
  const RT plx = pl * px, ply = pl * py;

  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;

  const Orientation oxy_pqr = Orientation (sign_of_determinant2x2 (qx_,qy_, rx_,ry_));
  if (oxy_pqr != ZERO) return oxy_pqr;

  const Orientation oyz_pqr = Orientation (sign_of_determinant2x2 (qy_,ql_, ry_,rl_));
  if (oyz_pqr != ZERO) return oyz_pqr;

  return Orientation (sign_of_determinant2x2 (qx_,ql_, rx_,rl_));
}

template <class RT>
Comparison_result
moebius_compare_xyzC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		      const RT &qx, const RT &qy, const RT &ql, const RT &qm)
{
  Comparison_result c = CGAL_NTS compare (pl * px, ql * qx);
  if (c != EQUAL) return c;

  c = CGAL_NTS compare (pl * py, ql * qy);
  if (c != EQUAL) return c;

  return CGAL_NTS compare (ql, pl);
}



template <class RT>
Oriented_side
moebius_power_testC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		     const RT &qx, const RT &qy, const RT &ql, const RT &qm,
		     const RT &rx, const RT &ry, const RT &rl, const RT &rm,
		     const RT &sx, const RT &sy, const RT &sl, const RT &sm,
		     const RT &tx, const RT &ty, const RT &tl, const RT &tm)
{
  static const RT _2 = 2, _4 = 4;
  const RT tlx = tl * tx, tly = tl * ty;

  const RT px_ = pl * px - tlx, py_ = pl * py - tly, pl_ = tl - pl;
  const RT qx_ = ql * qx - tlx, qy_ = ql * qy - tly, ql_ = tl - ql;
  const RT rx_ = rl * rx - tlx, ry_ = rl * ry - tly, rl_ = tl - rl;
  const RT sx_ = sl * sx - tlx, sy_ = sl * sy - tly, sl_ = tl - sl;

  const RT tl_t2 = tl * (tx * tx + ty * ty), T = tl_t2 - tm;
  const RT pw_ = (tl * (_4 * (tl_t2 - pl * (px * tx + py * ty)) + tl - pl) +
		  _2 * (pl * (px * px + py * py) - pm - T));
  const RT qw_ = (tl * (_4 * (tl_t2 - ql * (qx * tx + qy * ty)) + tl - ql) +
		  _2 * (ql * (qx * qx + qy * qy) - qm - T));
  const RT rw_ = (tl * (_4 * (tl_t2 - rl * (rx * tx + ry * ty)) + tl - rl) +
		  _2 * (rl * (rx * rx + ry * ry) - rm - T));
  const RT sw_ = (tl * (_4 * (tl_t2 - sl * (sx * tx + sy * ty)) + tl - sl) +
		  _2 * (sl * (sx * sx + sy * sy) - sm - T));

  return Oriented_side (- sign_of_determinant4x4 (px_, py_, pl_, pw_,
						  qx_, qy_, ql_, qw_,
						  rx_, ry_, rl_, rw_,
						  sx_, sy_, sl_, sw_));
}


template <class RT>
Oriented_side
moebius_power_testC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		     const RT &qx, const RT &qy, const RT &ql, const RT &qm,
		     const RT &rx, const RT &ry, const RT &rl, const RT &rm,
		     const RT &tx, const RT &ty, const RT &tl, const RT &tm)
{
  static const RT _2 = 2, _4 = 4;
  const RT plx = pl * px, ply = pl * py;
  const RT qlx = ql * qx, qly = ql * qy;
  const RT rlx = rl * rx, rly = rl * ry;
  const RT tlx = tl * tx, tly = tl * ty;

  const RT px_ = plx - tlx, py_ = ply - tly, pl_ = tl - pl;
  const RT qx_ = qlx - tlx, qy_ = qly - tly, ql_ = tl - ql;
  const RT rx_ = rlx - tlx, ry_ = rly - tly, rl_ = tl - rl;

  const RT tl_t2 = tl * (tx * tx + ty * ty), T = tl_t2 - tm;
  const RT pw_ = (tl * (_4 * (tl_t2 - pl * (px * tx + py * ty)) + tl - pl) +
		  _2 * (pl * (px * px + py * py) - pm - T));
  const RT qw_ = (tl * (_4 * (tl_t2 - ql * (qx * tx + qy * ty)) + tl - ql) +
		  _2 * (ql * (qx * qx + qy * qy) - qm - T));
  const RT rw_ = (tl * (_4 * (tl_t2 - rl * (rx * tx + ry * ty)) + tl - rl) +
		  _2 * (rl * (rx * rx + ry * ry) - rm - T));

  Sign cmp;
  cmp = sign_of_determinant3x3 (px_, py_, pw_,
				qx_, qy_, qw_,
				rx_, ry_, rw_);
  if (cmp != ZERO)
    return Oriented_side (cmp * sign_of_determinant2x2 (plx-rlx, ply-rly,
							qlx-rlx, qly-rly));

  cmp = sign_of_determinant3x3 (px_, pl_, pw_,
				qx_, ql_, qw_,
				rx_, rl_, rw_);
  if (cmp != ZERO)
    return Oriented_side (cmp * sign_of_determinant2x2 (plx-rlx, rl-pl,
							qlx-rlx, rl-ql));
  
  cmp = sign_of_determinant3x3 (py_, pl_, pw_,
				qy_, ql_, qw_,
				ry_, rl_, rw_);
  return Oriented_side (cmp * sign_of_determinant2x2 (ply-rly, rl-pl,
						      qly-rly, rl-ql));
}


template <class RT>
Oriented_side
moebius_power_testC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		     const RT &qx, const RT &qy, const RT &ql, const RT &qm,
		     const RT &tx, const RT &ty, const RT &tl, const RT &tm)
{
  static const RT _2 = 2, _4 = 4;
  const RT plx = pl * px, ply = pl * py;
  const RT qlx = ql * qx, qly = ql * qy;
  const RT tlx = tl * tx, tly = tl * ty;

  const RT px_ = plx - tlx, py_ = ply - tly, pl_ = tl - pl;
  const RT qx_ = qlx - tlx, qy_ = qly - tly, ql_ = tl - ql;

  const RT tl_t2 = tl * (tx * tx + ty * ty), T = tl_t2 - tm;
  const RT pw_ = (tl * (_4 * (tl_t2 - pl * (px * tx + py * ty)) + tl - pl) +
		  _2 * (pl * (px * px + py * py) - pm - T));
  const RT qw_ = (tl * (_4 * (tl_t2 - ql * (qx * tx + qy * ty)) + tl - ql) +
		  _2 * (ql * (qx * qx + qy * qy) - qm - T));

  Comparison_result cmp;

  cmp = CGAL_NTS compare (plx, qlx);
  if (cmp != EQUAL)
    return Oriented_side (cmp * sign_of_determinant2x2 (px_, pw_, qx_, qw_));

  cmp = CGAL_NTS compare (ply, qly);
  if (cmp != EQUAL)
    return Oriented_side (cmp * sign_of_determinant2x2 (py_, pw_, qy_, qw_));

  cmp = CGAL_NTS compare (ql, pl);
  return Oriented_side (cmp * sign_of_determinant2x2 (pl_, pw_, ql_, qw_));
}

template <class RT>
Oriented_side
moebius_power_testC3 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		     const RT &qx, const RT &qy, const RT &ql, const RT &qm)
{
  static const RT _4 = 4;
  const RT ql2 = ql * ql, pl2 = pl * pl;
  const RT q2 = qx * qx + qy * qy, p2 = px * px + py * py;
  return Oriented_side (CGAL_NTS compare (ql2 - _4 * ((ql - ql2) * q2 - qm),
					  pl2 - _4 * ((pl - pl2) * p2 - pm)));
}


CGAL_END_NAMESPACE

#endif//CHR_MOEBIUS_DIAGRAM_REGULAR_FTC2_H
