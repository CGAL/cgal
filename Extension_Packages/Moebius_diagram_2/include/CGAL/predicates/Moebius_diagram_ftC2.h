#ifndef CHR_MOEBIUS_DIAGRAM_FTC2_H
#define CHR_MOEBIUS_DIAGRAM_FTC2_H
CGAL_BEGIN_NAMESPACE


template <class RT>
bool
compare_on_oriented_lineC2 (const RT &a, const RT &b,
			    const RT &x1, const RT &y1, 
			    const RT &x2, const RT &y2)
{
  const RT _0 = 0;
  switch (CGAL_NTS sign (b)) {
  case CGAL::ZERO: break;
  case CGAL::POSITIVE: if (x1 > x2) return false; break;
  case CGAL::NEGATIVE: if (x1 < x2) return false; break;
  }
  switch (CGAL_NTS sign (a)) {
  case CGAL::ZERO: break;
  case CGAL::POSITIVE: if (y1 < y2) return false; break;
  case CGAL::NEGATIVE: if (y1 > y2) return false; break;
  }
  return true;

  return 
    (b == _0 || (b > _0 && x1 >= x2) || x1 <= x2) &&
    (a == _0 || (a > _0 && y1 <= y2) || y1 >= y2);
}


template <class RT>
CGAL::Orientation
moebius_orientationC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		      const RT &qx, const RT &qy, const RT &ql, const RT &qm)
{
  return CGAL::Orientation (CGAL_NTS compare (pl, ql));
}

template <class RT>
bool
moebius_has_circleC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
		     const RT &qx, const RT &qy, const RT &ql, const RT &qm)
{
  const RT qpx = qx - px, qpy = qy - py;
  return pl * ql * (qpx * qpx + qpy * qpy) > (pl - ql) * (qm - pm);
}

template <class RT>
bool
moebius_circle_cross_lineC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
			    const RT &qx, const RT &qy, const RT &ql, const RT &qm,
			    const RT &rx, const RT &ry, const RT &rl, const RT &rm)
{
  static const RT _2 = 2, _4 = 4;
  
  const RT plx = pl * px, ply = pl * py;
  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  
  const RT qpx = qx - px, qpy = qy - py;
  const RT r_pq = pl * ql * (qpx * qpx + qpy * qpy) - ql_ * (qm - pm);

  const RT P = pl * (px * px + py * py) - pm;
  const RT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const RT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const RT a_ = det2x2_by_formula (ql_, qx_, rl_, rx_);
  const RT b_ = det2x2_by_formula (ql_, qy_, rl_, ry_);
  const RT c_ = det2x2_by_formula (Q_, ql_, R_, rl_);

  return CGAL_NTS square (_2 * (a_ * qx_ + b_ * qy_) - ql_ * c_) < _4 * r_pq * (a_ * a_ + b_ * b_);
}


template <class RT>
CGAL::Oriented_side
moebius_circle_side_of_centerC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
				const RT &qx, const RT &qy, const RT &ql, const RT &qm,
				const RT &rx, const RT &ry, const RT &rl, const RT &rm)
{
  static const RT _2 = 2;
  
  const RT plx = pl * px, ply = pl * py;
  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  
  const RT P = pl * (px * px + py * py) - pm;
  const RT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const RT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const RT a_ = det2x2_by_formula (ql_, qx_, rl_, rx_);
  const RT b_ = det2x2_by_formula (ql_, qy_, rl_, ry_);
  const RT c_ = det2x2_by_formula (Q_, ql_, R_, rl_);

  return CGAL::Oriented_side (CGAL_NTS compare (_2 * (a_ * qx_ + b_ * qy_), ql_ * c_));
}

template <class RT>
CGAL::Bounded_side
moebius_side_of_vertexC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
			 const RT &qx, const RT &qy, const RT &ql, const RT &qm,
			 const RT &rx, const RT &ry, const RT &rl, const RT &rm,
			 const RT &sx, const RT &sy, const RT &sl, const RT &sm)
{  
  static const RT _4 = 4;

  const RT plx = pl * px, ply = pl * py;

  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const RT sx_ = sl * sx - plx, sy_ = sl * sy - ply, sl_ = pl - sl;

  const RT P = pl * (px * px + py * py) - pm;
  const RT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const RT R_ = rl * (rx * rx + ry * ry) - rm - P;
  const RT S_ = sl * (sx * sx + sy * sy) - sm - P;

  const RT Dx = det3x3_by_formula (Q_,qy_,ql_, R_,ry_,rl_, S_,sy_,sl_);
  const RT Dy = det3x3_by_formula (qx_,Q_,ql_, rx_,R_,rl_, sx_,S_,sl_);
  const RT Dz = det3x3_by_formula (qx_,qy_,Q_, rx_,ry_,R_, sx_,sy_,S_);
  const RT D = det3x3_by_formula (qx_,qy_,ql_, rx_,ry_,rl_, sx_,sy_,sl_);

  CGAL_assertion (CGAL_NTS sign (D) == CGAL::POSITIVE);

  CGAL::Bounded_side result =
    CGAL::Bounded_side (CGAL_NTS compare (_4 * D * Dz, Dx * Dx + Dy * Dy));

  return result == CGAL::ON_BOUNDARY ? CGAL::ON_UNBOUNDED_SIDE : result;
}

template <class RT>
CGAL::Oriented_side
moebius_circle_side_of_vertexC2 (const RT &px, const RT &py, const RT &pl, const RT &pm,
				const RT &qx, const RT &qy, const RT &ql, const RT &qm,
				const RT &rx, const RT &ry, const RT &rl, const RT &rm,
				const RT &sx, const RT &sy, const RT &sl, const RT &sm)
{
  static const RT _2 = 2;
  const RT plx = pl * px, ply = pl * py;

  const RT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const RT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const RT sx_ = sl * sx - plx, sy_ = sl * sy - ply, sl_ = pl - sl;

  const RT P = pl * (px * px + py * py) - pm;
  const RT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const RT R_ = rl * (rx * rx + ry * ry) - rm - P;
  const RT S_ = sl * (sx * sx + sy * sy) - sm - P;

  const RT Dx = det3x3_by_formula (Q_,qy_,ql_, R_,ry_,rl_, S_,sy_,sl_);
  const RT Dy = det3x3_by_formula (qx_,Q_,ql_, rx_,R_,rl_, sx_,S_,sl_);
  //  const RT Dz = det3x3_by_formula (qx_,qy_,Q_, rx_,ry_,R_, sx_,sy_,S_);
  const RT D = det3x3_by_formula (qx_,qy_,ql_, rx_,ry_,rl_, sx_,sy_,sl_);

  CGAL_assertion (CGAL_NTS sign (D) != CGAL::ZERO);

  const RT a_ = det2x2_by_formula (qx_, ql_, rx_, rl_);
  const RT b_ = det2x2_by_formula (qy_, ql_, ry_, rl_);
  //  const RT c_ = det2x2_by_formula (Q_, ql_, R_, rl_);

  return CGAL::Oriented_side (CGAL_NTS compare (ql_ * (a_ * Dy - b_ * Dx),
						_2 * D * (qx_ * b_ - qy_ * a_)));
}

CGAL_END_NAMESPACE
#endif
