#ifndef CHR_MOEBIUS_CONSTRUCTION_DIAGRAM_FTC2_H
#define CHR_MOEBIUS_CONSTRUCTION_DIAGRAM_FTC2_H

#include <CGAL/Moebius_utils.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
void
moebius_regular_circumcenterC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			       const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			       const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			       const FT &sx, const FT &sy, const FT &sl, const FT &sm,
			       FT &x, FT &y)
{
  static const FT _two = 2;

/*   std::cout << "regular_center "  */
/* 	    <<"("<<px<<", "<<py<<", "<<pl<<", "<<pm<<") " */
/* 	    <<"("<<qx<<", "<<qy<<", "<<ql<<", "<<qm<<") " */
/* 	    <<"("<<rx<<", "<<ry<<", "<<rl<<", "<<rm<<") " */
/* 	    <<"("<<sx<<", "<<sy<<", "<<sl<<", "<<sm<<")\n"; */

  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const FT sx_ = sl * sx - plx, sy_ = sl * sy - ply, sl_ = pl - sl;

  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;
  const FT S_ = sl * (sx * sx + sy * sy) - sm - P;


  const FT Rx = det3x3_by_formula (Q_,qy_,ql_, R_,ry_,rl_, S_,sy_,sl_);
  const FT Ry = det3x3_by_formula (qx_,Q_,ql_, rx_,R_,rl_, sx_,S_,sl_);
  const FT R = _two * det3x3_by_formula (qx_,qy_,ql_, rx_,ry_,rl_, sx_,sy_,sl_);

  CGAL_assertion (CGAL_NTS sign (R) != CGAL::ZERO);

  x = Rx / R;
  y = Ry / R;
  //  std::cout << "               ("<<x<<", "<<y<<")\n";
}

template <class FT>
void
moebius_lineC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
	       const FT &qx, const FT &qy, const FT &ql, const FT &qm,
	       FT &a, FT &b, FT &c)
{
  static const FT _2 = 2;
  CGAL_precondition (pl == ql);
  
  const FT p = px * px + py * py;
  const FT q = qx * qx + qy * qy;
  a = _2 * pl * (qx - px);
  b = _2 * pl * (qy - py);
  c = pl * (p - q) - pm + qm;
}

template <class FT>
CGAL::Orientation
moebius_circleC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
		 const FT &qx, const FT &qy, const FT &ql, const FT &qm,
		 FT &cx, FT &cy, FT &r)
{
  CGAL_precondition (pl != ql);

  const FT l_ = ql - pl;

  cx = (ql * qx - pl * px) / l_;
  cy = (ql * qy - pl * py) / l_;

  const FT qpx = qx - px;
  const FT qpy = qy - py;

  const FT r_1 = pl * ql * (qpx*qpx + qpy*qpy);
  const FT r_2 = l_ * (pm - qm);

  r = (r_1 - r_2) / (l_ * l_);
  return CGAL::Orientation (CGAL_NTS sign (l_));
}

#if 1

template <class FT>
void
moebius_vertexC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
		 const FT &qx, const FT &qy, const FT &ql, const FT &qm,
		 const FT &rx, const FT &ry, const FT &rl, const FT &rm,
		 FT &x, FT &y)
{
  static const FT _2 = 2;
  CGAL_assertion (pl == ql);
  CGAL_assertion (pl == rl);
  
  const FT qx_ = qx - px, qy_ = qy - py;
  const FT rx_ = rx - px, ry_ = ry - py;

  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;
  
  const FT dx = det2x2_by_formula (Q_, qy_, R_, ry_);
  const FT dy = det2x2_by_formula (qx_, Q_, rx_, R_);
  const FT d = _2 * pl * det2x2_by_formula (qx_, qy_, rx_, ry_);

  x = dx / d;
  y = dy / d;
}

// solves 4 * a * x^2 - 4 * b * x + c = 0
// x1 = b + sqrt()... x2 = b - sqrt()...
template <class FT>
void
moebius_solveC2 (const FT &a, const FT &b, const FT &c,
		FT &x1, FT &x2)
{
  static const FT _2 = 2;
  FT delta = b * b - a * c;

  // we are assuming that there is solutions
  if (delta < 0) {
    // hoping delta is not much negative...
    delta = 0;
  }

  const FT delta_ = CGAL_NTS sqrt (delta);

  x1 = (b - delta_) / (_2 * a);
  x2 = (b + delta_) / (_2 * a);
}


// actually solving { l*(x^2+y2) + 2*x0*x + 2*y0*y - r = 0
//                  { 2*a*x + 2*b*y + c = 0  
template <class FT>
void
moebius_compute_verticesC2 (const FT &l, const FT &x0, const FT &y0, const FT &r,
			   const FT &a, const FT &b, const FT &c,
			   FT &x1, FT &y1, FT &x2, FT &y2)
{
  static const FT _2 = 2, _4 = 4;

  

  const FT a2 = a * a;
  const FT b2 = b * b;
  const FT a_ = l * (a2 + b2);

  TRACE ("      compute, 2*"<<a<<"*x + 2*"<<b<<"*y + "<<c<<" = 0\n"
	 << "               "<<l<<"*(x^2+y^2) + 2*"<<x0<<"*x + +2*"<<y0<<"*y - "<<r<<" = 0\n");

  if (CGAL_NTS abs (a) < CGAL_NTS abs (b)) {
    TRACE ("      |a|<|b|\n");
    const FT b_ = _2 * b * (a * y0 - b * x0) - l * a * c;
    const FT c_ = l * c * c - _4 * b * (c * y0 + b * r);

    if (CGAL_NTS sign (b) == CGAL::POSITIVE) {
      moebius_solveC2 (a_, b_, c_, x1, x2);
    } else {
      moebius_solveC2 (a_, b_, c_, x2, x1);
    }
    y1 = (- c - _2 * a * x1) / (_2 * b);
    y2 = (- c - _2 * a * x2) / (_2 * b);

  } else {
    TRACE ("      |a|>=|b|\n");
    const FT b_ = _2 * a * (b * x0 - a * y0) - l * b * c;
    const FT c_ = l * c * c - _4 * a * (c * x0 + a * r);

    if (CGAL_NTS sign (a) == CGAL::POSITIVE) {
      moebius_solveC2 (a_, b_, c_, y2, y1);
    } else {
      moebius_solveC2 (a_, b_, c_, y1, y2);
    }
    x1 = (- c - _2 * b * y1) / (_2 * a);
    x2 = (- c - _2 * b * y2) / (_2 * a);
  }
}

template <class FT>
void
moebius_line_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			FT &x1, FT &y1, FT &x2, FT &y2)
{
  CGAL_assertion (pl == ql);
  
  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;

  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  moebius_compute_verticesC2 (rl_, rx_, ry_, R_,
			     qx_, qy_, - Q_,
			     x1, y1, x2, y2);
}

template <class FT>
void
moebius_circle_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			  const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			  const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			  FT &x1, FT &y1, FT &x2, FT &y2)
{
  CGAL_assertion (pl != ql);
  
  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;

  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const FT a = det2x2_by_formula (qx_, ql_, rx_, rl_);
  const FT b = det2x2_by_formula (qy_, ql_, ry_, rl_);
  const FT c = det2x2_by_formula (ql_, Q_, rl_, R_);

  if (pl > ql)
    moebius_compute_verticesC2 (ql_, qx_, qy_, Q_,
			       a, b, c,
			       x1, y1, x2, y2);
  else
    moebius_compute_verticesC2 (ql_, qx_, qy_, Q_,
			       a, b, c,
			       x2, y2, x1, y1);
}

#else
// solve a * x^2 - 2 * b * x + c = 0
template <class FT>
bool
moebius_solve (const FT &a, const FT &b, const FT &c,
	      FT &x1, FT &x2)
{
  //  static const FT _2 = 2, _4 = 4;
  // std::cout << "solve (" << a << "*x^2 -2*" << b << "*x + " << c << ")\n";

  const FT delta = b * b - a * c;
  // std::cout << "  delta = " << delta << "\n";
  if (delta < 0) return true;

  //  std::cout << " delta = " << delta << "\n";

  //  CGAL_assertion (delta >= 0);

  const FT delta_ = CGAL_NTS sqrt (delta);
  x1 = (b + delta_) / a;
  x2 = (b - delta_) / a;
  // std::cout << " " << x1 << ", " << x2 << "\n";
  return false;
}


template <class FT>
int
moebius_circle_2_line_2_intersectionC2 (const FT &cx, const FT &cy, const FT &r,
				       const FT &a, const FT &b, const FT &c,
				       FT &x1, FT &y1, FT &x2, FT &y2)
{
  static const FT _2 = 2;

  CGAL_precondition (CGAL_NTS sign (a) != CGAL::ZERO ||
		     CGAL_NTS sign (b) != CGAL::ZERO);

  // std::cout << "inter_clC2 (("<<cx<<", "<<cy<<"), "<<r<<") ("<<a<<", "<<b<<", "<<c<<")\n";

  const FT a2 = a * a;
  const FT b2 = b * b;
  const FT a_ = a2 + b2;

  if (CGAL_NTS abs (a) < CGAL_NTS abs (b)) {
    std::cout << "inter_cl: |a|<|b|\n";
    const FT b_ = b * (b * cx - a * cy) - a * c;
    const FT c_ = b2 * (cx * cx + cy * cy - r) + _2 * b * c * cy + c * c;

    bool no_solution = (b < 0 ?
			moebius_solve (a_, b_, c_, x1, x2) :
			moebius_solve (a_, b_, c_, x2, x1));
    if (no_solution) return 0;

    y1 = (- c - a * x1) / b; 
    y2 = (- c - a * x2) / b; 
  } else {
    std::cout << "inter_cl: |a|>=|b|\n";
    const FT b_ = a * (a * cy - b * cx) - b * c;
    const FT c_ = a2 * (cx * cx + cy * cy - r) + _2 * a * c * cx + c * c;
    bool no_solution = (a > 0 ?
			moebius_solve (a_, b_, c_, y1, y2) :
			moebius_solve (a_, b_, c_, y2, y1));
    if (no_solution) return 0;

    x1 = (- c - b * y1) / a;
    x2 = (- c - b * y2) / a;
  }
  // std::cout << "inter_cl2C2 = ("<<x1<<", "<<y1<<") ("<<x2<<", "<<y2<<")\n";
  return 2;
}

template <class FT>
int
moebius_circle_2_circle_2_intersectionC2 (const FT &cx1, const FT &cy1, const FT &r1,
					 const FT &cx2, const FT &cy2, const FT &r2,
					 FT &x1, FT &y1, FT &x2, FT &y2)
{
  static const FT _2 = 2;
  const FT a = _2 * (cx2 - cx1);
  const FT b = _2 * (cy2 - cy1);
  const FT c = cx1*cx1 + cy1*cy1 - r1 - (cx2*cx2 + cy2*cy2 - r2);
  std::cout << "inter_cc -> line " << a << " " << b << " " << c << "\n";
  return moebius_circle_2_line_2_intersectionC2 (cx1, cy1, r1, a, b, c, x1, y1, x2, y2);
}


template <class FT>
void
moebius_solve (const FT &a, const FT &b, const FT &c,
	      FT &x1, FT &x2)
{
  static const FT _2 = 2;

  //  std::cout << "solve (" << a << ", " << b << ", " << c << ")\n";
  const FT delta = b * b - a * c;
  //  std::cout << " delta = " << delta << "\n";
  CGAL_assertion (delta >= 0);
  const FT delta_ = CGAL_NTS sqrt (delta);
  x1 = (b + delta_) / (_2 * a);
  x2 = (b - delta_) / (_2 * a);
  //  std::cout << x1 << ", " << x2 << "\n";
}

template <class FT>
void
moebius_quadratic (const FT &qx_, const FT &qy_, const FT &ql_, const FT &Q_,
		  const FT &a, const FT &b, const FT &c,
		  FT &x1, FT &y1, FT &x2, FT &y2)
{
  static const FT _2 (2), __2 (-2), _4 (4);

  if (CGAL_NTS sign (b) == CGAL::ZERO) {
    CGAL_assertion (CGAL_NTS sign (a) != CGAL::ZERO);
    
    const FT a_ = ql_ * a * a;
    const FT b_ = __2 * qy_ * a * a;
    const FT c_ = _4 * a * (c * qx_ - a * Q_) + ql_ * c * c;
 
    x1 = x2 = c / (_2 * a);
    moebius_solve (a_, b_, c_, y1, y2);
    
  } else if (CGAL_NTS sign (a) == CGAL::ZERO) {
    const FT a_ = ql_ * b * b;
    const FT b_ = __2 * qx_ * b * b;
    const FT c_ = _4 * b * (c * qy_ - b * Q_) + ql_ * c * c;
 
    x1 = x2 = c / (_2 * b);
    moebius_solve (a_, b_, c_, y1, y2);
    
  } else if (CGAL_NTS abs (a) < CGAL_NTS abs (b)) {
    const FT a_ = ql_ * (a * a + b * b);
    const FT b_ = _2 * b * (b * qx_ + a * qy_) - ql_ * a * c;
    const FT c_ = _4 * b * (c * qy_ - a * Q_) + ql_ * c * c;

    moebius_solve (a_, b_, c_, x1, x2);
    y1 = (c - _2 * a * x1) / (_2 * b);
    y2 = (c - _2 * a * x2) / (_2 * b);

  } else {
    const FT a_ = ql_ * (a * a + b * b);
    const FT b_ = _2 * a * (a * qy_ + b * qx_) - ql_ * b * c;
    const FT c_ = _4 * a * (c * qx_ - a * Q_) + ql_ * c * c;

    moebius_solve (a_, b_, c_, y1, y2);
    x1 = (c - _2 * b * y1) / (_2 * a);
    x2 = (c - _2 * b * y2) / (_2 * a);

  }
}

template <class FT>
void
moebius_circle_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			  const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			  const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			  FT &x1, FT &y1,
			  FT &x2, FT &y2)
{
  CGAL_precondition (pl != ql);

  if (pl == rl) {
    //    std::cout << "pif\n";
    return moebius_parabola_verticesC2 (px,py,pl,pm, rx,ry,rl,rm, qx,qy,ql,qm, x1,y1, x2,y2);
  }

  if (ql == rl) {
    //    std::cout << "paf\n";
    return moebius_parabola_verticesC2 (qx,qy,ql,qm, rx,ry,rl,rm, px,py,pl,pm, x1,y1, x2,y2);
  }

  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const FT a = ql_ * rx_ - rl_ * qx_;
  const FT b = ql_ * ry_ - rl_ * qy_;
  const FT c = ql_ * R_ - rl_ * Q_;

  moebius_quadratic (qx_, qy_, ql_, Q_, a, b, c,
		    x1, y1, x2, y2);
}


template <class FT>
void
moebius_parabola_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			    const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			    const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			    FT &x1, FT &y1,
			    FT &x2, FT &y2)
{
  CGAL_precondition (pl == ql);
  CGAL_precondition (pl != rl);

  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const FT a = qx_;
  const FT b = qy_;
  const FT c = Q_;

  moebius_quadratic (rx_, ry_, rl_, R_, a, b, c,
		    x1, y1, x2, y2);
}

template <class FT>
void
moebius_line_vertexC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
		      const FT &qx, const FT &qy, const FT &ql, const FT &qm,
		      const FT &rx, const FT &ry, const FT &rl, const FT &rm,
		      FT &x, FT &y)
{
  static const FT _2 = 2;

  CGAL_assertion (pl == ql);
  CGAL_assertion (pl == rl);

  const FT qx_ = qx - px, rx_ = rx - px;
  const FT qy_ = qy - py, ry_ = ry - py;
  const FT qm_ = qm - pm, rm_ = rm - pm;

  const FT den = _2 * pl * det2x2_by_formula (qx_, qy_, rx_, ry_);

  x = det2x2_by_formula (qm_, qy_, rm_, ry_) / den;
  y = det2x2_by_formula (qx_, qm_, rx_, rm_) / den;
}
#endif //0
#if 0

// (x1,y1) sees the regions of p, q and r in ccw order
// (x2,y2) sees the regions of p, q and r in cw order
template <class FT>
void
moebius_parabola_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			    const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			    const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			    FT &x1, FT &y1,
			    FT &x2, FT &y2)
{
  static const FT _2 (2), _4 (4);

  CGAL_precondition (pl == ql);
  CGAL_precondition (pl != rl);

  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  const FT a = qx_;
  const FT b = qy_;
  const FT c = Q_;

  const Sign sign_a = CGAL_NTS sign (a);
  const Sign sign_b = CGAL_NTS sign (b);
  const FT cd = rl_;
  if (sign_b != ZERO) {
    // -8*b*a*yc+8*b^2*xc+4*dc*a*c
    const FT tmp1 = a * ry_ - b * rx_;
    const FT X = a * c * cd - _2 * b * tmp1;
    
    const FT tmp2 = c * cd;
    const FT a2b2 = a*a + b*b;
    const FT delta = _4 * (tmp1*tmp1 + (a*rx_ + b*ry_)*c*cd + R_ * a2b2) - tmp2 * tmp2;
    const FT delta_ = b * CGAL_NTS sqrt (delta);

    const FT D = _2 * a2b2;
    
    const FT _2a = _2 * a;
    const FT _2b = _2 * b;

   if ((CGAL::Sign) (((int) sign_b) * (int) CGAL_NTS sign (rl_)) == CGAL::POSITIVE) {
      x1 = (X - delta_) / D;
      x2 = (X + delta_) / D;
    } else {
      x1 = (X + delta_) / D;
      x2 = (X - delta_) / D;
    }

    y1 = (c - _2a * x1) / _2b;
    y2 = (c - _2a * x2) / _2b;

  } else {
    CGAL_assertion (CGAL_NTS sign (a) != ZERO);
    
    // 4*dc*b*c-8*a*b*xc+8*a^2*yc
    //const FT X = b * c * cd - _2 * a * (b * cx - a * cy);
    const FT X = _2 * a * a * ry_;
    
    // 
    const FT tmp1 = a * ry_;
    const FT tmp2 = c * cd;
    const FT a2b2 = a*a;
    const FT delta = _4 * (tmp1*tmp1 + (a*rx_)*c*cd + R_ * a2b2) - tmp2 * tmp2;
    const FT delta_ = a * CGAL_NTS sqrt (delta);
    
    const FT D = _2 * a2b2;
    
    const FT _2a = _2 * a;
    
    if ((CGAL::Sign) (((int) sign_a) * (int) CGAL_NTS sign (rl_)) == CGAL::POSITIVE) {
      y1 = (X + delta_) / D;
      y2 = (X - delta_) / D;
    } else {
      y1 = (X - delta_) / D;
      y2 = (X + delta_) / D;
    }

    x1 = x2 = c / _2a;
  }
}

template <class FT>
void
moebius_circle_verticesC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
			  const FT &qx, const FT &qy, const FT &ql, const FT &qm,
			  const FT &rx, const FT &ry, const FT &rl, const FT &rm,
			  FT &x1, FT &y1,
			  FT &x2, FT &y2)
{
  static const FT _2 = 2, _4 = 4;

  CGAL_precondition (pl != ql);

  const FT plx = pl * px, ply = pl * py;
  const FT qx_ = ql * qx - plx, qy_ = ql * qy - ply, ql_ = pl - ql;
  const FT rx_ = rl * rx - plx, ry_ = rl * ry - ply, rl_ = pl - rl;
  const FT P = pl * (px * px + py * py) - pm;
  const FT Q_ = ql * (qx * qx + qy * qy) - qm - P;
  const FT R_ = rl * (rx * rx + ry * ry) - rm - P;

  FT a, b, c;
  if (pl > ql) {
    a = ql_ * rx_ - rl_ * qx_;
    b = ql_ * ry_ - rl_ * qy_;
    c = ql_ * R_ - rl_ * Q_;
  } else {
    a = qx_ * rl_ - rx_ * ql_;
    b = qy_ * rl_ - ry_ * ql_;
    c = Q_ * rl_ - R_ * ql_;
  }

  const FT cy = qy_;
  const FT cd = ql_;

  const Sign sign_a = CGAL_NTS sign (a);
  const Sign sign_b = CGAL_NTS sign (b);
  
  if (sign_b != CGAL::ZERO) {
    const FT tmp1 = a * qy_ - b * qx_;
    const FT X = a * c * cd - _2 * b * tmp1;
    
    const FT tmp2 = c * cd;
    const FT a2b2 = a*a + b*b;
    const FT delta = _4 * (tmp1*tmp1 + (a*qx_ + b*qy_)*c*cd + Q_ * a2b2) - tmp2 * tmp2;
    const FT delta_ = b * CGAL_NTS sqrt (delta);

    const FT D = _2 * a2b2;
    
    const FT _2a = _2 * a;
    const FT _2b = _2 * b;

    if (sign_b == CGAL::POSITIVE) {
      x1 = (X - delta_) / D;
      x2 = (X + delta_) / D;
    } else {
      x1 = (X + delta_) / D;
      x2 = (X - delta_) / D;
    }

    y1 = (c - _2a * x1) / _2b;
    y2 = (c - _2a * x2) / _2b;

  } else { // sign_b == CGAL::ZERO
    CGAL_assertion (sign_a != CGAL::ZERO);

    const FT X = _2 * a * a * qy_;
    const FT tmp1 = a * cy;
    const FT tmp2 = c * cd;
    const FT a2b2 = a*a;
    const FT delta = _4 * (tmp1*tmp1 + (a*qx_)*c*cd + R_ * a2b2) - tmp2 * tmp2;
    const FT delta_ = a * CGAL_NTS sqrt (delta);
    
    const FT D = _2 * a2b2;
    
    const FT _2a = _2 * a;
 

    if (sign_a == CGAL::POSITIVE) {
      y1 = (X + delta_) / D;
      y2 = (X - delta_) / D;
    } else {
      y1 = (X - delta_) / D;
      y2 = (X + delta_) / D;
    }
    x1 = x2 = c / _2a;
  }
}

template <class FT>
void
moebius_line_vertexC2 (const FT &px, const FT &py, const FT &pl, const FT &pm,
		      const FT &qx, const FT &qy, const FT &ql, const FT &qm,
		      const FT &rx, const FT &ry, const FT &rl, const FT &rm,
		      FT &x, FT &y)
{
  static const FT _2 = 2;

  CGAL_assertion (pl == ql);
  CGAL_assertion (pl == rl);

  const FT qx_ = qx - px, rx_ = rx - px;
  const FT qy_ = qy - py, ry_ = ry - py;
  const FT qm_ = qm - pm, rm_ = rm - pm;

  const FT den = _2 * pl * det2x2_by_formula (qx_, qy_, rx_, ry_);

  x = det2x2_by_formula (qm_, qy_, rm_, ry_) / den;
  y = det2x2_by_formula (qx_, qm_, rx_, rm_) / den;
}

#endif // 0

CGAL_END_NAMESPACE

#endif//CHR_CONSTRUCTION_MOEBIUS_DIAGRAM_FTC2_H
