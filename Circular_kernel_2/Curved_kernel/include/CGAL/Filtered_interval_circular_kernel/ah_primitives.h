#ifndef AH_PRIMITIVES_H_
#define AH_PRIMITIVES_H_

#include <CGAL/Filtered_interval_circular_kernel/filtered_interval_circular_kernel_table.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Bbox_2.h>

namespace CGAL{
 namespace CGALi{
 
typedef CGAL::Interval_nt<false>::Protector IntervalProtector;
typedef CGAL::Interval_nt<false> Interval;
typedef CGAL::Bbox_2 Bbox;
typedef CGAL::Filtered_interval_circular_kernel_table Table;
typedef Table::SEARCH_RESULT Res;
const Res INTERSECTS = Table::INTERSECTS; 
const Res DONT_INTERSECT = Table::DONT_INTERSECT;
const Res DONT_KNOW = Table::DONT_KNOW;

inline Res arc_passes_through(
  const Interval &sx, const Interval &sy,
  const Interval &tx, const Interval &ty,
  const Interval &cx, const Interval &cy,
  const Interval &x, const Interval &y) {
  // This function, with some more comparisons can be more accurate
  if(sy.inf() > cy.sup()) {
    if(ty.sup() < cy.inf()) {
      if(y.inf() > cy.sup()) {
        if(x.sup() < sx.inf()) return INTERSECTS;
        else if(x.inf() > sx.sup()) return DONT_INTERSECT;
        else return DONT_KNOW;
      } else if(y.sup() < cy.inf()) {
        if(x.sup() < tx.inf()) return INTERSECTS;
        else if(x.inf() > tx.sup()) return DONT_INTERSECT;
        else return DONT_KNOW;
      } else {
        if(sy.inf() > y.sup() && ty.sup() < y.inf()) return INTERSECTS;
        else return DONT_KNOW;
      }
    } else if(ty.inf() > cy.sup()) {
      if(sx.sup() < tx.inf()) {
        if(y.sup() < cy.inf()) return INTERSECTS;
        else if(y.inf() > cy.sup()) {
          if(x.inf() > sx.sup() && x.sup() < tx.inf()) 
            return DONT_INTERSECT;
          else if((x.sup() < sx.inf()) || (x.inf() > tx.sup())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if(x.inf() > cx.sup()) {
            if(y.sup() < ty.inf()) return INTERSECTS;
            else return DONT_KNOW;
          } else {
            if(y.sup() < sy.inf()) return INTERSECTS;
            else return DONT_KNOW;
          } 
        }
      } else if(sx.inf() > tx.sup()) {
        if(y.sup() < cy.inf()) return DONT_INTERSECT;
        else if(y.inf() > cy.sup()) {
          if(x.inf() > tx.sup() && x.sup() < sx.inf()) 
            return INTERSECTS;
          else if((x.sup() < tx.inf()) || (x.inf() > sx.sup())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if(x.inf() > cx.sup()) {
            if(y.sup() < sy.inf()) return DONT_INTERSECT;
            else return DONT_KNOW;
          } else {
            if(y.sup() < ty.inf()) return DONT_INTERSECT;
            else return DONT_KNOW;
          }
        }
      } else return DONT_KNOW;
    } else {
      if(tx.sup() < sx.inf()) {
        if(y.sup() < cy.inf()) {
          if((x.inf() > tx.sup()) || (y.sup() < ty.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else if(y.inf() > cy.sup()) {
          if(x.inf() > sx.sup()) return DONT_INTERSECT;
          else if(x.sup() < sx.inf() && (y.inf() > ty.sup())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if((x.inf() > tx.sup()) && (y.sup() < sy.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW; 
        }
      } else if(tx.inf() > sx.sup()) {
        if(y.inf() > cy.sup()) {
          if(x.sup() < sx.inf()) return INTERSECTS;
          else return DONT_KNOW;
        } else if(y.sup() < cy.inf()) {
          if((x.sup() < tx.inf()) || (y.sup() < ty.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if((x.sup() < tx.inf()) && (y.sup() < sy.inf())) 
            return INTERSECTS;
          else return DONT_KNOW;
        }
      } else return DONT_KNOW;
    }
  } else if(sy.sup() < cy.inf()) {
    if(ty.inf() > cy.sup()) {
      if(y.inf() > cy.sup()) {
        if(x.inf() > tx.sup()) return INTERSECTS;
        else if(x.sup() < tx.inf()) return DONT_INTERSECT;
        else return DONT_KNOW;
      } else if(y.sup() < cy.inf()) {
        if(x.inf() > sx.sup()) return INTERSECTS;
        else if(x.sup() < sx.inf()) return DONT_INTERSECT;
        else return DONT_KNOW;
      } else {
        if(sy.inf() > y.sup() && ty.sup() < y.inf()) return INTERSECTS;
        else return DONT_KNOW;
      }
    } else if(ty.sup() < cy.inf()) {
      if(tx.inf() > sx.sup()) {
        if(y.inf() > cy.sup()) return DONT_INTERSECT;
        else if(y.sup() < cy.inf()) {
          if((x.inf() > tx.sup()) || (x.sup() < sx.inf())) 
            return DONT_INTERSECT;
          else if((sx.sup() < x.inf()) && (x.sup() < tx.inf())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if(x.inf() > cx.sup()) {
            if(y.inf() > ty.sup()) return DONT_INTERSECT;
            else return DONT_KNOW;
          } else {
            if(y.inf() > sy.sup()) return DONT_INTERSECT;
            else return DONT_KNOW;
          }  
        }
      } else if(tx.sup() < sx.inf()) {
        if(y.inf() > cy.sup()) return INTERSECTS;
        else if(y.sup() < cy.inf()) {
          if((x.inf() > sx.sup()) || (x.sup() < tx.inf())) 
            return INTERSECTS;
          else if((tx.sup() < x.inf()) && (x.sup() < sx.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if(x.inf() > cx.sup()) {
            if(y.inf() > sy.sup()) return INTERSECTS;
            else return DONT_KNOW;
          } else {
            if(y.inf() > ty.sup()) return INTERSECTS;
            else return DONT_KNOW;
          }  
        }
      } else return DONT_KNOW;
    } else {
      if(tx.inf() > sx.sup()) {
        if(y.inf() > cy.sup()) {
          if((x.sup() < tx.inf()) || (y.inf() > ty.sup())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else if(y.sup() < cy.inf()) {
          if((x.sup() < sx.inf())) return DONT_INTERSECT;
          else if((sx.sup() < x.inf()) && (x.sup() < tx.inf())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if(x.sup() < tx.inf()) {
            if(x.sup() < sx.inf()) return DONT_INTERSECT;
            else return DONT_KNOW; 
          } else return DONT_KNOW;
        }
      } else if(tx.sup() < sx.inf()) {
        if(y.inf() > cy.sup()) {
          if((x.inf() > tx.sup()) || (y.inf() > ty.sup())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else if(y.sup() < cy.inf()) {
          if((x.inf() > sx.sup())) return INTERSECTS;
          else if((tx.sup() < x.inf()) && (x.sup() < sx.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if(x.inf() > tx.sup()) {
            if(x.inf() > sx.sup()) return INTERSECTS;
            else return DONT_KNOW; 
          } else return DONT_KNOW; 
        } 
      } else return DONT_KNOW;
    }
  } else {
    if(sx.sup() < cx.inf()) {
      if(ty.sup() < cy.inf()) {
        if(y.inf() > cy.sup()) {
          if(x.inf() > sx.sup()) return DONT_INTERSECT;
          else return DONT_KNOW;
        } else if(y.sup() < cy.inf()) {
          if(x.inf() > tx.sup()) return DONT_INTERSECT;
          else if((sx.sup() < x.inf()) && (x.sup() < tx.inf())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if(x.inf() > sx.sup()) {
            if(x.inf() > tx.sup()) return DONT_INTERSECT;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        }
      } else if(ty.inf() > cy.sup()) {
        if(y.sup() < cy.inf()) {
          if(x.inf() > sx.sup()) return INTERSECTS;
          else return DONT_KNOW;
        } else if(y.inf() > cy.sup()) {
          if(x.inf() > tx.sup()) return INTERSECTS;
          else if((sx.sup() < x.inf()) && (x.sup() < tx.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if(x.inf() > sx.sup()) {
            if(x.inf() > tx.sup()) return INTERSECTS;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        }
      } else {
        if(tx.inf() > sx.sup()) {
          if(y.inf() > cy.sup()) {
            if((x.inf() > sx.sup()) && (x.sup() < tx.inf())) 
              return DONT_INTERSECT;
            else return DONT_KNOW;
          } else if(y.sup() < cy.inf()) {
            if((x.inf() > sx.sup()) && (x.sup() < tx.inf())) 
              return INTERSECTS;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        } else return DONT_KNOW;
      }
    } else {
      if(ty.sup() < cy.inf()) {
        if(y.inf() > cy.sup()) {
          if(x.sup() < sx.inf()) return INTERSECTS;
          else return DONT_KNOW;
        } else if(y.sup() < cy.inf()) {
          if(x.sup() < tx.inf()) return INTERSECTS;
          else if((tx.sup() < x.inf()) && (x.sup() < sx.inf())) 
            return DONT_INTERSECT;
          else return DONT_KNOW;
        } else {
          if(x.sup() < sx.inf()) {
            if(x.sup() < tx.inf()) return INTERSECTS;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        }
      } else if(ty.inf() > cy.sup()) {
        if(y.sup() < cy.inf()) {
          if(x.sup() < sx.inf()) return DONT_INTERSECT;
          else return DONT_KNOW;
        } else if(y.inf() > cy.sup()) {
          if(x.sup() < tx.inf()) return DONT_INTERSECT;
          else if((tx.sup() < x.inf()) && (x.sup() < sx.inf())) 
            return INTERSECTS;
          else return DONT_KNOW;
        } else {
          if(x.sup() < sx.inf()) {
            if(x.sup() < tx.inf()) return DONT_INTERSECT;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        }
      } else {
        if(tx.sup() < sx.inf()) {
          if(y.inf() > cy.sup()) {
            if((x.inf() > tx.sup()) && (x.sup() < sx.inf())) 
              return INTERSECTS;
            else return DONT_KNOW;
          } else if(y.sup() < cy.inf()) {
            if((x.inf() > tx.sup()) && (x.sup() < sx.inf())) 
              return DONT_INTERSECT;
            else return DONT_KNOW;
          } else return DONT_KNOW;
        } else return DONT_KNOW;
      }
    }
  }
}

template <typename CK>
inline Res intersect(const typename CK::Circular_arc_2 &arc0,
                     const typename CK::Circular_arc_2 &arc1) {	
  bool maybe = false;	
  IntervalProtector ip;
  Interval r0_2 = CGAL::to_interval(arc0.squared_radius());
  Interval r1_2 = CGAL::to_interval(arc1.squared_radius());
  Interval a0c_x = CGAL::to_interval(arc0.center().x());
  Interval a0c_y = CGAL::to_interval(arc0.center().y());
  Interval a1c_x = CGAL::to_interval(arc1.center().x());
  Interval a1c_y = CGAL::to_interval(arc1.center().y());
  Interval x1_m_x0 = a1c_x - a0c_x;
  Interval y1_m_y0 = a1c_y - a0c_y;
  Interval d2 = CGAL::square(x1_m_x0) + CGAL::square(y1_m_y0);
  if(d2.inf() == 0.0) { 
  	if(d2.sup() != 0.0) return DONT_KNOW;
	if(!r0_2.do_overlap(r1_2)) return DONT_INTERSECT;
	return DONT_KNOW;
  }
  Interval r0_2_d2 = r0_2 / d2;
  Interval r1_2_d2 = r1_2 / d2;
  Interval a_d = 0.5 * (r0_2_d2 - r1_2_d2 + 1);
  Interval h_d;
  Res rs1 = Table::get_h_div_d_with_excellent_percision(r0_2_d2, a_d, h_d);
  if(rs1 == DONT_INTERSECT) return rs1;
  if(rs1 == DONT_KNOW) {
  	maybe = true;
  }
  Interval x0_p_ad_x_x1mx0 = a0c_x + a_d * x1_m_x0;
  Interval y0_p_ad_x_y1my0 = a0c_y + a_d * y1_m_y0;
  Interval hd_y1_m_y0 = h_d * y1_m_y0;
  Interval hd_x1_m_x0 = h_d * x1_m_x0;
  Interval x1_res = x0_p_ad_x_x1mx0 + hd_y1_m_y0;
  Interval x2_res = x0_p_ad_x_x1mx0 - hd_y1_m_y0;
  Interval y1_res = y0_p_ad_x_y1my0 - hd_x1_m_x0;
  Interval y2_res = y0_p_ad_x_y1my0 + hd_x1_m_x0;
  Res rs2;  
  
  // compare end_points with roots             
  // we assume that the bboxes (x1_res, y1_res) (x2_res,y2_res) 
  // dont include the center of the circle 
  if(arc0.is_full() && arc1.is_full()) {
   	return INTERSECTS;
  } else if(arc0.is_full()) {
   	// those stuffs are declared here for efficiency issues
   	// since we may not need them outside
   	Bbox source1 = arc1.source().bbox();
   	Bbox target1 = arc1.target().bbox();
   	Interval arc1_sx = Interval(source1.xmin(),source1.xmax());
   	Interval arc1_sy = Interval(source1.ymin(),source1.ymax());
   	Interval arc1_tx = Interval(target1.xmin(),target1.xmax());
   	Interval arc1_ty = Interval(target1.ymin(),target1.ymax());

  	if((rs1 = arc_passes_through(arc1_sx, arc1_sy,
   	                             arc1_tx, arc1_ty,
   	                             a1c_x, a1c_y,
   	                             x1_res, y1_res)) == 
   	                             INTERSECTS) 
  	  return maybe ? DONT_KNOW : rs1; 
    	
    if((rs2 = arc_passes_through(arc1_sx, arc1_sy,
                                 arc1_tx, arc1_ty,
                                 a1c_x, a1c_y,
                                 x2_res, y2_res)) == 
                                 INTERSECTS) 
      return maybe ? DONT_KNOW : rs2;
    	
   	if(rs1 == DONT_INTERSECT && rs2 == DONT_INTERSECT)                             
   	  return DONT_INTERSECT;
    return DONT_KNOW;
  } else if(arc1.is_full()) {
    // those stuffs are declared here for efficiency issues
    // since we may not need them outside
    Bbox source0 = arc0.source().bbox();
   	Bbox target0 = arc0.target().bbox();
   	Interval arc0_sx = Interval(source0.xmin(),source0.xmax());
   	Interval arc0_sy = Interval(source0.ymin(),source0.ymax());
   	Interval arc0_tx = Interval(target0.xmin(),target0.xmax());
   	Interval arc0_ty = Interval(target0.ymin(),target0.ymax());
    	
    if((rs1 = arc_passes_through(arc0_sx, arc0_sy,
                                 arc0_tx, arc0_ty,
                                 a0c_x, a0c_y,
                                 x1_res, y1_res)) == 
                                 INTERSECTS) 
      return maybe ? DONT_KNOW : rs1; 
    	
    if((rs2 = arc_passes_through(arc0_sx, arc0_sy,
                                 arc0_tx, arc0_ty,
                                 a0c_x, a0c_y,
                                 x2_res, y2_res)) == 
                                 INTERSECTS) 
      return maybe ? DONT_KNOW : rs2;
    	
      if(rs1 == DONT_INTERSECT && rs2 == DONT_INTERSECT)                             
        return DONT_INTERSECT;
      return DONT_KNOW;
  } else {
    // those stuffs are declared here for efficiency issues
    // since we may not need them outside
    Bbox source0 = arc0.source().bbox();
   	Bbox target0 = arc0.target().bbox();
   	Interval arc0_sx = Interval(source0.xmin(),source0.xmax());
   	Interval arc0_sy = Interval(source0.ymin(),source0.ymax());
   	Interval arc0_tx = Interval(target0.xmin(),target0.xmax());
   	Interval arc0_ty = Interval(target0.ymin(),target0.ymax());
   	Bbox source1 = arc1.source().bbox();
   	Bbox target1 = arc1.target().bbox();
   	Interval arc1_sx = Interval(source1.xmin(),source1.xmax());
   	Interval arc1_sy = Interval(source1.ymin(),source1.ymax());
   	Interval arc1_tx = Interval(target1.xmin(),target1.xmax());
   	Interval arc1_ty = Interval(target1.ymin(),target1.ymax());
        
    Res rs3, rs4;

    if(((rs1 = arc_passes_through(arc0_sx, arc0_sy,
    	                          arc0_tx, arc0_ty,
    	                          a0c_x, a0c_y,
    	                          x1_res, y1_res)) == 
    	                          INTERSECTS) &&
       ((rs2 = arc_passes_through(arc1_sx, arc1_sy,
    	                          arc1_tx, arc1_ty,
    	                          a1c_x, a1c_y,
    	                          x1_res, y1_res)) == 
    	                          INTERSECTS)) 
      return  maybe ? DONT_KNOW : rs1; 
    	
    if(((rs3 = arc_passes_through(arc0_sx, arc0_sy,
                                  arc0_tx, arc0_ty,
    	                          a0c_x, a0c_y,
    	                          x2_res, y2_res)) == 
    	                          INTERSECTS) &&
       ((rs4 = arc_passes_through(arc1_sx, arc1_sy,
    	                          arc1_tx, arc1_ty,
    	                          a1c_x, a1c_y,
    	                          x2_res, y2_res)) == 
    	                          INTERSECTS)) 
      return  maybe ? DONT_KNOW : rs3; 
    	
    if(((rs1 == DONT_INTERSECT) || 
        (rs2 == DONT_INTERSECT)) &&
       ((rs3 == DONT_INTERSECT) || 
    	(rs4 == DONT_INTERSECT)))                             
      return DONT_INTERSECT;                               	                           
    return DONT_KNOW;
  }
}

template <typename CK>
inline Res intersect(const typename CK::Circle_2 &c0,
                     const typename CK::Circle_2 &c1) {		
  IntervalProtector ip;
  Interval r0_2 = CGAL::to_interval(c0.squared_radius());
  Interval r1_2 = CGAL::to_interval(c1.squared_radius());
  Interval a0c_x = CGAL::to_interval(c0.center().x());
  Interval a0c_y = CGAL::to_interval(c0.center().y());
  Interval a1c_x = CGAL::to_interval(c1.center().x());
  Interval a1c_y = CGAL::to_interval(c1.center().y());
  Interval d2 = CGAL::square(a1c_x - a0c_x) + CGAL::square(a1c_y - a0c_y);
  if(d2.inf() == 0.0) { 
	if(!r0_2.do_overlap(r1_2)) return DONT_INTERSECT;
	else return DONT_KNOW;
  }
  Interval r0_2_d2 = r0_2 / d2;
  Interval r1_2_d2 = r1_2 / d2;
  Interval a_d = 0.5 * (r0_2_d2 - r1_2_d2 + 1);
  Interval a2_d2 = CGAL::square(a_d);
  if(a2_d2.inf() > r0_2_d2.sup()) return DONT_INTERSECT;
  if(a2_d2.sup() >= r0_2_d2.inf()) return DONT_KNOW;
  return INTERSECTS;
}

inline Res seg_passes_through(
  const Interval &sx, const Interval &sy,
  const Interval &tx, const Interval &ty,
  const Interval &x, const Interval &y) {
  if(sx.inf() > tx.sup()) {
    if(x.inf() > sx.sup()) return DONT_INTERSECT;
    else if(x.sup() < tx.inf()) return DONT_INTERSECT;
    else if((tx.sup() < x.inf()) && (x.sup() < sx.inf())) return INTERSECTS;
    else {
      // Y_comparison - same block
      if(sy.inf() > ty.sup()) {
        if(y.inf() > sy.sup()) return DONT_INTERSECT;
        else if(y.sup() < ty.inf()) return DONT_INTERSECT;
        else if((ty.sup() < y.inf()) && (y.sup() < sy.inf())) return INTERSECTS;
        else return DONT_KNOW;
      } else if(sy.sup() < ty.inf()) {
        if(y.inf() > ty.sup()) return DONT_INTERSECT;
        else if(y.sup() < sy.inf()) return DONT_INTERSECT;
        else if((sy.sup() < y.inf()) && (y.sup() < ty.inf())) return INTERSECTS;
        else return DONT_KNOW;
      } else {
        if(y.inf() > std::max(sy.sup(), ty.sup())) return DONT_INTERSECT;
        else if(y.sup() < std::min(sy.inf(), ty.inf())) return DONT_INTERSECT;
        else return DONT_KNOW;
      }
    }
  } else if(sx.sup() < tx.inf()) {
    if(x.inf() > tx.sup()) return DONT_INTERSECT;
    else if(x.sup() < sx.inf()) return DONT_INTERSECT;
    else if((sx.sup() < x.inf()) && (x.sup() < tx.inf())) return INTERSECTS;
    else {
      // Y_comparison - same block
      if(sy.inf() > ty.sup()) {
        if(y.inf() > sy.sup()) return DONT_INTERSECT;
        else if(y.sup() < ty.inf()) return DONT_INTERSECT;
        else if((ty.sup() < y.inf()) && (y.sup() < sy.inf())) return INTERSECTS;
        else return DONT_KNOW;
      } else if(sy.sup() < ty.inf()) {
        if(y.inf() > ty.sup()) return DONT_INTERSECT;
        else if(y.sup() < sy.inf()) return DONT_INTERSECT;
        else if((sy.sup() < y.inf()) && (y.sup() < ty.inf())) return INTERSECTS;
        else return DONT_KNOW;
      } else {
        if(y.inf() > std::max(sy.sup(), ty.sup())) return DONT_INTERSECT;
        else if(y.sup() < std::min(sy.inf(), ty.inf())) return DONT_INTERSECT;
        else return DONT_KNOW;
      }
    }
  } else {
    // Y_comparison - same block
    if(sy.inf() > ty.sup()) {
      if(y.inf() > sy.sup()) return DONT_INTERSECT;
      else if(y.sup() < ty.inf()) return DONT_INTERSECT;
      else if((ty.sup() < y.inf()) && (y.sup() < sy.inf())) return INTERSECTS;
      else return DONT_KNOW;
    } else if(sy.sup() < ty.inf()) {
      if(y.inf() > ty.sup()) return DONT_INTERSECT;
      else if(y.sup() < sy.inf()) return DONT_INTERSECT;
      else if((sy.sup() < y.inf()) && (y.sup() < ty.inf())) return INTERSECTS;
      else return DONT_KNOW;
    } else {
      if(y.inf() > std::max(sy.sup(), ty.sup())) return DONT_INTERSECT;
      else if(y.sup() < std::min(sy.inf(), ty.inf())) return DONT_INTERSECT;
      else return DONT_KNOW;
    }
  }	
}

template <typename CK>
inline Res intersect(const typename CK::Line_arc_2 &arc0,
                     const typename CK::Circular_arc_2 &arc1) {	
  bool maybe = false;
  IntervalProtector ip;
  Interval r_2 = CGAL::to_interval(arc1.squared_radius());
  Interval cx = CGAL::to_interval(arc1.center().x());
  Interval cy = CGAL::to_interval(arc1.center().y()); 
  Bbox source0 = arc0.source().bbox();
  Bbox target0 = arc0.target().bbox();
  Interval x1 = Interval(source0.xmin(),source0.xmax());
  Interval y1 = Interval(source0.ymin(),source0.ymax());
  Interval x2 = Interval(target0.xmin(),target0.xmax());
  Interval y2 = Interval(target0.ymin(),target0.ymax());
  Interval dx = x2 - x1;
  Interval dy = y2 - y1;
  Interval dr_2 = square(dx) + square(dy);
  Interval D = (x1-cx)*(y2-cy) - (x2-cx)*(y1-cy);
  Interval delta = r_2 * dr_2 - square(D);
  if(delta.sup() < 0) return DONT_INTERSECT;
  if(delta.inf() <= 0) {maybe = true;}//return DONT_KNOW;
  Interval sqrt_delta = sqrt(delta); 
  Interval cx_p_Ddy_div_dr2 = cx + (D * dy) / dr_2;
  Interval cy_m_Ddx_div_dr2 = cy - (D * dx) / dr_2;
  Interval dx_sqrt_delta_dr2 = (dx * sqrt_delta) / dr_2;
  Interval dy_sqrt_delta_dr2 = (dy * sqrt_delta) / dr_2;
  Interval x1_res, y1_res, x2_res, y2_res;
  x1_res = cx_p_Ddy_div_dr2 + dx_sqrt_delta_dr2;
  y1_res = cy_m_Ddx_div_dr2 + dy_sqrt_delta_dr2;
  x2_res = cx_p_Ddy_div_dr2 - dx_sqrt_delta_dr2;
  y2_res = cy_m_Ddx_div_dr2 - dy_sqrt_delta_dr2;
  Res rs1, rs2;
  if(arc1.is_full()) {
    if((rs1 = seg_passes_through(x1,y1,x2,y2,x1_res,y1_res)) == INTERSECTS) 
      return maybe ? DONT_KNOW : rs1;
    if((rs2 = seg_passes_through(x1,y1,x2,y2,x1_res,y1_res)) == INTERSECTS) 
      return maybe ? DONT_KNOW : rs2;
    if((rs1 == DONT_INTERSECT) && (rs2 == DONT_INTERSECT)) return DONT_INTERSECT;
    return DONT_KNOW;
  } else {
    // those stuffs are declared here for efficiency issues
    // since we may not need them outside
    Bbox source1 = arc1.source().bbox();
    Bbox target1 = arc1.target().bbox();
    Interval arc1_sx = Interval(source1.xmin(),source1.xmax());
    Interval arc1_sy = Interval(source1.ymin(),source1.ymax());
    Interval arc1_tx = Interval(target1.xmin(),target1.xmax());
    Interval arc1_ty = Interval(target1.ymin(),target1.ymax());        
  
    Res rs3, rs4;

    if(((rs1 = seg_passes_through(x1,y1,x2,y2,x1_res,y1_res)) == INTERSECTS) &&
       ((rs2 = arc_passes_through(arc1_sx, arc1_sy,
                                  arc1_tx, arc1_ty,
                                  cx, cy,
                                  x1_res, y1_res)) == 
    	                          INTERSECTS)) return maybe ? DONT_KNOW : rs1; 
    	
    if(((rs3 = seg_passes_through(x1,y1,x2,y2,x2_res,y2_res)) == INTERSECTS) &&
       ((rs4 = arc_passes_through(arc1_sx, arc1_sy,
    	                          arc1_tx, arc1_ty,
    	                          cx, cy,
    	                          x2_res, y2_res)) == 
    	                          INTERSECTS)) return maybe ? DONT_KNOW : rs3; 
    	
    if(((rs1 == DONT_INTERSECT) || 
        (rs2 == DONT_INTERSECT)) &&
       ((rs3 == DONT_INTERSECT) || 
        (rs4 == DONT_INTERSECT)))
       return DONT_INTERSECT;
    return DONT_KNOW;
  }
}

template <typename CK>
inline Res intersect(const typename CK::Circular_arc_2 &arc0,
                     const typename CK::Line_arc_2 &arc1) {	
  return intersect(arc1,arc0);
}

template <typename CK>
inline Res intersect(const typename CK::Line_arc_2 &l0,
                     const typename CK::Line_arc_2 &l1) {	
  IntervalProtector ip;
  Bbox source0 = l0.source().bbox();
  Bbox target0 = l0.target().bbox();
  Interval l0sx = Interval(source0.xmin(),source0.xmax());
  Interval l0sy = Interval(source0.ymin(),source0.ymax());
  Interval l0tx = Interval(target0.xmin(),target0.xmax());
  Interval l0ty = Interval(target0.ymin(),target0.ymax());
  Bbox source1 = l1.source().bbox();
  Bbox target1 = l1.target().bbox();
  Interval l1sx = Interval(source1.xmin(),source1.xmax());
  Interval l1sy = Interval(source1.ymin(),source1.ymax());
  Interval l1tx = Interval(target1.xmin(),target1.xmax());
  Interval l1ty = Interval(target1.ymin(),target1.ymax());
   
  if(std::min(l0sx.inf(), l0tx.inf()) > 
     std::max(l1sx.sup(), l1tx.sup())) return DONT_INTERSECT;
  if(std::min(l1sx.inf(), l1tx.inf()) > 
     std::max(l0sx.sup(), l0tx.sup())) return DONT_INTERSECT;
  if(std::min(l0sy.inf(), l0ty.inf()) > 
     std::max(l1sy.sup(), l1ty.sup())) return DONT_INTERSECT;
  if(std::min(l1sy.inf(), l1ty.inf()) > 
     std::max(l0sy.sup(), l0ty.sup())) return DONT_INTERSECT;
   
  if(std::min(l0sx.sup(), l0tx.sup()) > 
     std::max(l1sx.inf(), l1tx.inf())) return DONT_KNOW;
  if(std::min(l1sx.sup(), l1tx.sup()) > 
     std::max(l0sx.inf(), l0tx.inf())) return DONT_KNOW;
  if(std::min(l0sy.sup(), l0ty.sup()) > 
     std::max(l1sy.inf(), l1ty.inf())) return DONT_KNOW;
  if(std::min(l1sy.sup(), l1ty.sup()) > 
     std::max(l0sy.inf(), l0ty.inf())) return DONT_KNOW;

  Interval vs20s10x = l1sx - l0sx;
  Interval vs20s10y = l1sy - l0sy;
  Interval vs11s10x = l0tx - l0sx;
  Interval vs11s10y = l0ty - l0sy;

  Interval cross1 = vs20s10x * vs11s10y - vs20s10y * vs11s10x;

  int sign1;
  if(cross1.inf() > 0) {
    sign1 = 1;
  } else if(cross1.sup() < 0) {
    sign1 = -1;
  } else return DONT_KNOW;

  Interval vs21s10x = l1tx - l0sx;
  Interval vs21s10y = l1ty - l0sy;

  Interval cross2 = vs21s10x * vs11s10y - vs21s10y * vs11s10x;

  int sign2;
  if(cross2.inf() > 0) {
    sign2 = 1;
  } else if(cross2.sup() < 0) {
    sign2 = -1;
  } else return DONT_KNOW; 
  
  if(sign1 == sign2) return DONT_INTERSECT;

  Interval vs21s20x = l1tx - l1sx;
  Interval vs21s20y = l1ty - l1sy;

  Interval cross3 = vs20s10y * vs21s20x - vs20s10x * vs21s20y;
 
  int sign3;
  if(cross3.inf() > 0) {
    sign3 = 1;
  } else if(cross3.sup() < 0) {
    sign3 = -1;
  } else return DONT_KNOW; 

  Interval vs11s20x = l0tx - l1sx;
  Interval vs11s20y = l0ty - l1sy;
  Interval cross4 = vs11s20x * vs21s20y - vs11s20y * vs21s20x;

  int sign4;
  if(cross4.inf() > 0) {
    sign4 = 1;
  } else if(cross4.sup() < 0) {
    sign4 = -1;
  } else return DONT_KNOW; 

  if(sign3 == sign4) return DONT_INTERSECT;

  return INTERSECTS;
}
   
 	
 }
}

#endif /*AH_PRIMITIVES_H_*/
