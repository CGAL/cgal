// Copyright (c) 2006,2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Kerber    <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ARR_ALGEBRAIC_SEGMENT_TRAITS
#define CGAL_ARR_ALGEBRAIC_SEGMENT_TRAITS

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>

namespace CGAL {

template< class Coefficient_ > 
class Arr_algebraic_segment_traits_2 {
    
public:

    enum Site_of_point {
        POINT_IN_INTERIOR = 0,
        MIN_ENDPOINT = -1,
        MAX_ENDPOINT = 1 };
    
    // Template argument
    
    typedef Coefficient_ Coefficient;
    
    // 'derivation'
    typedef CGAL::Algebraic_kernel_d_1< Coefficient > Algebraic_kernel_d_1;

    typedef CGAL::Algebraic_curve_kernel_2< Algebraic_kernel_d_1 > 
        Algebraic_kernel_d_2;
    
    typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_kernel_d_2 > CKvA_2;

    typedef  CGAL::Arr_algebraic_segment_traits_2<Coefficient> Self;

    // Default constructor
    Arr_algebraic_segment_traits_2 () {}


    // Copy constructor
    Arr_algebraic_segment_traits_2 (const  Self& /* s */) { /* No state...*/}

    // Assignement operator
    const Self& operator= (const Self& s)
        {return s;}
    
    // public types
    
    typedef typename Algebraic_kernel_d_2::Algebraic_real_1 Algebraic_real_1;
    typedef typename Algebraic_kernel_d_2::Bound Bound;

    typedef typename Algebraic_kernel_d_2::Polynomial_2 Polynomial_2;

    /// public types for ArrangementTraits_2
    
    typedef typename Algebraic_kernel_d_2::Curve_analysis_2 Curve_2;
    
    typedef typename CKvA_2::Arc_2 X_monotone_curve_2;

    typedef typename CKvA_2::Point_2 Point_2;


    typedef typename CKvA_2::Has_left_category Has_left_category;
    typedef typename CKvA_2::Has_merge_category Has_merge_category;
    typedef typename CKvA_2::Has_do_intersect_category 
      Has_do_intersect_category;

    typedef typename CKvA_2::Left_side_category Left_side_category;
    typedef typename CKvA_2::Bottom_side_category Bottom_side_category;
    typedef typename CKvA_2::Top_side_category Top_side_category;
    typedef typename CKvA_2::Right_side_category Right_side_category;

    typedef typename CKvA_2::Multiplicity Multiplicity;

    typedef typename CKvA_2::Compare_x_2 Compare_x_2;
    Compare_x_2 compare_x_2_object() const {
        return CKvA_2::instance().compare_x_2_object();
    }

    typedef typename CKvA_2::Compare_xy_2 Compare_xy_2;
    Compare_xy_2 compare_xy_2_object() const {
        return CKvA_2::instance().compare_xy_2_object();
    }

    typedef typename CKvA_2::Compare_endpoints_xy_2 Compare_endpoints_xy_2;
    Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const {
        return CKvA_2::instance().compare_endpoints_xy_2_object();
    }

    typedef typename CKvA_2::Equal_2 Equal_2;
    Equal_2 equal_2_object() const {
        return CKvA_2::instance().equal_2_object();
    }

    typedef typename CKvA_2::Parameter_space_in_y_2 Parameter_space_in_y_2;
    Parameter_space_in_y_2 parameter_space_in_y_2_object() const {
        return CKvA_2::instance().parameter_space_in_y_2_object();
    }

    typedef typename CKvA_2::Compare_y_near_boundary_2 
       Compare_y_near_boundary_2;
    Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const {
        return CKvA_2::instance().compare_y_near_boundary_2_object();
    }

    typedef typename CKvA_2::Parameter_space_in_x_2 Parameter_space_in_x_2;
    Parameter_space_in_x_2 parameter_space_in_x_2_object() const {
        return CKvA_2::instance().parameter_space_in_x_2_object();
    }

    typedef typename CKvA_2::Compare_x_at_limit_2 Compare_x_at_limit_2;
    Compare_x_at_limit_2 compare_x_at_limit_2_object() const {
        return CKvA_2::instance().compare_x_at_limit_2_object();
    }

    typedef typename CKvA_2::Compare_x_near_limit_2  Compare_x_near_limit_2;
    Compare_x_near_limit_2 compare_x_near_limit_2_object() const {
        return CKvA_2::instance().compare_x_near_limit_2_object();
    }

    typedef typename CKvA_2::Construct_min_vertex_2 Construct_min_vertex_2;
    Construct_min_vertex_2 construct_min_vertex_2_object() const {
        return CKvA_2::instance().construct_min_vertex_2_object();
    }

    typedef typename CKvA_2::Construct_max_vertex_2 Construct_max_vertex_2;
    Construct_max_vertex_2 construct_max_vertex_2_object() const {
        return CKvA_2::instance().construct_max_vertex_2_object();
    }

    typedef typename CKvA_2::Construct_opposite_2 Construct_opposite_2;
    Construct_opposite_2 construct_opposite_2_object() const {
        return CKvA_2::instance().construct_opposite_2_object();
    }
  
    typedef typename CKvA_2::Is_vertical_2 Is_vertical_2;
    Is_vertical_2 is_vertical_2_object() const {
        return CKvA_2::instance().is_vertical_2_object();
    }

    typedef typename CKvA_2::Compare_y_at_x_2 Compare_y_at_x_2;
    Compare_y_at_x_2 compare_y_at_x_2_object() const {
        return CKvA_2::instance().compare_y_at_x_2_object();
    }

    typedef typename CKvA_2::Compare_y_at_x_left_2 Compare_y_at_x_left_2;
    Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const {
        return CKvA_2::instance().compare_y_at_x_left_2_object();
    }

    typedef typename CKvA_2::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
    Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const {
        return CKvA_2::instance().compare_y_at_x_right_2_object();
    }
    
    typedef typename CKvA_2::Intersect_2 Intersect_2;
    Intersect_2 intersect_2_object() const {
        return CKvA_2::instance().intersect_2_object();
    }

    typedef typename CKvA_2::Split_2 Split_2;
    Split_2 split_2_object() const {
        return CKvA_2::instance().split_2_object();
    }

    typedef typename CKvA_2::Are_mergeable_2 Are_mergeable_2;
    Are_mergeable_2 are_mergeable_2_object() const {
        return CKvA_2::instance().are_mergeable_2_object();
    }

    typedef typename CKvA_2::Merge_2 Merge_2;
    Merge_2 merge_2_object() const {
        return CKvA_2::instance().merge_2_object();
    }

    typedef typename CKvA_2::Make_x_monotone_2 Make_x_monotone_2;
    Make_x_monotone_2 make_x_monotone_2_object() const {
        return Make_x_monotone_2(&CKvA_2::instance());
    }

    typedef typename CKvA_2::Is_on_2 Is_on_2;

    Is_on_2 is_on_2_object() const {
        return Is_on_2(&CKvA_2::instance());
    }

    typedef typename CKvA_2::Construct_point_2 Construct_point_2;

    Construct_point_2 construct_point_2_object() const {
        return Construct_point_2(&CKvA_2::instance());
    }

    
    class Construct_x_monotone_segment_2 : public 
    CGAL::internal::Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< CKvA_2 > {
        
    public:

        typedef CGAL::internal::Curved_kernel_via_analysis_2_Functors::
        Curved_kernel_via_analysis_2_functor_base< CKvA_2 > Base;

        Construct_x_monotone_segment_2(CKvA_2* kernel) : Base(kernel) {}


    private:

        // helper function, returning the number of incident branches
        // of a point wrt to a curve
        std::pair<std::pair<int,int>, bool>
            branch_numbers_and_vertical(Curve_2 cv, Point_2 p) const {

            Equal_2 equal = this->_ckva()->equal_2_object();
            
            Algebraic_real_1 x = p.x();
            int no;
            bool is_event;
            cv.x_to_index(x,no,is_event);
            if(! is_event) {
                return std::make_pair(std::make_pair(1,1),false);
            }
            typename Curve_2::Status_line_1 status_line
                = cv.status_line_at_event(no);
            bool vertical = status_line.covers_line();
            for(int i = 0; i < status_line.number_of_events(); i++) {
                Point_2 curr_point(x,cv,i);
                if(equal(p,curr_point)) {
                    return std::make_pair(
                            status_line.number_of_incident_branches(i),
                            vertical);
                }
            }
            return std::make_pair(std::make_pair(0,0),vertical);
        }

        // abbrevation for convenience
        bool is_one_one(Curve_2 cv, Point_2 p) const {

            std::pair<std::pair<int,int>,bool> branches 
                = branch_numbers_and_vertical(cv,p);
            return branches.first.first==1 && branches.first.second==1
                && !branches.second;

        }

        template<typename OutputIterator> OutputIterator
        x_monotone_segment (Curve_2 cv, 
                            Point_2 p, 
                            boost::optional<Point_2> start,
                            boost::optional<Point_2> end,
                            OutputIterator out) const {
            
            //CGAL_assertion(is_one_one(cv,p));

            std::list<X_monotone_curve_2> segs;

            Is_on_2 on_arc 
                = this->_ckva()->is_on_2_object();
            Construct_min_vertex_2 left 
                = this->_ckva()->construct_min_vertex_2_object();
            Construct_max_vertex_2 right 
                = this->_ckva()->construct_max_vertex_2_object();
            Equal_2 equal = this->_ckva()->equal_2_object();

            std::vector<CGAL::Object> arcs;
            this->_ckva()->make_x_monotone_2_object() 
                (cv,std::back_inserter(arcs));
            typename std::vector<CGAL::Object>::const_iterator 
                it = arcs.begin(),helper;
            X_monotone_curve_2 it_seg;
            CGAL::assign(it_seg,*it);
            while(it!=arcs.end()) {
                if( on_arc(p,it_seg) ) {
                    break;
                }
                CGAL_assertion(it!=arcs.end());
                it++;
                CGAL::assign(it_seg,*it);
            }
            
            bool left_on_arc = start && on_arc(start.get(),it_seg);
            bool right_on_arc = end && on_arc(end.get(),it_seg);
            
            if( left_on_arc && right_on_arc ) {
                segs.push_back(it_seg.trim(start.get(),end.get()));
            }
            if(left_on_arc && (!right_on_arc)) {
                if(!it_seg.is_finite(CGAL::ARR_MAX_END) ||
                     !equal(start.get(),right(it_seg))) {
                  if(it_seg.is_finite(CGAL::ARR_MIN_END) && equal(start.get(),left(it_seg))) {
                        segs.push_back(it_seg);
                    } else {
                        X_monotone_curve_2 split1,split2;
                        it_seg.split(start.get(),split1,split2);
                        segs.push_back(split2);
                    }
                }
            }
            if((!left_on_arc) && right_on_arc) {
                if(!it_seg.is_finite(CGAL::ARR_MIN_END) ||
                   ! equal(left(it_seg),end.get())) {
                    if(it_seg.is_finite(CGAL::ARR_MAX_END) && equal(end.get(),right(it_seg))) {
                        segs.push_back(it_seg);
                    } else {
                        X_monotone_curve_2 split1,split2;
                        it_seg.split(end.get(),split1,split2);
                        segs.push_back(split1);
                    }
                }
            } 
            if( (!left_on_arc) && (!right_on_arc)) {
                segs.push_back(it_seg);
            }
            helper=it; // for later usage
            
            if(! left_on_arc) {
                
                Point_2 point_it;
                while(true) {
                    if(it_seg.is_finite(CGAL::ARR_MIN_END) &&
                       is_one_one(cv,left(it_seg) ) ) {
                        point_it=left(it_seg);
                    } else {
                        CGAL_assertion(! start);
                        break;
                    }
                    CGAL_assertion(it!=arcs.begin());
                    it--;
                    CGAL::assign(it_seg,*it);
                    while(! on_arc(point_it,it_seg)) {
                        CGAL_assertion(it!=arcs.begin());
                        it--;
                        CGAL::assign(it_seg,*it);
                    }
                    if(start && on_arc(start.get(),it_seg)) {
                        segs.push_front(it_seg.trim(start.get(),
                                                    right(it_seg)));
                        break;
                    } else {
                        segs.push_front(it_seg);
                    }
                }
            }
            if(! right_on_arc) {
                it=helper; // reset
                CGAL::assign(it_seg,*it);
                Point_2 point_it;
                while(true) {
                    if(it_seg.is_finite(CGAL::ARR_MAX_END) &&
                       is_one_one(cv,right(it_seg) ) ) {
                        point_it=right(it_seg);
                    } else {
                        CGAL_assertion(! end);
                        break;
                    }
                    it++;
                    CGAL_assertion(it!=arcs.end());
                    CGAL::assign(it_seg,*it);
                    while(! on_arc(point_it,it_seg)) {
                        it++;
                        CGAL_assertion(it!=arcs.end());
                        CGAL::assign(it_seg,*it);
                    }
                    if(end && on_arc(end.get(),it_seg)) {
                        segs.push_back(it_seg.trim(left(it_seg),end.get()));
                        break;
                    } else {
                        segs.push_back(it_seg);
                    }
                }
            }
                    
            std::copy(segs.begin(),segs.end(),out);
            return out;
            
        }

    public:
        
        template<typename OutputIterator> OutputIterator
            operator() (Curve_2 cv, 
                        Point_2 p,
                        Site_of_point site_of_p,
                        OutputIterator out) const {
            if(site_of_p==POINT_IN_INTERIOR) {
                return x_monotone_segment(cv,p,boost::none, boost::none,out);
            } else if(site_of_p==MIN_ENDPOINT) {
                return x_monotone_segment(cv,
                                       p,
                                       boost::optional<Point_2>(p), 
                                       boost::none,
                                       out);
            }
            CGAL_assertion(site_of_p==MAX_ENDPOINT);
            return x_monotone_segment(cv,
                                   p,
                                   boost::none,
                                   boost::optional<Point_2>(p), 
                                   out);
        }


        template<typename OutputIterator> OutputIterator
            operator() (Curve_2 cv, 
                        Point_2 end_left,
                        Point_2 end_right,
                        OutputIterator out) const {
                    
            

            CGAL_assertion((branch_numbers_and_vertical(cv,end_left).
                            first.second==1 &&
                            !branch_numbers_and_vertical(cv,end_left).second)||
                           (branch_numbers_and_vertical(cv,end_right).
                            first.first==1 &&
                            !branch_numbers_and_vertical(cv,end_right).second)
                           ||
                           (branch_numbers_and_vertical(cv,end_left).
                            first.second==0 &&
                            branch_numbers_and_vertical(cv,end_left).second) ||
                           (branch_numbers_and_vertical(cv,end_right).
                            first.first==0 &&
                            branch_numbers_and_vertical(cv,end_right).second));
            CGAL_assertion_code(
                    if((branch_numbers_and_vertical(cv,end_left).
                        first.second==0 &&
                        branch_numbers_and_vertical(cv,end_left).second) ||
                       (branch_numbers_and_vertical(cv,end_right).
                        first.first==0 &&
                        branch_numbers_and_vertical(cv,end_right).second)) {
                        Compare_x_2 equal_x 
                            = this->_ckva()->compare_x_2_object(); )
                    CGAL_assertion(equal_x(end_left,end_right)==CGAL::EQUAL);
            CGAL_assertion_code(});
            
            
            if( (branch_numbers_and_vertical(cv,end_left).
                 first.second==0 &&
                 branch_numbers_and_vertical(cv,end_left).second) ||
                (branch_numbers_and_vertical(cv,end_left).
                 first.second==1 &&
                 !branch_numbers_and_vertical(cv,end_left).second)) {
                
                return x_monotone_segment
                    (cv,
                     end_left,
                     boost::optional<Point_2>(end_left),
                     boost::optional<Point_2>(end_right),
                     out);
            } else {
                return x_monotone_segment
                    (cv,
                     end_right,
                     boost::optional<Point_2>(end_left),
                     boost::optional<Point_2>(end_right),
                     out);
            }
        }

        template<typename OutputIterator> OutputIterator
	  operator() (Point_2 p, Point_2 q, OutputIterator out) {
	  bool same_x=(this->_ckva()->compare_x_2_object()(p,q)==CGAL::EQUAL);
	  if(same_x) {
            Algebraic_real_1 x = p.x();
            Polynomial_2 f = Polynomial_2(x.polynomial());
            Curve_2 cv 
	      = this->_ckva()->kernel().construct_curve_2_object() (f);
            *out++=typename CKvA_2::Arc_2(p,q,cv);
	    return out;
	  }
	  Algebraic_real_1 px
	    =this->_ckva()->kernel().compute_x_2_object()(p.xy()),
	    py=this->_ckva()->kernel().compute_y_2_object()(p.xy()), 
	    qx=this->_ckva()->kernel().compute_x_2_object()(q.xy()), 
	    qy=this->_ckva()->kernel().compute_y_2_object()(q.xy());
	  bool p_rat=px.is_rational() && py.is_rational();
	  bool q_rat=qx.is_rational() && qy.is_rational();
	  if(p_rat && q_rat) {
	    typedef typename Algebraic_real_1::Rational Rational;
	    Rational pxr=px.rational(), pyr=py.rational(),
	      qxr=qx.rational(), qyr=qy.rational();
	    typedef CGAL::Fraction_traits<Rational> FT;
	    typename FT::Decompose decompose;
	    typedef typename FT::Numerator_type Numerator;
	    typedef typename FT::Denominator_type Denominator;
	    Rational term_at_y=qxr-pxr, term_at_x=-qyr+pyr,
	      term_at_1=pxr*qyr-pyr*qxr;
	    Denominator denom_curr;
	    Numerator term_at_y_int, term_at_x_int, term_at_1_int;
	    decompose(term_at_y, term_at_y_int, denom_curr);
	    term_at_x=term_at_x*denom_curr;
	    term_at_1=term_at_1*denom_curr;
	    decompose(term_at_x,term_at_x_int,denom_curr);
	    term_at_y_int=term_at_y_int*denom_curr;
	    term_at_1=term_at_1*denom_curr;
	    decompose(term_at_1,term_at_1_int,denom_curr);
	    term_at_y_int=term_at_y_int*denom_curr;
	    term_at_x_int=term_at_x_int*denom_curr;
	    typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
	      ::Coefficient_type Polynomial_1;
	    Polynomial_2 pol(Polynomial_1(term_at_1_int,term_at_x_int),
			     Polynomial_1(term_at_y_int));
	    Curve_2 curve=this->_ckva()->kernel().construct_curve_2_object()
	      (pol);
            
	    CGAL_assertion(this->_ckva()->is_on_2_object()(p,curve));
	    CGAL_assertion(this->_ckva()->is_on_2_object()(q,curve));
	    return this->operator()(curve,p,q,out);
	  }
	  CGAL_precondition(same_x || (p_rat && q_rat));
	  return out;
	}
    };
    
    Construct_x_monotone_segment_2 construct_x_monotone_segment_2_object() 
        const {
        return Construct_x_monotone_segment_2(&CKvA_2::instance());
    }

/*
    
    class Construct_vertical_segment_2 : public 
    CGAL::internal::Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< CKvA_2 > {
        
    public:

        typedef CGAL::internal::Curved_kernel_via_analysis_2_Functors::
        Curved_kernel_via_analysis_2_functor_base< CKvA_2 > Base;

        Construct_vertical_segment_2(CKvA_2* kernel) : Base(kernel) {}


    public:

        X_monotone_curve_2 operator() (Point_2 p, Point_2 q) {
            Algebraic_real_1 x = p.x();
            Polynomial_2 f = Polynomial_2(x.polynomial());
            Curve_2 cv 
                = this->_ckva()->kernel().construct_curve_2_object() (f);
            return typename CKvA_2::Arc_2(p,q,cv);
        }

        
        template<typename OutputIterator> OutputIterator
            operator() (Point_2 p,
                        Site_of_point site_of_p) const {
            Algebraic_real_1 x = p.x();
            Polynomial_2 f = Polynomial_2(x.polynomial());
            Curve_2 cv = this->_ckva().kernel().construct_curve_2_object() (f);

            if(site_of_p==POINT_IN_INTERIOR) {
                return typename CKvA_2::Arc_2(x,cv);
            } else if(site_of_p==MIN_ENDPOINT) {
                return typename CKvA_2::Arc_2(p,ARR_MIN_END,cv);
            }
            CGAL_assertion(site_of_p==MAX_ENDPOINT);
            return typename CKvA_2::Arc_2(p,ARR_MAX_END,cv);
                            
        }
    };

    Construct_vertical_segment_2 construct_vertical_segment_2_object() const {
        return Construct_vertical_segment_2(&CKvA_2::instance());
    }
  
*/      

    class Construct_curve_2 : public 
    CGAL::internal::Curved_kernel_via_analysis_2_Functors::
    Curved_kernel_via_analysis_2_functor_base< CKvA_2 > {

    public:

        typedef CGAL::internal::Curved_kernel_via_analysis_2_Functors::
        Curved_kernel_via_analysis_2_functor_base< CKvA_2 > Base;

        Construct_curve_2(CKvA_2* kernel) : Base(kernel) {}

        Curve_2 operator() (Polynomial_2 p) const {
            return this->_ckva()->kernel().construct_curve_2_object() (p);
        }
    };

    Construct_curve_2 construct_curve_2_object() const {
        return Construct_curve_2(&CKvA_2::instance());
    }


    

    
/*

  // additional functionality (for not introducing a "general" arc)

  class Connect_points_2 {

    
    template< class OutputIterator >
    OutputIterator operator() (Curve_2, Point_2, Point_2, ..., 
                               OutputIterator) {

      
      while (false || true) {
        *oi++ = X_monotone_curve_2(...);
      }
      
      return oi;
    }

    template< class OutputIterator >
    OutputIterator operator() (Curve_2, Point_2, Point_2, int, ...,
                               OutputIterator) {
      
      
      while (false || true) {
        *oi++ = X_monotone_curve_2(...);
      }
      
      return oi;
    }

  };

*/
 
};

} //namespace CGAL

#endif // CGAL_ARR_ALGEBRAIC_SEGMENT_TRAITS
