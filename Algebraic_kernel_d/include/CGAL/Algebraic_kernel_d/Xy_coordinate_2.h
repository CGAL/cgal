// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H

#include <CGAL/basic.h>
#include <boost/numeric/interval.hpp>

#include <CGAL/Bbox_2.h>

#include <CGAL/Arithmetic_kernel.h>

namespace CGAL {

namespace internal {

template < class AlgebraicCurveKernel_2, class Rep_, 
      class HandlePolicy_ ,
      class Allocator_>
        //::boost::fast_pool_allocator<Rep_> >
class Xy_coordinate_2;


template < class AlgebraicCurveKernel_2 >
class Xy_coordinate_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Xy_coordinate_2_rep<Algebraic_curve_kernel_2> Self;

    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
        Curve_analysis_2;

    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 
        Algebraic_real_1;

    typedef CGAL::Bbox_2 Bbox_2;

    typedef CGAL::Handle_with_policy<Self>
        Xy_coordinate_2_inst;

    // constructors
public:
    // default constructor ()
    Xy_coordinate_2_rep() : _m_arcno(-1) {
    }
    
    // standard constructor
    Xy_coordinate_2_rep(const Algebraic_real_1& x,
                        const Curve_analysis_2& curve, int arcno) 
      : _m_kernel(curve.kernel()),_m_x(x), _m_curve(curve), _m_arcno(arcno) {
    }

    // data

    const Algebraic_curve_kernel_2* _m_kernel;

    // x-coordinate
    Algebraic_real_1 _m_x;
    
    // supporting curve
    mutable Curve_analysis_2 _m_curve;
    
    // arc number on curve
    mutable int _m_arcno;

    // y-coordinate
    mutable boost::optional< Algebraic_real_1 > _m_y;

    //! A bounding box for the given point
    mutable boost::optional< std::pair<double,Bbox_2> > _m_bbox_2_pair;

};

//! \brief class \c Xy_coordinate_2 represents a single root of a system of 
//! two polynomial equations in two variables that are models 
//! \c AlgebraicCurveKernel_2::Polynomial_2
//!
//! \c Xy_coordinate_2 coordinate is represented by an \c Algebraic_real_1,
//! a supporting curve and an arcno and is valid only for finite solutions,
//! i.e., it cannot represent points at infinity 
template <class AlgebraicCurveKernel_2, 
          class Rep_ = internal::Xy_coordinate_2_rep<AlgebraicCurveKernel_2>,
          class HandlePolicy_= CGAL::Handle_policy_union, 
          class Allocator_ = CGAL_ALLOCATOR(Rep_) >
class Xy_coordinate_2 : 
    public ::CGAL::Handle_with_policy<Rep_, HandlePolicy_, Allocator_> 
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;
    
    //! this instance's third template parameter
    typedef HandlePolicy_ Handle_policy;
    
    //! this instance's fourth template parameter
    typedef Allocator_ Allocator;

    //! this instance itself
    typedef Xy_coordinate_2<Algebraic_curve_kernel_2, Rep, Handle_policy,
        Allocator> Self;
        
    //! an instance of AlgebraicKernel_1
    typedef typename Algebraic_curve_kernel_2::Algebraic_kernel_d_1 
        Algebraic_kernel_d_1;
    
    typedef typename Algebraic_curve_kernel_2::Polynomial_1 Polynomial_1;
    typedef CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;

    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Polynomial_2;
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    //! type of (explicit) x- and y-coordinates
    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 
        Algebraic_real_1;

    //! Coefficient type
    typedef typename Algebraic_curve_kernel_2::Coefficient Coefficient;

    //! type of curve pair analysis
    typedef typename Algebraic_curve_kernel_2::Curve_pair_analysis_2
                Curve_pair_analysis_2;
    
    //! type of curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
                Curve_analysis_2;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep, Handle_policy, Allocator> Base;

    //! type for approximation boundaries
    typedef typename Algebraic_curve_kernel_2::Bound Bound;

    //! type for bound intervals
    typedef boost::numeric::interval<Bound> Bound_interval;

    //! Type for the bounding box
    typedef typename Rep::Bbox_2 Bbox_2;
    
    //!@}
private:
    //! \name private methods
    //!@{

    /*!\brief
     * Simplifies the representation of two points whose supporting curves
     * share a common part.
     */
    bool _simplify(const Xy_coordinate_2& p, const Xy_coordinate_2& q) const
    {
        std::vector<Curve_analysis_2> parts_of_f, parts_of_g, common;

        if(kernel()->decompose_2_object()(p.curve(), q.curve(), 
            std::back_inserter(parts_of_f), std::back_inserter(parts_of_g),
                std::back_inserter(common))) {

            CGAL_assertion((parts_of_f.size() == 1 ||
                       parts_of_g.size() == 1) && common.size() == 1);
            if(parts_of_f.size() == 1) {
                p.simplify_by(kernel()->construct_curve_pair_2_object()(
                    parts_of_f[0], common[0]));
            } 
            
            if(parts_of_g.size() == 1) {
                q.simplify_by(kernel()->construct_curve_pair_2_object()(
                    parts_of_g[0], common[0]));
            } 
            return true;
        }
        return false;
    }
    
    //!@}
public:
    //!\name Constructors
    //!@{

    /*!\brief 
     * default constructor
     *
     * A default-constructed point supports no operation other than
     * having \c CGAL::degree(curve()) return \c -1. 
     */
    Xy_coordinate_2() : 
        Base(Rep()) { 
    }

    /*!\brief
     * copy constructor
     */
    Xy_coordinate_2(const Self& p) : 
        Base(static_cast<const Base&>(p)) {  
    }

    /*!\brief
     * Point at \c x, on \c curve with \c arcno. Finite points on vertical arcs
     * are also constructed in this way
     */
    Xy_coordinate_2(const Algebraic_real_1& x, const Curve_analysis_2& curve,
                    int arcno) :
          Base(Rep(x, curve, arcno)) {
            
        CGAL_precondition(arcno >= 0);
        CGAL_precondition_code(
            typename Curve_analysis_2::Status_line_1 v =
                curve.status_line_for_x(x);
        );
        CGAL_precondition(arcno >= 0 && arcno < v.number_of_events());
    }
    
    /*!\brief
     * constructs a point from a given represenation
     */
    Xy_coordinate_2(Rep rep) : 
        Base(rep) {  
    }
   
    //!@}
public:
    //!\name Access functions
    //!@{
    
    /*!\brief 
     * x-coordinate of the point
     */
    const Algebraic_real_1& x() const { 
        return this->ptr()->_m_x; 
    }

    /*!
     * \brief y-coordinate of this point
     *
     * Note: In general, this method results in a extremly large polynomial
     * for the y-coordinate. It is recommended to use it carefully,
     * and using get_approximation_y() instead whenever approximations suffice.
     */
    Algebraic_real_1 y() const {

        typedef std::vector< Algebraic_real_1 > Roots;
        // EBEB 2012-07-05 deactivated map for y-roots for not being used
        // typedef typename Curve_analysis_2::Status_line_1 Key;
        // EBEB 2012-07-05 deactivated map for y-roots for not being used
        // typedef Roots Data;
        // EBEB 2012-07-05 deactivated map for y-roots for not being used
        //        typedef std::map< Key, Data, CGAL::Handle_id_less_than< Key > > 
        //    Y_root_map;
        
        // EBEB 2012-07-05 deactivated map for y-roots for not being used
        //static Y_root_map y_root_map;

        if (!this->ptr()->_m_y) {
            
            Polynomial_2 f = curve().primitive_polynomial_2();
            // This will be the defining polynomial of y
            Polynomial_1 y_pol;

            // Filter: If we know that the point is critical, we can use
            // the resultant of f and f_y with respect to x as polynomial
            bool point_is_certainly_critical = false;
            typename Curve_analysis_2::Status_line_1 line =
                curve().status_line_at_exact_x(x());
            
            // EBEB 2012-07-05 deactivated map for y-roots for not being used
            //typename Y_root_map::iterator yit = 
            //    y_root_map.find(line);

            // TODO: Cache resultant computation
            // exacus-related code shouldn't be used here
            //curve().x_to_index(x(),i,is_event);
            if (line.is_event()) {
                //typename Internal_curve_2::Event1_info ev_info =
                //   curve().event_info(i);
                typename Curve_analysis_2::Status_line_1::Arc_pair ipair =
                    line.number_of_incident_branches(arcno());
                
                if (ipair.first != 1 || ipair.second != 1) {
                    point_is_certainly_critical = true;
                    y_pol = CGAL::make_square_free(
                              CGAL::resultant
                                (typename Polynomial_traits_2::Swap() (f,0,1),
                                 typename Polynomial_traits_2::Swap() 
                                   (CGAL::differentiate(f),0,1))
                        );
                    // BUGFIX: y_pol might be zero:
                    if(y_pol.is_zero()) {
                        // force re-computation with bigger resultant
                        point_is_certainly_critical=false;
                    }                             
                    
                }
            }
            
            if (!point_is_certainly_critical) {
                
                Polynomial_2 r(x().polynomial());
                y_pol = CGAL::make_square_free(
			  CGAL::resultant
                            (typename Polynomial_traits_2::Swap() (f,0,1),
                             typename Polynomial_traits_2::Swap() (r,0,1))
                );
                
            }
            typename Algebraic_kernel_d_1::Solve_1 real_roots;
            
            Roots y_roots;
            real_roots(y_pol, std::back_inserter(y_roots), false ); 
            
            long prec = 16;
	    
	    typename Algebraic_curve_kernel_2::Approximate_absolute_y_2
	      approx_y=kernel()->approximate_absolute_y_2_object();
	    
	    std::pair<Bound,Bound> y_pair = approx_y(*this,prec);
	    
	    Bound_interval y_iv(y_pair.first,y_pair.second);
            
            typedef typename std::vector<Algebraic_real_1>::const_iterator
                Iterator;
            
            std::list< Iterator > candidates;
            
            for (Iterator it = y_roots.begin(); it != y_roots.end(); it++) {
                Bound_interval it_interval(it->low(), it->high());
                if (boost::numeric::overlap(it_interval, y_iv)) {
                    candidates.push_back(it);
                }
            }
            CGAL_assertion(!candidates.empty());

            while (candidates.size() > 1) {
	        prec*=2;
	        y_pair = approx_y(*this,prec);
	    
                y_iv = Bound_interval(y_pair.first,y_pair.second);

                for (typename std::list< Iterator >::iterator dit, cit =
                         candidates.begin(); cit != candidates.end(); ) {
                    bool remove = false;
                    Bound_interval 
                        cit_interval((*cit)->low(), (*cit)->high());
                    if (!boost::numeric::overlap(cit_interval, y_iv)) {
                        dit = cit;
                        remove = true;
                    }
                    cit++;
                    if (remove) {
                        candidates.erase(dit);
                    }
                }
            }
            CGAL_assertion(static_cast< int >(candidates.size()) == 1);
            this->ptr()->_m_y = 
                Algebraic_real_1(
                        (*candidates.begin())->polynomial(), 
                        (*candidates.begin())->low(), 
                        (*candidates.begin())->high()
                );
        }
        CGAL_postcondition(bool(this->ptr()->_m_y));
        return *this->ptr()->_m_y;
    }
    
    /*!\brief
     * supporting curve of the point
     */
    Curve_analysis_2 curve() const {
        return this->ptr()->_m_curve; 
    }
    
    /*!\brief
     * arc number of point
     *
     */
    int arcno() const { 
        return this->ptr()->_m_arcno; 
    }

    //!@}
public:
    //!\name comparison predicates
    //!@{

    /*!\brief
     * compares x-coordinates of \c *this with \c q
     * 
     * do we need this method or one should use Algebraic_curve_kernel_2
     * directly ?
     */
    CGAL::Comparison_result compare_x(const Self& q) const {
    
        if(this->is_identical(q)) {
            return CGAL::EQUAL;
        }
        return kernel()->compare_1_object()(this->x(), q.x());
    }

    /*!\brief
     * compares \c *this with \c q lexicographically
     */
    CGAL::Comparison_result compare_xy(const Self& q, 
        bool equal_x = false) const {
        
        if(this->is_identical(q)) 
            return CGAL::EQUAL;

        CGAL::Comparison_result res = (equal_x ? CGAL::EQUAL : compare_x(q)); 
        if(res == CGAL::EQUAL) {
            res = _compare_y_at_x(q);
        }
        return res;
    }
    
    //! equality
    bool operator == (const Self& q) const {return q.compare_xy(*this)== 0;}
    
    //! inequality
    bool operator != (const Self& q) const {return q.compare_xy(*this)!= 0;}

    //! less than in (x,y) lexicographic order
    bool operator <  (const Self& q) const {return q.compare_xy(*this)> 0;}

    //! less-equal in (x,y) lexicographic order
    bool operator <= (const Self& q) const {return q.compare_xy(*this)>= 0;}

    //! greater than in (x,y) lexicographic order
    bool operator >  (const Self& q) const {return q.compare_xy(*this)< 0;}

    //! greater-equal in (x,y) lexicographic order
    bool operator >= (const Self& q) const {return q.compare_xy(*this)<= 0;}
    
    //!@}

    //!@{
    //! \name 

    const Algebraic_curve_kernel_2* kernel() const {
        return this->ptr()->_m_kernel;
    }

private:

    /*!\brief
     * compares y-coordinates for covertical points \c *this and \c q
     *
     * \pre x() == q.x()
     */
    CGAL::Comparison_result _compare_y_at_x(const Self& q) const 
    {
        CGAL_precondition(this->compare_x(q) == CGAL::EQUAL);
    
        Curve_analysis_2 f = curve(), g = q.curve();
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - q.arcno());
        if(Self::_simplify(*this, q)) 
            // restart since supporting curves might be equal now
            return _compare_y_at_x(q);
                        
        Curve_pair_analysis_2 cpa_2 =
            kernel()->construct_curve_pair_2_object()(f, g);
            
            
        typename Curve_pair_analysis_2::Status_line_1 vline =
            cpa_2.status_line_for_x(x());
        return CGAL::sign(vline.event_of_curve(arcno(), f) -
                    vline.event_of_curve(q.arcno(), g));
    }
    
    //!@}
public:
    //!\name Reconstructing functions
    //!@{
    
    /*!\brief
     * Simplifies the representation of a point.
     * 
     * Given a decomposition of the point's supporting \c curve() into 
     * a pair of two curves \c pair, this function searches this point
     * in the curve pair and resets the curve and the arcno to this
     * found arc. It can happen, that both curves of the pair fit this 
     * condition (intersection of the two curves at this point), then it
     * chooses the simpler one (less total degree).
     *
     * \pre pair must be a decomposition of curve()
     */
    void simplify_by(const Curve_pair_analysis_2& cpa_2) const {
    
        CGAL_precondition_code(
            Polynomial_2 mult =
                    cpa_2.curve_analysis(0).polynomial_2() *
                    cpa_2.curve_analysis(1).polynomial_2();
        );
        // common parts
        CGAL_precondition(CGAL::resultant(mult, 
                                          curve().polynomial_2()).is_zero());
        // full parts
        CGAL_precondition(CGAL::degree(mult) == 
                          CGAL::degree(curve().polynomial_2()));
        CGAL_precondition(CGAL::total_degree(mult) ==
                          CGAL::total_degree(curve().polynomial_2()));

        typename Curve_pair_analysis_2::Status_line_1 cpv_line =
                cpa_2.status_line_for_x(x());
        // # of arcs must match
        CGAL_precondition_code(
            typename Curve_analysis_2::Status_line_1 cv_line =
                curve().status_line_for_x(x());
        );
        CGAL_precondition(cpv_line.number_of_events() == 
            cv_line.number_of_events());

        bool cid = false;
        std::pair<int, int> p = cpv_line.curves_at_event(arcno());
        if(p.first != -1 && p.second != -1) {
            // both curves involved: choose simpler one
            // Remark: In this case, a vertical line in the curves can be
            // ignored, since it has not been considered when constructing
            // the point from the composed curved (also including this vertical
            // line). Therefore, the old arc number is also valid in the curve
            // pair.
            Polynomial_2 ff = cpa_2.curve_analysis(0).polynomial_2(),
	                 gg = cpa_2.curve_analysis(1).polynomial_2();
            if(total_degree(ff) > total_degree(gg)) 
                cid = true;
        } else 
          cid = (p.first == -1);
        // overwrite data
        this->ptr()->_m_curve = cpa_2.curve_analysis(cid);
        this->ptr()->_m_arcno = (cid == 0 ? p.first : p.second);
    }
    
    //! befriending output iterator
   // friend std::ostream& operator << <>(std::ostream& os, const Self& pt);

    //!@}
public:
    
    //! Returns whether the x-coordinate equals zero
    bool is_x_zero() const {
      return CGAL::is_zero(this->ptr()->_m_x);
    }

    //! Returns whether the y-coordinate equals zero
    bool is_y_zero() const {
      CGAL::Sign lower_sign = CGAL::sign(this->lower_bound_y()),
                 upper_sign = CGAL::sign(this->upper_bound_y());
      if( lower_sign == CGAL::ZERO ||upper_sign == CGAL::ZERO) {
	if(lower_sign==upper_sign) { //both zero
	  return true;
	} else { // one zero, one not...isol interval is OPEN
	  return false;
	}
      } else if( lower_sign==upper_sign) { // zero not in isol interval
	return false;
      } else { // zero in interval, need to check
	Polynomial_1 constant_pol =
 	  CGAL::get_coefficient(curve().primitive_polynomial_2(),0);
	bool zero_is_root_of_local_pol 
  	  = kernel()->is_zero_at_1_object()(constant_pol,this->ptr()->_m_x);
        // Since we know that y_iv is an _isolating_ interval,
        // we can immediately return
        return zero_is_root_of_local_pol;
      }
    }

    // returns a double approximation of the point
    std::pair<double, double> to_double() const {

      typedef typename CGAL::Get_arithmetic_kernel<Bound>::Arithmetic_kernel 
        AT;
        typedef typename AT::Bigfloat_interval BFI; 

        long old_prec = get_precision(BFI());
        
        set_precision (BFI(), 53);

	// rely on double conversion of the x-type
        double double_x = CGAL::to_double(this->ptr()->_m_x);
        double double_y;


        if (this->lower_bound_y()==this->upper_bound_y()) {
            double_y = CGAL::to_double(convert_to_bfi(this->lower_bound_y()));
        } else if(is_y_zero()) {
            double_y = 0.;
        } else {
            while(CGAL::sign(this->lower_bound_y()) != 
                  CGAL::sign(this->upper_bound_y()) ) {
                this->refine_y();
            }
            long final_prec = set_precision(BFI(),get_precision(BFI())+4);
            
            BFI bfi = CGAL::hull(convert_to_bfi(this->lower_bound_y()), 
				 convert_to_bfi(this->upper_bound_y()));
            
            while( !singleton(bfi) &&  
                   get_significant_bits(bfi) < final_prec  ){
                this->refine_y();
                bfi = CGAL::hull(
                        convert_to_bfi(this->lower_bound_y()), 
                        convert_to_bfi(this->upper_bound_y()));
            }
            double_y 
                = CGAL::to_double((CGAL::lower(bfi)+ CGAL::upper(bfi)) / 2);
        }
        set_precision(BFI(),old_prec);
        return std::make_pair(double_x, double_y); 
    }

    public:

    
    void refine_y() const {
        this->curve().status_line_at_exact_x(this->x()).refine(this->arcno());
    }

    Bound lower_bound_y() const {
      return this->curve().status_line_at_exact_x(this->x()).
          lower_bound(this->arcno());
    }

    Bound upper_bound_y() const {
      return this->curve().status_line_at_exact_x(this->x()).
          upper_bound(this->arcno());
    }

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE

    void refine_x() const {
      this->x().refine();
    }

    Bound lower_bound_x() const {
      return this->x().low();
    }

    Bound upper_bound_x() const {
      return this->x().high();
    }

#endif


     // friend function to provide a fast hashing
    friend std::size_t hash_value(const Self& x) {
        return static_cast<std::size_t>(x.id());
    }

    //!@}

}; // class Xy_coordinate_2

template < class AlgebraicCurveKernel_2, class Rep> 
std::ostream& operator<< (std::ostream& os, 
    const Xy_coordinate_2<AlgebraicCurveKernel_2, Rep>& pt)
{
  switch (::CGAL::get_mode(os)) {
  case ::CGAL::IO::PRETTY: {
    os << "[x-coord: " << CGAL::to_double(pt.x()) << "; curve: " <<
      pt.curve().polynomial_2() << 
      "; arcno: " << pt.arcno() << "]\n";
    break;
  } 
  case ::CGAL::IO::BINARY:
    std::cerr << "BINARY format not yet implemented" << std::endl;
    break;
  default:
    // ASCII
    os << "Algebraic_real_xca_2(";
    os << pt.x();
    os << ",";
    os << pt.curve();
    os << ",";
    os << pt.arcno();
    os << ")";
  }
  return os;    
}

template < class AlgebraicCurveKernel_2, class Rep_ > 
std::istream& operator >> (
    std::istream& is,
    Xy_coordinate_2< AlgebraicCurveKernel_2, Rep_>& pt) {

  CGAL_precondition(CGAL::is_ascii(is));
  
  // this instance's first template argument
  typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
  
  // this instance's second template argument
  typedef Rep_ Rep;

  // myself
  typedef Xy_coordinate_2< Algebraic_curve_kernel_2, Rep > Xy_coordinate_2;
  
  typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 
    Algebraic_real_1;

  typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
    Curve_analysis_2;
  
  // x-coordinate
  Algebraic_real_1 x;
    
  // supporting curve
  Curve_analysis_2 curve;
    
  // arc number on curve
  int arcno;
  
  // read "Algebraic_real_xca_2("
  swallow(is, 'A');
  swallow(is, 'l');
  swallow(is, 'g');
  swallow(is, 'e');
  swallow(is, 'b');
  swallow(is, 'r');
  swallow(is, 'a');
  swallow(is, 'i');
  swallow(is, 'c');
  swallow(is, '_');
  swallow(is, 'r');
  swallow(is, 'e');
  swallow(is, 'a');
  swallow(is, 'l');
  swallow(is, '_');
  swallow(is, 'x');
  swallow(is, 'c');
  swallow(is, 'a');
  swallow(is, '_');
  swallow(is, '2');
  swallow(is, '(');
  
  
  // read values
  is >> x;
  swallow(is, ',');
  
  is >> curve;
  swallow(is, ',');
  
  is >> arcno;
  
  // read the ")
  swallow(is, ')'); 
  
  pt = Xy_coordinate_2(x, curve, arcno);
  
  return is;
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H
