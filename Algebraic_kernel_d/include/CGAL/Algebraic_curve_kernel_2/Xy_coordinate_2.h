// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//				   Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurveKernel_2, class Rep_ > 
class Xy_coordinate_2;

namespace CGALi {

template < class AlgebraicCurveKernel_2 >
class Xy_coordinate_2_rep {

public:
	// this first template argument
	typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

	// myself
    typedef Xy_coordinate_2_rep<Algebraic_curve_kernel_2> Self;

	typedef typename Algebraic_curve_kernel_2::Curve_pair_2 Curve_pair_2;

	typedef typename Curve_pair_2::Algebraic_curve_2 Curve_2; 
	
	typedef typename Curve_2::X_coordinate X_coordinate_1; 

	// constructors
public:
    // default constructor ()
    Xy_coordinate_2_rep() : arcno_(-1)
	{   }
    
    // standard constructor
    Xy_coordinate_2_rep(const X_coordinate_1& x, const Curve_2& curve, 
				int arcno) : x_(x), curve_(curve), arcno_(arcno)
	{   }

    // data
    // x-coordinate
    X_coordinate_1 x_;
    
    // supporting curve
    mutable Curve_2 curve_;
    
    // arc number on curve
    mutable int arcno_;

    // befriending the handle
    friend class Xy_coordinate_2<Algebraic_curve_kernel_2, Self>;
};

//! \brief class \c Xy_coordinate_2 represents a single root of a system of 
//! two polynomial equations in two variables that are models 
//! \c AlgebraicCurveKernel_2::Polynomial_2
//!
//! \c Xy_coordinate_2 coordinate is represented by an \c X_coordinate_1,
//! a supporting curve and an arcno and is valid only for finite solutions,
//! i.e., it cannot represent points at infinity 
template <class AlgebraicCurveKernel_2, 
		  class Rep_ = CGALi::Xy_coordinate_2_rep<AlgebraicCurveKernel_2> >
class Xy_coordinate_2 : public ::CGAL::Handle_with_policy< Rep_ > 
{
public:

    //! this instance's first template parameter
	typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
	typedef Rep_ Rep;

    //! this instance itself
    typedef Xy_coordinate_2<Algebraic_curve_kernel_2, Rep> Self;
	
	//! type of a curve pair 
	typedef typename Algebraic_curve_kernel_2::Curve_pair_2 Curve_pair_2;

	//! type of a curve analysis (replacement to CCPA_2)
	typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
		Curve_analysis_2;

	//! type of an algabraic curve
	typedef typename Curve_pair_2::Algebraic_curve_2 Curve_2; 

	//! type of X_coordinate
	typedef typename Curve_2::X_coordinate X_coordinate_1;

	//! type of curve pair analysis
	typedef typename Algebraic_curve_kernel_2::Curve_pair_analysis_2
				Curve_pair_analysis_2;
	
	//! type of pair vertical line
	typedef typename Curve_pair_analysis_2::Curve_pair_vertical_line_1
				Curve_pair_vertical_line_1;

	//! type of curve analysis
	typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
				Curve_analysis_2;
	
	//! type of vertical line
	typedef typename Curve_analysis_2::Curve_vertical_line_1
				Curve_vertical_line_1;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
    
public:
    // Rebind
    /*!\brief
     * An auxiliary structure for rebinding the point with a new rep
     */
    template < typename NewRep >
    class Rebind_curve_point_2
    {
    public:
        typedef CGALi::Xy_coordinate_2<Algebraic_curve_kernel_2, NewRep> 
			Other;
        
        /*!\brief
         * constructs a point of type \c Other from the point \c pt 
         * of type \c Self.
         *
         * All known items of the base class rep will be copied.
         */
        Other operator()(const Self& pt) {
            NewRep newrep;
            newrep.x_ = pt.ptr()->x_;
            newrep.curve_ = pt.ptr()->curve_;
            newrep.arcno_ = pt.ptr()->arcno_;
            return Other(newrep);
        }
    };
    
private:
    /*!\brief
     * Simplifies the representation of two points whose supporting curves
     * share a common part.
     */
    inline 
    static bool simplify(const Xy_coordinate_2& p, const Xy_coordinate_2& q) 
    {
		std::vector<Curve_2> parts_of_f, parts_of_g, common;
		Algebraic_curve_kernel_2 ak_2;

        if(ak_2.decompose_2_object()(p.curve(), q.curve(), 
			std::back_inserter(parts_of_f), std::back_inserter(parts_of_g),
			std::back_inserter(common))) {
            CGAL_assertion(static_cast<int>(parts_of_f.size()) == 2 ||
                       static_cast<int>(parts_of_g.size()) == 2);
			
			// ATTENTION: here the cache must be used !!
            if (static_cast< int >(parts_of_f.size()) == 2) {
                p.simplify_by(Curve_pair_analysis_2(
					Curve_pair_2(parts_of_f[0], parts_of_f[1])));
			} // else nothing to replace
            
            if (static_cast< int >(parts_of_g.size()) == 2) {
                q.simplify_by(Curve_pair_analysis_2(
					Curve_pair_2(parts_of_g[0], parts_of_g[1])));
            } // else nothing to replace
            return true;
        }
        return false;
    }
    
public:
    //!\name Constructors
    //!@{

    /*!\brief 
     * default constructor
     *
     * A default-constructed point supports no operation other than
     * having \c curve().degree() return \c -1. 
     */
    Xy_coordinate_2() : Base(Rep()) 
    {  }

    /*!\brief
     * copy constructor
     */
    Xy_coordinate_2(const Self& p) : Base(static_cast<const Base&>(p)) 
	{  }
    
    /*!\brief
     * Point at \c x, on \c curve with \c arcno. Finite points on vertical arcs
     * are also constructed in this way
	 * 
	 * \pre y-coordinate of this point must be finite
     */

	// from outside Curve_analysis_2 you must not directly use the methods
	// of Curve_2 class - only through available Curve_analysis_2 interface
    Xy_coordinate_2(const X_coordinate_1& x, const Curve_2& curve, int arcno) :
			Base(Rep(x, curve, arcno)) 
    {
        CGAL_precondition(arcno >= 0);
        CGAL_precondition_code(
			int i = -1;
			Curve_analysis_2 ca(curve);
			typename Curve_analysis_2::Curve_vertical_line_1 v = 
				ca.vertical_line_for_x(x);
		);
        CGAL_precondition(arcno < v.number_of_events());
    }
    
    /*!\brief
     * constructs a point from a given represenation
     */
    Xy_coordinate_2(Rep rep) : Base(rep) 
    {  }
    
    //!@}
    
public:
    //!\name Destructors
    //!@{
    
    /*!\brief 
     * Empty desctructor
     */
    virtual ~Xy_coordinate_2() {
        // empty
    }
    
    //!@}
    
public:
    //!\name Access functions
    //!@{
    
    /*!\brief 
     * x-coordinate of the point
     */
    X_coordinate_1 x() const { 
        return this->ptr()->x_; 
    }
    
    /*!\brief
     * supporting curve of the point
     */
    Curve_2 curve() const { 
        return this->ptr()->curve_; 
    }
    
    /*!\brief
     * arc number of point
     *
     */
    int arcno() const { 
        return this->ptr()->arcno_; 
    }

    //!@}
    //!\name Traits modifiers
    //!@{
protected:
    
    /*!\brief
     * test wether comparison of y_order of two covertical points is known
     * in advance
     */
    inline
#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    virtual 
#endif
    bool knows_y_order(const Self& p) const {
        //std::cout << "ACP: knows_y_order" << std::endl;
        return false;
    }

#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    /*!\brief
     * tests wether comparison of y-order of two covertical points is known
     * in advance
     */
    template < class Point_2 >
    inline bool knows_y_order(const Point_2& p) const {
        //std::cout << "ACPf: knows_y_order" << std::endl;
        if (dynamic_cast<const Point_2*>(this)) {
            return p.knows_y_order(*(dynamic_cast<const Point_2*>(this)));
        }
        return knows_y_order(static_cast<Self>(p));
    }
#endif
    
    /*!\brief
     * computes known y-order of covertical points
     */
    inline
#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    virtual 
#endif
    CGAL::Comparison_result known_y_order(const Self& p) const {
        //std::cout << "ACP: known_y_order" << std::endl;
        return CGAL::EQUAL;
    }

#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    /*!\brief
     * computes known y-order of covertical points
     */
    template < class Point_2 >
    inline CGAL::Comparison_result known_y_order(const Point_2& t) const {
        //std::cout << "ACPf: known_y_order" << std::endl;
        if (dynamic_cast<const Point_2*>(this)) {
            return -t.known_y_order(*(dynamic_cast<const Point_2*>(this)));
        }
        return known_y_order(static_cast<Self>(t));
    }
#endif
    
    /*!\brief
     * tests wether y_order of covertical points need to be reversed
     */
    inline
#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    virtual 
#endif
    bool reverse_y_order(const Self& p) const {
        //std::cout << "ACP: reverse_y_order" << std::endl;
        return false;
    }
    
#ifndef SoX_GAPS_NO_VIRTUAL_DISPATCH
    /*!\brief
     * tests wether y_order of covertical points needs to be reversed
     */
    template < class Point_2 >
    inline bool reverse_y_order(const Point_2& t) const {
        //std::cout << "ACPf: reverse_y_order" << std::endl;
        if (dynamic_cast<const Point_2*>(this)) {
            return t.reverse_y_order(*(dynamic_cast<const Point_2*>(this)));
        }
        return reverse_y_order(static_cast<Self>(t));
    }
#endif
    
    //!@}

    //!\name Comparisons
    //!@{
public:
    /*!\brief
     * compares x-coordinates of \c *this with \c q
     * 
     * do we need this method or one should use Algebraic_curve_kernel_2
	 * directly ?
     */
    template < class Point_2 >
    CGAL::Comparison_result compare_x(const Point_2& q) const {
        if (this->is_identical(q)) {
            return CGAL::EQUAL;
        }
		Algebraic_curve_kernel_2 ak_2;
        return ak_2.compare_x_2_object()(this->x(), q.x());
    }

    /*!\brief
     * compares \c *this with \c q lexicographically
     */
    template < class Point_2 >
    CGAL::Comparison_result compare_xy(const Point_2& q, 
                                 bool equal_x = false) const 
    {
        if (this->is_identical(q)) {
            return CGAL::EQUAL;
        }
        if (!equal_x) {
            CGAL::Comparison_result c = this->compare_x(q);
            if (c != CGAL::EQUAL) {
                CGAL_assertion(c == CGAL::SMALLER || c == CGAL::LARGER);
                return c;
            } 
        }
        // else
        return this->compare_y_at_x(q);
    }

    //! equality
    template < class Point_2 >
    bool operator == (const Point_2& q) const {return q.compare_xy(*this)== 0;}
    
    //! inequality
    template < class Point_2 >
    bool operator != (const Point_2& q) const {return q.compare_xy(*this)!= 0;}

    //! less than in (x,y) lexicographic order
    template < class Point_2 >
    bool operator <  (const Point_2& q) const {return q.compare_xy(*this)> 0;}

    //! less-equal in (x,y) lexicographic order
    template < class Point_2 >
    bool operator <= (const Point_2& q) const {return q.compare_xy(*this)>= 0;}

    //! greater than in (x,y) lexicographic order
    template < class Point_2 >
    bool operator >  (const Point_2& q) const {return q.compare_xy(*this)< 0;}

    //! greater-equal in (x,y) lexicographic order
    template < class Point_2 >
    bool operator >= (const Point_2& q) const {return q.compare_xy(*this)<= 0;}
    
private:
    /*!\brief
     * compares y-coordinates for covertical points \c *this and \c q
     *
     * \pre x() == q.x()
     */
    template < class Point_2 >
    CGAL::Comparison_result compare_y_at_x(const Point_2& q) const 
	{
        CGAL_precondition(q.compare_x(*this) == CGAL::EQUAL);

        // known filter
        if (this->knows_y_order(q)) {
            return this->known_y_order(q);
        }

        Curve_2 f = this->curve();
        Curve_2 g = q.curve();
        if (!f.is_identical(g)) {  // common parts of curves?
            if (this->simplify(*this,q)) {
                // ask for predicate again
                // since this->curve() and this->curve() can be the 
                // equal now
                return this->compare_y_at_x(q);
            }
        } 
        if (f.is_identical(g)) {
            CGAL::Sign result = CGAL::sign(this->arcno() - q.arcno());
            if (reverse_y_order(q))  
				return -result;
        	return result;
        }
		// attention: here the cache should be used
		Curve_pair_vertical_line_1& vline = 
			Curve_pair_analysis_2(
				Curve_pair_2(f, g)).vertical_line_for_x(x());
        CGAL::Sign result = 
            CGAL::sign(vline.get_event_of_curve(0, this->arcno()) - 
					vline.get_event_of_curve(1, q.arcno()));
		if(reverse_y_order(q)) 
				return -result;
        return result;
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
    void simplify_by(const Curve_pair_analysis_2& pair) const 
	{ 
        CGAL_precondition_code(
			typename Curve_2::Poly_d mult =
					pair.get_curve_analysis(0).get_polynomial_2().f() *
					pair.get_curve_analysis(1).get_polynomial_2().f();
			typename CGAL::Polynomial_traits_d<typename Curve_2::Poly_d>::
				Total_degree total_degree;
        );
        // common parts
        CGAL_precondition(NiX::resultant(mult, curve().f()).is_zero());
        // full parts
        CGAL_precondition(mult.degree() == curve().f().degree());
        CGAL_precondition(total_degree(mult) == total_degree(curve().f()));

		Curve_pair_vertical_line_1& vpair = pair.vertical_line_for_x(x());
        // # of arcs must match
		CGAL_precondition_code(
			Curve_vertical_line_1& vline = Curve_analysis_2(curve()).
				vertical_line_for_x(x);
		);
        CGAL_precondition(vpair.number_of_events() == 
			vline.number_of_events());

        int cid = 0;
		std::pair<int, int> p = vpair.get_curves_at_event(arcno());

        if(p.first != -1 && p.second != -1) {
            // both curves involved: choose simpler one
            // Remark: In this case, a vertical line in the curves can be
            // ignored, since it has not been considered when constructing
            // the point from the composed curved (also including this vertical
            // line). Therefore, the old arc number is also valid in the curve
            // pair.
			typename CGAL::Polynomial_traits_d<typename Curve_2::Poly_d>::
				Total_degree total_degree;
			if(total_degree(pair.get_curve_analysis(0).get_polynomial_2().f())
				 > 
			total_degree(pair.get_curve_analysis(1).get_polynomial_2().f())) 
                cid = 1;
        } else 
			cid = (p.first == -1 ? 1 : 0);
        // overwrite data
        this->ptr()->curve_ = pair.get_curve_analysis(cid).get_polynomial_2();
		this->ptr()->arcno_ = (cid == 0 ? p.first : p.second);
    }

    //!@}

}; // class Xy_coordinate_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H
