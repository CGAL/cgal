// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_ALCIX_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_ALCIX_H



#include <vector>
#include <set>
#include <map>

#include <CGAL/basic.h>
#include <CGAL/Cache.h>
#include <CGAL/function_objects.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Arr_enums.h>

#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/alg_real_utils.h>
#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Bitstream_descartes_bfs.h>
#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Bitstream_descartes_traits_on_vert_line.h>

#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Non_generic_position_exception.h>

#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_curve_kernel_2/analyses/macros.h>



#include <CGAL/Algebraic_curve_kernel_2/Status_line_CA_1.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Event_line_builder.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/subresultants.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Shear_controller.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Shear_transformation.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Zero_resultant_exception.h>
#include <CGAL/Polynomial/sturm_habicht_sequence.h>

#include <CGAL/Algebraic_curve_kernel_2/analyses/shear.h>

CGAL_BEGIN_NAMESPACE

template<typename AlgebraicKernel_2, 
         typename Rep_>
class Curve_analysis_2;

namespace CGALi {

// \brief Representation class for algebraic curves.
template< typename AlgebraicKernel_2>
class Curve_analysis_2_rep {
    
public:
    //! this instance's template parameter
    typedef AlgebraicKernel_2 Algebraic_kernel_2;

    //! the class itself
    typedef Curve_analysis_2_rep Self;
    
    //! The handle class
    typedef Curve_analysis_2<Algebraic_kernel_2,Self> Handle;
    
    //protected:
public:

    typedef int size_type;
    
    typedef Handle Curve_analysis_2;

    CGAL_ACK_SNAP_ALGEBRAIC_CURVE_TYPEDEFS;

    typedef std::map< Boundary, Status_line_1 > 
    Vert_line_at_rational_map;
    
#if 1 // TODO using x() results in much slower running times (pre-precision)
    typedef 
    std::map< X_coordinate_1, 
              Status_line_1, 
              CGAL::Handle_id_less_than< X_coordinate_1 > >
    Vert_line_map;
#else
    typedef std::map< X_coordinate_1, Status_line_1 >
    Vert_line_map;
#endif
    
    //!\name Constructors
    //!@{

    //! Default constructor
    Curve_analysis_2_rep()
    {
    }
    
    //! Constructor with polynomial
    Curve_analysis_2_rep(Polynomial_2 poly, 
                         CGAL::CGALi::Degeneracy_strategy strategy) :  
        f(poly), degeneracy_strategy(strategy)
    {
    }
    
    //!@}
    
private:

    typedef CGALi::Event_line_builder<Handle> Event_line_builder;
    

    // Internal information struct about x-coordinates
    struct Event_coordinate_1 {
        X_coordinate_1 val;
        size_type mult_of_prim_res_root;
        size_type index_of_prim_res_root;
        size_type mult_of_content_root;
        size_type index_of_content_root;
        size_type mult_of_prim_lcoeff_root;
        size_type index_of_prim_lcoeff_root;
        boost::optional<Status_line_1> stack;
    };
    
    // Functor to get the X_coordinate of an Event_coordinate
    struct Val_functor {
        typedef Event_coordinate_1 argument_type;
        typedef X_coordinate_1 result_type;
        result_type operator() (argument_type event) const {
            return event.val;
        }
    };


    //! The object holding the information about events, as an optional
    mutable boost::optional<std::vector<Event_coordinate_1> > 
        event_coordinates;

    //! The polynomial defining the curve
    boost::optional<Polynomial_2> f;

    //! How degenerate situations are handled
    size_type degeneracy_strategy;

    /*!
     * \brief The polynomial without its content (the gcd of the coeffs).
     *
     * The content is the greatest common divisor of the coefficients of \c f
     * considered as polynomial <tt>y</tt>. \c The polynomial f_primitive is
     * \c f/cont(f). The corresponding curve is equal to the curve of \c f,
     * only without vertical line components.
     */
    mutable boost::optional<Polynomial_2> f_primitive;
    
    //! the polynomial containing all roots of the resultant of the primitive
    //! part of f and its y-derivative
    mutable boost::optional<Polynomial_1> resultant_primitive_f_fy;

    //! the polynomial containing all roots of the resultant of the primitive
    //! part of f and its x-derivative
    mutable boost::optional<Polynomial_1> resultant_primitive_f_fx;

    //! The Sturm-Habicht polynomials of f
    mutable boost::optional<std::vector<Polynomial_2> > sturm_habicht_primitive_f;

    //! The content of f
    mutable boost::optional<Polynomial_1> content;

    //! The object for building event lines
    //mutable boost::optional<Event_line_builder> event_line_builder;
    
    //! The non-working shear factors, as far as known
    mutable std::set<Integer> bad_shears;

    //! The already known shear factors
    mutable std::map<Integer,Handle> sheared_curves;

    //! Has the curve vertical line components
    mutable boost::optional<bool> has_vertical_component;

    //! The intermediate values
    mutable boost::optional<std::vector<boost::optional<Boundary> > > 
    intermediate_values;

    //! stores Y_values at rational coordinate
    mutable Vert_line_at_rational_map vert_line_at_rational_map;
    
    //! stores vert_lines
    mutable Vert_line_map vert_line_map;

    /**! \brief Information about whether arcs at +/- infty 
     *   are asymptotic to y=beta,
     *   or go to +/- infty also in y-direction
     */
    mutable boost::optional<std::vector<CGAL::Object> >
    horizontal_asymptotes_left, horizontal_asymptotes_right;

    //! friends
    friend class ::CGAL::Curve_analysis_2<Algebraic_kernel_2,Self>;

}; // class Curve_analysis_2_rep
} // namespace CGALi


/*!
 * \brief Algebraic curves of arbitrary degree. 
 *
 * TODO: Update doc
 * A realisation of SoX::AlgebraicCurve_2, defined in the \ref SoX_intro 
 * SweepX-library.\n
 * The object contains two vector of objects that describe
 * the curve at given x-values. One vector (called \c event_lines) stores
 * the data for x-values where singularities and x-extreme points can occur.
 * The second one (calles \c intermediate_lines) contains the information for
 * some rational value between two critical values. In other words, the
 * following invariant is satisfied:
 * \code intermediate_lines[i].x() <= event_lines[i].x() 
 * <= intermediate_lines[i+1].x() \endcode
 * for all valid values of \c i.\n 
 * The data at the x-values is represented by specialisations
 * of the AcX::Vert_line class. This class is itself derived from 
 * SoX::Event1_info. So, this framework is compatible with 
 * \ref SoX_GAPS "GAPS".
 */
template<typename AlgebraicKernel_2, 
  typename Rep_ 
   = CGALi::Curve_analysis_2_rep< AlgebraicKernel_2> 
>
class Curve_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ > {
  
public:
    //! this instance' first template parameter
    typedef AlgebraicKernel_2 Algebraic_kernel_2;
  
    //! this instance' second template parameter
    typedef Rep_ Rep;

private:
  
    //! The internal type for event coordinates
    typedef typename Rep::Event_coordinate_1 Event_coordinate_1;

    typedef typename Rep::Event_line_builder Event_line_builder;

    // Base class
    typedef ::CGAL::Handle_with_policy<Rep> Base;
    
    // This type
    typedef Curve_analysis_2<Algebraic_kernel_2,Rep> Self;
    
public:

    typedef typename Rep::size_type size_type;
    
    //! Needed by the concept
    typedef CGAL::Handle_id_less_than< Self > Less_than;

    typedef Self Curve_2;
    CGAL_ACK_SNAP_ALGEBRAIC_CURVE_TYPEDEFS;

    //! Traits type for Polynomial_2
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    typedef CGAL::Coercion_traits<Boundary, Coefficient> Coercion;
    typedef typename Coercion::Type Coercion_type;

    typedef CGAL::Polynomial< Coercion_type > Poly_coer_1;

#if DOXYGEN_RUNNUNG
    //! type for x-coordinates;
    typedef Algebraic_real X_coordinate_1;

    //! type for x-coordinates;
    typedef Algebraic_real Y_coordinate_1;
#endif

    typedef typename Algebraic_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

    //! type of horizontal asymtote values
    typedef CGAL::Object Asymptote_y;

    //! type of Event1_info
    typedef Status_line_1 Event1_info;

private:

    // Helping structs
    
    struct Event_functor {
        Event_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Status_line_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->status_line_at_event(index);
        }
    };

    struct Intermediate_functor {
        Intermediate_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Status_line_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->status_line_of_interval(index);
        }
    };

    struct Stha_functor {
        Stha_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Polynomial_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->principal_sturm_habicht_primitive_f(index);
        }
    };

public:

    //! Type for an iterators over the Event-lines
    typedef boost::transform_iterator<Event_functor, 
                              boost::counting_iterator<size_type> > 
    Event_line_iterator;

    //! Type for an iterators over the Intermediate-lines
    typedef boost::transform_iterator<Intermediate_functor, 
                              boost::counting_iterator<size_type> > 
    Intermediate_line_iterator;

    typedef boost::transform_iterator<Stha_functor, 
                              boost::counting_iterator<size_type> > 
    Principal_sturm_habicht_iterator;

public:

    //!\name Constructors
    //!@{  
      
    //! \brief Default constructor
    Curve_analysis_2() :Base(Rep()) {
    }

    /*! 
     * \brief Constructor with polynomial.
     *
     * Analyses the curve given by the polynomial. If the parameter is 
     * ommitted, the polynomial <tt>f=1</tt> is taken as standard argument.
     */
    Curve_analysis_2(Polynomial_2 f) 
        throw(CGALi::Zero_resultant_exception<Polynomial_2>)
        : Base(Rep(f,CGAL::CGALi::SHEAR_STRATEGY))
    {

    }

    /*! 
     * \brief Constructor with polynomial and degeneracy strategie.
     *
     * TODO: doc
     */
    Curve_analysis_2(Polynomial_2 f, 
                     CGAL::CGALi::Degeneracy_strategy strategy) 
        throw(CGALi::Zero_resultant_exception<Polynomial_2>)
        : Base(Rep(f,strategy))
    {

    } 

    //! \brief Copy constructor
    Curve_analysis_2(const Self& alg_curve)
        : Base(static_cast<const Base&>(alg_curve)) 
    {
    }


    //!@}
      
public:

    /*
     * \brief Creates data using event-line-iterator
     *
     * Creates a Algebraic_curve_2 object with the given defining polynomial
     * and the event lines from the iterator range
     */
    template<typename InputIterator1,typename InputIterator2>
    void set_event_lines(InputIterator1 event_begin,
                         InputIterator1 event_end,
                         InputIterator2 intermediate_begin,
                         InputIterator2 intermediate_end) const {
        
        if(! this->ptr()->event_coordinates) {
            
            std::vector<Event_coordinate_1> event_coordinate_vector;
            
            for(InputIterator1 it = event_begin; it != event_end; it++) {
                Event_coordinate_1 curr_event;
                curr_event.val = it->x();
                event_coordinate_vector.push_back(curr_event);
            }
            this->ptr()->event_coordinates = event_coordinate_vector;

        }

        InputIterator1 it1 = event_begin;
        for(size_type i = 0; i < num_events() ; i++ ) {
            this->ptr()->vert_line_map[event_coordinates()[i].val] = *it1; 
            event_coordinates()[i].stack = *it1;

            it1++;
        }
        CGAL_assertion(it1 == event_end);

        if(! this->ptr()->intermediate_values) {
            this->ptr()->intermediate_values 
                = std::vector<boost::optional<Boundary> >(num_events()+1);
        }

        InputIterator2 it2 = intermediate_begin;
        for(size_type i = 0; 
            i < static_cast<int>(intermediate_values().size()); 
            i++,it2++) {
            
            CGAL_assertion(it2->x().is_rational());
            Boundary q = it2->x().rational();
            
            intermediate_values()[i] = q;
            
            this->ptr()->vert_line_map[it2->x()] = *it2;
            this->ptr()->vert_line_at_rational_map[q] = *it2;
            
        }
        CGAL_assertion(it2 == intermediate_end);
        
    }
        
public:
    
    /*! \brief Sets the polynomial.
     *
     * Only possible when the object has no initialised polynomial yet.
     */
    void set_f(Polynomial_2 f) {
        CGAL_precondition(! has_defining_equation());
        if((! this->ptr()->f) || f!=this->ptr()->f.get()) {
            this->copy_on_write();
            this->ptr()->f=f;
        }
    }

public:

    /*! \brief Returns whether the curve has a valid polynomial
     */
    bool has_defining_equation() const {
        return this->ptr()->f;
    }

public:

    //! \brief Returns whether the curve is y-regular
    bool is_y_regular() {
        return f().lcoeff().degree() == 0;
    }
    
public:

    /*!\brief
     * returns \c true iff curve has vertical component
     */
    bool has_vertical_component() const {
        if(is_y_regular()) {
	    this->ptr()->has_vertical_component = false;
        }
        if(! this->ptr()->has_vertical_component) {
            // This is computed as side effect 
            // when the event coordinates are computed
            event_coordinates();
            CGAL_assertion(this->ptr()->has_vertical_component);
        }
        return this->ptr()->has_vertical_component.get();
    }

public:

    //! \brief Defining polynomial
    Polynomial_2 polynomial_2() const {
        CGAL_precondition(this->ptr()->f);
        return this->ptr()->f.get();
    }

public:
    //! Shortcut for convenience
    Polynomial_2 f() const {
        return polynomial_2();
    }

public:

    //! \brief Returns the number of event lines
    size_type number_of_status_lines_with_event() const {
        CGAL_precondition(this->ptr()->f);
        return static_cast<size_type>(event_coordinates().size());
    }
      
public:
    //! Shortcut for convenience
    size_type num_events() const {
        return number_of_status_lines_with_event();
    }

/* TODO: Remove completely

public:

    //! \brief Returns the X_coorinate of the ith event
    X_coordinate_1 event_x(size_type i) const {
        return event_coordinates()[i].val;
    }

*/
      
public:

    /*! 
     * \brief Gives the index of an event line over a given 
     * \c Algebraic_real, or the interval it falls in.
     */
    void x_to_index(X_coordinate_1 x,size_type& i,bool& is_event) const {
        CGAL_precondition(has_defining_equation());
        typename Rep::Val_functor xval;
        i = std::lower_bound(
                ::boost::make_transform_iterator(event_coordinates().begin(), 
                                                 xval),
                ::boost::make_transform_iterator(event_coordinates().end(),
                                                 xval),
                x
        ) - ::boost::make_transform_iterator(event_coordinates().begin(), 
                                             xval);
        is_event = (i < static_cast<size_type>(event_coordinates().size()) && 
                    (event_coordinates()[i].val == x) );
    }

public:

    //! \brief Returns the ith event line.
    Status_line_1& status_line_at_event(size_type i) const {
        CGAL_precondition(has_defining_equation());
        size_type n = static_cast<size_type>(event_coordinates().size());
        (void)n;
        CGAL_precondition(i>=0 && i<n);
        if(! event_coordinates()[i].stack) {
            Status_line_1 event_line = create_status_line_at_event(i);
            this->ptr()->vert_line_map[event_coordinates()[i].val] 
                = event_line; 
            event_coordinates()[i].stack = event_line;
        }
        return event_coordinates()[i].stack.get();
    }
    
public:    

    /*!
     * \brief Returns a vert line for the rational <tt>x</tt>-coordinate b
     *
     * If a vert-line object exists for the specified value, it is
     * simply returned. Otherwise, it is newly created.
     */
    Status_line_1& status_line_at_exact_x(Boundary b) const {
        return status_line_at_exact_x(X_coordinate_1(b));
    }

private:

    Status_line_1& status_line_at_exact_non_event_x(X_coordinate_1 alpha) 
        const {

        if(alpha.is_rational()) {
            
            typename Rep::Vert_line_at_rational_map::iterator it =
                this->ptr()->vert_line_at_rational_map.find
                (alpha.rational());
            
            if (it != this->ptr()->vert_line_at_rational_map.end()) {
                return it->second;
            }
        }
        
        typename Rep::Vert_line_map::iterator it =
            this->ptr()->vert_line_map.find(alpha);
        
        if (it != this->ptr()->vert_line_map.end()) {
            return it->second;
        }
        
        
        // Not stored yet, so create it and store it
        Status_line_1 cvl 
            = create_status_line_at_non_event(alpha);
        this->ptr()->vert_line_map[alpha] = cvl;
        
        if(alpha.is_rational()) {
            this->ptr()->vert_line_at_rational_map[alpha.rational()] = cvl;
        }
        return this->ptr()->vert_line_map[alpha];
    }

public:

    /*!
     * \brief Returns a vert line for the <tt>x</tt>-coordinate alpha
     *
     * If a vert-line object exists for the specified value, it is
     * simply returned. Otherwise, it is newly created.
     */
    Status_line_1& status_line_at_exact_x(X_coordinate_1 alpha) const {
        bool is_event_value;
        size_type index;
        this->x_to_index(alpha,index,is_event_value);
        if(is_event_value) {
            return status_line_at_event(index);
        }
        else {
            return status_line_at_exact_non_event_x(alpha);
        }
    }

private:
    
    Status_line_1 create_status_line_at_event(size_type index) const 
        throw(CGAL::CGALi::Non_generic_position_exception) {

        try {
            
            Event_coordinate_1& event = event_coordinates()[index];
        
            X_coordinate_1 x = event.val;
            
#if AcX_SHEAR_ALL_NOT_Y_REGULAR_CURVES
            if(event.mult_of_prim_lcoeff_root > 0) {
                throw CGAL::CGALi::Non_generic_position_exception();
            }
#else
            if(event.mult_of_prim_lcoeff_root > 0) {
                if(event.mult_of_prim_lcoeff_root > 1 ||
                   event.mult_of_prim_res_root > 1) {
                    throw CGAL::CGALi::Non_generic_position_exception();
                }
            }
        
#endif
        
#if AcX_DEBUG_PRINT
            double ev_approx = CGAL::to_double(x);
            AcX_DSTREAM((index+1) << "th line: "
                        << std::setw(6) << std::setprecision(3)
                        << ev_approx
                        << ".."
                        << std::flush);
#endif	
            size_type left_arcs 
                = status_line_for_x(x,CGAL::NEGATIVE).number_of_events();
            size_type right_arcs 
                = status_line_for_x(x,CGAL::POSITIVE).number_of_events();
        
            bool root_of_resultant=(event.mult_of_prim_res_root>0);
            bool root_of_content=(event.mult_of_content_root>0);
        
            size_type mult_of_resultant  = event.mult_of_prim_res_root;
            
            //AcX_DSTREAM("Event line for " << index << " " << root_of_resultant << " " << root_of_content << " " << mult_of_resultant << " " << left_arcs << " " << right_arcs << std::endl);
            Status_line_1 ev_line 
                = event_line_builder().create_event_line(index,
                                                         x,
                                                         left_arcs,
                                                         right_arcs,
                                                         root_of_resultant,
                                                         root_of_content,
                                                         mult_of_resultant);
        
            event.stack = ev_line;
            AcX_DSTREAM("done" << std::endl);
            return ev_line;
        } catch(CGAL::CGALi::Non_generic_position_exception exc) {
            switch(this->ptr()->degeneracy_strategy) {
            case(CGAL::CGALi::EXCEPTION_STRATEGY): {
                throw CGAL::CGALi::Non_generic_position_exception();
                break;
            }
            case(CGAL::CGALi::SHEAR_STRATEGY): {
                return create_non_generic_event_with_shear(index);
                break;
            }
            }
        }
        // !!! Never reached
        return Status_line_1();
    }

private:

    /**! 
     * \brief Method for non-generic situations, using shear and backshear 
     *
     * Note that this methods creates <b>all</b> event lines of the object
     */
    Status_line_1 create_non_generic_event_with_shear(size_type index) const {

        AcX_DSTREAM("Use sheared technique..." << std::endl);
        CGALi::Shear_controller<Integer> shear_controller;
        Integer s(0);
        while(true) {
            try {
                s = shear_controller.get_shear_factor();

                AcX_DSTREAM("Trying shear factor " << s << std::endl);
                // TODO: Move shear somewhere else
                Self D(CGAL::CGALi::shear(f_primitive(),Coefficient(s)),
                       CGAL::CGALi::EXCEPTION_STRATEGY);
                Shear_transformation< Self > 
                    shear_transformation;
                shear_transformation.report_sheared_disc_roots
                    (boost::make_transform_iterator(
                             event_coordinates().begin(),
                             typename Rep::Val_functor()),
                     boost::make_transform_iterator(
                             event_coordinates().end(),
                             typename Rep::Val_functor()) 
                    );
              
                // Store the sheared curve for later use
                this->ptr()->sheared_curves.insert(std::make_pair(s,D));
                shear_transformation(D,-s,(Self&)*this,false);
                set_vertical_line_components();
                
                break;
            }
            catch(CGAL::CGALi::Non_generic_position_exception err) {

                shear_controller.report_failure(s);
                AcX_DSTREAM("Bad shear factor, retrying..." << std::endl);
            }
        }
        
        return status_line_at_event(index);
    }

public:

    Status_line_1 status_line_of_interval(size_type i) const
    {
        CGAL_precondition(i >= 0 && i <= number_of_status_lines_with_event());
        
        Boundary b = boundary_value_in_interval(i);

        return status_line_at_exact_non_event_x(X_coordinate_1(b));
    }
    

public:

    Status_line_1 status_line_for_x(X_coordinate_1 x,
                                    CGAL::Sign perturb = CGAL::ZERO) const
    {
        size_type i;
        bool is_evt;
        x_to_index(x, i, is_evt);
        if(is_evt) {
            if(perturb == CGAL::ZERO)
                return status_line_at_event(i);
            if(perturb == CGAL::POSITIVE)
                i++;
        } 
        return status_line_of_interval(i);
    }


private:

    /*!
     * \brief Creates an intermediate line at rational position \c b
     */
    Status_line_1
    create_status_line_at_non_event(Boundary b) const {
        return create_status_line_at_non_event(X_coordinate_1(b));
    }

private:

    /*! 
     * \brief Creates an intermediate line at position \c ar.
     *
     * It is required that none of the following situations occurs at position
     * <tt>ar</tt>: singularity, vertical tangent line, vertical asymptote.\n
     * Otherwise, the method might run into an infinite loop. 
     * 
     * Note that the returned object is fixed and has \c id -1.
     */
    Status_line_1
    create_status_line_at_non_event(X_coordinate_1 ar, int index = -1) const {

        if(index==-1) {
            bool event;
            x_to_index(ar,index,event);
            CGAL_assertion(!event);
        }
        CGAL_assertion(index>=0);

        // TODO .. delay creation of refinement object 
        // especially when ar is rational
        
        Bitstream_traits traits(ar);

        Bitstream_descartes 
            bitstream_descartes(CGAL::CGALi::Square_free_descartes_tag(),
                                f_primitive(),
                                traits);

        size_type root_number=bitstream_descartes.number_of_real_roots();

        Status_line_1 status_line(ar, index, *this, 
                                  root_number);
        status_line.set_isolator(bitstream_descartes);

        return status_line;
    }

private:

   /**!
    * \brief Returns an instance of an event line builder
    *
    * Note: So far, a new instance is created each time the function is called
    * TODO: Fix this 
    */
    Event_line_builder event_line_builder() const {
        
        /*
        if(! this->ptr()->event_line_builder) {
            this->ptr()->event_line_builder 
                = Event_line_builder(this, f_primitive());
        }
                
        return this->ptr()->event_line_builder.get();
        */
        return Event_line_builder(*this, f_primitive());
    }

public:

    //! \brief Number of arcs over the given interval
    size_type arcs_over_interval(size_type i) const {
        //std::cout << "Polynomial=" << this->ptr()->f.get() 
        //<< ".." << std::flush;
        CGAL_precondition(has_defining_equation());
        //std::cout << "okay!" << std::endl;
        size_type n = static_cast<size_type>(intermediate_values().size());
        (void)n;
        CGAL_precondition(i>=0 && i<=n);
        return status_line_of_interval(i).number_of_events();
    }

public:

    //! \brief Rational number in the <tt>i</tt>th interval between events
    Boundary boundary_value_in_interval(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>(intermediate_values().size()));
        if(! intermediate_values()[i]) {
          // Create it
            if(event_coordinates().size()==0) {
                CGAL_assertion(i==0);
                intermediate_values()[0]=Boundary(0);
            } else {
                if(i==0) {
                    intermediate_values()[i] 
                        = simple_rational_left_of(event_coordinates()[i].val);
                } else if(i == static_cast<size_type>(event_coordinates().size())) {
                    intermediate_values()[i] 
                        = simple_rational_right_of
                        (event_coordinates()[i-1].val);
                    
                } else {
                    intermediate_values()[i]
                        = simple_rational_between(event_coordinates()[i-1].val,
                                                  event_coordinates()[i].val);
                    //= event_coordinates()[i-1].rational_between
                    //(event_coordinates()[i]);
                }
            }
        }
        return intermediate_values()[i].get();
    }

public:

    Y_coordinate_1 y_at(const Boundary& r, size_type arcno) const {
        CGAL_assertion(arcno >= 0);

        typedef typename CGAL::Fraction_traits<Poly_coer_1> FT;
        
        CGAL_assertion
            (static_cast<bool>((boost::is_same
                                < typename FT::Numerator_type,
                                  Polynomial_1 >::value)));

        typename FT::Numerator_type p;
        typename FT::Denominator_type denom;
        
        Poly_coer_1 f_at_r_with_denom = f().evaluate(r);

        typename FT::Decompose()(f_at_r_with_denom,p,denom);

        Status_line_1 cvl = event_info_at_x(r);
        return Y_coordinate_1(p, cvl.lower_boundary(arcno),
                            cvl.upper_boundary(arcno));
    }
    

    typedef std::pair<Polynomial_2,Polynomial_2> Polynomial_pair;

    template<typename T> class Gcd {
    public:

        T operator() (std::pair<T,T> pair) {
            return CGAL::CGALi::gcd(pair.first,pair.second);
        }
    } ;     

    template<typename T> class Pair_cannonicalize {
    public:
        std::pair<T,T> operator() (std::pair<T,T> pair) {
            if(pair.first > pair.second) {
                return std::make_pair(pair.second,pair.first);
            }
            else {
                return pair;
            }
        }
    };

    typedef CGAL::Pair_lexicographical_less_than<Polynomial_2,Polynomial_2,
      std::less<Polynomial_2>,
      std::less<Polynomial_2> > Pair_compare;
    
    typedef CGAL::Cache<Polynomial_pair,Polynomial_2,Gcd<Polynomial_2>,
      Pair_cannonicalize<Polynomial_2>,
      Pair_compare> Gcd_cache;

    static Gcd_cache& get_gcd_cache() {
        static Gcd_cache cache;
        return cache;
    }

    // end of adds for static caching

public:

    //! Returns the content of the polynomial
    Polynomial_1 content() const {
        if(! this->ptr()->content) {
            compute_content_and_primitive_part();
        }
        return this->ptr()->content.get();
    }

public:

    //! Returns the primitive part of the polynomial
    Polynomial_2 f_primitive() const {
        if(! this->ptr()->f_primitive) {
            compute_content_and_primitive_part();
        }
        return this->ptr()->f_primitive.get();
    }

private:

    //! compute and set content and primitive part
    void compute_content_and_primitive_part() const {

        CGAL_assertion(has_defining_equation());

        AcX_DSTREAM("Computing the content..." << std::flush);
        this->ptr()->content = typename CGAL::Polynomial_traits_d< Polynomial_2 >::Univariate_content_up_to_constant_factor()( f() );
        if(content().degree()==0) {
            AcX_DSTREAM("no vertical lines as components" << std::endl);
            this->ptr()->f_primitive=f();
        }
        else {
            AcX_DSTREAM("non-trivial content found" << std::endl);
            // Content must be square free, because the curve is square free
            CGAL_assertion( typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free()(content()));
            this->ptr()->f_primitive=f() / content();
	    
        }

    }

private:

    //! Returns the Sturm-Habicht sequence of the primitive part of f
    std::vector<Polynomial_2>& sturm_habicht_primitive_f() const 
    throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->sturm_habicht_primitive_f) {
            compute_sturm_habicht_primitive_f();
        }  
        return this->ptr()->sturm_habicht_primitive_f.get();
    }

public: 

    //! Returns the ith Sturm-Habicht polynomial of the primitive part of f
    Polynomial_2 sturm_habicht_primitive_f(size_type i) const 
      throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        CGAL_assertion(i>=0 && 
                    i < static_cast<size_type>(sturm_habicht_primitive_f().size()));
        return sturm_habicht_primitive_f()[i];
    }

public:

    //! Returns the ith principal sturm habicht coefficient
    Polynomial_1 principal_sturm_habicht_primitive_f(size_type i) const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        CGAL_assertion(i>=0 && 
                    i < static_cast<size_type>
                       (sturm_habicht_primitive_f().size()));

        CGAL_assertion(sturm_habicht_primitive_f()[i].degree()<=i);
        if(sturm_habicht_primitive_f()[i].degree() < i) {
            return Polynomial_1(0);
        } // else:
        return sturm_habicht_primitive_f()[i][i];
    }

public:

    //! Returns the ith principal sturm habicht coefficient
    Polynomial_1 coprincipal_sturm_habicht_primitive_f(size_type i) const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        CGAL_assertion(i>=1 && 
                    i < static_cast<size_type>(sturm_habicht_primitive_f().size()));
        CGAL_assertion(sturm_habicht_primitive_f()[i].degree()<=i);
        if(sturm_habicht_primitive_f()[i].degree() < i-1) {
            return Polynomial_1(0);
        } // else:
        return sturm_habicht_primitive_f()[i][i-1];
    }

public:

    //! Returns an iterator to the principal Sturm-Habicht coefficients,
    //! starting with the resultant
    Principal_sturm_habicht_iterator principal_sturm_habicht_begin() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(0),
             Stha_functor(this));
    }

    //! Returns an iterator to the end of principal Sturm-Habicht coefficients
    Principal_sturm_habicht_iterator principal_sturm_habicht_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(sturm_habicht_primitive_f().size()),
             Stha_functor(this));
    }

private:

    void compute_sturm_habicht_primitive_f() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        
        AcX_DSTREAM("Compute Sturm-Habicht.." << std::flush);
        std::vector<Polynomial_2> stha;
        
        // Fix a problem for constant primitive part.
        // In this case, the St.-Ha. sequence is never needed
        if(f_primitive().degree() == 0) {
            // Set the resultant
            stha.push_back(f_primitive());
        } else {
            
#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
#warning USES BEZOUT MATRIX FOR SUBRESULTANTS
            CGAL::bezout_polynomial_subresultants
                (f_primitive(),
                 CGAL::diff(f_primitive()),
                 std::back_inserter(stha));
            stha.push_back(f_primitive());
            size_type p = f_primitive().degree();
            CGAL_assertion(static_cast<size_type>(stha.size()) == p+1);
            for(size_type i=0;i<p; i++) {
                if((p-i)%4==0 || (p-i)%4==1) {
                    stha[i] = stha[i];
                } else {
                    stha[i] = -stha[i];
                }
            }
            
#else
            typename Polynomial_traits_2::Sturm_habicht_sequence()
                (f_primitive(),std::back_inserter(stha));
#endif
        }
        // Also set the resultant, if not yet set
        if(! this->ptr()->resultant_primitive_f_fy) {
            this->ptr()->resultant_primitive_f_fy = stha[0][0];
            if(this->ptr()->resultant_primitive_f_fy.get().is_zero()) {
                throw CGALi::Zero_resultant_exception<Polynomial_2>(f());
            }
        }
        
        this->ptr()->sturm_habicht_primitive_f = stha;
        CGAL_assertion(resultant_primitive_f_fy() == 
                       principal_sturm_habicht_primitive_f(0) ||
                       resultant_primitive_f_fy() == 
                       -principal_sturm_habicht_primitive_f(0) );
        AcX_DSTREAM("done" << std::endl);
    }

public:

    //! Returns the resultant of the primitive part of f with its y-derivative
    Polynomial_1 resultant_primitive_f_fy() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->resultant_primitive_f_fy) {
            compute_resultant_primitive_f_fy();
        }
        return this->ptr()->resultant_primitive_f_fy.get();
    }

    //! Returns the resultant of the primitive part of f with its x-derivative
    Polynomial_1 resultant_primitive_f_fx() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->resultant_primitive_f_fx) {
            compute_resultant_primitive_f_fx();
        }
        return this->ptr()->resultant_primitive_f_fx.get();
    }

private:

    void compute_resultant_primitive_f_fy() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        
        AcX_DSTREAM("Compute resultant.." << std::flush);

        CGAL_assertion(has_defining_equation());

#if AcX_SPEED_UP_FOR_REGULAR_CURVES
        bool speed_up=true;
#else
#ifdef AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL
        bool speed_up = f().degree() > 
            AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL;
#else
        bool speed_up=false;
#endif
#endif
        
        if(! speed_up) {
            
            // Compute resultant using the Sturm-Habicht sequence
	  if(f().degree() == 0) {
	    this->ptr()->resultant_primitive_f_fy = Polynomial_1(1);
	  } else {
            this->ptr()->resultant_primitive_f_fy 
	      = principal_sturm_habicht_primitive_f(0);
	  }
            
        } else {
            
            this->ptr()->resultant_primitive_f_fy
                = CGAL::CGALi::resultant
                    (f_primitive(),CGAL::diff(f_primitive()));

        }

        AcX_DSTREAM("done" << std::endl);

        if(resultant_primitive_f_fy().is_zero()) {
            throw CGALi::Zero_resultant_exception<Polynomial_2>(f());
        }
    }
    

    void compute_resultant_primitive_f_fx() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        
        AcX_DSTREAM("Compute x-resultant.." << std::flush);

        CGAL_assertion(has_defining_equation());

        // Transpose the polynomial
        Polynomial_2 f_yx = typename Polynomial_traits_2::Swap() (f(),0,1);

        if( f_yx.degree() == 0 ) {
            // Polynomial only consists of horizontal lines
            // primitive resultant is set to 1
            this->ptr()->resultant_primitive_f_fx = Polynomial_1(1);
        } else {
            
            Polynomial_2 f_yx_primitive;
            
            Polynomial_1 content_yx = typename CGAL::Polynomial_traits_d< Polynomial_2 >::Univariate_content_up_to_constant_factor()( f_yx );
            if(content_yx.degree()==0) {
                f_yx_primitive=f_yx;
            }
            else {
                CGAL_assertion(typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free()(content_yx));
                f_yx_primitive=f_yx / content_yx;
                
            }
            
            this->ptr()->resultant_primitive_f_fx
                = CGAL::CGALi::resultant
                (typename Polynomial_traits_2::Swap() (f_yx_primitive,0,1),
                 typename Polynomial_traits_2::Swap() 
                     (CGAL::diff(f_yx_primitive),0,1) );
        }
        
        AcX_DSTREAM("done" << std::endl);

        if(resultant_primitive_f_fx().is_zero()) {
            throw CGALi::Zero_resultant_exception<Polynomial_2>(f());
        }
    }




private:

    std::vector<Event_coordinate_1>& event_coordinates() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->event_coordinates) {
            compute_event_coordinates();
        }
        return this->ptr()->event_coordinates.get();
    }

private:

    std::vector<boost::optional<Boundary> >& intermediate_values() const 
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
        
        if(! this->ptr()->intermediate_values) {
            // This is created during event_coordiantes()
            event_coordinates();
            CGAL_assertion(this->ptr()->intermediate_values);
        }
        return this->ptr()->intermediate_values.get();
    }


private:

    /*!
     * \brief Constructs the event values of the curve.
     *
     * Isolates the real roots of the resultant of the curve with its 
     * derivative. The result is stored in the event_values vector.
     * If <tt>res</tt> is specified, this polynomials' roots are isolated
     * So the function relies on the fact that you only pass the resultant
     * (up to a constant).
     * Also, if the curve has vertical lines as components, the positions
     * of them are inserted into the sequence
     */
    void compute_event_coordinates() const
        throw(CGALi::Zero_resultant_exception<Polynomial_2>) {
         
        AcX_DSTREAM("compute events..." << std::flush);
         
        Solve_1 solve_1;
         
        std::vector<X_coordinate_1> content_roots;
        std::vector<size_type> content_mults;
        solve_1(content(),
                std::back_inserter(content_roots),
                std::back_inserter(content_mults));
        
        // Set the vertical_line_components flag as side effect
        this->ptr()->has_vertical_component = (content_roots.size() > 0);

        std::vector<X_coordinate_1> res_roots;
        std::vector<size_type> res_mults;
        solve_1(resultant_primitive_f_fy(),
                std::back_inserter(res_roots),
                std::back_inserter(res_mults));
        
        std::vector<X_coordinate_1> lcoeff_roots;
        std::vector<size_type> lcoeff_mults;
        solve_1(f_primitive().lcoeff(),
                std::back_inserter(lcoeff_roots),
                std::back_inserter(lcoeff_mults));
        

        //Now, merge the vertical line positions with the resultant roots
        typename CGAL::Real_embeddable_traits<X_coordinate_1>::Compare compare;

        std::vector<X_coordinate_1> event_values;
        std::vector<CGAL::CGALi::Three_valued> event_values_info;

        CGAL::CGALi::set_union_with_source
            (res_roots.begin(),
             res_roots.end(),
             content_roots.begin(),
             content_roots.end(),
             std::back_inserter(event_values),
             std::back_inserter(event_values_info),
             compare);

        // Now, build the Event_coordinate_1 entries 
        // for each element of event_values
        size_type curr_res_index = 0, curr_content_index = 0, curr_lcoeff_index = 0;
        std::vector<Event_coordinate_1> event_coordinate_vector;

        for(size_type i = 0; i < static_cast<size_type>(event_values.size()); i++ ) {
            
            Event_coordinate_1 curr_event;
            curr_event.val = event_values[i];
            switch(event_values_info[i]) {
            
            case(CGAL::CGALi::ROOT_OF_FIRST_SET): {
                curr_event.index_of_prim_res_root = curr_res_index;
                CGAL_expensive_assertion(res_roots[curr_res_index] == 
                                         event_values[i]);
                curr_event.mult_of_prim_res_root 
                    = res_mults[curr_res_index];
                curr_res_index++;
                if(curr_lcoeff_index < static_cast<size_type>(lcoeff_roots.size()) &&
                   event_values[i]==lcoeff_roots[curr_lcoeff_index]) {
                    // We have a root of the leading coefficient
                    // of the primitve polynomial
                    curr_event.index_of_prim_lcoeff_root = curr_lcoeff_index;
                    curr_event.mult_of_prim_lcoeff_root 
                        = lcoeff_mults[curr_lcoeff_index];
                    curr_lcoeff_index++;
                } else {
                    curr_event.index_of_prim_lcoeff_root = -1;
                    curr_event.mult_of_prim_lcoeff_root = 0;
                }
                
                curr_event.index_of_content_root = -1;
                curr_event.mult_of_content_root = 0;
                break;
            }
            case(CGAL::CGALi::ROOT_OF_SECOND_SET): {
                curr_event.index_of_content_root = curr_content_index;
                CGAL_expensive_assertion(content_roots[curr_content_index] == 
                                         event_values[i]);
                curr_event.mult_of_content_root 
                    = content_mults[curr_content_index];
                curr_content_index++;
                curr_event.index_of_prim_res_root = -1;
                curr_event.mult_of_prim_res_root = 0;
                CGAL_expensive_assertion(event_values[i]!=
                                         lcoeff_roots[curr_lcoeff_index]);
                curr_event.index_of_prim_lcoeff_root = -1;
                curr_event.mult_of_prim_lcoeff_root = 0;
                break;
            }
            case(CGAL::CGALi::ROOT_OF_BOTH_SETS): {
                curr_event.index_of_prim_res_root = curr_res_index;
                CGAL_expensive_assertion(res_roots[curr_res_index] == 
                                         event_values[i]);
                curr_event.mult_of_prim_res_root 
                    = res_mults[curr_res_index];
                curr_res_index++;
                if(curr_lcoeff_index < static_cast<size_type>(lcoeff_roots.size()) &&
                   event_values[i]==lcoeff_roots[curr_lcoeff_index]) {
                    // We have a root of the leading coefficient
                    // of the primitve polynomial
                    curr_event.index_of_prim_lcoeff_root = curr_lcoeff_index;
                    curr_event.mult_of_prim_lcoeff_root 
                        = lcoeff_mults[curr_lcoeff_index];
                    curr_lcoeff_index++;
                } else {
                    curr_event.index_of_prim_lcoeff_root = -1;
                    curr_event.mult_of_prim_lcoeff_root = 0;
                }
                curr_event.index_of_content_root = curr_content_index;
                CGAL_expensive_assertion(content_roots[curr_content_index] == 
                                         event_values[i]);
                curr_event.mult_of_content_root 
                    = content_mults[curr_content_index];
                curr_content_index++;
                break;
            }
            } // of switch
            /*
              AcX_DSTREAM("Constructed event_coordinate: " 
              << CGAL::to_double(curr_event.val) << " " 
              << "\nmult_of_prim_res_root : "
              << curr_event.mult_of_prim_res_root
              << "\nindex_of_prim_res_root : "
              << curr_event.index_of_prim_res_root
              << "\nmult_of_content_root : "
              << curr_event.mult_of_content_root
              << "\nindex_of_content_root : "
              << curr_event.index_of_content_root
              << "\nmult_of_lcoeff_root : "
              << curr_event.mult_of_prim_lcoeff_root
              << "\nindex_of_lcoeff_root : "
              << curr_event.index_of_prim_lcoeff_root
              << std::endl);
            */
            event_coordinate_vector.push_back(curr_event);
        }
        

        CGAL_assertion(curr_lcoeff_index == 
                       static_cast<size_type>(lcoeff_roots.size()));
        CGAL_assertion(curr_res_index == 
                       static_cast<size_type>(res_roots.size()));
        CGAL_assertion(curr_content_index == 
                       static_cast<size_type>(content_roots.size()));

        this->ptr()->intermediate_values 
            = std::vector<boost::optional<Boundary> >
            (event_coordinate_vector.size()+1);
        this->ptr()->event_coordinates = event_coordinate_vector;
      
        AcX_DSTREAM("done" << std::endl);

    }

public:    

    /*! 
     * \brief Applies a shear of a curve.
     *
     * Uses internal caching to avoid repeated shears
     */
    Self& shear_primitive_part(Integer s) const
        throw(CGAL::CGALi::Non_generic_position_exception)
    {
        CGAL_assertion(s!=0);
        if(this->ptr()->bad_shears.find(s) !=
           this->ptr()->bad_shears.end()) {
            throw CGAL::CGALi::Non_generic_position_exception();
        }
        typedef typename std::map<Integer,Self>::iterator 
            Map_iterator;
        Map_iterator it = this->ptr()->sheared_curves.find(s);
        if(it != this->ptr()->sheared_curves.end()) {
            return it->second;
        }
        try {
            Shear_transformation<Self> shear_transformation;
            Self D=shear_transformation((Self&)*this, s);
            std::pair<Map_iterator,bool> insertion =
                this->ptr()->sheared_curves.insert(std::make_pair(s,D));
            CGAL_assertion(insertion.second);
            return insertion.first->second;
        }
        catch(CGAL::CGALi::Non_generic_position_exception err) {
            this->ptr()->bad_shears.insert(s);
            throw CGAL::CGALi::Non_generic_position_exception();
        }
    }
    
public:

    //! Begin of the sheared curves
    typename std::map<Coefficient,Self>::const_iterator shear_begin() {
        return this->ptr()->sheared_curves.begin();
    }

    //! End of the sheared curves
    typename std::map<Coefficient,Self>::const_iterator shear_end() {
        return this->ptr()->sheared_curves.end();
    }

private:	
  
    void set_vertical_line_components() const {
        
        for(size_type i = 0; 
            i < static_cast<size_type>(event_coordinates().size()); 
            i++ ) {
            
            if(event_coordinates()[i].mult_of_content_root > 0) {
                status_line_at_event(i)._set_v_line();
            }
        }
         
    }
    

public:

    /*!
     * \brief Increases the precision of all Event_lines and Intermediate_lines
     *
     * For each Event_line and Intermediate_line, every isolating interval of
     * an arc is refined until its size is smaller then \c precision.
     */

    void refine_all(Boundary precision) {

        for(size_type i=0;i<static_cast<size_type>(event_coordinates().size());i++) {
            //	AcX_DSTREAM(i << ": " << std::flush);
            Status_line_1& el = status_line_at_event(i);
            for(size_type j=0;j<el.number_of_events();j++) {
                //AcX_DSTREAM(j << " " << std::flush);
                el.refine_to(j,precision);
            }
        }
        for(size_type i=0;i<static_cast<size_type>(intermediate_values().size());i++) {
            Status_line_1 il = status_line_of_interval(i);
            for(size_type j=0;j<il.number_of_events();j++) {
                il.refine_to(j,precision);
            }
        }
    }

public:

    //! \brief Iterator for the \c Status_line_1s
    Event_line_iterator event_begin() const {
        return boost::make_transform_iterator(boost::counting_iterator<size_type>(0),
                                              Event_functor(this));
    }

    //! \brief Iterator for the \c Status_line_1s
    Event_line_iterator event_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(num_events()),
             Event_functor(this));
    }

public:
   
    //! \brief Iterator for the \c Intermediate_lines
    Intermediate_line_iterator intermediate_begin() const {
        return boost::make_transform_iterator(boost::counting_iterator<size_type>(0),
                                              Intermediate_functor(this));
    }

    //! \brief Iterator for the \c Intermediate_lines
    Intermediate_line_iterator intermediate_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(intermediate_values().size()),
             Intermediate_functor(this));
    }

public:

    /*!
     * \brief 
     * Returns asymptotic information about arcs to -infty in x-direction
     *
     * Returns beta, if the \c i th arc to -infty is asymptotic for y=beta,
     * or +/-infty, if it converges to +/-infinity also in y-direction
     */
    CGAL::Object horizontal_asymptote_for_arc_to_minus_infinity(size_type i) const {
        if(! this->ptr()->horizontal_asymptotes_left) {
            compute_horizontal_asymptotes();
        }
        std::vector<CGAL::Object>& asym_info 
            = this->ptr()->horizontal_asymptotes_left.get();
        CGAL_precondition(i>=0 && 
                          i<static_cast<size_type>(asym_info.size()));
        return asym_info[i];
    }

public:

    CGAL::Object asymptotic_value_of_arc(CGAL::Arr_parameter_space loc,
                                         size_type arcno) const {
        
        CGAL_precondition(loc == CGAL::ARR_LEFT_BOUNDARY ||
                          loc == CGAL::ARR_RIGHT_BOUNDARY);
        
        if(loc == CGAL::ARR_LEFT_BOUNDARY) {
            
            if(! this->ptr()->horizontal_asymptotes_left) {
                compute_horizontal_asymptotes();
            }
            std::vector<CGAL::Object>& asym_info 
                = this->ptr()->horizontal_asymptotes_left.get();
            CGAL_precondition(arcno>=0 && 
                              arcno<static_cast<size_type>(asym_info.size()));
            return asym_info[arcno];
        } // else loc == CGAL::ARR_RIGHT_BOUNDARY

        if(! this->ptr()->horizontal_asymptotes_right) {
            compute_horizontal_asymptotes();
        }
        std::vector<CGAL::Object>& asym_info 
            = this->ptr()->horizontal_asymptotes_right.get();
        CGAL_precondition(arcno>=0 && 
                          arcno<static_cast<size_type>(asym_info.size()));
        return asym_info[arcno];
        
    }


private:

    void compute_horizontal_asymptotes() const {
      
        // TODO: Filter out curves with no arc to +/- infty

        typename CGAL::Fraction_traits<Boundary>::Decompose decompose;

        Integer num,denom;
        
        Solve_1 solve_1;
        std::vector<size_type> dummy_multiplicities;

        Polynomial_1 leading_coefficient_in_x 
            = typename Polynomial_traits_2::Swap() (this->f(),0,1).lcoeff();
        std::vector<Y_coordinate_1> roots_of_lcoeff;
        
        solve_1(leading_coefficient_in_x,
                std::back_inserter(roots_of_lcoeff),
                std::back_inserter(dummy_multiplicities));
        std::vector<Boundary> stripe_bounds;
        find_intermediate_values(roots_of_lcoeff.begin(),
                                 roots_of_lcoeff.end(),
                                 std::back_inserter(stripe_bounds));
        Boundary leftmost_bound = boundary_value_in_interval(0),
            rightmost_bound = boundary_value_in_interval(this->num_events());
        for(size_type i=0;i<static_cast<size_type>(stripe_bounds.size());i++) {
            Boundary& beta = stripe_bounds[i];
            decompose(beta,num,denom);
            Polynomial_1 poly_at_beta 
                = this->f().evaluate_homogeneous(num,denom);
            std::vector<X_coordinate_1> x_coordinates_at_beta;
            solve_1(poly_at_beta,std::back_inserter(x_coordinates_at_beta),
                    std::back_inserter(dummy_multiplicities));
            size_type number_of_roots 
                = static_cast<size_type>(x_coordinates_at_beta.size());
            if(number_of_roots>0) {
                if(leftmost_bound > x_coordinates_at_beta[0].low()) {
                    leftmost_bound = x_coordinates_at_beta[0].low();
                }
                if(rightmost_bound 
                   < x_coordinates_at_beta[number_of_roots-1].high()) {
                    rightmost_bound 
                        = x_coordinates_at_beta[number_of_roots-1].high();
                }
            }     
        }
        
        // Just to be sure...
        leftmost_bound = leftmost_bound - 1;
        rightmost_bound = rightmost_bound + 1;

        decompose(leftmost_bound,num,denom);
        Polynomial_1 curve_at_left_end 
            = typename Polynomial_traits_2::Swap() (this->f(),0,1)
            .evaluate_homogeneous(num,denom);
        std::vector<Y_coordinate_1> roots_at_left_end;
        solve_1(curve_at_left_end,std::back_inserter(roots_at_left_end));
        size_type number_of_roots_at_left_end 
            = static_cast<size_type>(roots_at_left_end.size());
        std::vector<CGAL::Object> asym_left_info;
        size_type current_stripe=0,i=0;
        while(i<number_of_roots_at_left_end) {
            if(current_stripe==static_cast<size_type>(stripe_bounds.size())) {
                asym_left_info.push_back( CGAL::make_object
                                              (CGAL::ARR_TOP_BOUNDARY) );
                i++;
                continue;
            }
            if(roots_at_left_end[i].low() > stripe_bounds[current_stripe]) {
                current_stripe++;
                continue;
            }        
            if(roots_at_left_end[i].high() < stripe_bounds[current_stripe]) {
                if(current_stripe==0) {
                    asym_left_info.push_back(CGAL::make_object
                                                 (CGAL::ARR_BOTTOM_BOUNDARY));
                    i++;
                    continue;
                } else {
                    asym_left_info.push_back(CGAL::make_object
                                 (roots_of_lcoeff[current_stripe-1]));
                    i++;
                    continue;
                }
            }
            roots_at_left_end[i].refine();
        }
        this->ptr()->horizontal_asymptotes_left = asym_left_info;
         
        decompose(rightmost_bound,num,denom);
        Polynomial_1 curve_at_right_end 
            = typename Polynomial_traits_2::Swap() (this->f(),0,1)
            .evaluate_homogeneous(num,denom);
        std::vector<Y_coordinate_1> roots_at_right_end;
        solve_1(curve_at_right_end,std::back_inserter(roots_at_right_end));
        size_type number_of_roots_at_right_end 
            = static_cast<size_type>(roots_at_right_end.size());
        std::vector<CGAL::Object> asym_right_info;
        current_stripe=0;
        i=0;
        while(i<number_of_roots_at_right_end) {
            if(current_stripe==static_cast<size_type>(stripe_bounds.size())) {
                asym_right_info.push_back(CGAL::make_object
                                              (CGAL::ARR_TOP_BOUNDARY) );
                i++;
                continue;
            }
            if(roots_at_right_end[i].low() > stripe_bounds[current_stripe]) {
                current_stripe++;
                continue;
            }        
            if(roots_at_right_end[i].high() < stripe_bounds[current_stripe]) {
                if(current_stripe==0) {
                    asym_right_info.push_back(CGAL::make_object
                                                  (CGAL::ARR_BOTTOM_BOUNDARY));
                    i++;
                    continue;
                } else {
                    asym_right_info.push_back
                        (CGAL::make_object(roots_of_lcoeff[current_stripe-1]));
                    i++;
                    continue;
                }
            }
            roots_at_right_end[i].refine();
        }
        this->ptr()->horizontal_asymptotes_right = asym_right_info;
 
    }
  

    // friend function for id-based hashing
    friend std::size_t hash_value(const Self& x) {
        return static_cast<std::size_t>(x.id());
    }
    
    
}; // class Algebraic_curve_2_2


//! \brief Prints the objects.
template<typename AlgebraicKernel_1, typename Rep_>
std::ostream& operator<< (
        std::ostream& out, 
        const Curve_analysis_2< AlgebraicKernel_1, 
        Rep_ >& curve) {
    typedef Curve_analysis_2< AlgebraicKernel_1, 
        Rep_ > Curve;

    typedef typename Curve::size_type size_type;

    out << "--------------- Analysis results ---------------" << std::endl;
    out << "Number of constructed event lines: " << curve.num_events() 
	<< std::endl;
    out << "(Horizontal) asymptotes at -infty: " << std::flush;
    for(size_type i = 0;i < curve.arcs_over_interval(0);i++) {
        
        const CGAL::Object& curr_asym_info_obj 
            = curve.asymptotic_value_of_arc(CGAL::ARR_LEFT_BOUNDARY,i);
        typename Curve::Y_coordinate_1 curr_asym_info;
        bool is_finite = CGAL::assign(curr_asym_info,curr_asym_info_obj);
        if(! is_finite) {
            CGAL::Arr_parameter_space loc;
            CGAL_assertion_code(bool is_valid = )
                CGAL::assign(loc, curr_asym_info_obj);
            CGAL_assertion(is_valid);
            if(loc == CGAL::ARR_TOP_BOUNDARY) {
                out << "+infty " << std::flush;
            } else {
                CGAL_assertion(loc == CGAL::ARR_BOTTOM_BOUNDARY);
                out << "-infty " << std::flush;
            }
        } else { // is_finite
            out << CGAL::to_double(curr_asym_info) 
                << " " << std::flush;
        }
    }

    out << std::endl;

    out << "Intermediate line at " 
	<< CGAL::to_double(curve.boundary_value_in_interval(0))
	<< ": " << curve.arcs_over_interval(0) << " passing arcs" << std::endl 
	<< std::endl;
    for(size_type i = 0; i<curve.num_events();i++) {
        out << curve.status_line_at_event(i) << std::endl;
        out << "Intermediate line at " 
            << CGAL::to_double(curve.boundary_value_in_interval(i+1))
            << ": " << curve.arcs_over_interval(i+1) 
            << " passing arcs" << std::endl
            << std::endl;
    }
    out << "(Horizontal) asymptotes at +infty: " << std::flush;
    for(size_type i = 0;i < curve.arcs_over_interval(curve.num_events());i++) {

        const CGAL::Object& curr_asym_info_obj 
            = curve.asymptotic_value_of_arc(CGAL::ARR_RIGHT_BOUNDARY,i);
        typename Curve::Y_coordinate_1 curr_asym_info;
        bool is_finite = CGAL::assign(curr_asym_info,curr_asym_info_obj);
        if(! is_finite) {
            CGAL::Arr_parameter_space loc;
            CGAL_assertion_code(bool is_valid = )
                CGAL::assign(loc, curr_asym_info_obj);
            CGAL_assertion(is_valid);
            if(loc == CGAL::ARR_TOP_BOUNDARY) {
                out << "+infty " << std::flush;
            } else {
                CGAL_assertion(loc == CGAL::ARR_BOTTOM_BOUNDARY);
                out << "-infty " << std::flush;
            }
        } else { // is_finite
            out << CGAL::to_double(curr_asym_info) 
                << " " << std::flush;
        }
    }

    out << std::endl;
    
    out << "------------------------------------------------" << std::endl;
    
    return out;
}

//! \brief Reads the objects from stream.
template<typename AlgebraicKernel_1, typename Rep_>
std::istream& operator>> (
        std::istream& in, 
        Curve_analysis_2< AlgebraicKernel_1, Rep_ >& curve) {
    typename Curve_analysis_2< AlgebraicKernel_1, Rep_ >::Polynomial_2 f;
    in >> f;
    curve.set_f(f);
    return in;
}
  

CGAL_END_NAMESPACE



#endif // ALGEBRAIC_CURVE_2_H
