// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_BASE_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_BASE_H 1

/*!\file include/CGAL/Algebraic_kernel_d/Algebraic_surface_3_z_at_xy_isolator_traits_base.h
 * \brief Traits to work with restricted cad
 */

#include <CGAL/config.h>

#include <CGAL/Cache.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/iterator/counting_iterator.hpp>


#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_in_z_for_xy_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! type of model of CGAL::ArrangementTraits_2 concept
template < class CurvedKernelViaAnalysis_2 >
class Wrap_ckva_for_z_at_xy_isolator_traits_2 : 
        public CurvedKernelViaAnalysis_2 {
public:  

    //! this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! base type
    typedef Curved_kernel_via_analysis_2 Base;
    
    //! type of curve kernel
    typedef typename Base::Curve_kernel_2 Curve_kernel_2;
    
    //! the class itself
    typedef Wrap_ckva_for_z_at_xy_isolator_traits_2< 
        Curved_kernel_via_analysis_2 
    >
    Self;

    //! default constructor
    Wrap_ckva_for_z_at_xy_isolator_traits_2() :
        Base() {
    }
    
    //! constructor from given curve kernel instance
    Wrap_ckva_for_z_at_xy_isolator_traits_2(const Curve_kernel_2& kernel) :
        Base(kernel) {
    }
    
    //! do not allow merges
    typedef CGAL::Tag_false           Has_merge_category;
};

} // namespace CGALi

/*!\brief
 * Abstract base class for all model of ZAtXyIsolatorTraits
 */
template < class CurvedKernelViaAnalysis_2, class Surface_3_>
class Algebraic_surface_3_z_at_xy_isolator_traits_base {
    
public:
    //////////////////////////////////////////////////////////////////////////
    // Types
    
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef Surface_3_ Surface_3;
    
    //! the class itself
    typedef 
    Algebraic_surface_3_z_at_xy_isolator_traits_base< 
        Curved_kernel_via_analysis_2, Surface_3 
    > Self;

    
public:    
    //!\name Public types
    //!@{
    
    //! type of curve kernel
    typedef typename 
    Curved_kernel_via_analysis_2::Curve_kernel_2 Curve_kernel_2;
   
    //! type of arrangement traits
    typedef CGAL::CGALi::Wrap_ckva_for_z_at_xy_isolator_traits_2 <
        Curved_kernel_via_analysis_2 > Arrangement_traits_2;

    // ensure that there is no merge
    BOOST_STATIC_ASSERT(
            (boost::is_same<CGAL::Tag_false, 
             typename Arrangement_traits_2::Has_merge_category>::value)
    );
    
    //! type of projected point
    typedef typename Arrangement_traits_2::Point_2 Point_2;

    //! type of projected arc
    typedef typename Arrangement_traits_2::Arc_2 X_monotone_curve_2;

    //! type of projected curve
    typedef typename Surface_3::Polynomial_3 Polynomial_3;
    
    //! type of projected curve
    typedef typename Polynomial_3::NT Polynomial_2;
    
    //! type of Isolator traits
    typedef CGAL::Bitstream_in_z_for_xy_traits< Self > Isolator_traits;

    //! type of isolator
    typedef CGAL::CGALi::Bitstream_descartes< Isolator_traits > 
    Z_at_xy_isolator;

    //! type of Rational
    typedef typename Isolator_traits::Rational Rational;

    //! type of Interval
    typedef typename Isolator_traits::Interval Interval;
    
    // these functors must be refined
    //! type of multiplicities of curve events
    typedef CGAL::Null_functor Multiplicities_of_curve_events;
    
    //! type of construct isolator
    typedef CGAL::Null_functor Construct_isolator;
    
    //! type of equal_z functor
    typedef CGAL::Null_functor Equal_z;

    //! type of adjacency functor
    typedef CGAL::Null_functor Adjacency;
    
    //!@}

    // TODO Implicit requirements:
    // * Surface_3::surface_cache() and cached surface construction,
    //   i.e., for each trivariate polynomial there exists exactly one
    //   Surface_3 instance!
    // * Surface_pair_3::surface_pair_cache() 
    //   -> but actually only needed in specializations

public:
    /*!\brief
     * functor to compute square-free factorization that also deals with
     * non-trivial contents
     */
    class Square_free_factorization_3 {
    public:
        template < class SurfaceFactorIterator, class MultIterator >
        int operator()(Polynomial_3 f, 
                       SurfaceFactorIterator cit, 
                       MultIterator mit) {
            
            int n = 0;

            if(CGAL::is_zero(f)) {
                return n;
            }
            
            std::vector< Polynomial_3 > factors;
            
            typedef typename Polynomial_3::NT Polynomial_2;
            
            // Remark: The latter analysis (mults) may rely on the fact 
            // that vertical components are surfaces on its own. Therefore
            // we split them off
            Polynomial_2 content = f.content();
            
            std::vector< Polynomial_2 > cfactors;
            std::vector< int > cmultiplicities;
            // TODO what if content itself has a non-trivial content
            CGAL::CGALi::filtered_square_free_factorize(
                    content, 
                    std::back_inserter(cfactors),
                    std::back_inserter(cmultiplicities)
            );
            CGAL_assertion(cfactors.size() == cmultiplicities.size());
            
            typename std::vector< Polynomial_2 >::iterator cfit = 
                cfactors.begin();
            
            typename std::vector< int >::iterator cmit = 
                cmultiplicities.begin();
            
            for (;cfit != cfactors.end(); cfit++, cmit++) {
                if (cfit->degree() > 0) { 
                    factors.push_back(Polynomial_3(*cfit));
                    *mit++ = *cmit;
                    n++;
                }
            }
            
            f /= content;
            
            n += CGAL::CGALi::filtered_square_free_factorize(
                    f, std::back_inserter(factors), mit
            );
            
            // create surface
            int count = 0;
            for (typename std::vector< Polynomial_3 >::const_iterator fit =
                     factors.begin(); fit != factors.end(); fit++) {
                if (CGAL::total_degree(*fit) > 0) {
                    Surface_3 surface(Surface_3::surface_cache()(*fit));
                    *cit++ = surface;
                    count++;
                }
            }
            
            return count;
        }
    };

    /*!\brief
      * returns instance of Square_free_factorization_3
     */
    Square_free_factorization_3 square_free_factorization_3_object() const { 
        return Square_free_factorization_3(); 
    }
    
public:

    //////////////////////////////////////////////////////////////////////////
    // Curves

protected:
    /*!\brief
     * functor to computes square-free facortization that also deals with
     * non-trivial contents
     */
    class Square_free_factorization_2 {
    protected:
        //! key type
        typedef Polynomial_2 Key;
        
        typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
        
        //! data type
        typedef std::list< std::pair< Curve_analysis_2, int > > Data;
        
        //! functor to canonicalize keys
        class Canonicalizer {
        public:
            /*!\brief
             * Returns canonicalized version of \c p
             */
            Polynomial_2 operator()(const Polynomial_2& p) {
                return CGAL::CGALi::canonicalize_polynomial(p);
            }
        };
        
        //! create data from canonicalized key
        class Creator {
        public:
            Data operator()(Polynomial_2 f) const {
                
                std::vector< Polynomial_2 > factors;
                std::vector< int > multiplicities;
                
                typedef typename Polynomial_2::NT Polynomial_1;
                
                // first for content
                // Remark: The latter analysis (mults) may rely on the fact 
                // that vertical components are curves on its own. Therefore
                // we split them off
                Polynomial_1 content = f.content();
                
                std::vector< Polynomial_1 > cfactors;
                std::vector< int > cmultiplicities;
                CGAL::CGALi::filtered_square_free_factorize(
                        content, 
                        std::back_inserter(cfactors),
                        std::back_inserter(cmultiplicities)
                );
                CGAL_assertion(cfactors.size() == cmultiplicities.size());
                
                typename std::vector< Polynomial_1 >::iterator cfit = 
                    cfactors.begin();
                
                typename std::vector< int >::iterator cmit = 
                    cmultiplicities.begin();
                
                for (;cfit != cfactors.end(); cfit++, cmit++) {
                    if (cfit->degree() > 0) { 
                        // TASK set to true and debug segfault
                        factors.push_back(Polynomial_2(*cfit));
                        multiplicities.push_back(*cmit);
                    }
                }
                
                f /= content;
                
                // then for content-free polynomial
                CGAL::CGALi::filtered_square_free_factorize(
                        f, 
                        std::back_inserter(factors), 
                        std::back_inserter(multiplicities)
                );
                

                // output
                typename std::vector< Polynomial_2 >::iterator fit = 
                    factors.begin();
                
                typename std::vector< int >::iterator mit = 
                    multiplicities.begin();
                
                Data outlist;

                // create curve
                for (; fit != factors.end(); fit++, mit++) {
                    CGAL_assertion(mit != multiplicities.end());
                    if (CGAL::total_degree(*fit) > 0) {
                        typename 
                            Curve_kernel_2::Construct_curve_2
                            construct_curve = 
                            Arrangement_traits_2::instance().kernel().
                            construct_curve_2_object();
                        Curve_analysis_2 curve =
                            construct_curve(*fit);
                        outlist.push_back(std::make_pair(curve, *mit));
                    }
                }
                
                return outlist;
            }
        };
        
        //! Cache type
        typedef 
        CGAL::Cache< Polynomial_2, Data, Creator, Canonicalizer > Curves_cache;
        
        //! returns static instance of curves cache
        static Curves_cache& curves_cache() {
            static Curves_cache cache;
            return cache;
        }
        
    public:
        
        /*!\brief
         * returns pairs of square-free factors of \c f with multiplicity
         */
        template < class OutputIterator >
        int operator()(Polynomial_2 f, OutputIterator oi) const {

            // Catch the zero polynomial
            if(CGAL::is_zero(f)) {
                return 0;
            }
            
            const Data& clist = curves_cache()(f);
            
            std::copy(clist.begin(), clist.end(), oi);

            return clist.size();
        }
    };
    
    /*!\brief
      * returns instance of Square_free_factorization_2
     */
    Square_free_factorization_2 square_free_factorization_2_object() const { 
        return Square_free_factorization_2(); 
    }

public:
    /*!\brief
     * Compute projected silhouette curves of given surface
     */
    class Construct_projected_surface_curves_2 {
    public:

        /*!\brief
         * returns through \c oi projected
         * boundary curve(s) of \c s with multiplicities 
         *
         * value_type of OutputIterator is std::pair< P_curve_2, int>
         */
        template < class OutputIterator >
        int operator()(const Surface_3& surface, OutputIterator oi) const {
            
            Polynomial_3 surface_f = surface.f();
            CGAL_precondition(typename CGAL::Polynomial_traits_d< Polynomial_3 >::Is_square_free()(surface_f));
            
            // compute correct polynomial
            CGAL_assertion(surface_f.degree() > 0 ||
                           CGAL::total_degree(surface_f.content()) > 0);
            
            // TODO: remove factor a_n from resultant if a_n != c.
            Polynomial_2 poly = (surface_f.degree() == 0 ? 
                                 surface_f.content() : 
                                 surface.resultant_f_fz(true, true)
            );

            typename Self::Square_free_factorization_2 sqf;
            // TODO square_free_factorization_2_object()
            
            return sqf(poly, oi);
        }


        /*!\brief
         * returns through \c oi projected
         * boundary curve(s) of \c s with multiplicities 
         *
         * value_type of OutputIterator is std::pair< P_curve_2, int>
         */
        template < class OutputIterator >
        int operator()(const Surface_3& surface, int i,
                       OutputIterator oi) const {
            
            Polynomial_3 surface_f = surface.f();
            CGAL_precondition(typename CGAL::Polynomial_traits_d< Polynomial_3 >::Is_square_free()(surface_f));
            
            CGAL_precondition(0 <= i);
            CGAL_precondition(i <= surface_f.degree());
            
            Polynomial_2 poly = surface_f[i];
            
            typename Self::Square_free_factorization_2 sqf;
            // TODO square_free_factorization_2_object()
            
            return sqf(poly, oi);
        }
        
        /*!\brief
         * indicates whether the meant "curve" "covers"  the whole plane
         */
        bool operator()(const Surface_3& surface, int i) {

            Polynomial_3 surface_f = surface.f();
            CGAL_precondition(typename CGAL::Polynomial_traits_d< Polynomial_3 >::Is_square_free()(surface_f));
            
            CGAL_precondition(0 <= i);
            CGAL_precondition(i <= surface_f.degree());
            
            Polynomial_2 poly = surface_f[i];
            
            return poly == Polynomial_2(0);
        }


        /*!\brief
         * returns through \c oi projected
         * boundary curve(s) of \c s with multiplicities 
         */
        template < class OutputIterator >
        int operator()(const Surface_3& surface, int k, int n,
                       OutputIterator oi) const {
            
            Polynomial_3 surface_f = surface.f();
            CGAL_precondition(typename CGAL::Polynomial_traits_d< Polynomial_3 >::Is_square_free()(surface_f));

            CGAL_precondition(0 <= k);
            CGAL_precondition(k <= surface_f.degree());

            CGAL_precondition(0 <= n);
            CGAL_precondition(n <= surface_f.degree());
            
            Polynomial_2 poly = surface.principal_sturm_habicht_coefficient(
                    k, n, true
            );
            
            typename Self::Square_free_factorization_2 sqf;
            // TODO square_free_factorization_2_object()
            
            return sqf(poly, oi);
        } 
        
        /*!\brief
         * indicates whether the meant "curve" "covers"  the whole plane
         */
        bool operator()(const Surface_3& surface, int k, int n) {
            
            Polynomial_3 surface_f = surface.f();
            CGAL_precondition(typename CGAL::Polynomial_traits_d< Polynomial_3 >::Is_square_free()(surface_f));
            
            CGAL_precondition(0 <= k);
            CGAL_precondition(k <= surface_f.degree());

            CGAL_precondition(0 <= n);
            CGAL_precondition(n <= surface_f.degree());
            
            Polynomial_2 poly = surface.principal_sturm_habicht_coefficient(
                    k, n, true
            );
            return poly == Polynomial_2(0);
        }
        

    };
    
    /*!\brief
     * returns instance of Construct_projected_silhouettes_2
     */
    Construct_projected_surface_curves_2 
    construct_projected_surface_curves_2_object() 
        const { 
        return Construct_projected_surface_curves_2(); 
    }
    
public:
   
    

    /*!\brief 
     * Computes projected cut curves for two two given surfaces
     */
    class Construct_projected_cuts_2 {
    public:
        /*!\brief
         * returns through \c oi projected
         * intersection curve(s) of \c s1 and \c s2
         */
        template < class OutputIterator >
        int operator()(const Surface_3& surface1, const Surface_3& surface2,
                       OutputIterator oi) const {
            
            Polynomial_3 surface_f1 = surface1.f();
            Polynomial_3 surface_f2 = surface2.f();
            // TODO replace with coprimality check
            CGAL_precondition(CGAL::CGALi::gcd(surface_f1, surface_f2).degree() < 1);
            
            // compute correct polynomial // TODO is CGAL::CGALi::resultant best?
            Polynomial_2 f = CGAL::CGALi::resultant(surface_f1, surface_f2);

            // compute cached squarefree factorization of intersection curves
            typename Self::Square_free_factorization_2 sqf;
            // TODO square_free_factorization_2_object()

            return sqf(f, oi);
        }
    };

    /*!\brief
     * returns instance of Construct_projected_cuts_2
     */
    Construct_projected_cuts_2
    construct_projected_cuts_2_object() const { 
        return Construct_projected_cuts_2(); 
    }


    /*!\brief
     * Constructs isolator traits for given point
     */
    class Construct_isolator_traits {
    public:
        
    private:
        // point less 
        struct Point_less {
            bool operator()(const Point_2& p1, const Point_2& p2) const {
                // TASK cache for points-xy or id?
                return p1.id() < p2.id();
                //return p1 < p2;
            }
        };
        
        //! type of traits map for a given point
        typedef std::map< Point_2, 
                          std::pair< 
                              boost::optional <Isolator_traits >,
                              boost::optional <Isolator_traits >
                          >, 
                          Point_less > 
        Isolator_traits_map;
        
        //! returns static member
        static 
        Isolator_traits_map& isolator_traits_map() {
            static Isolator_traits_map map;
            return map;
        }
        
    public:

#if 0
        //!\name Constructors
        //!@{
        
        /*!\brief
         * default constructor
         */
        Construct_isolator_traits() {
        }
        
        /*!\brief
         * Stores traits
         */
        Construct_isolator_traits(const Self& traits) :
            _m_traits(traits) {
        }
        
        //!@}
#endif
        
        /*!\brief
         * returns a cached isolator traits for given point \c pt
         */
        Isolator_traits operator()(const Point_2& point, 
                                   bool use_artificial_x_interval) {
            
            typename Isolator_traits_map::iterator itit = 
                isolator_traits_map().find(point);
            
            if (itit == isolator_traits_map().end()) {
                // construct one
                // if feature is a vertex, than we have to activate
                // the artificial x_interval mode
                itit = isolator_traits_map().insert(
                            itit, 
                            std::make_pair(
                                    point, std::make_pair(
                                            boost::none, boost::none
                                    )
                            )
                );
            }
            
            // if artificial interval is forced
            if (use_artificial_x_interval) {
                // create one
                if (!itit->second.second) {
                    itit->second.second = 
                        Isolator_traits(point, true);
                }
                // return cached
                return *(itit->second.second);
            }
            // if non-special available return it
            if (itit->second.first) {
                return *(itit->second.first);
            }
            // else if a special one is available return it
            if (itit->second.second) {
                return *(itit->second.second);
            }
            // else a normal one does the job as well
            CGAL_assertion(!itit->second.first);
            itit->second.first = 
                Isolator_traits(point, false);
            
            // finally return it
            return *(itit->second.first);
        }
        
#if 0
    private:
        // members
        //! instances of traits class
        mutable Self _m_traits;
#endif
    };
    
    /*!\brief
     * returns instance of Is_singular_2
     */
    Construct_isolator_traits construct_isolator_traits_2_object() const { 
        return Construct_isolator_traits(); 
    }

public:    
    /*!\brief
     * Functor to compute whether a given point lies on a given curve
     */
    class Point_on_curve_2 {
    public:
        
#if 0
        //!\name Constructors
        //!@{
        
        /*!\brief
         * default constructor
         */
        Point_on_curve_2() {
        }
        
        /*!\brief
         * Stores traits
         */
        Point_on_curve_2(const Self& traits) :
            _m_traits(traits) {
        }
        
        //!@}
#endif

    public:
        /*!\brief
         * returns \c true of pt lies on curve defined by f
         */
        bool operator()(const Point_2& pt, 
                        const Polynomial_2& non_square_free_f,
                        bool no_gcd = false) const {
            
            typename Self::Construct_isolator_traits 
                construct_isolator_traits; // TODO single instance
            Isolator_traits traits = 
                construct_isolator_traits(pt, false);
            
            return this->operator()(traits, non_square_free_f, no_gcd);
        }

        /*!\brief
         * returns \c true of pt lies on curve defined by f
         */
        bool operator()(const Isolator_traits& traits,
                        const Polynomial_2& non_square_free_f,
                        bool no_gcd = false) const {
            
            Polynomial_2 p2 = CGAL::CGALi::make_square_free(non_square_free_f);
            
            typename 
                Arrangement_traits_2::Curve_kernel_2::Construct_curve_2
                construct_curve = 
                Arrangement_traits_2::instance().kernel().
                construct_curve_2_object();
            typename 
                Arrangement_traits_2::Curve_kernel_2::Curve_analysis_2 curve =
                construct_curve(p2);
            
            return traits.point().is_on(curve);
        }
        
#if 0
    private:
        // members
        //! instances of traits class
        mutable Self _m_traits;
#endif
    };


    /*!\brief
     * returns instance of Point_on_curve_2
     */
    Point_on_curve_2 point_on_curve_2_object() const { 
        return Point_on_curve_2(); 
    }

public:

    class Intermediate_values_for_stack {
        
    public:
        
        typedef std::vector<Rational> result_type;
        
        template<typename OutputIterator>
        OutputIterator operator() (Z_at_xy_isolator isolator,
        OutputIterator out) const {
  
            int n = isolator.number_of_real_roots();

            if( n == 0 ) {

                *out++ = Rational(0);

            } else {
                
                *out++ = isolator.left_boundary(0) - Rational(1);

                for( int i = 1; i < n; i++ ) {
                    
                    *out++ = (isolator.left_boundary(i) + 
                              isolator.right_boundary(i-1)) / 2 ;
                    
                    // Note that even if both are equal, this value is always
                    // an intermediate one
                }
                
            *out++ = isolator.right_boundary(n-1)+Rational(1);

            }
            
            
            return out;
        }
    };

    class Assign_roots_to_buckets {

    public:
        
        template<typename BoundIterator, 
                 typename IndexIterator, 
                 typename OutputIterator>
        OutputIterator operator() (Z_at_xy_isolator isolator,
                                   IndexIterator index_begin,
                                   IndexIterator index_end,
                                   BoundIterator bound_begin,
                                   BoundIterator bound_end,
                                   OutputIterator out) {

            int bucket_roots = 0; //number of roots in current bucket

            BoundIterator it = bound_begin;
            IndexIterator i = index_begin;

            while( (it != bound_end) || i!=index_end) {

                if(i == index_end) {
                    *out++ = bucket_roots;
                    bucket_roots=0;
                    it++;
                    continue;
                }
                if(it == bound_end) {
                    bucket_roots++;
                    i++;
                    continue;
                }
                // A this point, we have a bucket with upper bound *it, and some
                // root with index i. Question: Is that root above or below *it?
                while( isolator.left_boundary(*i) <= *it &&
                       *it <= isolator.right_boundary(*i) ) {
                
                    isolator.refine_interval(*i);

                }
                
                if( isolator.left_boundary(*i) > *it ) {
                    *out++ = bucket_roots;
                    bucket_roots=0;
                    it++;
                    continue;
                } else {
                    bucket_roots++;
                    i++;
                    continue;
                } 
            }
            
            // At the end, there are those roots left which are in the highest bucket
            *out++ = bucket_roots;

            return out;

        }
       
        template<typename BoundIterator, typename OutputIterator>
        OutputIterator operator() (Z_at_xy_isolator isolator,
                                   BoundIterator bound_begin,
                                   BoundIterator bound_end,
                                   OutputIterator out) {
            
            int n = isolator.number_of_real_roots();

            return this->operator() (isolator,
                                     boost::counting_iterator<int>(0),
                                     boost::counting_iterator<int>(n),
                                     bound_begin,
                                     bound_end,
                                     out);

        }
    };

};

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_BASE_H
// EOF
