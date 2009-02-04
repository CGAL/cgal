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

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_H 1

/*!\file include/CGAL/Algebraic_kernel_d/Algebraic_surface_3_z_at_xy_isolator_traits.h
 * \brief Traits to work with restricted cad
 */

#include <CGAL/config.h>

#include <CGAL/Handle_with_policy.h>
#include <CGAL/function_objects.h>
#include <CGAL/Cache.h>

#include <CGAL/number_utils.h>

#include <CGAL/Polynomial/sturm_habicht_sequence.h>
#include <CGAL/Algebraic_curve_kernel_2/alg_real_utils.h>

#include <CGAL/Arrangement_2l/Z_stack_helpers.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_enums.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_accessor.h>
#include <CGAL/Arrangement_2l/Adjacencies_3.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_surface_3_z_at_xy_isolator_traits_base.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_surface_pair_3.h>

CGAL_BEGIN_NAMESPACE

/*!\brief
 * Model of ZAtXyIsolatorTraits for surfaces
 */
template < class CurvedKernelViaAnalysis_2, class Surface_3_ >
class Algebraic_surface_3_z_at_xy_isolator_traits : public 
CGAL::Algebraic_surface_3_z_at_xy_isolator_traits_base< 
    CurvedKernelViaAnalysis_2, Surface_3_
> {

public:
    //////////////////////////////////////////////////////////////////////////
    // Types
    
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef Surface_3_ Surface_3;
    
    //! the class itself
    typedef Algebraic_surface_3_z_at_xy_isolator_traits< 
        Curved_kernel_via_analysis_2, Surface_3 
    > Self;

    //! type of base
    typedef CGAL::Algebraic_surface_3_z_at_xy_isolator_traits_base< 
        Curved_kernel_via_analysis_2, Surface_3
    >
    Base;

    // TODO repeat all these types?
    //! type of arranagement traits
    typedef typename Base::Arrangement_traits_2 Arrangement_traits_2;

    //! type of projected point
    typedef typename Base::Point_2 Point_2;

    //! type of projected arc
    typedef typename Base::X_monotone_curve_2 X_monotone_curve_2;

    //! type of curve kernel
    typedef typename Base::Curve_kernel_2 Curve_kernel_2;

    //! type of X_coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;

    //! type of projected curve
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //! type of trivariate polynomial
    typedef typename Base::Polynomial_3 Polynomial_3;
    
    //! type of bivariate polynomial
    typedef typename Base::Polynomial_2 Polynomial_2;

    typedef typename Polynomial_2::NT Polynomial_1;

    //! type of Isolator traits
    typedef typename Base::Isolator_traits  Isolator_traits;
    
    //! type of isolator
    typedef typename Base::Z_at_xy_isolator Z_at_xy_isolator;

    //! type of rational number
    typedef typename Base:: Rational Rational;

    //! type of rational interval
    typedef typename Base:: Interval Interval;

    
    //////////////////////////////////////////////////////////////////////////
    // Isolators
public:
    /*!\brief
     * Constructs isolator for a surface at a given point
     */
    class Construct_isolator {
    private:
        
        struct Point_less {
            bool operator()(const Point_2& p1, const Point_2& p2) const {
                // TASK cache for points-xy or id?
                return p1.id() < p2.id();
                //return p1 < p2;
            }
        };

        //! type of isolator map for a given surface
        typedef std::map< Point_2, Z_at_xy_isolator, Point_less > Isolator_map;
        
        typedef CGAL::Polynomial<Polynomial_3> Polynomial_4;

        //! type of map of isolators for all surfaces
        typedef std::map< Surface_3, Isolator_map > Isolators_map;
        
        // TODO move to traits!
        //! returns static member
        static 
        Isolators_map& isolators_map() {
            static Isolators_map map;
            return map;
        }


        Z_at_xy_isolator _isolation_for_vertical_line_on_surface
        (const Surface_3 surface, Isolator_traits traits) const {

#if !NDEBUG
            std::cout << "Create pseudo-stack for vertical line.." 
                      << std::endl;
#endif 

            typedef CGAL::Polynomial_traits_d<Polynomial_3> 
                Polynomial_traits_3;

            Polynomial_3 f = surface.f();

            Polynomial_3 f_y = typename Polynomial_traits_3::Differentiate() 
                (f,1);
            
            Polynomial_3 disc_z( surface.resultant_f_fz() );

            
            std::vector< Polynomial_2 > sres_f_fy;

#if !NDEBUG
            std::cout << "Swapping f and f_y.." << std::flush;
#endif

            Polynomial_3 f_xzy = typename Polynomial_traits_3::Swap() (f,1,2);
            Polynomial_3 f_y_xzy = typename Polynomial_traits_3::Swap()
                (f_y,1,2);

#if !NDEBUG
            std::cout << "done" << std::flush;
#endif

#if !NDEBUG
            std::cout << "PRS of f and f_y.." << std::flush;
#endif
            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Principal_subresultants()(f_xzy, f_y_xzy, 
                                            std::back_inserter(sres_f_fy));
            sres_f_fy.push_back(f_y_xzy.lcoeff());

#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
#if !NDEBUG
            std::cout << "Res sq_free.." << std::flush;
#endif
            Polynomial_2 f1 = CGAL::make_square_free(sres_f_fy[0]);
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif            
#if !NDEBUG
            std::cout << "Res of f and disc(f)_z.." << std::flush;
#endif
            Polynomial_2 f2_non_sq 
                = typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Resultant()(f_xzy,disc_z);
#if !NDEBUG
            std::cout << "sq_free.." << std::flush;
#endif
            Polynomial_2 f2 
                = CGAL::make_square_free(f2_non_sq);
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
            
#if !NDEBUG
            std::cout << "Look up non-vanishing sres.." << std::flush;
#endif

            int i = 0;

            //typename Base::Point_on_curve_2 point_on_curve; 

            X_coordinate_1 x_val = traits.point().x();

            //Polynomial_3 curr = f;

/*            while( point_on_curve(traits,curr.content_utcf()) ) {
                i++;
                curr=typename Polynomial_traits_3::Derivative() (curr,1);
            }
*/

            while( (sres_f_fy[i].is_zero()) ||
                   x_val.is_root_of( typename CGAL::Polynomial_traits_d< Polynomial_2 >::Univariate_content_up_to_constant_factor()( sres_f_fy[i] ) ) ) {
                i++;
            }

#if !NDEBUG
            std::cout << "sq-free.." << std::flush;
#endif
            //Polynomial_3 f3 = CGAL::make_square_free(curr);
            Polynomial_2 f3 = CGAL::make_square_free(sres_f_fy[i]);

#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
            
            // Make coprime
#if !NDEBUG
            std::cout << "f1: " << f1 << std::endl;
            std::cout << "f2: " << f2 << std::endl;
            std::cout << "f3: " << f3 << std::endl;
            std::cout << "Coprime 1.." << std::flush;
#endif      
            Polynomial_2 gcd12 = CGAL::CGALi::gcd(f1,f2);
            f2 = CGAL::integral_division(f2,gcd12);
#if !NDEBUG
            std::cout << "2.." << std::flush;
#endif
            Polynomial_2 gcd13 = CGAL::CGALi::gcd(f1,f3);
            f3 = CGAL::integral_division(f3,gcd13);
#if !NDEBUG
            std::cout << "3.." << std::flush;
#endif            
            Polynomial_2 gcd23 = CGAL::CGALi::gcd(f2,f3);
            f3 = CGAL::integral_division(f3,gcd23);
#if !NDEBUG
            std::cout << "done" << std::endl;
            std::cout << "f1: " << f1 << std::endl;
            std::cout << "f2: " << f2 << std::endl;
            std::cout << "f3: " << f3 << std::endl;
#endif
#if !NDEBUG
            std::cout << "Curve cache.." << std::flush;
#endif
            typename 
                Curve_kernel_2::Construct_curve_2 construct_curve = 
                Arrangement_traits_2::instance().kernel().
                construct_curve_2_object();
            Curve_analysis_2 c1 = construct_curve(f1);
            Curve_analysis_2 c2 = construct_curve(f2);
            Curve_analysis_2 c3 = construct_curve(f3);
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif

#if !NDEBUG
            std::cout << "Vert lines.." << std::flush;
#endif            

            typedef typename Curve_analysis_2::Status_line_1 Vert_line;
            Vert_line vl1 = c1.status_line_at_exact_x( x_val ),
                vl2 = c2.status_line_at_exact_x( x_val ),
                vl3 = c3.status_line_at_exact_x( x_val );
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif            

            typedef typename 
                Curve_kernel_2::Curve_pair_analysis_2::Status_line_1 
                CPA_Status_line_1;

           
#if !NDEBUG
            std::cout << "Curve pairs.." << std::flush;
#endif
            typedef typename Curve_kernel_2::Curve_pair_analysis_2 
                Curve_pair_analysis_2;

            typename 
                Curve_kernel_2::Construct_curve_pair_2 construct_curve_pair = 
                Arrangement_traits_2::instance().kernel().
                construct_curve_pair_2_object();

            Curve_pair_analysis_2 p12 = construct_curve_pair( c1, c2 );
            Curve_pair_analysis_2 p13 = construct_curve_pair( c1, c3 );
            Curve_pair_analysis_2 p23 = construct_curve_pair( c2, c3 );
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif

            // We encode the slice information in a int-vector
            // We do so because the interface distinguishes between
            // slices at intervals and event, but we don't want to do so
            // 0 stands for an intersection
#if !NDEBUG
            std::cout << "Slice to int-vec.." << std::flush;
#endif

            std::vector<int> sl12, sl13, sl23;
            
            CPA_Status_line_1 slice_12 = p12.status_line_for_x(x_val);
            for (int i=0; i < slice_12.number_of_events(); i++) {
                std::pair< int, int > arcs 
                    = slice_12.curves_at_event(i,c1,c2);
                if (arcs.first != -1 && arcs.second != -1) {
                    sl12.push_back(0);
                } else if (arcs.first != -1) {
                    sl12.push_back(1);
                } else {
                    sl12.push_back(2);
                }
            }
            CPA_Status_line_1 slice_13 = p13.status_line_for_x(x_val);
            for (int i=0; i < slice_13.number_of_events(); i++) {
                std::pair< int, int > arcs 
                    = slice_13.curves_at_event(i,c1,c3);
                if (arcs.first != -1 && arcs.second != -1) {
                    sl13.push_back(0);
                } else if (arcs.first != -1) {
                    sl13.push_back(1);
                } else {
                    sl13.push_back(3);
                }
            }
            CPA_Status_line_1 slice_23 = p23.status_line_for_x(x_val);
            for (int i=0; i < slice_23.number_of_events(); i++) {
                std::pair< int, int > arcs 
                    = slice_23.curves_at_event(i,c2,c3);
                if (arcs.first != -1 && arcs.second != -1) {
                    sl23.push_back(0);
                } else if (arcs.first != -1) {
                    sl23.push_back(2);
                } else {
                    sl23.push_back(3);
                }
            }
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
#if !NDEBUG
            std::cout << "Merge.." << std::flush;
#endif
            std::vector<std::pair<Vert_line, int> > pseudo_stack;

            // To arrange the vert-line entries, there is a big case distinction
            int vl_id1 = 0, vl_id2 = 0, vl_id3 = 0;
            int vl_n1 = vl1.number_of_events(), 
                vl_n2 = vl2.number_of_events(), 
                vl_n3 = vl3.number_of_events();
            int sl_id12 = 0, sl_id13 = 0, sl_id23 = 0;

            while(vl_id1 < vl_n1 || vl_id2 < vl_n2 || vl_id3 < vl_n3) {

                // Determine which one comes next

                if(vl_id1 == vl_n1) {
                    
                    // 2 or 3
                    if(vl_id2 == vl_n2) {
                        pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                        vl_id3++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                    if(vl_id3 == vl_n3) {
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        sl_id12++;
                        sl_id23++;
                        continue;
                    }
                    CGAL_assertion( sl_id23 < static_cast<int>(sl23.size()) );
                    if(sl23[sl_id23] == 2) {
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        sl_id12++;
                        sl_id23++;
                        continue;
                    } else if(sl23[sl_id23] == 3) {
                        pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                        vl_id3++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    } else {
                        CGAL_assertion( sl23[sl_id23] == 0 );
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        vl_id3++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                }
                if(vl_id2 == vl_n2) {
                    
                    // 1 or 3
                    if(vl_id3 == vl_n3) {
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        sl_id12++;
                        sl_id13++;
                        continue;
                    }
                    CGAL_assertion( sl_id13 < static_cast<int>(sl13.size()) );
                    if(sl13[sl_id13] == 1) {
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        sl_id12++;
                        sl_id13++;
                        continue;
                    } else if(sl13[sl_id13] == 3) {
                        pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                        vl_id3++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    } else {
                        CGAL_assertion( sl13[sl_id13] == 0 );
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        vl_id3++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                }
                if(vl_id3 == vl_n3) {
                    
                    // 1 or 2
                    CGAL_assertion( sl_id12 < static_cast<int>(sl12.size()) );
                    if(sl12[sl_id12] == 1) {
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        sl_id12++;
                        sl_id13++;
                        continue;
                    } else if(sl12[sl_id12] == 2) {
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        sl_id12++;
                        sl_id23++;
                        continue;
                    } else {
                        CGAL_assertion( sl12[sl_id12] == 0 );
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        vl_id2++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                }
                CGAL_assertion( vl_id1 != vl_n1 && 
                                vl_id2 != vl_n2 && 
                                vl_id3 != vl_n3 );

                if(sl12[sl_id12] == 0) {

                    if(sl13[sl_id13] == 0) {
                        // intersection of all three
                        CGAL_assertion(sl23[sl_id23] == 0);
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        vl_id2++;
                        vl_id3++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    } else if(sl13[sl_id13] == 1 ) {
                        CGAL_assertion(sl23[sl_id23] == 2);
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        vl_id2++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    } else {
                        CGAL_assertion(sl13[sl_id13] == 3);
                        CGAL_assertion(sl23[sl_id23] == 3);
                        pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                        vl_id3++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                } else if(sl12[sl_id12] == 1) {
                    if(sl13[sl_id13] == 0) {
                        CGAL_assertion(sl23[sl_id23] == 3);
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        vl_id3++;
                        sl_id12++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    } else if(sl13[sl_id13] == 1 ) {
                        pseudo_stack.push_back(std::make_pair(vl1, vl_id1));
                        vl_id1++;
                        sl_id12++;
                        sl_id13++;
                        continue;
                    } else {
                        CGAL_assertion(sl13[sl_id13] == 3);
                        CGAL_assertion(sl23[sl_id23] == 3);
                        pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                        vl_id3++;
                        sl_id13++;
                        sl_id23++;
                        continue;
                    }
                } else {
                    CGAL_assertion(sl12[sl_id12] == 2);
                    if(sl13[sl_id13] == 0) {
                        CGAL_assertion(sl23[sl_id23] == 2);
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        sl_id12++;
                        sl_id23++;
                        continue;
                    } else if(sl13[sl_id13] == 1 ) {
                        CGAL_assertion(sl23[sl_id23] == 2);
                        pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                        vl_id2++;
                        sl_id12++;
                        sl_id23++;
                        continue;
                    } else {
                        CGAL_assertion(sl13[sl_id13] == 3);
                        if(sl23[sl_id23] == 0) {
                            pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                            vl_id2++;
                            vl_id3++;
                            sl_id12++;
                            sl_id13++;
                            sl_id23++;
                            continue;
                        } else if(sl23[sl_id23] == 2) {
                            pseudo_stack.push_back(std::make_pair(vl2, vl_id2));
                            vl_id2++;
                            sl_id12++;
                            sl_id23++;
                            continue;
                        } else {
                            CGAL_assertion(sl23[sl_id23] == 3);
                            pseudo_stack.push_back(std::make_pair(vl3, vl_id3));
                            vl_id3++;
                            sl_id13++;
                            sl_id23++;
                            continue;
                        }
                    }
                }
            }
            CGAL_assertion( sl_id12 == static_cast<int>(sl12.size()) );
            CGAL_assertion( sl_id13 == static_cast<int>(sl13.size()) );
            CGAL_assertion( sl_id23 == static_cast<int>(sl23.size()) );
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif

#if !NDEBUG
            std::cout << "Isolator.." << std::flush;
#endif
            Z_at_xy_isolator isol(CGAL::CGALi::Vert_line_adapter_descartes_tag(),
                                  pseudo_stack.begin(),
                                  pseudo_stack.end(),
                                  traits);
#if !NDEBUG
            std::cout << "done" << std::flush;
#endif            

            return isol;
            
        }

    public:
        /*!\brief 
         * constructs isolator for \c surface over \c point. 
         * \c feature indicates whether point lies on a projected
         * on silhouette curve of \c surface, i.e., CGAL::FACE means no.
         */
        Z_at_xy_isolator operator()(const Surface_3& surface, 
                                    const Point_2& point,
                                    const CGAL::Nk& nk,
                                    CGAL::Dcel_feature feature) {
#if !NDEBUG
            std::cout << "Construct_isolator for surface " << surface.id()
                      << " at " << point 
                      << " lying on " << feature 
                      << " with mult=" << nk.mult() << ", n=" << nk.n()
                      << ", k=" << nk.k() << " ... " << std::flush;
#endif
            
            // Nk allows access to 
            // * .mult(), degree of silhouette, if feature == EDGE
            // * .n(), local degree: deg f(point,z) in R[z]
            // * .k(), local gcd degree: deg gcd (f(point,z),f'(point,z))
            
#if CGAL_CAD_BENCHMARK_TIMERS
            isol_timer.start();
#endif
            
            
            typename Isolators_map::iterator it = 
                isolators_map().find(surface);
            
            if (it == isolators_map().end()) {
                // if no map exists, insert new one
                it = isolators_map().insert(
                        it, std::make_pair(surface, Isolator_map())
                );
            }
            
            // search for isolator
            typename Isolator_map::iterator pit = it->second.find(point);
            
            // if none is found
            if (pit == it->second.end()) {
                
                // create a new one

                Z_at_xy_isolator isol;

                Polynomial_3 surface_f = surface.f();

                typename Base::Construct_isolator_traits 
                    construct_isolator_traits; // TODO single instance
                Isolator_traits traits = construct_isolator_traits(
                        point, feature == CGAL::VERTEX
                );

                CGAL_assertion(traits.point() == point);
                
                typename Base::Point_on_curve_2 point_on_curve; 
                // TODO single instance                        

#if !NDEBUG
                std::cout << "Compute n.." << std::endl;
#endif
                int n = nk.n();
                // Compute n, if necessary
                if(n == -2) {
                    n = surface_f.degree();
                    while(point_on_curve(traits,surface_f[n])) {
                        n--;
                    }
                }
#if !NDEBUG
                std::cout << "done. n = " << n << std::endl;
#endif
                
                // for backward compatibility
                bool has_vertical_line = (feature == CGAL::VERTEX && n == -1);

                CGAL_assertion(!has_vertical_line || feature == CGAL::VERTEX);

                
                // If there is a vertical line: Special treatment
                if(has_vertical_line) {
#if !NDEBUG
                    std::cout << "Vertical line detected" << std::endl;
#endif
                    isol = _isolation_for_vertical_line_on_surface
                        ( surface, traits );
                    
                } else {
                    
                    Polynomial_3 local_f = surface.f(n);
#if !NDEBUG
                    std::cout << "Local f = " << local_f << std::endl;
#endif

                    int k = nk.k();

                    if (k == 0) {
                        
                        isol = Z_at_xy_isolator
                            ( CGAL::CGALi::Square_free_descartes_tag(),
                              local_f,
                              traits );

                    } else {
                        
                        // Compute 
                        // m = #roots of F_point
                        // k = deg gcd (F_point,F_point')
                        // for F the surface (using subresultants)
                        
#if !NDEBUG
                        std::cout << "Compute m and k..." << std::flush;
#endif
                        
                        
                        bool k_fixed;

                        if(k==-1) {
                            k = ( n == surface.f().degree() && 
                                  feature != CGAL::FACE ) ? 1 : 0;
                            k_fixed = false;
                        } else {
                            k_fixed = true;
                        }

                        std::vector< int > signs;
                        
                        for( int i = 0; i < k; i++) {
                            signs.push_back(0);
                        }

                        for( int i = k; i <= local_f.degree(); i++ ) {
                            Polynomial_2 curr_pol = 
                                surface.
                                principal_sturm_habicht_coefficient(
                                        i, n, true);
                            // Omit point on curve test in first iteration,
                            // if k was already known
                            if(! (i==k && k_fixed) ) {
                                if (point_on_curve(traits, curr_pol)) {
                                    if(!k_fixed) {
                                        k++;
                                    }
                                    signs.push_back(0);
                                    continue;
                                }
                            
                            }                                
                            
                            k_fixed = true;
                            Interval approx 
                                = traits.approximation(curr_pol);
                            while(CGAL::sign(approx.lower()) != 
                                  CGAL::sign(approx.upper())) {
                                traits.refine();
                                approx = traits.approximation(curr_pol);
                            }
                            signs.push_back(
                                    CGAL::sign(approx.lower()) 
                                    == CGAL::NEGATIVE
                                    ? -1 : 1);
                            
                        }                                          
                        
                        // Function ignores leading zeroes
#if CGAL_USE_M_K_DESCARTES
                        int m = 
                            CGAL::CGALi::stha_count_number_of_real_roots<int>(
                                    signs.begin(), signs.end()
                            );


                        
#if !NDEBUG
                        std::cout << "done, m=" 
                                  << m << ", k=" 
                                  << k 
                                  << std::endl;
#endif

#endif
                        
                        // Try m-k-Descartes first
                        try {
                            
#if CGAL_USE_M_K_DESCARTES
#if !NDEBUG
                            std::cout << "Try m-k-Descartes.." << std::flush;
#endif
                            isol = Z_at_xy_isolator(CGAL::CGALi::M_k_descartes_tag(),
                                                    local_f,
                                                    m,k,traits);
#if !NDEBUG
                            std::cout << "success" << std::endl;
#endif
#else
                            throw CGAL::CGALi::Non_generic_position_exception();
#endif
                        } catch (CGAL::CGALi::Non_generic_position_exception 
                                 Ex) {
#if !NDEBUG
                            std::cout << "get cofactor" << std::endl;
#endif
                            // The cofactor of k-1 determines the 
                            // square free polynomial
                            Polynomial_3 sq_free_f_a_b = 
                                surface.cofactor_for_fz(k-1,n);
                            
                            CGAL::CGALi::Square_free_descartes_tag sqfr;
                            
#if !NDEBUG
                            std::cout << "isolate sq-free part" << std::endl;
#endif
                            isol = Z_at_xy_isolator(sqfr,
                                                    sq_free_f_a_b,
                                                    traits);
                        }

                    }
                    
                }
             
                
#if !NDEBUG
                std::cout << "done." << std::endl;
#endif
                // Insert isolator into map
                pit = it->second.insert(pit, std::make_pair(point, isol));       
                
            }
            
#if CGAL_CAD_BENCHMARK_TIMERS
            isol_timer.stop();
#endif
            // return stored value
            return pit->second;
            
            
        }
    };
    
        
    template<typename Poly_2>
    static Interval evaluate_polynomial_2_at_approximated_point(Poly_2 f, 
                                                         Interval x_iv, 
                                                         Interval y_iv) {
        
        CGAL_assertion(f.degree()>=0);
        int n=f.degree();
        Interval ret(0,0);
        for(int i=n;i>=0;i--) {
            ret *= y_iv;
            ret += CGAL::CGALi::evaluate_iv(f[i],x_iv);
        }
        return ret;
    }


    //////////////////////////////////////////////////////////////////////////
    // Intervals

public:    
    
    /*!\brief
     * determines whether two z-intervals are equal special case for quadrics
     */
    class Equal_z {
    public:
        
        Polynomial_3 local_gcd(const Surface_3& surface1,
                               const Surface_3& surface2,
                               const Isolator_traits& traits,
                               int& k) const {
            CGAL_precondition(k >= 1);
#if !NDEBUG
            std::cout << "Compute the gcd..." << std::flush;
#endif
            //!  type of surface pair
            typedef CGAL::Algebraic_surface_pair_3< Self > 
                Algebraic_surface_pair_3;
            
            Algebraic_surface_pair_3 pair = 
                Algebraic_surface_pair_3::surface_pair_cache()(
                        std::make_pair(surface1, surface2)
                );

            // IMPORTANT: pair is used as a cache, 
            //            DO NOT try to access RS_3 members, as there are
            //            computed using equal_z, i.e. this operator
            
            
            // Compute gcd using subresultant sequence (like in AS_3)
            typename Base::Point_on_curve_2 point_on_curve; 
            
            while(true) {
#if !NDEBUG
                std::cout << "sres_" << k << "=" 
                          << pair.get_principal_subresultant(k) << std::endl;
#endif
                if (point_on_curve(traits,
                                   pair.get_principal_subresultant(k))) {
                    k++;
                } else {
                    break;
                }
            }
            
            Polynomial_3 gcd = pair.get_polynomial_subresultant(k);
#if !NDEBUG
            std::cout << "done" << std::endl;
            std::cout << "Local gcd: " << gcd << std::endl;
#endif            
            return gcd;
        }
        
        /*!\brief
         * compute whether \c i1-th interval of \c isolator1 defined by 
         * \c surface1 is equal to \c i2-th interval of \c isolator2 defined
         * by \c surface2, if \c mult is multiplicity of projected point.
         */
        bool operator()(
                const Surface_3& surface1, 
                const Z_at_xy_isolator& isolator1, 
                int i1, 
                bool on_silhouette1,
                const Surface_3& surface2, 
                const Z_at_xy_isolator& isolator2, 
                int i2,
                bool on_silhouette2) const {

            // TODO: be careful if an isolator is "faked" at vertical line
#if !NDEBUG
            std::cout << "Start equal_z.." << std::endl;
#endif

            // REMARK: Only called for overlapping interval pair 
            //         if mult > 1 (edges), or if
            //         adjacence test of vertex failed and at least ONE
            //         silhouette exists!
                  
            CGAL_precondition(surface1 != surface2);

            CGAL_precondition(isolator1.polynomial().degree() > 0);
            CGAL_precondition(isolator2.polynomial().degree() > 0);
            
            // test equalities
            CGAL_assertion(isolator1.traits().point() == 
                           isolator2.traits().point());

            Isolator_traits traits = isolator1.traits();

            int k = 1;
            
            Point_2 point = traits.point();
            

#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
            
            Polynomial_3 local_gcd_p = local_gcd(surface1, surface2, traits,k);
            
            // Question now: Is local_gcd_p square free ?
            
            // if the gcd is of degree one or
            // if at most one silhouette exists at point -> gcd is square-free
            
#if !NDEBUG
            std::cout << "k=" << k 
                      << "\non silhouette 1: " << on_silhouette1 
                      << "\non silhouette 2: " << on_silhouette2 << std::endl;
#endif 
            if (k == 1 || !on_silhouette1 || !on_silhouette2) {
                
#if !NDEBUG
                std::cout << "Evaluate signs of gcd.." << std::endl;
                
                std::cout << "#1..." << std::flush;
#endif
                
                // Does the gcd change its sign inside the isolating intervals?
                
                CGAL::Sign sign_left1, sign_right1;
                
                typedef typename Isolator_traits::Integer Integer;
                
                typedef typename CGAL::Coercion_traits<Polynomial_2,Rational>
                    ::Type Polynomial_rat_2;

                typedef typename CGAL::Fraction_traits<Polynomial_rat_2>
                    FT;

                CGAL_assertion
                    (static_cast<bool>((boost::is_same
                                        < typename FT::Numerator_type,
                                        Polynomial_2 >::value)));

                typename FT::Denominator_type denom;
                
                Polynomial_rat_2 local_gcd_at_left1_with_denom
                    = local_gcd_p.evaluate(isolator1.left_boundary(i1));

                typename FT::Numerator_type local_gcd_at_left1;

                typename FT::Decompose() (local_gcd_at_left1_with_denom,
                                          local_gcd_at_left1,
                                          denom);
                
                Interval approx = traits.approximation(local_gcd_at_left1);
                while (CGAL::sign(approx.lower()) != 
                       CGAL::sign(approx.upper())) {
                    traits.refine();
                    approx = traits.approximation(local_gcd_at_left1);
                }
                
                sign_left1 = CGAL::sign(approx.lower());
                
#if !NDEBUG
                std::cout << "done\n" << "#2..." << std::flush;
#endif
                
                Polynomial_rat_2 local_gcd_at_right1_with_denom
                    = local_gcd_p.evaluate(isolator1.right_boundary(i1));

                typename FT::Numerator_type local_gcd_at_right1;
                
                typename FT::Decompose() (local_gcd_at_right1_with_denom,
                                          local_gcd_at_right1,
                                          denom);
                                
                approx = traits.approximation(local_gcd_at_right1);
                while(CGAL::sign(approx.lower()) != 
                      CGAL::sign(approx.upper())) {
                    traits.refine();
                    approx = traits.approximation(local_gcd_at_right1);
                }
                
                sign_right1 = CGAL::sign(approx.lower());
                
#if !NDEBUG
                std::cout << "done, exiting" << std::endl;
#endif
                
                return sign_left1 != sign_right1;
                
            } else {
                
                // if both silhouettes exists -> gcd may be square-full
                
#if !NDEBUG
                std::cout << "Construct isolator for gcd..." << std::flush;
#endif
                
                // TODO: construct_isolator_object()?
                Construct_isolator construct_isolator;
                
                // no need to activate artificial_x_interval_mode for
                // possible rational x of point
                // TODO: is EDGE correct?
                // both requirements imply EDGE as choice
                
                CGAL::Nk mynk;
                // Set mult to "I don't know"
                mynk.set_mult(0);
                mynk.set_n(local_gcd_p.degree());
                // implicitly set: mynk.set_k(-2);
                
                Z_at_xy_isolator gcd_isol =
                    construct_isolator(local_gcd_p, traits.point(), 
                                       mynk, CGAL::EDGE);
                // false as isolators define finite number of roots
                // -> gcd must also have finite number of roots
                
#if !NDEBUG
                std::cout << "done" << std::endl;
#endif
                
                Rational left_bound = isolator1.left_boundary(i1);
                Rational right_bound = isolator1.right_boundary(i1);
                
                for(int i = 0; i < gcd_isol.number_of_real_roots(); i++) {
                    
                    if (gcd_isol.left_boundary(i) > right_bound ||
                        gcd_isol.right_boundary(i) < left_bound ) {
                        continue;
                    } else {
                        if (gcd_isol.left_boundary(i) >= left_bound &&
                            gcd_isol.right_boundary(i) <= right_bound) {
                            // Completely contained, so its a root of gcd
                            return true;
                        } else {
                            gcd_isol.refine_interval(i);
                        }
                    }
                }
                
                return false;
            }
            
            /* NOT REACHED */ return false; 
        }
    };
    
    /*!\brief
     * returns instance of Equal_z
     */
    Equal_z equal_z_object() const { 
        return Equal_z(); 
    }
    
    /*!\brief
     * Determines the adjacency between two incident cells of the silhouette
     * for a given surface.
     */
    class Adjacency {
    public:
        typedef typename CGAL::Adjacencies_3::Adjacency_pair Adjacency_pair;

        typedef typename CGAL::Adjacencies_3::Adjacency_interval 
        Adjacency_interval;
        
        typedef CGAL::Restricted_cad_3< Self > Restricted_cad_3;
        typedef CGAL::Restricted_cad_3_accessor< Restricted_cad_3 > 
        Accessor;
        
        typedef typename Restricted_cad_3::Rep Arrangement_2;

        typedef typename Restricted_cad_3::Face_const_handle Face_const_handle;
        typedef typename Accessor::Halfedge_const_handle Halfedge_const_handle;
        typedef typename Restricted_cad_3::Vertex_const_handle 
        Vertex_const_handle;

        typedef typename Restricted_cad_3::X_coordinate_1 Algebraic_real_1;
        
        typedef typename 
            CGAL::Fraction_traits<Rational>::Numerator_type Integer;

        typedef typename Restricted_cad_3::Arrangement_traits_2::Curve_kernel_2
            ::Algebraic_kernel_1::Solve_1 Real_roots;

        typedef std::pair<Polynomial_2, int> Root_info;


    private:

        template < class Handle_ >
        struct Handle_less {
            typedef Handle_ Handle;
            bool operator()( Handle h1, Handle h2) {
                //return std::distance(h1, h2) > 0;
                CGAL::Handle_id_less_than< typename Restricted_cad_3::Data > 
                    less;
                return less(*(h1->data()), *(h2->data()));
            }
        };

        // Data about a higher-dimension-object 
        // wrt to a vertex
        struct High_dim_cell_info_for_vertex {
            
            boost::optional<Point_2> sample_point;
            boost::optional<CGAL::Adjacencies_3> adjacencies;

        };

        typedef std::map< Face_const_handle, 
                          High_dim_cell_info_for_vertex, 
                          Handle_less< Face_const_handle > > 
        Facepoints_vertex;

        typedef std::map< Halfedge_const_handle, 
                          High_dim_cell_info_for_vertex,
                          Handle_less< Halfedge_const_handle > > 
        Halfedgepoints;

        struct Vertex_cell_info {
            
            boost::optional<std::vector<Rational> > intermediate_values;
            Facepoints_vertex facepoints;
            Halfedgepoints halfedgepoints;
            
        };

        typedef typename std::map< Vertex_const_handle, 
                                   Vertex_cell_info,
                                   Handle_less< Vertex_const_handle > > 
        Vertex_info;


        // Data about a face wrt to an edge
        struct Face_cell_info_for_edge {
             
            // Really needed? boost::optional<Rational> sample_value;
            boost::optional<CGAL::Adjacencies_3> adjacencies;

        };

        typedef std::map< Face_const_handle, 
                          Face_cell_info_for_edge, 
                          Handle_less< Face_const_handle > > 
        Facepoints_edge;

        struct Edge_cell_info {
            
            boost::optional< Polynomial_2 > planar_curve;
            boost::optional< Algebraic_real_1> sample_value;
            Facepoints_edge facepoints;
            
        };
        

        typedef typename std::map< Halfedge_const_handle, 
                                   Edge_cell_info,
                                   Handle_less< Halfedge_const_handle > > 
        Halfedge_info;

        
        
        static
        Vertex_info& vertex_info() {
            static Vertex_info vinfo;
            return vinfo;
        }

        static
        Halfedge_info& halfedge_info() {
            static Halfedge_info hinfo;
            return hinfo;
        }

        //! Computes adjacencies in case that isolator1 is a m-k-instance
        CGAL::Adjacencies_3 _adjacency_by_m_k_descartes
            (const Z_at_xy_isolator& isolator1,
             const Z_at_xy_isolator& isolator2) const {
            
            CGAL_precondition(isolator1.type() == CGAL::CGALi::M_K_DESCARTES);
            CGAL_precondition(isolator2.type() == CGAL::CGALi::SQUARE_FREE_DESCARTES);
            
            CGAL_precondition( isolator1.polynomial() ==
                               isolator2.polynomial() );

            typename CGAL::Adjacencies_3::Adjacency_vector adj_vec;

            const int 
                n1 = isolator1.number_of_real_roots(),
                n2 = isolator2.number_of_real_roots();
            
            int index1 = 0, index2 = 0;
            
            while (index1 < n1) {
                
                if (isolator1.is_certainly_simple_root(index1)) {
                    adj_vec.push_back(std::make_pair(index1, index2));
                    index1++;
                    index2++;
                } else {
                    for( int i = 0; i < n2 - n1 + 1; i++ ) {
                        adj_vec.push_back(std::make_pair(index1, index2));
                        index2++;
                    }
                    index1++;
                }
            }
            
            CGAL_postcondition(index1 == n1);
            CGAL_postcondition(index2 == n2);
        
            return CGAL::Adjacencies_3(adj_vec);
            
        }

        class Evaluate_homogeneous_x {
            
        public:

            Evaluate_homogeneous_x(Rational r) {
                typename CGAL::Fraction_traits<Rational>::Decompose decompose;
                decompose(r,num,denom);  
            }

            typedef Polynomial_2 argument_type;
            typedef Polynomial_1 result_type;

            Polynomial_1 operator() (Polynomial_2 p) const {
                
                typename CGAL::Polynomial_traits_d<Polynomial_2>::
                    Evaluate_homogeneous evh;
                typename CGAL::Polynomial_traits_d<Polynomial_2>::
                    Move move;
                
                return evh(move(p,0,1),num,denom);
            }
            
        private:
            
            Integer num,denom;

        };

        class Evaluate_homogeneous_y {
            
        public:

            Evaluate_homogeneous_y(Rational r) {
                typename CGAL::Fraction_traits<Rational>::Decompose decompose;
                Integer x_num, x_denom;
                decompose(r,num,denom);  
            }

            typedef Polynomial_2 argument_type;
            typedef Polynomial_1 result_type;

            Polynomial_1 operator() (Polynomial_2 p) const {
                
                typename CGAL::Polynomial_traits_d<Polynomial_2>::
                    Evaluate_homogeneous evh;

                return evh(p,num,denom);
            }
            
        private:
            
            Integer num,denom;

        };


        // needed to transform Curves into Polynomials
        struct Curve_to_pol {
            
            typedef std::pair< Curve_analysis_2, int > argument_type;
            typedef Polynomial_2 result_type;

            Polynomial_2 operator() (argument_type c) const {
                return c.first.polynomial_2();
            }
        };

        // Just for testing
        void print_point(const Point_2 p) const {
            Isolator_traits traits(p);
            typename Isolator_traits::Box box = traits.approximation_square(52);
            std::cout << "(" <<  CGAL::to_double(box.first.upper())
                      << "," <<  CGAL::to_double(box.second.upper())
                      << std::endl;
        }


        class Box_boundary_square_full_exception {
            
            // Empty, just indicates that a polynomial on the the boudary
            // of a box had a multiple root

        };

        template<typename A, typename B>
        class Pair_compare_first {
            
            typedef std::pair<A,B> Pair;

        public:

            typedef Pair first_argument_type;
            typedef Pair second_argument_type;
            typedef bool result_type;
            
            bool operator() (const Pair& a, const Pair& b) const {
                return (a.first < b.first);
            }
            

        };

        template<typename A, typename B>
        class Pair_first_argument {
            
            typedef std::pair<A,B> Pair;

            int arg;

        public:
            
            typedef Pair argument_type;
            typedef A result_type;

            A operator() (const Pair& a) const {
                return a.first;
            }
            

        };


        template<typename A, typename B>
        class Pair_second_argument {
            
            typedef std::pair<A,B> Pair;

            int arg;

        public:
            
            typedef Pair argument_type;
            typedef B result_type;

            B operator() (const Pair& a) const {
                return a.second;
            }
            

        };

        enum Box_side {
            LEFT = 0,
            RIGHT = 1,
            BOTTOM = 2,
            TOP = 3 };


        template<typename InputIterator, 
                 typename RootOutputIterator,
                 typename RootInfoOutputIterator,
                 typename Transformation>
        void _merge_real_roots( InputIterator curves_begin,
                                InputIterator curves_end,
                                RootOutputIterator roots_out,
                                RootInfoOutputIterator root_infos_out,
                                Transformation transform,
                                Box_side box_side,
                                Interval iv) const 
            throw(Box_boundary_square_full_exception)  {
            
            typedef std::vector<std::pair<Algebraic_real_1, Root_info> > 
                Root_vector_with_info;
            
            Root_vector_with_info roots, roots_help, curr_roots_with_info;
            Real_roots real_roots;
            for(InputIterator it = curves_begin; it!=curves_end; it++) {

                Polynomial_1 pol_1 = transform(*it);

                if(! pol_1.is_zero()) {

                    if(! CGAL::CGALi::is_square_free(pol_1) ) {
                        throw Box_boundary_square_full_exception();
                    }

                    std::vector<Algebraic_real_1> curr_roots;
                    
                    real_roots( pol_1, 
                                std::back_inserter(curr_roots) );

                    curr_roots_with_info.clear();
                    for(int i = 0; 
                        i < static_cast<int>(curr_roots.size()); 
                        i++) {
                        Root_info root_info;
                        if(box_side == LEFT || box_side == RIGHT) {
                            root_info = std::make_pair(
                                    CGAL::canonicalize(*it), i);
                        } else {
                            root_info = std::make_pair(
                                    CGAL::canonicalize(*it), -1);
                        }
                        curr_roots_with_info.push_back
                            (std::make_pair(curr_roots[i], root_info));
                    }

                    // Throw away roots not in iv
                    typename Root_vector_with_info::iterator root_it
                        = curr_roots_with_info.begin();
                    while( root_it != curr_roots_with_info.end() ) {

                        if(root_it->first < iv.lower() || 
                           root_it->first > iv.upper() ) {
                            root_it = curr_roots_with_info.erase(root_it);
                        } else {
                            root_it++;
                        }
                    }
                    
                    
                    // Now merge the roots into the "big" container
                    
                    roots_help.clear();
                    
                    std::set_union( curr_roots_with_info.begin(),
                                    curr_roots_with_info.end(),
                                    roots.begin(),
                                    roots.end(),
                                    std::back_inserter(roots_help),
                                    Pair_compare_first<Algebraic_real_1,
                                                       Root_info >() );
                    roots.clear();
                    std::copy( roots_help.begin(), roots_help.end(), 
                               std::back_inserter(roots) );
                    
                    
                }
            }
            
            std::copy( boost::make_transform_iterator
                       ( roots.begin(), 
                         Pair_first_argument
                         <Algebraic_real_1,Root_info>()),
                       boost::make_transform_iterator
                       ( roots.end(), 
                         Pair_first_argument
                         <Algebraic_real_1,Root_info>()),
                       roots_out );
            std::copy( boost::make_transform_iterator
                       ( roots.begin(), 
                         Pair_second_argument
                         <Algebraic_real_1,Root_info>()),
                       boost::make_transform_iterator
                       ( roots.end(), 
                         Pair_second_argument
                         <Algebraic_real_1,Root_info>()),
                       root_infos_out );
            // TODO add postcondition: OIs got same number of elements
        }

        void _edge_adjacency_with_curve_analysis 
            ( const Surface_3& surface,
              const Z_at_xy_isolator& isolator1, 
              const Halfedge_const_handle& he_handle ) const {
            
            Halfedge_info& hinfo = halfedge_info();
            typename Halfedge_info::iterator he_it = hinfo.find(he_handle);
            CGAL_precondition(he_it != hinfo.end());
            

            // Look for a "half-rational" point on the edge.

            Isolator_traits traits = isolator1.traits();
            
            Rational rat_val;


            bool is_vertical = he_handle->curve().is_vertical();

            if(! is_vertical) {

                CGAL_assertion(traits.point().x().is_rational());
                
                rat_val 
                    = traits.point().x().rational();
            } else {

                // Use the defining equation of the point
                Point_2 sample_point 
                    = isolator1.traits().point();
                Curve_analysis_2 sample_curve = sample_point.curve();
                CGAL_assertion(sample_curve.polynomial_2().degree() == 1);
                typename CGAL::Fraction_traits<Rational>::Compose compose;
                CGAL_assertion(sample_curve.polynomial_2()[0].degree()<=0);
                CGAL_assertion(sample_curve.polynomial_2()[1].degree()<=0);
                rat_val = compose(-sample_curve.polynomial_2()[0][0],
                                  sample_curve.polynomial_2()[1][0] );

            }

            if(is_vertical) {            
                
                // easy, take the x-coordinate of the point
                he_it->second.sample_value = traits.point().x();
                

            } else {

                typedef typename CGAL::Coercion_traits<Polynomial_1,Rational>
                    ::Type Polynomial_rat_1;

                typedef typename CGAL::Fraction_traits<Polynomial_rat_1>
                    FT;

                CGAL_assertion
                    (static_cast<bool>((boost::is_same
                                        < typename FT::Numerator_type,
                                        Polynomial_1 >::value)));

                typename FT::Denominator_type denom;
                
                Polynomial_rat_1 p_at_rat_val_with_denom 
                    = typename CGAL::Polynomial_traits_d<Polynomial_2>::Swap()
                        (traits.point().curve().polynomial_2(),0,1)
                            .evaluate(rat_val);

                typename FT::Numerator_type local_sil;

                typename FT::Decompose() (p_at_rat_val_with_denom,
                                          local_sil,
                                          denom);


                Real_roots real_roots;

                std::vector<Algebraic_real_1> y_roots;

                real_roots( local_sil, std::back_inserter(y_roots) );
                
                // Find the "right" y-root for the edge - easy with arcno
                int arcno = traits.point().arcno();

                CGAL_assertion(arcno >=0 && 
                               arcno < static_cast<int> ( y_roots.size() ) );

                he_it->second.sample_value = y_roots[arcno];
            }

            Polynomial_2 surface_section;
            
            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Evaluate_homogeneous evh;
            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Move move;
            
            if(! is_vertical) {
                
                typename CGAL::Fraction_traits<Rational>::Decompose decompose;
                
                Integer x_num, x_denom;
                
                decompose(rat_val,x_num,x_denom); 

                // rat_val is an x-coordinate
                surface_section = 
                    evh(move(surface.f(),0,2),x_num,x_denom);
                

            } else {

                typename CGAL::Fraction_traits<Rational>::Decompose decompose;
                
                Integer x_num, x_denom;
                
                decompose(rat_val,x_num,x_denom); 

                // rat_val is an x-coordinate
                surface_section =
                    evh(move(surface.f(),1,2),x_num,x_denom);
            }
                
            
            he_it->second.planar_curve = 
                CGAL::canonicalize(CGAL::make_square_free(surface_section));

            // Now, use Curve analysis for adjacencies
            typename 
                Curve_kernel_2::Construct_curve_2 construct_curve = 
                Arrangement_traits_2::instance().kernel().
                construct_curve_2_object();

            Curve_analysis_2 section = 
                construct_curve(he_it->second.planar_curve.get());
            
            typename Curve_analysis_2::Status_line_1 vert_line = 
                section.status_line_at_exact_x(
                        he_it->second.sample_value.get()
                );
            
            // Read off the adjacencies from the vert-line 
            
            typename CGAL::Adjacencies_3::Adjacency_vector 
                adj_vec_left, adj_vec_right;

            int face_left_index = 0, face_right_index = 0;

            // First vertical asymptotes from minus infty
            
            int to_minus_from_left, to_minus_from_right, 
                to_plus_from_left, to_plus_from_right;

            std::pair< int, int > minus_inf = 
                vert_line.number_of_branches_approaching_minus_infinity();

            std::pair< int, int > plus_inf = 
                vert_line.number_of_branches_approaching_plus_infinity();
            
            to_minus_from_left = minus_inf.first;
            to_minus_from_right = minus_inf.second;

            to_plus_from_left = plus_inf.first;
            to_plus_from_right = plus_inf.second;
            
            for(int i=0; i<to_minus_from_left; i++) {

                adj_vec_left.push_back(std::make_pair(-1,face_left_index));
                face_left_index++;

            }

            for(int i=0; i<to_minus_from_right; i++) {

                adj_vec_right.push_back(std::make_pair(-1,face_right_index));
                face_right_index++;

            }

            // Now the "regular arcs"

            int num_arcs = vert_line.number_of_events();
            
            for(int i = 0; i < num_arcs; i++) {
                
                std::pair< int, int > arcs = 
                    vert_line.number_of_incident_branches(i);

                for( int j = 0; j < arcs.first; j++ ) {

                    adj_vec_left.push_back
                        (std::make_pair(i,face_left_index));
                    face_left_index++;
                }
                
                for( int j = 0; j < arcs.second; j++ ) {

                    adj_vec_right.push_back
                        (std::make_pair(i,face_right_index));
                    face_right_index++;

                }
                
            }
            
            // Finally, aymptotic arcs again

            for(int i=0; i<to_plus_from_left; i++) {

                adj_vec_left.push_back(std::make_pair(num_arcs,
                                                      face_left_index));
                face_left_index++;

            }

            for(int i=0; i<to_plus_from_right; i++) {

                adj_vec_right.push_back(std::make_pair(num_arcs,
                                                    face_right_index));
                face_right_index++;

            }
            
            if (is_vertical) {
                // important!
                std::swap(adj_vec_right, adj_vec_left);
            }
            
            CGAL::Adjacencies_3 
                adj_this_face(adj_vec_right), 
                adj_twin_face(adj_vec_left);
            
            Face_cell_info_for_edge this_face, twin_face;
            this_face.adjacencies = adj_this_face;
            twin_face.adjacencies = adj_twin_face;
            
            // Now, we have to adjacency vectors, and we have to decide
            // which belongs to which face
            
            Facepoints_edge& facepoints = he_it->second.facepoints;

            if(he_handle->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                // edge goes from left to right, so its *this* face 
                facepoints[he_handle->face()] = this_face;
                facepoints[he_handle->twin()->face()] = twin_face;
            } else {
                // edge goes from right to left, so it is *twin* face
                facepoints[he_handle->face()] = twin_face;
                facepoints[he_handle->twin()->face()] = this_face;
            }
        }
              
        
        CGAL::Adjacencies_3 _edge_face_adjacency
            (const Surface_3& surface,
             const Z_at_xy_isolator& isolator1, 
             CGAL::Object dcel_handle1,
             const Z_at_xy_isolator& isolator2, 
             CGAL::Object dcel_handle2) const {

#if !NDEBUG
            std::cout << "EDGE-FACE-ADJACENCY" << std::endl;
#endif

            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);
            Accessor acc(cad); 
            
            // Look in the cache
            Halfedge_const_handle he_handle;
            CGAL_assertion_code(bool check_edge = )
                    CGAL::assign(he_handle, dcel_handle1);
            CGAL_assertion(check_edge);
            CGAL_assertion(
                    acc.point_on_dcel_handle(isolator1.traits().point(), 
                                             he_handle)
            );
            
            Halfedge_info& hinfo = halfedge_info();
            typename Halfedge_info::iterator he_it = hinfo.find(he_handle);
            if( he_it == hinfo.end() ) {
                // Create entry for Halfedge
                hinfo[he_handle] = Edge_cell_info();
                he_it = hinfo.find(he_handle);
            }
            CGAL_assertion(he_it != hinfo.end());

            Facepoints_edge& facepoints = he_it->second.facepoints;
            Face_const_handle face_handle;
            CGAL_assertion_code(bool check_face = )
                    CGAL::assign(face_handle, dcel_handle2);
            CGAL_assertion(check_face);
            CGAL_assertion(
                     acc.point_on_dcel_handle(isolator2.traits().point(), 
                                              face_handle)
            );
            typename Facepoints_edge::iterator fp_it 
                = facepoints.find(face_handle);
            if( fp_it == facepoints.end() ) {
                // Create entry for Face
                facepoints[face_handle] = Face_cell_info_for_edge();
                fp_it = facepoints.find(face_handle);
            }
            CGAL_assertion(fp_it != facepoints.end());
            
            Face_cell_info_for_edge& face_info = fp_it->second;
            
            if(! face_info.adjacencies) {
                // Compute adjacencies

                // 1) Use m-k-filter:

                if(isolator1.type() == CGAL::CGALi::M_K_DESCARTES &&
                   isolator1.polynomial() == isolator2.polynomial() ) {

                    face_info.adjacencies = _adjacency_by_m_k_descartes
                        (isolator1, isolator2);
                } else {

                    // 2a) Use the AlciX-idea:
                    _edge_adjacency_with_curve_analysis(surface, 
                                                        isolator1,
                                                        he_handle);
                }
            }
            CGAL_assertion(face_info.adjacencies);
            
            CGAL_assertion
                (static_cast<int>(face_info.adjacencies.get().size()) ==
                 isolator2.number_of_real_roots());
            return face_info.adjacencies.get();

        }
        
        template<typename InputIterator>
        int _box_around_point_away_from
        ( Point_2 p,
          InputIterator points_begin,
          InputIterator points_end,
          int initial_bound ) const {
            
            typedef typename Isolator_traits::Box Box;

            // First, create IsolatorTraits for all points

            Isolator_traits vertex_traits(p);

            std::vector<Isolator_traits> traits_points;
            for( InputIterator it = points_begin; it!= points_end; it++ ) {
                traits_points.push_back( Isolator_traits(*it) );
            }

            // Set precision to initial precision
            int bound = initial_bound;

            Box pbox = vertex_traits.approximation_square(bound);
            
            // Iterate through traits with variable it
            typename std::vector<Isolator_traits>::iterator it 
                = traits_points.begin();

            // Now, refine boxes until all boxes are away from p's box
            while(! traits_points.empty() ) {
                CGAL_assertion( it != traits_points.end() );
                Box cbox = it->approximation_square(bound);
                // check for an overlap
                
                bool overlap_x 
                    = ! ( cbox.first.lower() > pbox.first.upper() ||
                          cbox.first.upper() < pbox.first.lower() );
                bool overlap_y 
                    = ! ( cbox.second.lower() > pbox.second.upper() ||
                          cbox.second.upper() < pbox.second.lower() );

                if(! ( overlap_x && overlap_y ) ) {
                    // this point is outside the pbox
                    it = traits_points.erase(it);
                } else {
                    it++;
                }

                if( it == traits_points.end() && ! (traits_points.empty() ) ) {
                    
                    bound++;
                    pbox = vertex_traits.approximation_square(bound);
                    it = traits_points.begin();

                }
                
            }
            
            CGAL_postcondition( bound >= initial_bound );

            return bound;

        }
  

        int _find_arcno(Algebraic_real_1 ar) const {

            Real_roots real_roots;
            std::vector<Algebraic_real_1> roots;
            real_roots(ar.polynomial(), std::back_inserter(roots));

            int n = static_cast<int>(roots.size());
            
            int i = 0;
            
            while(i < n) {

                if( roots[i].high() < ar.low() ) {
                    i++;
                } else if( roots[i].high() >= ar.high() ) {
                    break;
                } else {
                    ar.refine();
                }

            }

            CGAL_postcondition(i<n);
            return i;
        }

        int _exclude_nk_boundary_points_around_vertex
            (const Surface_3& surface,
             Vertex_const_handle v_handle,
             int bound=0) const {

            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            std::vector<Point_2> edge_endpoints;

            // use Halfedge_around_vertex_iterator
            
            typedef typename 
                Restricted_cad_3::Halfedge_around_vertex_const_circulator 
                Halfedges_around_vertex_const_iterator;
            
            Halfedges_around_vertex_const_iterator edge_it 
                = v_handle->incident_halfedges();
            
            Accessor acc(cad);

            for( int i=0; 
                 i < static_cast<int>(v_handle->degree()); 
                 i++, edge_it++ ) {
                
                CGAL_assertion(edge_it->source() != v_handle);

                CGAL_assertion(edge_it->data());
                CGAL::Nk edge_nk = acc.nk(edge_it,surface);
                int n = edge_nk.n();
                int k = edge_nk.k();

                Halfedge_const_handle he_handle = edge_it;

                int number_of_edges = 1; // Needed to handle self loops
                while( (!(he_handle->source()->is_at_infinity())) && 
                       (!(acc.nk(he_handle,surface).n() !=n )) &&
                       (!(acc.nk(he_handle,surface).k() != k)) &&
                       (!(he_handle->source()== v_handle) ) &&
                       (!(he_handle->source()->degree() > 2) ) ) {

                    CGAL_assertion( he_handle->source()->degree() == 2 );
                    he_handle = he_handle->prev();
                    number_of_edges++;

                }

                if( he_handle->source()->is_at_infinity() ) {
                    continue;
                }
                if(he_handle->source() == v_handle) {
                    // Self loop
                    CGAL_assertion(number_of_edges >= 2);
                    he_handle = edge_it;
                    for( int i = 0; i < ( number_of_edges - 1 ) / 2; i++ ) {
                        he_handle = he_handle->prev();
                    }
                    CGAL_assertion(he_handle->source() != v_handle);
                }


                edge_endpoints.push_back(he_handle->source()->point());
                
            }
            
            
            return _box_around_point_away_from( v_handle->point(),
                                                edge_endpoints.begin(),
                                                edge_endpoints.end(),
                                                bound );
        }

        template<typename Arr_2>
        CGAL::Object _locate_point_on_vertex_incident_path(
                const Point_2 p, 
                const Arr_2& arr) const {
#if 0 // search edges
#if !NDEBUG
            std::cout << "Locate point in edges (vertices) .." << std::flush;
#endif
            for (typename Arr_2::Halfedge_const_iterator eit =
                     arr.halfedges_begin(); 
                 eit != arr.halfedges_end(); eit++) {
                if (eit->direction() == CGAL::RIGHT_TO_LEFT) {
                    continue;
                }

                typename Base::Point_on_curve_2
                    point_on_curve;  // TODO single instance
                    
                if (!point_on_curve(p, eit->curve().support().f())) {
                    continue;
                }
                    
                if (eit->curve().equal_y_at_x(p)) {
                    if (!eit->target()->is_at_infinity() &&
                        eit->target()->point() == p) {
#if !NDEBUG
                        std::cout << " done. Found target." << std::endl;
#endif
                        return (CGAL::make_object(eit->target()));
                    }
                    // else 
                    if (!eit->source()->is_at_infinity() &&
                        eit->source()->point() == p) {
#if !NDEBUG
                        std::cout << " done. Found source." << std::endl;
#endif
                        return (CGAL::make_object(eit->source()));
                    }
                    // else 
#if !NDEBUG
                    std::cout << " done. Took edge" << std::endl;
#endif
                    return (CGAL::make_object(eit));
                }
            }

            // NOT REACHED
            CGAL_error_msg("Not allowed to reach here");
            // return empty object
            return CGAL::Object();
#else // use point locations
#if !NDEBUG
            std::cout << "Simple point location .." << std::flush;
#endif
            CGAL::Object obj = 
                CGAL::Arr_naive_point_location<Arr_2>(arr).
                locate( p );
#if !NDEBUG
            std::cout << "done" << std::endl;
#endif
            return obj;
#endif
        }
        
        template<typename Arr_2, 
                 typename Map1, 
                 typename Map2,
                 typename RootInputIterator,
                 typename RootInfoInputIterator >
        void _fill_map_points_for_box_boundary
        (const Arr_2& arr,
         Map1& sample_point_map,
         Map2& sample_point_to_box_boundary_map,
         RootInputIterator box_roots_begin,
         RootInputIterator box_roots_end,
         RootInfoInputIterator box_info_begin,
         RootInfoInputIterator box_info_end,
         Rational rat_val,
         Box_side box_side) const {

            CGAL_assertion(
                    std::distance(box_roots_begin, box_roots_end) ==
                    std::distance(box_info_begin, box_info_end)
            );
            
            RootInfoInputIterator info_it = box_info_begin;

            for (RootInputIterator it = box_roots_begin; 
                 it != box_roots_end; 
                 it++, info_it++) {
                
                Point_2 p;

                if (box_side == LEFT || box_side==RIGHT) {
                    
                    Algebraic_real_1 p_x( rat_val );
                    typename 
                        Curve_kernel_2::Construct_curve_2 construct_curve = 
                        Arrangement_traits_2::instance().kernel().
                        construct_curve_2_object();

                    Curve_analysis_2 p_curve = construct_curve(info_it->first);
                    p = Point_2(p_x, p_curve, info_it->second);
                } else { //box_side == BOTTOM || TOP
                    // TODO construct point on p_curve(info_it->first)
                    p = Accessor::construct_point_with_rational_y
                        ( *it, rat_val );
                }
                
#if 0 // TODO implement point location prevention 
                
                // First filter: Perform a linear search on the edges, 
                // and check whether they have the correct
                // arcno, supporting curve and the x-value in their interior

                typename 
                    Curve_kernel_2::Construct_curve_2 construct_curve = 
                    Arrangement_traits_2::instance().kernel().
                    construct_curve_2_object();
                Curve_analysis_2 support = construct_curve(info_it->first);
                int arcno = info_it->second;

                typename Arr_2::Halfedge_const_iterator edge_it;



                if(arcno == -1) {
                    
                    edge_it = arr.halfedges_end();

                } else {
#if !NDEBUG
                    std::cout << "Search by arcno.." << std::flush;
#endif 
#if 0
                    std::cout << "Search for: "<< std::endl;
                    std::cout << support.id() << std::endl;
                    std::cout << arcno << std::endl;
                    std::cout << CGAL::to_double(p.x()) 
                              << std::endl;
#endif
                    for(edge_it = arr.halfedges_begin();
                        edge_it != arr.halfedges_end();
                        edge_it++) {
                        
                        typename Arr_2::Traits_2::X_monotone_curve_2 segment
                            = edge_it->curve();
#if 0                  
                        std::cout << "Now at: " << std::endl;
                        std::cout << segment << std::endl;
#endif
                        if( segment.arcno() == arcno && 
                            segment.support().id()
                            == support.id() &&
                            segment.is_in_x_range_interior(p.x()) ) {
                            break;
                        }
                    }
                }

                if (edge_it != arr.halfedges_end()) {
#if! NDEBUG
                    std::cout << "ok!" << std::endl;
#endif
                    sample_point_map[edge_it] = p;
#if 1
                    CGAL_assertion_code
                        (
                                CGAL::Object obj 
                                = CGAL::Arr_naive_point_location<Arr_2>(arr).
                                locate( p );      
                                typename Arr_2::Halfedge_const_handle he_handle;
                                if(CGAL::assign(he_handle, obj) ) {
                                    CGAL_assertion(he_handle==edge_it || 
                                                   he_handle==edge_it->twin());
                                }
                        );
#endif
                    
                } else {
#if! NDEBUG
                    if(arcno != -1) {
                        std::cout << "failed" << std::endl;
                    }
#endif
                    
                    // Point location
#if !NDEBUG
                    std::cout << "Point location.." << std::flush;
#endif
                    CGAL::Object obj 
                        = CGAL::Arr_naive_point_location<Arr_2>(arr).
                        locate( p );
#if !NDEBUG
                    std::cout << "done" << std::endl;
#endif
                    CGAL_assertion_code
                        (typename Arr_2::Face_const_handle face_handle);
                    CGAL_assertion(! CGAL::assign(face_handle, obj));
                    
                    typename Arr_2::Vertex_const_handle v_handle;
                    if(CGAL::assign(v_handle, obj) ) {
                        
                        // A vertex, we are only interested in it if it has
                        // degree 2 (otherwise, it can not be part of the
                        // star-arrangement of the vertex
                        if(v_handle->degree()==2) {
                            // Take one of the adjacent edges for the map
                            // (it doesn't matter which one)
                            sample_point_map[v_handle->incident_halfedges()] = p;
                        }
                    } else {
                        typename Arr_2::Halfedge_const_handle he_handle;
                        CGAL_assertion_code(bool check_edge = )
                            CGAL::assign(he_handle, obj);
                        CGAL_assertion(check_edge);
                        sample_point_map[he_handle] = p;
                    }
                }
#endif

                // Point location
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[2].start();
#endif
                CGAL::Object obj = 
                    _locate_point_on_vertex_incident_path(p, arr);
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[2].stop();
#endif
                CGAL_assertion_code
                    (typename Arr_2::Face_const_handle face_handle);
                CGAL_assertion(! CGAL::assign(face_handle, obj));

                typename Arr_2::Vertex_const_handle v_handle;
                if(CGAL::assign(v_handle, obj) ) {

                    // A vertex, we are only interested in it if it has
                    // degree 2 (otherwise, it can not be part of the
                    // star-arrangement of the vertex
                    if(v_handle->degree()==2) {
                        // Take one of the adjacent edges for the map
                        // (it doesn't matter which one)
                        sample_point_map[v_handle->incident_halfedges()] = p;

                    }
                } else {
                    typename Arr_2::Halfedge_const_handle he_handle;
                    CGAL_assertion_code(bool check_edge = )
                        CGAL::assign(he_handle, obj);
                    CGAL_assertion(check_edge);
                    sample_point_map[he_handle] = p;
                }
                
                sample_point_to_box_boundary_map[p] 
                    = std::make_pair(box_side, *it);
                
            }
            
            CGAL_assertion(info_it == box_info_end);
        }
        
        

        Point_2 _face_point_next_to_edge_sample_point
        ( std::pair<Box_side, Algebraic_real_1> box_info,
          const Interval& x_iv,
          const Interval& y_iv,
          const std::vector<Algebraic_real_1>& box_roots_left,
          const std::vector<Algebraic_real_1>& box_roots_right,
          const std::vector<Algebraic_real_1>& box_roots_bottom,
          const std::vector<Algebraic_real_1>& box_roots_top ) const {

            Point_2 sample_point_for_face;
            
            switch(box_info.first) {
                
            case(LEFT): {
                
                if( box_info.second == y_iv.upper() ) {
                    // Upper left corner
                    int n = box_roots_top.size();
                    
                    CGAL_assertion( n >=1);
                    CGAL_assertion( box_roots_top[0] ==
                                    x_iv.upper() );
                    
                    Rational rat_x;
                    
                    if(n == 1) {
                        rat_x = ( x_iv.upper() + x_iv.lower() ) / 2;
                    } else {
                        rat_x = CGAL::CGALi::simple_rational_between
                            ( box_roots_top[0],
                              box_roots_top[1] );
                    }
                    
                    sample_point_for_face 
                        = Accessor::construct_point_with_rational_y
                        ( Algebraic_real_1(rat_x), y_iv.upper() );
                    break;
                }
                
                // Now, the non special case...
                
                typename 
                    std::vector< Algebraic_real_1 >
                    ::const_iterator root_it, root_it_succ;
                root_it = box_roots_left.begin();
                while(root_it != box_roots_left.end()) {
                    if( *root_it == box_info.second ) {
                        break;
                    } else {
                        root_it++;
                    }
                }
                CGAL_assertion( root_it != box_roots_left.end() );
                        
                Rational rat_y; 
                        
                root_it_succ = root_it;
                root_it_succ++;                        

                if(root_it_succ == box_roots_left.end() ) {
                    // Choose the corner
                    rat_y = y_iv.upper();
                } else {
#if !NDEBUG
                    std::cout << "root_it: " << CGAL::to_double(*root_it) 
                              << std::endl;
                    std::cout << "root_it_succ: " 
                              << CGAL::to_double(*root_it_succ) 
                              << std::endl;
#endif
                    rat_y = CGAL::CGALi::simple_rational_between
                        ( *root_it,
                          *root_it_succ );
                }

                sample_point_for_face 
                    = Accessor::construct_point_with_rational_y
                    ( Algebraic_real_1(x_iv.lower()), rat_y );
                break;
                        
            }

            case(RIGHT): {
                if( box_info.second == y_iv.lower() ) {
                    // Lower right corner
                    int n = static_cast<int>(box_roots_bottom.size());
                    CGAL_assertion( n  >= 1 );
                    CGAL_assertion( 
                            box_roots_bottom[n-1] == x_iv.upper() );

                    Rational rat_x;

                    if(n == 1) {
                        rat_x = ( x_iv.upper() + x_iv.lower() ) / 2;
                    } else {
                        rat_x = CGAL::CGALi::simple_rational_between
                            ( box_roots_bottom[n-2],
                              box_roots_bottom[n-1] );
                    }
                            
                    sample_point_for_face 
                        = Accessor::construct_point_with_rational_y
                        ( Algebraic_real_1(rat_x), y_iv.lower() );
                    break;
                }
                        
                // Now, the non special case...
                        
                typename 
                    std::vector< Algebraic_real_1 >
                    ::const_iterator root_it, root_it_prec;
                root_it = box_roots_right.begin();
                while(root_it != box_roots_right.end()) {
                    if( *root_it == box_info.second ) {
                        break;
                    } else {
                        root_it++;
                    }
                }
                CGAL_assertion( root_it != box_roots_right.end() );
                        
                Rational rat_y; 
                        
                if(root_it == box_roots_right.begin() ) {
                    // Choose the corner
                    rat_y = y_iv.lower();
                } else {

                    root_it_prec = root_it;
                    root_it_prec--;                        
                    rat_y = CGAL::CGALi::simple_rational_between
                        ( *root_it_prec,
                          *root_it );
                }

                sample_point_for_face 
                    = Accessor::construct_point_with_rational_y
                    ( Algebraic_real_1(x_iv.upper()), rat_y );
                break;
                        
            }
                     
            case(BOTTOM) : {
                if( box_info.second == x_iv.lower() ) {
                    // Lower left corner
                    int n = static_cast<int>(box_roots_left.size());
                    CGAL_assertion( n  >= 1 );
                    CGAL_assertion( 
                            box_roots_left[0] == y_iv.lower() );

                    Rational rat_y;

                    if(n == 1) {
                        rat_y = ( y_iv.upper() + y_iv.lower() ) / 2;
                    } else {
                        rat_y = CGAL::CGALi::simple_rational_between
                            ( box_roots_left[0],
                              box_roots_left[1] );
                    }
                            
                    sample_point_for_face 
                        = Accessor::construct_point_with_rational_y
                        ( Algebraic_real_1(x_iv.lower()), rat_y );
                    break;
                }
                        
                // Now, the non special case...
                        
                typename 
                    std::vector< Algebraic_real_1 >
                    ::const_iterator root_it, root_it_prec;
                root_it = box_roots_bottom.begin();
                while(root_it != box_roots_bottom.end()) {
                    if( *root_it == box_info.second ) {
                        break;
                    } else {
                        root_it++;
                    }
                }
                CGAL_assertion( root_it != box_roots_bottom.end() );
                        
                Rational rat_x; 
                        
                if(root_it == box_roots_bottom.begin() ) {
                    // Choose the corner
                    rat_x = x_iv.lower();
                } else {

                    root_it_prec = root_it;
                    root_it_prec--;                        
                    rat_x = CGAL::CGALi::simple_rational_between
                        ( *root_it_prec,
                          *root_it );
                }

                sample_point_for_face 
                    = Accessor::construct_point_with_rational_y
                    ( Algebraic_real_1(rat_x), y_iv.lower() );
                break;
                        
            }

            case(TOP) : {
                if( box_info.second == x_iv.upper() ) {
                    // Upper right corner
                    int n = static_cast<int>(box_roots_right.size());
                    CGAL_assertion( n  >= 1 );
                    CGAL_assertion( 
                            box_roots_right[n-1] == y_iv.upper() );

                    Rational rat_y;

                    if(n == 1) {
                        rat_y = ( y_iv.upper() + y_iv.lower() ) / 2;
                    } else {
                        rat_y = CGAL::CGALi::simple_rational_between
                            ( box_roots_right[n-2],
                              box_roots_right[n-1] );
                    }
                            
                    sample_point_for_face 
                        = Accessor::construct_point_with_rational_y
                        ( Algebraic_real_1(x_iv.upper()), rat_y );
                    break;
                }
                        
                // Now, the non special case...
                        
                typename 
                    std::vector< Algebraic_real_1 >
                    ::const_iterator root_it, root_it_succ;
                root_it = box_roots_top.begin();
                while(root_it != box_roots_top.end()) {
                    if( *root_it == box_info.second ) {
                        break;
                    } else {
                        root_it++;
                    }
                }
                CGAL_assertion( root_it != box_roots_top.end() );
                        
                Rational rat_x; 
                        
                root_it_succ = root_it;
                root_it_succ++;

                if(root_it_succ == box_roots_top.end() ) {
                    // Choose the corner
                    rat_x = x_iv.upper();
                } else {

                    rat_x = CGAL::CGALi::simple_rational_between
                        ( *root_it,
                          *root_it_succ );
                }

                sample_point_for_face 
                    = Accessor::construct_point_with_rational_y
                    ( Algebraic_real_1(rat_x), y_iv.upper() );
                break;
            }
   
            } // of switch
                 
            return sample_point_for_face;
        }

        void _create_sample_points_around_non_vertical_vertex
            ( const Surface_3& surface,
              const Z_at_xy_isolator& isolator1, 
              Vertex_const_handle v_handle ) const {

            // We need the cad later for point location and n-k-queries
            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            Vertex_info& v_info = vertex_info();
            typename Vertex_info::iterator v_it = v_info.find(v_handle);

            CGAL_precondition(v_it != v_info.end());
            
            // Intermediate_values known?
            if(! v_it->second.intermediate_values) {
                
                std::vector<Rational> intermediates;
                
                typename Base::Intermediate_values_for_stack 
                    intermediate_func;
                
                intermediate_func( isolator1, 
                                   std::back_inserter(intermediates) );

                v_it->second.intermediate_values = intermediates;
                
            }
            CGAL_assertion(v_it->second.intermediate_values);

            typedef typename 
                CGAL::Coercion_traits<Rational, Polynomial_2>::Type 
                Rat_poly_2;
                        
            int bound = 0;
   
            for( typename std::vector<Rational>::iterator it 
                     = v_it->second.intermediate_values.get().begin();
                 it != v_it->second.intermediate_values.get().end();
                 it++ ) {

                Rat_poly_2 f_at_z = surface.f().evaluate(*it);

                while(true) {
                    std::pair<Interval,Interval> box 
                        = isolator1.traits().approximation_square(bound);
                
                    Interval approx 
                        = evaluate_polynomial_2_at_approximated_point
                        (f_at_z, box.first, box.second);
                    if(approx.lower() <=0 && approx.upper() >=0) {
                        bound++;
                    } else {
                        break;
                    }
                }
                
            }

            if( v_handle->is_isolated() ) {
                
                // We only need a sample point in the unique adjacent face

                // Iterate over ALL silhouette curves
                std::vector< std::pair< Curve_analysis_2, int> > 
                    silhouette_curves;
                typename Base::Construct_projected_surface_curves_2()
                    ( surface, 
                      std::back_inserter(silhouette_curves)
                    );
                Interval y_iv 
                    = isolator1.traits().approximation_square(bound).second;

                Rational y_val = y_iv.upper();

                Algebraic_real_1 x_val = isolator1.traits().point().x();
                    
                typename Base::Point_on_curve_2 point_on_curve; 

                typename Curve_kernel_2::Lower_boundary_y_2
                    lower_boundary_y = 
                    Arrangement_traits_2::instance().kernel().
                    lower_boundary_y_2_object();
                typename Curve_kernel_2::Upper_boundary_y_2
                    upper_boundary_y = 
                    Arrangement_traits_2::instance().kernel().
                    upper_boundary_y_2_object();

                typename Curve_kernel_2::Refine_y_2
                    refine_y = 
                    Arrangement_traits_2::instance().kernel().
                    refine_y_2_object();
                
                for( typename std::vector<std::pair< Curve_analysis_2, int > >
                         ::iterator it
                         = silhouette_curves.begin();
                     it != silhouette_curves.end();
                     it++ ) {
                    
                    typedef typename Curve_analysis_2::Status_line_1 Vert_line;

                    Vert_line vl = it->first.status_line_at_exact_x(x_val);

                    if( point_on_curve(isolator1.traits(),
                                       it->first.polynomial_2()) ) {
                        
                        // Search the arcno of the point
                        int i_down = 0, i_up = vl.number_of_events() - 1;
                        while (i_down != i_up) {
                            
                            while(upper_boundary_y(vl.algebraic_real_2(i_down))
                                  < y_iv.lower() ) {
                                i_down++;
                            } 
                            while(lower_boundary_y(vl.algebraic_real_2(i_up)) 
                                  > y_iv.upper() ) {
                                i_up--;
                            }
                            if(i_down != i_up) {
                                y_iv = isolator1.traits().
                                    approximation_square(++bound).second; 
                            }
                            
                        }
                        CGAL_assertion(i_up == i_down);
                        if( i_up < vl.number_of_events() - 1 ) {
                            
                            y_val = std::min(y_val, 
                                             lower_boundary_y(
                                                     vl.algebraic_real_2(
                                                             i_up+1
                                                     ) 
                                             )
                            );
                            
                        }
                        
                    } else { // !point_on_curve(isolator1.traits(), it->first)
                        
                        // Search the arcno of the point above
                        int i = 0;
                        while (i < vl.number_of_events()) {
                            
                            if (upper_boundary_y(vl.algebraic_real_2(i)) 
                                < y_iv.lower()) {
                                i++;
                            } else if (lower_boundary_y(
                                               vl.algebraic_real_2(i)
                                       ) > y_iv.upper()) {
                                break;
                            } else {
                                
                                if (upper_boundary_y(
                                            vl.algebraic_real_2(i) 
                                    ) - 
                                    lower_boundary_y(
                                            vl.algebraic_real_2(i) 
                                    ) >=
                                    (y_iv.upper() - y_iv.lower())) {
                                    y_iv = isolator1.traits().
                                        approximation_square(++bound).second; 
                                    
                                } else {
                                    refine_y(vl.algebraic_real_2(i));
                                }
                            }
                        }
                        if( i < vl.number_of_events() ) {
                            
                            y_val = std::min(y_val, 
                                             lower_boundary_y(
                                                     vl.algebraic_real_2(i) 
                                             )
                            );
                        }

                    } // end of !point_on_curve

                } // end of iteration over all curves
                
                Point_2 sample_point 
                    = Accessor::construct_point_with_rational_y( x_val,y_val );
                
                CGAL_assertion_code(CGAL::Object obj=cad.locate(sample_point));
                CGAL_assertion_code(Face_const_handle face_handle_2);
                CGAL_assertion_code(bool check_face 
                                    = CGAL::assign(face_handle_2,obj));
                CGAL_assertion(check_face);
                Face_const_handle face_handle = v_handle->face();
                CGAL_assertion( face_handle == face_handle_2 );

                Facepoints_vertex& facepoints = v_it->second.facepoints;
                
                High_dim_cell_info_for_vertex face_cell;
                face_cell.sample_point = sample_point;

                facepoints[face_handle] = face_cell;
                
                return;

                
            } // end of if(v_handle_is_isolated())


            // Now, look for points adjacent to the vertex where the
            // nk-component changes

            bound = _exclude_nk_boundary_points_around_vertex( surface, 
                                                               v_handle, 
                                                               bound);


            
            typename Isolator_traits::Box v_box 
                = isolator1.traits().approximation_square(bound);

            typename Isolator_traits::Interval x_iv = v_box.first,
                y_iv = v_box.second;

#if !NDEBUG            
            std::cout << "x_iv: " << x_iv.lower() << " " << x_iv.upper()
                      << std::endl
                      << "y_iv: " << y_iv.lower() << " " << y_iv.upper()
                      << std::endl;
            

            std::cout << "start root isolation at the box boundaries.." 
                      << std::endl;
#endif

            // Iterate over ALL silhouette curves
            std::vector< std::pair< Curve_analysis_2, int > > 
                silhouette_curves;
            typename Base::Construct_projected_surface_curves_2()
                (surface, 
                 std::back_inserter(silhouette_curves)
                );
            
            std::vector< Algebraic_real_1 > 
                box_roots_left, box_roots_right,
                box_roots_bottom, box_roots_top;

            std::vector< Root_info > box_info_left, box_info_right,
                box_info_bottom, box_info_top;

            while(true) {
                box_roots_left.clear();
                box_roots_right.clear();
                box_roots_bottom.clear();
                box_roots_top.clear();

                box_info_left.clear();
                box_info_right.clear();
                box_info_bottom.clear();
                box_info_top.clear();
                
                try {
                
                    _merge_real_roots 
                        ( boost::make_transform_iterator
                          ( silhouette_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator
                          ( silhouette_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_left),
                          std::back_inserter(box_info_left),
                          Evaluate_homogeneous_x(x_iv.lower()),
                          LEFT,
                          y_iv );
                    CGAL_assertion(
                            box_roots_left.size() == box_info_left.size()
                    );

                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( silhouette_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( silhouette_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_right),
                          std::back_inserter(box_info_right),
                          Evaluate_homogeneous_x(x_iv.upper()),
                          RIGHT,
                          y_iv );
                    CGAL_assertion(
                            box_roots_right.size() == box_info_right.size()
                    );
                    
                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( silhouette_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( silhouette_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_bottom),
                          std::back_inserter(box_info_bottom),
                          Evaluate_homogeneous_y(y_iv.lower()),
                          BOTTOM,
                          x_iv );
                    CGAL_assertion(
                            box_roots_bottom.size() == box_info_bottom.size()
                    );
                    
                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( silhouette_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( silhouette_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_top),
                          std::back_inserter(box_info_top),
                          Evaluate_homogeneous_y(y_iv.upper()),
                          TOP,
                          x_iv );
                    CGAL_assertion(
                            box_roots_top.size() == box_info_top.size()
                    );

#if !NDEBUG
                    std::cout << "done" << std::endl;
#endif
                    break;
                    
                } catch (Box_boundary_square_full_exception ex) {
                    box_roots_left.clear();
                    box_roots_right.clear();
                    box_roots_bottom.clear();
                    box_roots_top.clear();
                    
                    bound++;

                    typename Isolator_traits::Box v_box 
                        = isolator1.traits().approximation_square(bound);
                    
                    x_iv = v_box.first;
                    y_iv = v_box.second;

                }
            }
            
            // Now, the "dirty part" starts

            // A map from halfedge handles to Points, storing the results
            // of the point location
            std::map<
                Halfedge_const_handle, 
                Point_2, 
                Handle_less< Halfedge_const_handle > > sample_point_map;
            
            // Another map to link points with their "origin"
            
            std::map<Point_2, std::pair<Box_side, Algebraic_real_1> >
                sample_point_to_box_boundary_map;

            const Arrangement_2* arr = Accessor(cad).rep();

            _fill_map_points_for_box_boundary
                (*arr,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_left.begin(),
                 box_roots_left.end(),
                 box_info_left.begin(),
                 box_info_left.end(),
                 x_iv.lower(),
                 LEFT);

            _fill_map_points_for_box_boundary
                (*arr,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_right.begin(),
                 box_roots_right.end(),
                 box_info_right.begin(),
                 box_info_right.end(),
                 x_iv.upper(),
                 RIGHT);

            _fill_map_points_for_box_boundary
                (*arr,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_bottom.begin(),
                 box_roots_bottom.end(),
                 box_info_bottom.begin(),
                 box_info_bottom.end(),
                 y_iv.lower(),
                 BOTTOM);

            _fill_map_points_for_box_boundary
                (*arr,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_top.begin(),
                 box_roots_top.end(),
                 box_info_top.begin(),
                 box_info_top.end(),
                 y_iv.upper(),
                 TOP);
            

            // end of map creation

            // Now, go through the edge-paths again, and find for each halfedge
            // the sample-point
            typedef typename 
                Restricted_cad_3::Halfedge_around_vertex_const_circulator 
                Halfedges_around_vertex_const_iterator;
            
            Halfedges_around_vertex_const_iterator edge_it      
                = v_handle->incident_halfedges();

            Vertex_cell_info& cell_info = vertex_info() [v_handle];
            Facepoints_vertex& facepoints = cell_info.facepoints;
            Halfedgepoints& halfedgepoints = cell_info.halfedgepoints;
            
            for( int i=0; 
                 i < static_cast<int>(v_handle->degree()); 
                 i++, edge_it++ ) {

                Halfedge_const_handle he_handle = edge_it;
                while( sample_point_map.find(he_handle) == 
                       sample_point_map.end() ) {
                    
                    CGAL_assertion(he_handle->source()->degree() == 2);
                    he_handle = he_handle->prev();
                    CGAL_assertion( he_handle->face() == edge_it->face() );
                    
                }
                
                Point_2 sample_point = sample_point_map[he_handle];
                
                if(halfedgepoints.find(edge_it) == halfedgepoints.end()) {
                    halfedgepoints[edge_it] = High_dim_cell_info_for_vertex();
                }
                High_dim_cell_info_for_vertex& edge_info 
                    = halfedgepoints[edge_it];
                if(! edge_info.sample_point) {
                    edge_info.sample_point = sample_point;
                }

                // Now for the face left of the edge
                Face_const_handle face_handle = edge_it->face();

                if(facepoints.find(face_handle) == facepoints.end()) {
                    facepoints[face_handle] = High_dim_cell_info_for_vertex();
                }
                High_dim_cell_info_for_vertex& face_info 
                    = facepoints[face_handle];

                if(! face_info.sample_point) {
                    
                    // Where does the edge sample point come from?
                    CGAL_assertion( sample_point_to_box_boundary_map.find
                                       (sample_point) !=
                                    sample_point_to_box_boundary_map.end() );
                    std::pair<Box_side,Algebraic_real_1> box_info
                        = sample_point_to_box_boundary_map[sample_point];
                    
                    face_info.sample_point 
                        = _face_point_next_to_edge_sample_point
                        (box_info,
                         x_iv,
                         y_iv,
                         box_roots_left,
                         box_roots_right,
                         box_roots_bottom,
                         box_roots_top);
   
                }
                
            }
            
        }

        class Non_coprime_exception {

        public:
            
            Non_coprime_exception(int i) : i(i) {}

            int i;
        };

        // Data object for the overlay
        
        class Vertical_overlay_data;


        class Vertical_overlay_data_rep {

        public:

            Vertical_overlay_data_rep() {}

            Vertical_overlay_data_rep( bool on_silhouette,
                                       boost::optional<Rational> bucket_value )
                : _on_silhouette(on_silhouette),
                  _bucket_value(bucket_value)
            {}

        private:

            bool _on_silhouette;
            
            boost::optional<Rational> _bucket_value;

            boost::optional<Halfedge_const_handle> _halfedge_const_handle;

            friend class Vertical_overlay_data;
        };

        
        class Vertical_overlay_data 
            : public CGAL::Handle_with_policy< Vertical_overlay_data_rep >

        {

        public:
            
            typedef CGAL::Handle_with_policy< Vertical_overlay_data_rep > Base;

            Vertical_overlay_data( bool on_silhouette = false,
                                   boost::optional<Rational> bucket_value 
                                   = boost::none ) 
                : Base(on_silhouette, bucket_value)
            {}

            Vertical_overlay_data( const Vertical_overlay_data d1,
                                   const Vertical_overlay_data d2 ) {
                this->ptr()->_on_silhouette
                    = ( d1.on_silhouette() || d2.on_silhouette() );
                this->ptr()->_bucket_value
                    = ( d1.on_bucket_curve() ? d1.ptr()->_bucket_value
                                             : d2.ptr()->_bucket_value );
                this->ptr()->_halfedge_const_handle
                    = (d1.has_halfedge_const_handle() ? 
                       d1.ptr()->_halfedge_const_handle :
                       d2.ptr()->_halfedge_const_handle );
            }

            bool on_silhouette() const {
                return this->ptr()->_on_silhouette;
            }
            bool on_bucket_curve() const {
                return this->ptr()->_bucket_value;
            }
            bool on_both() const {
                return on_silhouette() && on_bucket_curve();
            }
            Rational bucket_value() const {
                CGAL_precondition(on_bucket_curve());
                return this->ptr()->_bucket_value.get();
            }

            bool has_halfedge_const_handle() const {
                return this->ptr()->_halfedge_const_handle;
            }
            
            Halfedge_const_handle halfedge_const_handle() const {
                CGAL_assertion(has_halfedge_const_handle());
                return this->ptr()->_halfedge_const_handle.get();
            }

            void set_halfedge_const_handle(Halfedge_const_handle he_handle) {
                this->ptr()->_halfedge_const_handle = he_handle;
            }


        };


        typedef typename Arrangement_2::Traits_2 Traits_2;

        typedef CGAL::Arr_extended_dcel<Traits_2,Vertical_overlay_data,
                Vertical_overlay_data,Vertical_overlay_data>
        Extended_dcel;
        
        typedef CGAL::Arrangement_2<Traits_2, Extended_dcel> Overlayed_arr_2;
        
        typedef typename Overlayed_arr_2::Vertex_const_handle 
        Overlay_vertex_const_handle;
        
        typedef typename Overlayed_arr_2::Halfedge_const_handle 
        Overlay_halfedge_const_handle;
        
        typedef typename Overlayed_arr_2::Face_const_handle 
        Overlay_face_const_handle;

        typedef typename 
                Overlayed_arr_2::Halfedge_around_vertex_const_circulator 
                Overlay_halfedges_around_vertex_const_iterator;

        // Needed for comparison of halfedeges according to their data
        struct Handle_less_for_overlay_halfedges {
            typedef Overlay_halfedge_const_handle Handle;
            bool operator()( Handle h1, Handle h2) {
                //return std::distance(h1, h2) > 0;
                CGAL::Handle_id_less_than< Vertical_overlay_data > 
                    less;
                return less(h1->data(), h2->data());
            }
        };


        class Vertical_overlay_traits {
            
        public:
            
            typedef typename 
            Arrangement_2::Vertex_const_iterator Vertex_handle_A;
            typedef typename 
            Overlayed_arr_2::Vertex_const_iterator Vertex_handle_B;
            typedef typename 
            Overlayed_arr_2::Vertex_iterator Vertex_handle_R;
            
            typedef typename 
            Arrangement_2::Halfedge_const_iterator Halfedge_handle_A;
            typedef typename 
            Overlayed_arr_2::Halfedge_const_iterator Halfedge_handle_B;
            typedef typename 
            Overlayed_arr_2::Halfedge_iterator Halfedge_handle_R;
            
            typedef typename 
            Arrangement_2::Face_const_iterator Face_handle_A;
            typedef typename 
            Overlayed_arr_2::Face_const_iterator Face_handle_B;
            typedef typename 
            Overlayed_arr_2::Face_iterator Face_handle_R;

            
            /*!\brief
             * Constructs a new overlay traits
             */
            Vertical_overlay_traits
            (boost::optional<Rational> bucket_val,
             boost::optional<typename 
                             Arrangement_2::Vertex_const_handle> opt_v_handle
             = boost::none)
                : _data(! bucket_val, bucket_val), _v_handle(opt_v_handle),
                  _ov_handle(boost::none)
            { }

            //! Destructor
            virtual ~Vertical_overlay_traits() {
            }

        public:    
            
            void reset_bucket_val(boost::optional<Rational> bucket_val) {
                _data = Vertical_overlay_data(!bucket_val, bucket_val);
            }

            typename Overlayed_arr_2::Vertex_const_iterator ov_handle() const {
                CGAL_precondition(_ov_handle);
                return _ov_handle.get();
            }

            /*!
             * Create a vertex v that corresponds to the coinciding vertices 
             * v1 and v2.
             */
            virtual void create_vertex (Vertex_handle_A v1,
                                        Vertex_handle_B v2,
                                        Vertex_handle_R v)
            {
                CGAL_precondition(v1->point() == v->point());
                CGAL_precondition(v2->point() == v->point());
                
                if( (_v_handle && v1==_v_handle.get() ) || 
                    (_ov_handle && v2==_ov_handle.get() ) ) {
                    _ov_handle = v;
                }

                v->set_data(Vertical_overlay_data( _data,v2->data() ));
                
            }
            
            /*!
             * Create a vertex v that mathces v1, which lies of the edge e2.
             */
            virtual void create_vertex (Vertex_handle_A v1,
                                        Halfedge_handle_B e2,
                                        Vertex_handle_R v)
            {
                CGAL_precondition(v1->point() == v->point());
                if(_v_handle && v1==_v_handle.get()) {
                    _ov_handle = v;
                } 
                v->set_data(Vertical_overlay_data( _data,e2->data() ));
            }
            
            /*!
             * Create a vertex v that mathces v1, contained in the face f2.
             */
            virtual void create_vertex (Vertex_handle_A v1,
                                        Face_handle_B f2,
                                        Vertex_handle_R v)
            {
                CGAL_precondition(v1->point() == v->point());
                if(_v_handle && v1==_v_handle.get()) {
                    _ov_handle = v;
                } 
                v->set_data(Vertical_overlay_data( _data,f2->data() ));

            }
    
            /*!
             * Create a vertex v that mathces v2, which lies of the edge e1.
             */
            virtual void create_vertex (Halfedge_handle_A e1,
                                        Vertex_handle_B v2,
                                        Vertex_handle_R v)
            {
                CGAL_precondition(v2->point() == v->point());
                if(_ov_handle && v2==_ov_handle.get()) {
                    _ov_handle = v;
                } 
                Vertical_overlay_data new_data( _data,v2->data() );
                new_data.set_halfedge_const_handle(e1);
                v->set_data(new_data);

            }
            
            /*!
             * Create a vertex v that mathces v2, contained in the face f1.
             */
            virtual void create_vertex (Face_handle_A f1,
                                        Vertex_handle_B v2,
                                        Vertex_handle_R v)
            {
                if(_ov_handle && v2==_ov_handle.get()) {
                    _ov_handle = v;
                } 
                CGAL_precondition(v2->point() == v->point());
                
                Vertical_overlay_data empty_data;
                v->set_data(Vertical_overlay_data( empty_data,
                                                   v2->data() ));

            }
    
            /*!
             * Create a vertex v that mathces the intersection 
             * of the edges e1 and e2.
             */
            virtual void create_vertex (Halfedge_handle_A e1,
                                        Halfedge_handle_B e2,
                                        Vertex_handle_R v)
            {
                Vertical_overlay_data new_data( _data,e2->data() );
                new_data.set_halfedge_const_handle(e1);
                v->set_data(new_data);

            }
    
            /*!
             * Create an edge e that matches the overlap between e1 and e2.
             */
            virtual void create_edge (Halfedge_handle_A e1,
                                      Halfedge_handle_B e2,
                                      Halfedge_handle_R e)
            {
                Vertical_overlay_data new_data( _data, e2->data() );
                new_data.set_halfedge_const_handle(e1);
                e->set_data(new_data);
                e->twin()->set_data(new_data);

            }
    
            /*!
             * Create an edge e that matches the edge e1, 
             * contained in the face f2.
             */
            virtual void create_edge (Halfedge_handle_A e1,
                                      Face_handle_B f2,
                                      Halfedge_handle_R e)
            {
                Vertical_overlay_data new_data( _data, f2->data() );
                new_data.set_halfedge_const_handle(e1);
                e->set_data(new_data);
                e->twin()->set_data(new_data);

            }
    
            /*!
             * Create an edge e that matches the edge e2, 
             * contained in the face f1.
             */
            virtual void create_edge (Face_handle_A f1,
                                      Halfedge_handle_B e2,
                                      Halfedge_handle_R e)
            {
                Vertical_overlay_data empty_data;
                Vertical_overlay_data new_data( empty_data, e2->data() );
                e->set_data(new_data);
                e->twin()->set_data(new_data);

            }
    
            /*!
             * Create a face f that matches the overlapping region 
             * between f1 and f2.
             */
            virtual void create_face (Face_handle_A f1,
                                      Face_handle_B f2,
                                      Face_handle_R f)
            {
                Vertical_overlay_data empty_data;
                f->set_data(Vertical_overlay_data( empty_data,
                                                   f2->data() ));
            }

        private:
            
            Vertical_overlay_data _data;

            boost::optional<typename 
                            Arrangement_2::Vertex_const_iterator> _v_handle;
            boost::optional<typename 
                            Overlayed_arr_2::Vertex_const_iterator> _ov_handle;
            
        };


        
        template<typename InputIterator>
        Overlayed_arr_2 _compute_overlay_for_vertex
            ( const Surface_3& surface,
              Vertex_const_handle v_handle,
              Overlay_vertex_const_handle& ov_handle,
              InputIterator bucket_curves_begin,
              InputIterator bucket_curves_end) const {
            
#if !NDEBUG
            std::cout << "Start computing the overlay" << std::endl;
#endif

            Overlayed_arr_2 arr, arr2;
            
            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            Accessor acc(cad);
            
            const Arrangement_2* sil_arr = acc.rep();

            Vertical_overlay_traits vo_traits
                (boost::none, boost::optional<Vertex_const_handle>(v_handle));

            CGAL::overlay ( *sil_arr, 
                            arr, 
                            arr2, 
                            vo_traits);
            
            arr = arr2;

            typename Vertex_info::iterator v_it = vertex_info().find(v_handle);

            CGAL_precondition(v_it != vertex_info().end());
            
            for( InputIterator it = bucket_curves_begin;
                 it != bucket_curves_end;
                 it++ ) {

                int i = std::distance(bucket_curves_begin,it);
                Rational& z = v_it->second.intermediate_values.get()[i]; 
                Restricted_cad_3 bucket_cad 
                    = Accessor::construct_for_curve(*it);
                Accessor bucket_acc(bucket_cad);
                const Arrangement_2* bucket_arr = bucket_acc.rep();
                boost::optional<Rational> bucket_overlay_traits(z);
                vo_traits.reset_bucket_val(bucket_overlay_traits);

                CGAL::overlay ( *bucket_arr, 
                                arr, 
                                arr2, 
                                vo_traits);
                arr = arr2;
            }
            
#if !NDEBUG
            std::cout << "overlay computed" << std::endl;
#endif
            ov_handle = vo_traits.ov_handle();

            return arr;
        }
        
        template<typename Map1>
        void _compute_adjacencies_for_face_with_bucket_arcs
        ( const Surface_3& surface,
          Vertex_const_handle v_handle,
          Face_const_handle face_handle,
          High_dim_cell_info_for_vertex& face_info,
          Overlay_halfedges_around_vertex_const_iterator oedge_it,
          int no_bucket_arcs,
          Map1& sample_point_map
          
        ) const {

            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            std::vector<std::pair<int,Rational> > 
                patches_to_intervals;
            
            // Check which patches run into an interval
            Overlay_halfedges_around_vertex_const_iterator 
                oedge_it2 = oedge_it;
            oedge_it2++;
            for(int i = 0; i < no_bucket_arcs; i++,oedge_it2++) {
                CGAL_assertion
                    (oedge_it2->data().on_bucket_curve());
                Rational z = oedge_it2->data().bucket_value();
                Overlay_halfedge_const_handle he_handle 
                    = oedge_it2;
                while( sample_point_map.find(he_handle) == 
                       sample_point_map.end() ) {
                    CGAL_assertion
                        ( sample_point_map.find(he_handle->twin())
                          == sample_point_map.end() );
                    CGAL_assertion
                        (he_handle->source()->degree() == 2);
                    he_handle = he_handle->prev();
                    CGAL_assertion(! he_handle->data().on_both() );
                    CGAL_assertion( he_handle->face() == 
                                    oedge_it2->face() );
                }
                
                
                Point_2 bucket_sample_point 
                    = sample_point_map[he_handle];
                
                Z_at_xy_isolator bucket_isol 
                    = Construct_isolator()
                    (surface, 
                     bucket_sample_point,
                     cad.nk(face_handle,surface),
                     CGAL::FACE);
                // z must be in the z-stack
                int n = bucket_isol.number_of_real_roots();
                int j;
                for( j = 0; j < n; j++ ) {

                    if(bucket_isol.left_boundary(j) <= z &&
                       bucket_isol.right_boundary(j) >= z) {
                        break;
                    }

                }
                CGAL_assertion(j < n);
                patches_to_intervals.push_back
                    (std::make_pair(j,z));
                            
                if(i == no_bucket_arcs/2) {
                    face_info.sample_point = bucket_sample_point;
                    // each one would be ok
                }
            }
            CGAL_assertion(face_info.sample_point);
            
            Accessor acc(cad);

            CGAL_assertion(
                    acc.point_on_dcel_handle(*(face_info.sample_point), 
                                             face_handle)
            );

            std::sort(patches_to_intervals.begin(),
                      patches_to_intervals.end());
            patches_to_intervals.erase
                (std::unique( patches_to_intervals.begin(),
                              patches_to_intervals.end() ),
                 patches_to_intervals.end() );

            Z_at_xy_isolator face_isol 
                = Construct_isolator()
                (surface, 
                 face_info.sample_point.get(),
                 cad.nk(face_handle,surface),
                 CGAL::FACE);
            
            typename CGAL::Adjacencies_3::Adjacency_vector adj_vec;
            int nf = face_isol.number_of_real_roots();

            typename 
                std::vector<std::pair<int,Rational> >::iterator
                pti_it = patches_to_intervals.begin();
            int interm_id = 0;
            std::vector<Rational>& interm
                = vertex_info().find(v_handle)->
                second.intermediate_values.get();
            int interm_n = interm.size();
            
            for( int i = 0; i < nf; i++ ) {
                if(interm_id == interm_n) {
                    CGAL_assertion
                        (pti_it == patches_to_intervals.end());
                    adj_vec.push_back(std::make_pair(interm_id-1,i));
                }

                if(pti_it != patches_to_intervals.end() && 
                   pti_it->first == i) {
                    while(pti_it != patches_to_intervals.end() && 
                          pti_it->first == i) {
                        while(interm_id != interm_n &&
                              interm[interm_id]!=pti_it->second) {
                            CGAL_assertion(interm[interm_id]<pti_it->second);
                            interm_id++;
                        }
                        adj_vec.push_back
                            (std::make_pair(interm_id-1,i));
                        adj_vec.push_back
                            (std::make_pair(interm_id,i));
                        pti_it++;
                    }                                
                } else {
                                
                    while(true) {
                                    
                        while(interm_id < interm_n &&
                              face_isol.left_boundary(i) >
                              interm[interm_id]) {
                            interm_id++;
                        }
                        if(interm_id == interm_n) {
                            adj_vec.push_back
                                (std::make_pair(interm_id-1,i));
                            break;
                        }
                        if(face_isol.right_boundary(i) <=
                           interm[interm_id]) {
                            adj_vec.push_back
                                (std::make_pair(interm_id-1,i));
                            break;
                                    
                        } else {
                            face_isol.refine_interval(i);
                        }

                    }
                }
            } // of i-loop
            std::sort(adj_vec.begin(), adj_vec.end());
            adj_vec.erase( std::unique(adj_vec.begin(),adj_vec.end()),
                           adj_vec.end() );
            face_info.adjacencies = CGAL::Adjacencies_3(adj_vec);

            return;
        }


        void _compute_adjacencies_around_vertical_vertex
        (const Surface_3& surface,
         const Z_at_xy_isolator& isolator1,
         Vertex_const_handle v_handle) const {

            // We need the cad later for point location and n-k-queries
            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            Vertex_info& v_info = vertex_info();
            typename Vertex_info::iterator v_it = v_info.find(v_handle);

            CGAL_precondition(v_it != v_info.end());
            
            // Intermediate_values known?
            if(! v_it->second.intermediate_values) {
                
                std::vector<Rational> intermediates;
                
                typename Base::Intermediate_values_for_stack 
                    intermediate_func;

             
                intermediate_func( isolator1, 
                                   std::back_inserter(intermediates) );

                v_it->second.intermediate_values = intermediates;
                
            }
            CGAL_assertion(v_it->second.intermediate_values);

            typedef typename 
                CGAL::Coercion_traits<Rational, Polynomial_2>::Type 
                Rat_poly_2;
            
            bool all_coprime = false;
            
            while(! all_coprime) {
                
                try {
                    
                    for( int i = 0; 
                         i < static_cast<int>
                             (v_it->second.intermediate_values.get().size());
                         i++ ) {
                        Rational& z = v_it->second.intermediate_values.get()[i];
                        Rat_poly_2 f_at_z = surface.f().evaluate(z);
                        int j = i+1;
                        while(j < static_cast<int>
                              (v_it->second.intermediate_values.get().size())) {
                            Rational& z2
                                = v_it->second.intermediate_values.get()[j];
                            CGAL_assertion( z != z2 );
                            Rat_poly_2 f2_at_z = surface.f().evaluate(z2);
                            if(typename CGAL::Polynomial_traits_d<Rat_poly_2>
                               ::Resultant()(f_at_z, f2_at_z).is_zero()) {
                                throw Non_coprime_exception(j);
                            }
                            j++;
                        }
                    }
                    // If we reach this point, all curves are coprime
                    all_coprime=true;

                } catch ( Non_coprime_exception ex) {
                    int i = ex.i;
                    CGAL_assertion(i>0);
                    // Rechoose intermediate value
                    Rational& z = v_it->second.intermediate_values.get()[i];
                    while( z <= isolator1.right_boundary(i-1) ) {
                        isolator1.refine_interval(i-1);
                    }
                    z = ( z + isolator1.right_boundary(i-1) ) / 2;
                }
            }                    
                    

            //Create all bucket curves

            std::vector< Curve_analysis_2 > bucket_curves;

            typedef typename 
                CGAL::Coercion_traits<Rational, Polynomial_2>::Type 
                Rat_poly_2;
            typedef CGAL::Fraction_traits<Rat_poly_2> FT;
            typename FT::Decompose ft_decompose;
            typedef typename FT::Numerator_type Numerator;
            typedef typename FT::Denominator_type Denominator;
            
            int n = static_cast<int>
                (v_it->second.intermediate_values.get().size());

#if !NDEBUG
            std::cout << "Bucket curves: " << std::endl;
#endif
            
            for(int i = 0; i < n; i++) {
                Rational& z = v_it->second.intermediate_values.get()[i];
                Rat_poly_2 rat_poly_2 = surface.f().evaluate(z);
                Polynomial_2 bucket_int_poly;
                Denominator dummy;
                
                ft_decompose(rat_poly_2, bucket_int_poly, dummy);
                typename 
                    Curve_kernel_2::Construct_curve_2 construct_curve = 
                    Arrangement_traits_2::instance().kernel().
                    construct_curve_2_object();
                bucket_curves.push_back(
                        construct_curve(
                                CGAL::make_square_free(bucket_int_poly) 
                        )
                );
#if !NDEBUG
                std::cout << "," << CGAL::make_square_free(bucket_int_poly) 
                          << std::flush;
#endif
            }
#if !NDEBUG
            std::cout << std::endl;
#endif

                        
            int bound = _exclude_nk_boundary_points_around_vertex
                (surface, v_handle );

            // Now, we need the overlay of the silhouette with all 
            // bucket curves
            Overlay_vertex_const_handle ov_handle;

            Overlayed_arr_2 overlay = _compute_overlay_for_vertex
                (surface,v_handle, ov_handle,
                 bucket_curves.begin(), bucket_curves.end());

            // TODO Repair the overlay traits s.t. ov-handle==ov_handle2
            
            // Need to find the vertex in the overlayed arrangement
            // uses point location
#if !NDEBUG
            std::cout << "Locate vertex .." << std::flush;
#endif
            // TODO do not search vertex!
            CGAL::Object obj 
                = CGAL::Arr_naive_point_location<Overlayed_arr_2>(overlay).
                locate( v_handle->point() );
            Overlay_vertex_const_handle ov_handle2;
            CGAL_assertion_code(bool check_vertex = )
                CGAL::assign(ov_handle2, obj);
            CGAL_assertion(check_vertex);
            //CGAL_assertion(ov_handle == ov_handle2);
            ov_handle = ov_handle2;
            
            // TODO Handle isolated vertex
            if( ov_handle->is_isolated() ) {
                CGAL_error_msg("Isol. vertex with vert. line not handled yet");
            }

            std::vector<Point_2> edge_endpoints;

            // use Halfedge_around_vertex_iterator
            
            Overlay_halfedges_around_vertex_const_iterator oedge_it 
                = ov_handle->incident_halfedges();
            


            for( int i=0; 
                 i < static_cast<int>(ov_handle->degree()); 
                 i++, oedge_it++ ) {
                
                CGAL_assertion(oedge_it->source() != ov_handle);

                Overlay_halfedge_const_handle he_handle = oedge_it;

                int number_of_edges = 1; // Needed to handle self loops
                while(  (! he_handle->source()->is_at_infinity()) &&
                        ( he_handle->source()->degree() == 2 )  &&
                        (he_handle->source() != ov_handle) ) {
                    he_handle = he_handle->prev();
                    number_of_edges++;
                    
                }

                if( he_handle->source()->is_at_infinity() ) {
                    continue;
                }
                if(he_handle->source() == ov_handle) {
                    // Self loop
                    CGAL_assertion(number_of_edges >= 2);
                    he_handle = oedge_it;
                    for( int i = 0; i < ( number_of_edges - 1 ) / 2; i++ ) {
                        he_handle = he_handle->prev();
                    }
                    CGAL_assertion(he_handle->source() != ov_handle);
                }


                edge_endpoints.push_back(he_handle->source()->point());
                
            }
            
            
            bound = _box_around_point_away_from( v_handle->point(),
                                                 edge_endpoints.begin(),
                                                 edge_endpoints.end(),
                                                 bound );
            
            typename Isolator_traits::Box v_box 
                = isolator1.traits().approximation_square(bound);

            typename Isolator_traits::Interval x_iv = v_box.first,
                y_iv = v_box.second;
#if !NDEBUG            

            std::cout << "x_iv: " << x_iv.lower() << " " << x_iv.upper()
                      << std::endl
                      << "y_iv: " << y_iv.lower() << " " << y_iv.upper()
                      << std::endl;
            

            std::cout << "start root isolation at the box boundaries.." 
                      << std::endl;
#endif

            // Iterate over ALL silhouette curves
            std::vector< std::pair< Curve_analysis_2, int > > overlay_curves;
            typename Base::Construct_projected_surface_curves_2()
                (surface, 
                 std::back_inserter(overlay_curves)
                );
            // ... and add all bucket curves (with multiplicity one)
            for( typename std::vector< Curve_analysis_2 >::iterator
                     it = bucket_curves.begin();
                 it != bucket_curves.end();
                 it++ ) {
                overlay_curves.push_back(std::make_pair(*it,1));
            }

            std::vector< Algebraic_real_1 > 
                box_roots_left, box_roots_right,
                box_roots_bottom, box_roots_top;

            std::vector< Root_info > box_info_left, box_info_right,
                box_info_bottom, box_info_top;

            while(true) {
                
                try {
                
                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( overlay_curves.begin(),
                              Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( overlay_curves.end(),
                              Curve_to_pol() ),
                          std::back_inserter(box_roots_left),
                          std::back_inserter(box_info_left),
                          Evaluate_homogeneous_x(x_iv.lower()),
                          LEFT,
                          y_iv );
                
                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( overlay_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( overlay_curves.end(),
                              Curve_to_pol() ),
                          std::back_inserter(box_roots_right),
                          std::back_inserter(box_info_right),
                          Evaluate_homogeneous_x(x_iv.upper()),
                          RIGHT,
                          y_iv );
                
                    _merge_real_roots 
                        ( boost::make_transform_iterator
                          ( overlay_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( overlay_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_bottom),
                          std::back_inserter(box_info_bottom),
                          Evaluate_homogeneous_y(y_iv.lower()),
                          BOTTOM,
                          x_iv );
                
                    _merge_real_roots 
                        ( boost::make_transform_iterator 
                          ( overlay_curves.begin(),
                            Curve_to_pol() ),
                          boost::make_transform_iterator 
                          ( overlay_curves.end(),
                            Curve_to_pol() ),
                          std::back_inserter(box_roots_top),
                          std::back_inserter(box_info_top),
                          Evaluate_homogeneous_y(y_iv.upper()),
                          TOP,
                          x_iv );
#if !NDEBUG
                    std::cout << "done" << std::endl;
#endif
                    break;
                    
                } catch (Box_boundary_square_full_exception ex) {
                    box_roots_left.clear();
                    box_roots_right.clear();
                    box_roots_bottom.clear();
                    box_roots_top.clear();
                    
                    bound++;

                    typename Isolator_traits::Box v_box 
                        = isolator1.traits().approximation_square(bound);
                    
                    x_iv = v_box.first;
                    y_iv = v_box.second;

                }

                
            }
            
            
            
            // A map from halfedge handles to Points, storing the results
            // of the point location
            std::map<
                Overlay_halfedge_const_handle, 
                Point_2, 
                Handle_less_for_overlay_halfedges > sample_point_map;
            
            // Another map to link points with their "origin"
            
            std::map<Point_2, std::pair<Box_side, Algebraic_real_1> >
                sample_point_to_box_boundary_map;

            _fill_map_points_for_box_boundary
                (overlay,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_left.begin(),
                 box_roots_left.end(),
                 box_info_left.begin(),
                 box_info_left.end(),
                 x_iv.lower(),
                 LEFT);

            _fill_map_points_for_box_boundary
                (overlay,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_right.begin(),
                 box_roots_right.end(),
                 box_info_right.begin(),
                 box_info_right.end(),
                 x_iv.upper(),
                 RIGHT);

            _fill_map_points_for_box_boundary
                (overlay,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_bottom.begin(),
                 box_roots_bottom.end(),
                 box_info_bottom.begin(),
                 box_info_bottom.end(),
                 y_iv.lower(),
                 BOTTOM);

            _fill_map_points_for_box_boundary
                (overlay,
                 sample_point_map,
                 sample_point_to_box_boundary_map,
                 box_roots_top.begin(),
                 box_roots_top.end(),
                 box_info_top.begin(),
                 box_info_top.end(),
                 y_iv.upper(),
                 TOP);
            

            // end of map creation

            // Traverse the star arrangement and look for sample points

            oedge_it = ov_handle->incident_halfedges();
            
            int arcs_for_silhouette = 0;
            
            for( int i=0; i < static_cast<int>(ov_handle->degree()); i++ ) {
                
                CGAL_assertion(oedge_it->source() != ov_handle);
                CGAL_assertion(! oedge_it->data().on_both());
                if(oedge_it->data().on_silhouette()) {
                    arcs_for_silhouette++;
                }
                oedge_it++;
            }

#if !NDEBUG            
            std::cout << "There are " << arcs_for_silhouette  << " of " 
                      << ov_handle->degree()
                      << " vertex arcs from the silhouette" << std::endl;
#endif

            Halfedgepoints& halfedgepoints = v_it->second.halfedgepoints;
            Facepoints_vertex& facepoints = v_it->second.facepoints;
            
            // Special case: There are no silhouette arcs:
            if(arcs_for_silhouette == 0) {
                CGAL_assertion(v_handle->is_isolated());
                Face_const_handle face_handle = v_handle->face();
                if(facepoints.find(face_handle) == facepoints.end()) {
                    facepoints[face_handle] = High_dim_cell_info_for_vertex();
                }
                High_dim_cell_info_for_vertex& face_info 
                    = facepoints[face_handle];
                CGAL_assertion(! ov_handle->is_isolated());
                oedge_it = ov_handle->incident_halfedges();
                _compute_adjacencies_for_face_with_bucket_arcs
                    (surface,v_handle,face_handle,face_info,oedge_it,
                     ov_handle->degree(),sample_point_map);
                
                
            } else {
                
                
                // We compute the adjacencies for faces directly in this step,
                // if they contain bucket arcs running into the vertex
                // Otherwise, a quite complicated intermediate structure has to
                // be stored for the faces

                for( int i=0; 
                     i < static_cast<int>(ov_handle->degree()); 
                     i++, oedge_it++ ) {

                    CGAL_assertion(! oedge_it->data().on_both());

                    // Only look for silhouette arcs first
                    if(! oedge_it->data().on_silhouette()) {
                        continue;
                    }

                    Overlay_halfedge_const_handle he_handle = oedge_it;
                    while( sample_point_map.find(he_handle) == 
                           sample_point_map.end() ) {
                        CGAL_assertion(sample_point_map.find(he_handle->twin())
                                       == sample_point_map.end() );

                        CGAL_assertion(he_handle->source()->degree() == 2);
                        he_handle = he_handle->prev();
                        CGAL_assertion(! he_handle->data().on_both() );
                        CGAL_assertion( he_handle->face() == oedge_it->face() );
                    
                    }

                    Point_2 sample_point = sample_point_map[he_handle];

                    CGAL_assertion(oedge_it->data().has_halfedge_const_handle());

                    Halfedge_const_handle edge_handle 
                        = oedge_it->data().halfedge_const_handle();
                    CGAL_assertion(edge_handle->source() == v_handle ||
                                   edge_handle->target() == v_handle );
                    if(edge_handle->source() == v_handle) {
                        edge_handle = edge_handle->twin();
                    }
                
                    
                    if(halfedgepoints.find(edge_handle) == halfedgepoints.end()) {
                        halfedgepoints[edge_handle] 
                            = High_dim_cell_info_for_vertex();
                    }
                    High_dim_cell_info_for_vertex& edge_info 
                        = halfedgepoints[edge_handle];
                    if(! edge_info.sample_point) {
                        edge_info.sample_point = sample_point;
                    }

                    // Now for the face left of the edge
                    Face_const_handle face_handle = edge_handle->face();
                
                    if(facepoints.find(face_handle) == facepoints.end()) {
                        facepoints[face_handle] = High_dim_cell_info_for_vertex();
                    }
                    High_dim_cell_info_for_vertex& face_info 
                        = facepoints[face_handle];

                    if(! face_info.sample_point) {

                        // How may bucket arcs are in the face?
                    
                        Overlay_halfedges_around_vertex_const_iterator oedge_it2 
                            = oedge_it;
                        oedge_it2++;
                        int no_bucket_arcs = 0;
                        while(oedge_it2->data().on_bucket_curve()) {
                            CGAL_assertion(! oedge_it2->data().on_both());
                            oedge_it2++;
                            no_bucket_arcs++;
                        }

                        if( no_bucket_arcs == 0 ) {
                        
                            // Where does the edge sample point come from?
                            CGAL_assertion
                                ( sample_point_to_box_boundary_map.find
                                  (sample_point) !=
                                  sample_point_to_box_boundary_map.end() );

                            std::pair<Box_side,Algebraic_real_1> box_info
                                = sample_point_to_box_boundary_map[sample_point];
                    
                            face_info.sample_point 
                                = _face_point_next_to_edge_sample_point
                                (box_info,
                                 x_iv,
                                 y_iv,
                                 box_roots_left,
                                 box_roots_right,
                                 box_roots_bottom,
                                 box_roots_top);

                        } else { // if there are bucekt arcs
                            _compute_adjacencies_for_face_with_bucket_arcs
                                (surface,v_handle,face_handle,face_info,oedge_it,
                                 no_bucket_arcs,sample_point_map);

                        }
                    }
                
                    CGAL_assertion(face_info.sample_point);
                }
                
                
            }
            return;
        }
    
        struct Pair_equal_first_to_i {
            
            Pair_equal_first_to_i(int id) : _id(id) {}

            typedef std::pair<int,int> argument_type;
            typedef bool result_type;
            
            bool operator() (std::pair<int,int> p) {
                return p.first == _id;
            }

        private:
            int _id;
        };

        bool _adjacency_by_transitivity
        (const Surface_3& surface,
         Vertex_const_handle v_handle,
         Z_at_xy_isolator isolator1,
         Halfedge_const_handle he_handle,
         Z_at_xy_isolator isolator2,
         High_dim_cell_info_for_vertex* cell_info) const {
            
            int n = isolator2.number_of_real_roots();

            
            Restricted_cad_3  cad = Restricted_cad_3::cad_cache()(surface); 

            Face_const_handle f1_handle = he_handle->face(),
                f2_handle = he_handle->twin()->face();
            

            CGAL::Object v_obj = CGAL::make_object(v_handle);
            CGAL::Object he_obj = CGAL::make_object(he_handle);
            CGAL::Object f1_obj = CGAL::make_object(f1_handle);
            CGAL::Object f2_obj = CGAL::make_object(f2_handle);
            
            Z_at_xy_isolator isolator_f1 
                = cad.z_stack(f1_handle).isolator(surface),
                isolator_f2 = cad.z_stack(f2_handle).isolator(surface);

            // Get the two edge-to-face adjacencies
#if CGAL_CAD_BENCHMARK_TIMERS
            bool adj0 = adj_timers[0].is_running();
            bool adj1 = adj_timers[1].is_running();
            if (adj0) {
                adj_timers[0].stop();
            }
            if (adj1) {
                adj_timers[1].stop();
            }
#endif
            CGAL::Adjacencies_3 adj_ef1 = 
                this->operator() (surface,
                                  isolator2,CGAL::EDGE,
                                  he_obj, false,
                                  isolator_f1,CGAL::FACE,
                                  f1_obj, false);
            CGAL::Adjacencies_3 adj_ef2 = 
                this->operator() (surface,
                                  isolator2,CGAL::EDGE,
                                  he_obj, false,
                                  isolator_f2,CGAL::FACE,
                                  f2_obj, false);
#if CGAL_CAD_BENCHMARK_TIMERS
            if (adj0) {
                adj_timers[0].start();
            }
            if (adj1) {
                adj_timers[1].start();
            }
#endif

            // We need at least one adjacent patch for each element
            // of he_handle's z-stack

            for( int i = 0; i < n; i++ ) {
                
                Pair_equal_first_to_i comp(i);

                if( std::find_if(adj_ef1.begin(), adj_ef1.end(), comp) ==
                    adj_ef1.end() &&
                    std::find_if(adj_ef2.begin(), adj_ef2.end(), comp) ==
                    adj_ef2.end() ) {
                    // There is no adjacent patch for the i-th stack point
                    return false;
                }

            }

            // If we are here, we know that adjacency is computable
            // by transitivity
            typename CGAL::Adjacencies_3::Adjacency_vector adj_vec;
                           
            // Get the two face-to-vertex adjacencies

#if CGAL_CAD_BENCHMARK_TIMERS
            adj0 = adj_timers[0].is_running();
            adj1 = adj_timers[1].is_running();
            if (adj0) {
                adj_timers[0].stop();
            }
            if (adj1) {
                adj_timers[1].stop();
            }
#endif
            CGAL::Adjacencies_3 adj_f1v = 
                this->operator() (surface,
                                  isolator_f1,CGAL::FACE,
                                  f1_obj, false,
                                  isolator1,CGAL::VERTEX,
                                  v_obj, false);
            CGAL::Adjacencies_3 adj_f2v = 
                this->operator() (surface,
                                  isolator_f2,CGAL::FACE,
                                  f2_obj, false,
                                  isolator1,CGAL::VERTEX,
                                  v_obj, false);
#if CGAL_CAD_BENCHMARK_TIMERS
            if (adj0) {
                adj_timers[0].start();
            }
            if (adj1) {
                adj_timers[1].start();
            }
#endif      
            typedef typename 
                CGAL::Adjacencies_3::Const_adjacency_iterator Adj_iterator;

                               
            for( int i = 0; i < n; i++ ) {
                
                Pair_equal_first_to_i comp(i);
                Adj_iterator it 
                    = std::find_if(adj_ef1.begin(), adj_ef1.end(), comp);

                if( it != adj_ef1.end() ) {

                    CGAL_assertion(it->first == i);
                    Pair_equal_first_to_i fv_comp(it->second);
                    Adj_iterator it_fv 
                        = std::find_if(adj_f1v.begin(), adj_f1v.end(), fv_comp);
                    CGAL_assertion(it_fv != adj_f1v.end());
                    CGAL_assertion(it_fv->first == it->second);
                    adj_vec.push_back(std::make_pair(it_fv->second,i));

                } else {
                    it = std::find_if(adj_ef2.begin(), adj_ef2.end(), comp);
                    CGAL_assertion( it != adj_ef2.end() );
                    CGAL_assertion(it->first == i);
                    Pair_equal_first_to_i fv_comp(it->second);
                    Adj_iterator it_fv 
                        = std::find_if(adj_f2v.begin(), adj_f2v.end(), fv_comp);
                    CGAL_assertion(it_fv != adj_f2v.end());
                    CGAL_assertion(it_fv->first == it->second);
                    adj_vec.push_back(std::make_pair(it_fv->second,i));
                }                     
                
            }
            cell_info->adjacencies = CGAL::Adjacencies_3(adj_vec);
            
            return true;

        }
        

        CGAL::Adjacencies_3 _vertex_adjacency
            (const Surface_3& surface,
             const Z_at_xy_isolator& isolator1, 
             CGAL::Object dcel_handle1,
             bool has_vertical_line1,
             const Z_at_xy_isolator& isolator2, 
             CGAL::Dcel_feature feature2,
             CGAL::Object dcel_handle2) const {
#if !NDEBUG
            std::cout << "VERTEX-ADJACENCY" << std::endl;
#endif
            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);

            Accessor acc(cad);

            // Look in the cache
            Vertex_const_handle v_handle;
            CGAL_assertion_code(bool check_vertex = )
                    CGAL::assign(v_handle, dcel_handle1);
            CGAL_assertion(check_vertex);
            CGAL_assertion(
                     acc.point_on_dcel_handle(isolator1.traits().point(), 
                                              v_handle)
            );

            Vertex_info& vinfo = vertex_info();
            typename Vertex_info::iterator v_it = vinfo.find(v_handle);
            if( v_it == vinfo.end() ) {
                // Create entry for Vertex
                vinfo[v_handle] = Vertex_cell_info();
                v_it = vinfo.find(v_handle);
            }
            CGAL_assertion(v_it != vinfo.end());

            
            High_dim_cell_info_for_vertex* cell_info;

            CGAL::Nk cell_nk;

            Halfedge_const_handle he_handle;

            if( feature2 == CGAL::FACE) {
                Facepoints_vertex& facepoints = v_it->second.facepoints;
                Face_const_handle face_handle;
                CGAL_assertion_code(bool check_face = )
                    CGAL::assign(face_handle, dcel_handle2);
                CGAL_assertion(check_face);
                CGAL_assertion(
                        acc.point_on_dcel_handle(isolator2.traits().point(), 
                                                 face_handle)
                );
                
                cell_nk = cad.nk(face_handle, surface);

                typename Facepoints_vertex::iterator fp_it 
                    = facepoints.find(face_handle);
                if( fp_it == facepoints.end() ) {
                    // Create entry for Face
                    facepoints[face_handle] = High_dim_cell_info_for_vertex();
                    fp_it = facepoints.find(face_handle);
                }
                CGAL_assertion(fp_it != facepoints.end());
                cell_info = &fp_it->second;
            } else { // fesature2==CGAL::EDGE
                Halfedgepoints& halfedgepoints = v_it->second.halfedgepoints;
                
                CGAL_assertion_code(bool check_edge = )
                    CGAL::assign(he_handle, dcel_handle2);
                CGAL_assertion(check_edge);
                CGAL_assertion(
                        acc.point_on_dcel_handle(isolator2.traits().point(), 
                                                 he_handle)
                );
                
                cell_nk = acc.nk(he_handle, surface);

                typename Halfedgepoints::iterator he_it 
                    = halfedgepoints.find(he_handle);
                if( he_it == halfedgepoints.end() ) {
                    // Create entry for Edge
                    halfedgepoints[he_handle] 
                        = High_dim_cell_info_for_vertex();
                    he_it = halfedgepoints.find(he_handle);
                }
                CGAL_assertion(he_it != halfedgepoints.end());
                cell_info = &he_it->second;
                
            }
            
            
            // Adjacencies known?
            if(!cell_info->adjacencies) {

                // M-K-Filter
                //if(false) {
                if(isolator1.type() == CGAL::CGALi::M_K_DESCARTES &&
                      feature2 == CGAL::FACE &&
                      isolator1.polynomial() == isolator2.polynomial() ) {
                    cell_info->adjacencies = _adjacency_by_m_k_descartes
                        (isolator1, isolator2);
                } else { // No MK-Filter
                    
                    // Transitivity filter
                    if(has_vertical_line1 ||
                       feature2 != CGAL::EDGE ||
                       ! _adjacency_by_transitivity(surface,
                                                    v_handle,
                                                    isolator1,
                                                    he_handle,
                                                    isolator2,
                                                    cell_info) ) {
                        
                        //No transitivity filter
                        // Sample point known?
                        if(!cell_info->sample_point) {
#if !NDEBUG
                            std::cout << "Create sample point for vertex ";
                            print_point(v_handle->point()) ;
                            std::cout << std::endl;
#endif
                            if(! has_vertical_line1) {
                                _create_sample_points_around_non_vertical_vertex
                                    (surface, isolator1, v_handle);
                            } else {
                                _compute_adjacencies_around_vertical_vertex
                                    ( surface,
                                      isolator1,
                                      v_handle );
                            }
                        }
                        // Maybe, adjacencies have been computed as side effect
                        if(! cell_info->adjacencies) {
                            CGAL_assertion(cell_info->sample_point);
                            Construct_isolator construct_isolator;
                            
                            
                            Z_at_xy_isolator cell_isol 
                                = construct_isolator
                                (surface, 
                                 cell_info->sample_point.get(),
                                 cell_nk,
                                 feature2);
                            
                            CGAL_assertion(v_it->second.intermediate_values);
                            
                            std::vector<int> roots_in_bucket;
                            
                            typename Base::Assign_roots_to_buckets()
                                ( cell_isol,
                                  v_it->second.
                                  intermediate_values.get().begin(),
                                  v_it->second.
                                  intermediate_values.get().end(),
                                  std::back_inserter(roots_in_bucket)
                                );
                            
                            int s = static_cast<int>(roots_in_bucket.size());
                            CGAL_assertion(s >= 2);
                            
                            typename CGAL::Adjacencies_3::Adjacency_vector 
                                adj_vec;
                            
                            int cell_index = 0;
                            
                            for( int j = 0; j < s; j++ ) {
                                for( int i = 0; i < roots_in_bucket[j]; i++ ) {
                                    
                                    adj_vec.push_back( std::make_pair
                                                       ( j - 1, 
                                                         cell_index ) );
                                    cell_index++;
                                    
                                }
                                
                            }
                            if(! has_vertical_line1) {
                                CGAL_assertion
                                    ( static_cast<int>(adj_vec.size()) == 
                                      isolator2.number_of_real_roots() );
                            }
                            cell_info->adjacencies = 
                                CGAL::Adjacencies_3(adj_vec);
                        }
                    }
                }
                
            }         
            CGAL_assertion(cell_info->adjacencies);
#if !NDEBUG
            std::cout << "Computed adjacencies are: " << std::endl;
            cell_info->adjacencies.get().print();
            std::cout << std::endl;
#endif
            return cell_info->adjacencies.get();

        }

        
    public:
        /*!\brief
         * computes list of adjancy pairs when going on \c surface
         * from \c dcel_handle1 to \c dcel_handle2
         */
        CGAL::Adjacencies_3 operator()(const Surface_3& surface,
                                      const Z_at_xy_isolator& isolator1, 
                                      CGAL::Dcel_feature feature1,
                                      CGAL::Object dcel_handle1,
                                      bool has_vertical_line1,
                                      const Z_at_xy_isolator& isolator2, 
                                      CGAL::Dcel_feature feature2,
                                      CGAL::Object dcel_handle2,
                                      bool has_vertical_line2) const {
            
            CGAL_assertion(feature1 != feature2);
            if(feature2 == CGAL::VERTEX || feature1==CGAL::FACE) {
                CGAL::Adjacencies_3 adj 
                    = this->operator() (surface,
                                        isolator2,feature2,
                                        dcel_handle2,has_vertical_line2,
                                        isolator1,feature1,
                                        dcel_handle1,has_vertical_line1);
                
                return adj.swap();
            }
            
            CGAL_precondition(! has_vertical_line2);

            // Easiest case: One of the z-stacks is empty
            if( isolator1.number_of_real_roots() == 0 ||
                isolator2.number_of_real_roots() == 0 ) {

                return CGAL::Adjacencies_3();

            }
#if !NDEBUG
            std::cout << "Adjacency between " 
                      << ( (feature1 == CGAL::VERTEX) ? "Vertex" : "Edge" )
                      << " and " 
                      << ( (feature2 == CGAL::EDGE) ? "Edge" : "Face" )
                      << std::endl;
            std::cout << "Point 1: ";
            print_point(isolator1.traits().point());
            std::cout << "), Point 2: ";
            print_point(isolator2.traits().point());
            std::cout << ")" << std::endl;
#endif
            if(feature1 == CGAL::EDGE) {
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[0].start();
#endif
                CGAL::Adjacencies_3 adj = 
                    _edge_face_adjacency
                    ( surface, 
                      isolator1, dcel_handle1,
                      isolator2, dcel_handle2
                    );
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[0].stop();
#endif          
                return adj;

            } else {
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[1].start();

#endif
                CGAL::Adjacencies_3 adj = 
                    _vertex_adjacency
                    ( surface, 
                      isolator1, dcel_handle1, has_vertical_line1, 
                      isolator2, feature2, dcel_handle2 
                    );
#if CGAL_CAD_BENCHMARK_TIMERS
                adj_timers[1].stop();
#endif          
                return adj;
            }
             
        }

    };
    
    /*!\brief
     * returns instance of Adjacency
     */
    Adjacency adjacency_object() const { 
        return Adjacency(); 
    }


};

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_Z_AT_XY_ISOLATOR_TRAITS_H
// EOF
