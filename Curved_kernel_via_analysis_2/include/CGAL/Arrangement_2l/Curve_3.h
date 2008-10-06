// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
// File          : include/SoX/GAPS/Curve_3.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef SoX_GAPS_CURVE_3_H
#define SoX_GAPS_CURVE_3_H 1

/*!\file SoX/GAPS/Curve_3.h
 * \brief 
 * definition of \c Curve_3<>
 */

#include <CGAL/config.h>

#include <CGAL/Arrangement_2l/macros.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_functors.h>

#include <CGAL/Arrangement_2l/Surface_pair_3.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class SurfaceZAtXyIsolatorTraits >
class Curve_3_rep : 
        public CGAL::CGALi::Surface_pair_3_rep< SurfaceZAtXyIsolatorTraits > {
    
public:
    //! this instance's template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! this instance itself
    typedef Curve_3_rep< Surface_z_at_xy_isolator_traits > Self;

    //! the base classe
    typedef CGAL::CGALi::Surface_pair_3_rep< Surface_z_at_xy_isolator_traits > 
    Base;
    
    //! type of restricted cad 
    typedef typename Base::Restricted_cad_3 Restricted_cad_3;

    //! type of creator
    typedef typename Base::Creator Creator;

    //! type of creator
    typedef typename Base::Overlayer Overlayer;
    
public:
    //!\name Constructors
    //!@{
    
    Curve_3_rep(const Surface_3& surface1, const Surface_3& surface2) :
        Base(surface1, surface2) {
    }
    
    //!@}
    
    // add special members??

};

} // namespace CGALi

template < 
class SurfaceZAtXyIsolatorTraits, 
class Rep_ = CGAL::CGALi::Curve_3_rep< SurfaceZAtXyIsolatorTraits >
>
class Curve_3 : 
        public Surface_pair_3<SurfaceZAtXyIsolatorTraits, Rep_ > {
    
public:
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;

    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of Base
    typedef CGAL::Surface_pair_3< Surface_z_at_xy_isolator_traits, Rep > Base;

    //! type of restricted cad
    typedef typename Base::Restricted_cad_3 Restricted_cad_3;

    //!\name Constructors
    //!@{
    
    Curve_3(const Surface_3& surface1, const Surface_3& surface2) :
        Base(surface1, surface2) {
    }
    
    //!@}
    
private: // make public members inaccessible
    //!\name Cads
    //!@{

    //! rscad of first silhouette
    Restricted_cad_3 silhouette1() const {
        return Base::silhouette1();
    }

    //! rscad of second silhouette
    Restricted_cad_3 silhouette2() const {
        return Base::silhouette1();
    }

    //! rscad of cut
    Restricted_cad_3 cut() const {
        return Base::cut();
    }
    
    //! rscad of first silhouette with cut
    Restricted_cad_3 silhouettes_cut() const {
        return Base::silhouettes_cut();
    }
    
    //!@}
    
private: // make public members inaccessible
    //!\name Curves
    //!@{

    /*!\brief
     * returns silhouettes curves and points for a given \c surface in pair
     */
    template < class CurveOutputIterator, class PointOutputIterator >
    void silhouette_objects(
            const Surface_3& surface, 
            CurveOutputIterator coi, PointOutputIterator poi) const {
        return Base::silhouette_objects(surface, coi, poi);

    }

    /*!\brief
     * returns projected intersection curves and points of the pair
     */
    template < class CurveOutputIterator, class PointOutputIterator >
    void cut_objects(
            CurveOutputIterator coi, PointOutputIterator poi) const {
        return Base::cut_objects(coi, poi);
    }
    
    //!@}

public:

    //!\name Space Curve Objects
    //!@{
    /*!\brief
     * returns spatial curves/points as combination of curve/sheet or
     * point/sheet, where sheet is an int wrt to smaller-degree surface
     */
    template <
        class SurfaceArc_3,
        class CurveOutputIterator, 
        class PointOutputIterator
    >
    void objects(
            SurfaceArc_3 dummy,
            Surface_3& support,
            CurveOutputIterator coi, PointOutputIterator poi
    ) const {
        
        typedef SurfaceArc_3 Surface_arc_3;
        typedef typename Surface_arc_3::Surface_point_2l Surface_point_3;

        typedef typename Restricted_cad_3::Vertex_const_iterator
            Vertex_const_iterator;
        typedef typename Restricted_cad_3::Edge_const_iterator
            Edge_const_iterator;
        typedef typename Restricted_cad_3::Z_stack Z_stack;
        
        // Remark (1): 
        // The simple case is an edge with mult == 1, ask member for it
        // -> Refine until only one z-cell of f remains (f is simpler than g!)
        // Else: Multiple intersection can occur.
        // Either use local gcd (point_on_curve for SR_i, i > 0!!!, as point
        // already lies on SR_0,
        // or use silhouette of g (costly?!?!) to determine isolator :-(

        Surface_3 other;
        
        // TODO select smaller surface out of the two by traits?!?!
        if (this->surface2().f().degree() < this->surface1().f().degree()) {
            support = this->surface2();
            other = this->surface1();
        } else {
            support = this->surface1();
            other = this->surface2();
        }
        
        // TODO support is not allowed to be vertical!
        CGAL_assertion(support.f().degree() > 0);
        if (support.f().degree() != CGAL::total_degree(support.f())) {
            // TODO what if support is not z-regular
            // and what if other is not z-regular??!
            std::cerr << "Warning: Support is not z-regular, "
                      << "we might encounter problems with vertical "
                      << "objects, or asympotes!" 
                      << std::endl;
        }
        
        Restricted_cad_3 cad = 
            (support.id() == this->surface1().id() ? 
             this->_silhouette1_cut() :
             this->_silhouette2_cut());
        // we won't access z-stacks of cad directly, only for support!

        CGAL::CGALi::Refineable_interval_helper< Z_at_xy_isolator > 
            iv_helper;
        
        // iterate over all features of cut
        for (Edge_const_iterator eit = cad.edges_begin();
             eit != cad.edges_end(); eit++) {
            
            if (cad.has_cut(eit)) {
                
                boost::optional< int > mult =
                    eit->data()->multiplicity_of_cut(support, other);
                CGAL_assertion(mult);
                
                Point_2 point = cad.sample_point(eit);
                
                typename 
                    Surface_z_at_xy_isolator_traits::Construct_isolator
                    construct_isolator; // TODO single instance
                
                Z_at_xy_isolator support_isolator =
                    construct_isolator(support, point, 
                                       cad.nk(eit, support),
                                       (cad.has_silhouette(eit) ?
                                        CGAL::EDGE : CGAL::FACE)
                    );
                
                int n = support_isolator.number_of_real_roots();
                if (n == 0) {
                    continue;
                }
                
                std::list< int > sheets;

                if (*mult == 1) {
                    
                    Z_at_xy_isolator other_isolator =
                        construct_isolator(
                                other, point, 
                                cad.nk(cad.faces_begin(), support),
                                CGAL::FACE
                        );
                    std::cout << "RRs" 
                              << support_isolator.number_of_real_roots()
                              << std::endl;


                    std::cout << "RRo" 
                              << other_isolator.number_of_real_roots()
                              << std::endl;

                    CGAL_assertion(other_isolator.number_of_real_roots() > 0);
                    
                    CGAL_assertion(point == other_isolator.traits().point());

                    std::list< std::pair< int, int > > overlaps;
                    
                    iv_helper.refined_overlaps(
                            support_isolator, other_isolator,
                            std::back_inserter(overlaps)
                    );
                    
                    CGAL_assertion(!overlaps.empty());

                    std::pair< int, int > overlap =
                        iv_helper.unique_overlap(
                                support_isolator, other_isolator,
                                overlaps.begin(), overlaps.end()
                        );
                    
                    sheets.push_back(overlap.first);
                    
                } else {
                    
                    typename 
                        Surface_z_at_xy_isolator_traits::Equal_z
                        equal_z; // TODO single instance

                    typename 
                        Surface_z_at_xy_isolator_traits::Point_on_curve_2
                        point_on_curve; // TODO single instance
                    
                    // can start with 1 as it lies on k = 0
                    int k = 1;

                    // TODO what happens if surfaces are not z-regular??!
                    Polynomial_3 local_gcd = 
                        equal_z.local_gcd(support, other,
                                          support_isolator.traits(), k);
                    
                    // TODO make use of k
                    
                    // TODO be careful if a surface has a vertical line!

                    Surface_3 gcd_surface = 
                        Surface_3::surface_cache()(local_gcd);

                    Polynomial_2 sil_local_gcd = gcd_surface.resultant_f_fz();
                    
                    bool on_gcd_sil = point_on_curve(point, sil_local_gcd);
                    
                    Z_at_xy_isolator gcd_isolator =
                        construct_isolator(
                                local_gcd, point,
                                cad.nk(eit, support),
                                (on_gcd_sil ? CGAL::EDGE : CGAL::FACE)
                        );
                    
                    // compute overlaps
                    std::vector< std::pair< int, int > > overlaps;
                    
                    iv_helper.refined_overlaps(
                            support_isolator, gcd_isolator,
                            std::back_inserter(overlaps)
                    );
                    
                    // refine until a gcd-interval is fully include in
                    // support interval or has moved away
                    for (int k = 0; k < static_cast< int >(overlaps.size());
                         k++) {
                        
                        while (true) {

                            if (iv_helper.overlap_or_order(
                                        gcd_isolator, overlaps[k].second,
                                        support_isolator, overlaps[k].first
                                ) != CGAL::EQUAL) {
                                break;
                            }
                            // else
                            
                            if (iv_helper.is_included(
                                        gcd_isolator, overlaps[k].second,
                                        support_isolator, overlaps[k].first,
                                        true
                                )
                            ) {
                                sheets.push_back(overlaps[k].first);
                                break;
                            }
                            
                            // else 
                            gcd_isolator.refine_interval(overlaps[k].second);
                        }
                    }
                }
                
                for (std::list< int >::const_iterator sit = sheets.begin();
                     sit != sheets.end(); sit++) {
                    
                    //std::cout << "cv: " << eit->curve() << std::endl;
                    
                    bool min_finite = 
                        eit->curve().is_finite(CGAL::ARR_MIN_END);
                    bool max_finite = 
                        eit->curve().is_finite(CGAL::ARR_MAX_END);
                    
                    bool src_at_min = 
                        (eit->direction() == CGAL::ARR_LEFT_TO_RIGHT);
                    
                    Surface_point_3 min, max;

                    int sheet = *sit;
                    int sheet_at_min = -2;
                    int sheet_at_max = -2;
                    bool z_inf_at_min = false;
                    bool z_inf_at_max = false;

                    if (min_finite) {
                        Point_2 pmin = 
                            eit->curve().curve_end(CGAL::ARR_MIN_END);
                        
                        std::pair< Z_stack, std::pair< int, int > > adjmin =
                            cad.surface_adjacency(eit, support, sheet, 
                                                  (src_at_min ? 
                                                   eit->source() : 
                                                   eit->target()));
                        
                        sheet_at_min = adjmin.second.first;
                        z_inf_at_min = 
                            (sheet_at_min == -1 || sheet_at_min == n);
                        if (!z_inf_at_min) {
                            min = Surface_point_3(pmin, support, sheet_at_min);
                        }
                    }
                    
                    if (max_finite) {
                         Point_2 pmax = 
                             eit->curve().curve_end(CGAL::ARR_MAX_END);
                         
                         std::pair< Z_stack, std::pair< int, int > > adjmax =
                             cad.surface_adjacency(eit, support, sheet, 
                                                   (src_at_min ? 
                                                    eit->target() :
                                                    eit->source()));
                         
                         sheet_at_max = adjmax.second.first;
                         z_inf_at_max = 
                            (sheet_at_max == -1 || sheet_at_max == n);
                         if (!z_inf_at_max) {
                             max = Surface_point_3(pmax, support, sheet_at_max);
                         }
                    }
                    
                    if (min_finite && max_finite) {
                        // arc is bounded
                        
                        if (z_inf_at_min && z_inf_at_max) {
                            // no z-asymptote
                            *coi++ = Surface_arc_3(
                                    eit->curve(), min, max, support,
                                    sheet, sheet_at_min, sheet_at_max
                            );
                        } else if (!z_inf_at_min && !z_inf_at_max) {
                            // z-asymptotes at both sides
                            *coi++ = Surface_arc_3(
                                    eit->curve(),
                                    (sheet_at_min == -1 ? 
                                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END),
                                    (sheet_at_max == -1 ? 
                                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END),
                                    support,
                                    sheet
                            );
                        } else {
                            // z-asymptote at one-side
                            int sheet_at_point = sheet_at_min;
                            if (!z_inf_at_min) {
                                sheet_at_point = sheet_at_max;
                            }
                            *coi++ = Surface_arc_3(
                                    eit->curve(), 
                                    (z_inf_at_min ? max : min),
                                    (sheet_at_point == -1 ? 
                                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END),
                                    support,
                                    sheet,
                                    (z_inf_at_min ? sheet_at_max : sheet_at_min)
                            );
                        }
                    } else if (!min_finite && !max_finite) {
                        // planar arc is a branch
                        *coi++ = Surface_arc_3(
                                eit->curve(), 
                                support,
                                sheet
                        );
                    } else {
                        // arc is ray
                        if (!z_inf_at_min && !z_inf_at_max) {
                            // usual one
                            *coi++ = Surface_arc_3(
                                    eit->curve(), 
                                    (min_finite ? min : max), 
                                    support,
                                    sheet, 
                                    (min_finite ? sheet_at_min : sheet_at_max)
                            );
                            
                        } else {
                            int sheet_at_point = sheet_at_min;
                            if (!z_inf_at_min) {
                                sheet_at_point = sheet_at_max;
                            }
                            // z-asymptote at finite-side
                            *coi++ = Surface_arc_3 (
                                    eit->curve(), 
                                    (sheet_at_point == -1 ? 
                                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END),
                                    support,
                                    sheet
                            );
                        }
                    }
                }
            } 
        }
        
        // search only for isolated points in support and 
        // check only them with the help of local gcd!
        
        // iterate over all features of cut
        for (Vertex_const_iterator vit = cad.vertices_begin();
             vit != cad.vertices_end(); vit++) {
            if (vit->is_at_infinity()) {
                continue;
            }
            if (cad.has_cut(vit)) {
                
                // check whether vertex is isolated or whether support
                // has an isolated event here
                
                // actually we test whether each sheet k at vit is connected
                // to something
                
                Point_2 point = cad.sample_point(vit);
                
                typename 
                    Surface_z_at_xy_isolator_traits::Construct_isolator
                    construct_isolator; // TODO single instance
                
                Z_at_xy_isolator support_isolator =
                    // edge is ok, as we don't want to activate artifical
                    // x-interval, but other isolators (m-k..)
                    construct_isolator(
                            support, point, 
                            cad.nk(vit, support),
                            (cad.has_silhouette(vit) ? 
                             CGAL::EDGE : CGAL::FACE)
                    );
                
                int n = support_isolator.number_of_real_roots();
                
                if (cad.nk(vit, support).n() == -1) { // vertical line
                    // TODO sufficient condition? actually both surfaces
                    // need to be vertical, right?
                    // better criterion: both polynomials vanish?
                    
                    if (n == 0) {
                        // whole line
                        *coi++ = Surface_arc_3(point, support);
                    } else {
                        CGAL_assertion(n > 0);
                        Surface_point_3 min, max;
                        min = Surface_point_3(point, support, 0);
                        for (int sheet = 0; sheet < n; sheet++) {
                            if (sheet == 0) {
                                // ray from z=-oo
                                *coi++ = Surface_arc_3(min, 
                                                       CGAL::ARR_MIN_END,
                                                       support);
                            }
                            if (sheet + 1 == n) {
                                // ray to z==oo
                                *coi++ = Surface_arc_3(min, 
                                                       CGAL::ARR_MAX_END,
                                                       support);
                                break;
                            }
                            CGAL_assertion(sheet + 1 < n);
                            max = Surface_point_3(point, support, sheet + 1);
                            // bounded vertical arc
                            *coi++ = Surface_arc_3(min, max, support);
                            min = max;
                        }
                    }
                    return; // as no further isolated vertices will occur
                }

                // else

                if (n == 0) {
                    continue;
                }
                
                std::set< int > sheets;
                
                if (!vit->is_isolated()) {
                    for (int k = 0; k < n; k++) {
                        typename Restricted_cad_3::
                            Halfedge_around_vertex_const_circulator circ, 
                            start =
                            vit->incident_halfedges();
                        circ = start;
                        bool isolated = true;
                        do {
                            ++circ;
                            
                            std::pair< Z_stack, std::pair< int, int > > 
                                adjmin =
                                cad.surface_adjacency(vit, support, k, 
                                                      circ);
                            
                            if (adjmin.second.first <= adjmin.second.second) {
                                isolated = false;
                                break;
                            }
                            
                        } while (circ != start);
                        
                        if (isolated) {
                            sheets.insert(k);
                        }
                    }
                    
                    if (sheets.empty()) {
                        continue;
                    }
                } else {
                    for (int k = 0; k < n; k++) {
                        sheets.insert(k);
                    }
                }
                
                typename 
                    Surface_z_at_xy_isolator_traits::Equal_z
                    equal_z; // TODO single instance
                
                typename 
                    Surface_z_at_xy_isolator_traits::Point_on_curve_2
                    point_on_curve; // TODO single instance
                
                int k = 1;
                
                Polynomial_3 local_gcd = 
                    equal_z.local_gcd(support, other,
                                      support_isolator.traits(), k);
                
                // TODO make use of k
                
                // TODO be carefull if one of the surfaces has a vertical line

                Surface_3 gcd_surface = 
                    Surface_3::surface_cache()(local_gcd);
                
                Polynomial_2 sil_local_gcd = gcd_surface.resultant_f_fz();
                
                bool on_gcd_sil = point_on_curve(point, sil_local_gcd);
                
                Z_at_xy_isolator gcd_isolator =
                    construct_isolator(
                            local_gcd, point, 
                            cad.nk(vit, support),
                            (on_gcd_sil ? CGAL::EDGE : CGAL::FACE)
                    );
                
                // compute overlaps
                std::vector< std::pair< int, int > > overlaps;
                
                iv_helper.refined_overlaps(
                        support_isolator, gcd_isolator,
                        std::back_inserter(overlaps)
                );
                
                // refine until a gcd-interval is fully include in
                // support interval or has moved away
                for (int k = 0; k < static_cast< int >(overlaps.size());
                     k++) {
                    
                    if (sheets.find(overlaps[k].first) == sheets.end()) {
                        // the overlaps[k].first sheets is not isolated!
                        continue;
                    }

                    // here we are only faces with sheet levels
                    // that form an isolated point of support :-)
                    
                    while (true) {
                        
                        if (iv_helper.overlap_or_order(
                                    gcd_isolator, overlaps[k].second,
                                    support_isolator, overlaps[k].first
                            ) != CGAL::EQUAL) {
                            break;
                        }
                        // else
                        
                        if (iv_helper.is_included(
                                    gcd_isolator, overlaps[k].second,
                                    support_isolator, overlaps[k].first,
                                    true
                            )
                        ) {
                            *poi++ = Surface_point_3(vit->point(),support,
                                                     overlaps[k].first);
                            break;
                        }
                        
                        // else 
                        gcd_isolator.refine_interval(overlaps[k].second);
                    }
                }
            }
        }
    }
    
    //!@}

}; // Curve_3

CGAL_END_NAMESPACE

#endif // SoX_GAPS_CURVE_3_H
// EOF
