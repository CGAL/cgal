// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_FUNCTORS_H
#define CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_FUNCTORS_H 1

/*!\file include/CGAL/Arrangement_2l/Restricted_cad_3_functor.h
 * \brief Contains functor classes related to  
 * \link Restricted_cad_3 \endlink
 */

#include <CGAL/config.h>

#include <algorithm>
#include <vector>

// TODO remove Poly_traits_d
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Arr_enums.h>
#include <CGAL/Object.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Arr_overlay_2.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/macros.h>
#include <CGAL/Arrangement_2l/P_dcel_info.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief 
 * Computes projected cut curves for two two given surfaces
 */
template < class SurfaceZAtXyIsolatorTraits >
class Construct_nk_decomposition_2 {
    public:

    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

    typedef CGAL::Restricted_cad_3< Surface_z_at_xy_isolator_traits > 
    Restricted_cad_3;
    typedef CGAL::Restricted_cad_3_accessor< Restricted_cad_3 > 
    Accessor;
    
    /*!\brief
     * returns through \c oi projected
     * intersection curve(s) of \c s1 and \c s2
     */
    Restricted_cad_3 operator()(const Surface_3& surface) const {
        Restricted_cad_3 nkd;
        
        typename Surface_z_at_xy_isolator_traits::
            Construct_projected_surface_curves_2
            construct_projected_surface_curves
            // TODO (traits.construct_projected_surface_curves_2_object())
            ;
        
        // STEP 1: Construct Silhouette-Arr with mults
        
#if !NDEBUG
        std::cout << "Construct A(S) for surface: " 
                  << surface.id() << " ... " << std::flush;
#endif
        
        std::list< std::pair< Curve_analysis_2, int > > curves;

        construct_projected_surface_curves(
                surface, 
                std::back_inserter(curves)
        );
#if !NDEBUG
        std::cout << "#sil-curves: " << curves.size() << std::endl;
        
        for (typename std::list< std::pair< Curve_analysis_2, int > >::
                 const_iterator it = curves.begin(); it != curves.end();
             it++) {
            std::cout << "sil-curve: " 
                      << it->first.polynomial_2() << std::endl;
        }
#endif        
        Accessor acc(nkd);
        
        if (curves.empty()) {
            // empty arrangement 
            acc.init(surface);
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                nkd._set_nk_value(fit, surface, 
                                  CGAL::Nk::N, surface.f().degree());
                nkd._set_nk_value(fit, surface, 
                                  CGAL::Nk::K, 0);
            }            
            return nkd;
        }

        std::list< Restricted_cad_3 > nkds;
        
        // otherwise create for each part of silhouette curve
        for (typename std::list< std::pair< Curve_analysis_2, int > >::
                 const_iterator
                 cit = curves.begin(); cit != curves.end(); cit++) {
            // split curve and insert points/segments into the arrangement
            nkds.push_back(
                    _arrangement_of_curve(
                            surface, cit->first, 
                            CGAL::Nk::MULT, cit->second
                    )
            );
        }
        
        curves.clear();

        Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > 
            overlay; // TODO single instance
        
        nkd = overlay(nkds.begin(), nkds.end(), surface, CGAL::Nk::MULT);
        
#if !NDEBUG
        std::cout << "done." << std::endl;
#endif
      
        CGAL_assertion_code(
                unsigned int num_faces = nkd.number_of_faces();
                unsigned int num_unb_faces = nkd.number_of_unbounded_faces();
        );

        CGAL_assertion_code((
        {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                CGAL_assertion(nkd.nk(vit, surface).mult() != -1);
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                CGAL_assertion(nkd.nk(eit, surface).mult() != -1);
            }
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                CGAL_assertion(nkd.nk(fit, surface).mult() == -1);
            }
        })
        );
        
        // STEP 1-2: Set n and k for faces
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 nkd.faces_begin();
             fit != nkd.faces_end(); fit++) {
            nkd._set_nk_value(fit, surface, 
                              CGAL::Nk::N, surface.f().degree());
            nkd._set_nk_value(fit, surface, 
                              CGAL::Nk::K, 0);
        }            

        CGAL_assertion_code((
        {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                CGAL_assertion(nkd.nk(vit, surface).mult() != -1);
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                CGAL_assertion(nkd.nk(eit,surface).mult() != -1);
            }
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                CGAL_assertion(nkd.nk(fit,surface).mult() == -1);
            }
        })
        );
        
        // STEP 2: Refine edges of nkd with respect to a_i
        
        std::list< int > used_i;
        
        int max_i = surface.f().degree(); // == degree_z
        
        bool stop = false;
        
        // take care of a_i!
        for (int i = max_i; i >= 0; i--) {
            
#if !NDEBUG
            std::cout << "Refine A(S) wrt a_" 
                      << i
                      << " ... " << std::flush;
            bool skip = false;
#endif
            
            if (surface.f(i).degree() != i) {
#if !NDEBUG
                std::cout << "skipped." << std::endl;
                skip = true;
#else
                continue;
#endif
            }

            curves.clear();
            
            construct_projected_surface_curves(
                    surface, i,
                    std::back_inserter(curves)
            );

            bool whole_plane = 
                construct_projected_surface_curves(surface, i);
            
#if !NDEBUG
            std::cout << "#alpha-curves: " << curves.size() << std::endl;
            for (typename std::list< std::pair< Curve_analysis_2, int > >::
                     const_iterator it = curves.begin(); it != curves.end();
                 it++) {
                
                std::cout << "alpha-curve: " 
                          << it->first.polynomial_2() << std::endl;
            }
            if (skip) {
                continue;
            }
#endif
            bool set_i = false;
            stop = _refine_with_respect_to(
                    surface, nkd, curves.begin(), curves.end(), 
                    whole_plane, 
                    CGAL::Nk::N, set_i, i
            ); 
            
            if (set_i) {
                // store set i
                used_i.push_back(i);
            }
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
            if (stop) {
#if !NDEBUG
                std::cout << "All n set." << std::endl;
#endif
                break;
            }
        }

        if (!stop) {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                const CGAL::Nk& nk = nkd.nk(vit, surface);
                if (nk.n() == -2) {
                    nkd._set_nk_value(vit, surface, CGAL::Nk::N, -1);
                }
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                const CGAL::Nk& nk = nkd.nk(eit, surface);
                if (nk.n() == -2) {
                    nkd._set_nk_value(eit, surface, CGAL::Nk::N, -1);
                }
            }
        }

        CGAL_assertion(nkd.number_of_faces() == num_faces);
        CGAL_assertion(nkd.number_of_unbounded_faces() == num_unb_faces);

        CGAL_assertion_code((
        {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                const CGAL::Nk& nk = nkd.nk(vit, surface);
                CGAL_assertion(nk.mult() != -1);
                CGAL_assertion(nk.n() != -2);
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                const CGAL::Nk& nk = nkd.nk(eit, surface);
                CGAL_assertion(nk.mult() != -1);
                CGAL_assertion(nk.n() != -2);
            }
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                const CGAL::Nk& nk = nkd.nk(fit, surface);
                CGAL_assertion(nk.mult() == -1);
                CGAL_assertion(nk.n() != -2);
            }
        })
        );
        
        // STEP 3: Refine Sil-Arr with respect to stha_i^n
        // iterate over set a_i
        stop = false;
        
        
        for (typename std::list< int >::const_iterator nit = 
                 used_i.begin(); nit != used_i.end(); nit++) {


#if !NDEBUG            
            std::cout << "*nit: " << *nit << std::endl;
            std::cout << "f(n): " << surface.f(*nit).degree() << std::endl;
#endif      
      
            // start with 0, as this can happen for asymptotic faces
            for (int k = 0; k < surface.f(*nit).degree(); k++) {

#if !NDEBUG
                std::cout << "Refine A(S) wrt stha_" 
                          << k << "^F(" << *nit << ")"
                          << " ... " << std::flush;
                bool skip = false;
#endif
                
                // quick exit for silhouette itself
                if (k == 0 && 
                    *nit == typename CGAL::Polynomial_traits_d< Polynomial_3 >::
                    Total_degree()(surface.f())) {
#if !NDEBUG
                    std::cout << "skipped." << std::endl;
                    skip = true;
#else
                    continue;
#endif
                }
                
                curves.clear();
                
                construct_projected_surface_curves(
                        surface, k, *nit,
                        std::back_inserter(curves)
                );

                bool whole_plane = 
                    construct_projected_surface_curves(surface, k, *nit);
                
#if !NDEBUG
                std::cout << "#sigma-curves_" << k 
                          << ": " << curves.size() << std::endl;
                for (typename std::list< std::pair< Curve_analysis_2, int > >::
                         const_iterator it = curves.begin(); it != curves.end();
                     it++) {
                    
                    std::cout << "sigma-curve: " 
                              << it->first.polynomial_2() << std::endl;
                }
                if (skip) {
                    continue;
                }
#endif
                
                bool set_k = false;
                stop = _refine_with_respect_to(
                        surface, nkd, curves.begin(), curves.end(), 
                        whole_plane,
                        CGAL::Nk::K, set_k, k, *nit
                );
                
#if !NDEBUG
                std::cout << "done." << std::endl;
#endif
                
                if (stop) {
#if !NDEBUG
                    std::cout << "All k set." << std::endl;
#endif
                    break;
                }
            }
        }

        CGAL_assertion(nkd.number_of_faces() == num_faces);
        CGAL_assertion(nkd.number_of_unbounded_faces() == num_unb_faces);
        
        if (!stop) {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                const CGAL::Nk& nk = nkd.nk(vit, surface);
                if (nk.k() == -2) {
                    if (nk.n() <= 0) {
                        CGAL_assertion(nk.n() > -2);
                        nkd._set_nk_value(vit, surface, CGAL::Nk::K, nk.n());
                    } 
                }
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                const CGAL::Nk& nk = nkd.nk(eit, surface);
                if (nk.k() == -2) {
                    if (nk.n() <= 0) {
                        CGAL_assertion(nk.n() > -2);
                        nkd._set_nk_value(eit, surface, CGAL::Nk::K, nk.n());
                    }
                }
            }
        }

        // STEP 4: Delete unneccessary vertices
        
        for (typename Restricted_cad_3::Vertex_const_iterator vt, vit =
                 nkd.vertices_begin();
             vit != nkd.vertices_end();) {

            bool remove = false;
            
            if (vit->degree() == 2) {
                const CGAL::Nk& nk = nkd.nk(vit, surface);
                //std::cout << "vnk: " << nk << std::endl;
                
                
                typename 
                    Restricted_cad_3::Halfedge_around_vertex_const_circulator 
                    e1 = vit->incident_halfedges();
                
                const CGAL::Nk& nk1 = nkd.nk(e1, surface);
                
                if (nk1.same_nk(nk)) {
                    
                    typename 
                        Restricted_cad_3::
                        Halfedge_around_vertex_const_circulator e2 = 
                        e1->next();
                    CGAL_assertion(e1 != e2);
                    
                    const CGAL::Nk& nk2 = nkd.nk(e2, surface);
                    
                    if (nk2.same_nk(nk)) {
                        if (e2->curve().are_mergeable(e1->curve())) {
                            remove = true;
                        }
                    }
                }
            }
            
            if (remove) {
                vt = vit; 
                vit++;
                acc.remove_vertex(acc.non_const_handle(vt));
            } else {
                vit++;
            }
        }

        // FINISHED!
        CGAL_assertion_code((
        {
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                const CGAL::Nk& nk = nkd.nk(fit, surface);
                //std::cout << "fnk: " << nk << std::endl;
                CGAL_assertion(nk.mult() == -1);
                CGAL_assertion(nk.n() ==  surface.f().degree());
                CGAL_assertion(nk.k() == 0);
            }
            
            for (typename Restricted_cad_3::Edge_const_iterator 
                     eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                const CGAL::Nk& nk = nkd.nk(eit, surface);
                //std::cout << "hnk: " << nk << std::endl;
                CGAL_assertion(nk.mult() > 0);
                CGAL_assertion(nk.n() != -2);
                CGAL_assertion(nk.k() != -2);
                CGAL_assertion(
                        (nk.n() == -1 && nk.k() == -1) ||
                        (nk.n() == 0 && nk.k() == 0) || 
                        (nk.n() > 0 && nk.k() >= 0 && nk.k() < nk.n())
                );
            }
            
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                const CGAL::Nk& nk = nkd.nk(vit, surface);
                //std::cout << "vnk: " << nk << std::endl;
                CGAL_assertion(nk.mult() == 0);
                CGAL_assertion(nk.n() != -2);
                CGAL_assertion(nk.k() != -2);
                CGAL_assertion(
                        (nk.n() == 0 && nk.k() == 0) || 
                        (nk.n() > 0 && nk.k() >= 0 && nk.k() < nk.n()) ||
                        (nk.n() == -1 && nk.k() == -1)
                );
            }
        })
        );
        
        return nkd;
    }
    
private:
    
    /*!\brief
     * creates segments and points for a given curve and inserts them
     * into an arrangement
     */
    Restricted_cad_3 
    _arrangement_of_curve(const Surface_3& surface, 
                          const Curve_analysis_2& curve, 
                          CGAL::Nk::Value_type type, 
                          int value) const {
        
        Restricted_cad_3 nkd = Accessor::construct_for_curve(curve);

        Accessor acc(nkd);
        
        // construct data objects
        acc.init(surface);
        
        // set nk-data
        if (type == CGAL::Nk::MULT) {
            // set mult
            // TODO store multiplicity of singular point??!?!
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                nkd._set_nk_value(vit, surface, CGAL::Nk::MULT, 0);
                // set it to 0 to indicate
                // that it belongs to sil
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                nkd._set_nk_value(eit, surface, CGAL::Nk::MULT, value);
            }
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                CGAL_assertion(nkd.nk(fit, surface).mult() == -1);
            }            
        }

        CGAL_assertion_code((
        {
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     nkd.vertices_begin();
                 vit != nkd.vertices_end(); vit++) {
                CGAL_postcondition(vit->data());
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                 nkd.edges_begin();
                 eit != nkd.edges_end(); eit++) {
                CGAL_postcondition(eit->data());
                CGAL_postcondition(eit->twin()->data());
            }
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     nkd.faces_begin();
                 fit != nkd.faces_end(); fit++) {
                CGAL_postcondition(fit->data());
            }
        })
        );

        return nkd;
    }
    
    template < class InputIterator >
    bool _refine_with_respect_to(const Surface_3& surface, 
                                 Restricted_cad_3& nkd,
                                 InputIterator begin, InputIterator end,
                                 bool whole_plane,
                                 CGAL::Nk::Value_type type,
                                 bool& set_value,
                                 int value1,
                                 int value2 = -1) const {
        
        Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > 
            overlay; // TODO single instance
        
        Restricted_cad_3 nkrefine;

        set_value = false;

        if (begin == end) {

            Accessor acc(nkrefine);
            acc.init(surface);
            
        } else {
            
            CGAL_precondition(!whole_plane);
            
            std::list< Restricted_cad_3 > nkds;
            
            for (InputIterator it = begin; it != end; it++) {
                // split curve and insert points/segments 
                // into the arrangement
                nkds.push_back(
                        _arrangement_of_curve(
                                surface, it->first, type, -1
                        )
                );
            }
            
            nkrefine = overlay(nkds.begin(), nkds.end(), surface, type);
        }

        // overlay nkd with nkrefine
        Restricted_cad_3 tmp = overlay(nkd, nkrefine, surface, type);

        Accessor acc(tmp);

        if (begin != end) {

            // clean wrong halfedges within a face of nkd
            for (typename Restricted_cad_3::Edge_const_iterator 
                     ti, eit =
                     tmp.edges_begin();
                 eit != tmp.edges_end();) {
                if (tmp.nk(eit, surface).feature1() == CGAL::FACE) {
                    ti = eit; 
                    eit++;
                    acc.remove_edge(acc.non_const_handle(ti));
                } else {
                    eit++;
                }
            }

            // clean wrong vertices within a face of nkd
            for (typename Restricted_cad_3::Vertex_const_iterator 
                     vt, vit =
                     tmp.vertices_begin();
                 vit != tmp.vertices_end();) {
                if (tmp.nk(vit, surface).feature1() == CGAL::FACE) {
                    vt = vit; 
                    vit++;
                    acc.remove_vertex(acc.non_const_handle(vt));
                } else {
                    vit++;
                }
            }

            // clean vertices with already stored data
            for (typename Restricted_cad_3::Vertex_const_iterator 
                     vt, vit =
                     tmp.vertices_begin();
                 vit != tmp.vertices_end(); ) {
                const CGAL::Nk& nk = tmp.nk(vit, surface);
                bool remove = false;
                if (nk.feature1() == CGAL::EDGE) {
                    if (type == CGAL::Nk::N) {
                        if (nk.n() == -2) {
                            remove = true;
                        }
                    } else if (type == CGAL::Nk::K) {
                        if (nk.k() != -2 || 
                            nk.n() != value2) {
                            remove = true;
                        }
                    }
                }
                if (remove) {
                    vt = vit;
                    vit++;
                    acc.remove_vertex(acc.non_const_handle(vt));
                } else {
                    vit++;
                }
            }
        }
        
        bool stop = true;
        // set values for correct halfedges ...
        for (typename Restricted_cad_3::Edge_const_iterator 
                 eit =
                 tmp.edges_begin();
             eit != tmp.edges_end(); eit++) { 
            CGAL::Nk nk = tmp.nk(eit, surface);
            //std::cout << "H1: " << nk << std::endl;
            if (nk.feature2() == CGAL::FACE) {
                CGAL_assertion(nk.mult() != -1);
                if (type == CGAL::Nk::N) {
                    if (!whole_plane && 
                        nk.n() == -2) {
                        //std::cout << "HN: " << value1 << std::endl;
                        nkd._set_nk_value(eit, surface, CGAL::Nk::N, value1);
                        set_value = true;
                    } 
                } else if (type == CGAL::Nk::K) {
                    if (nk.n() == value2 &&
                        !whole_plane && 
                        nk.k() == -2) {
                        //std::cout << "HK: " << value1 << std::endl;
                        nkd._set_nk_value(eit, surface, CGAL::Nk::K, value1);
                        set_value = true;
                    } 
                }
            }
            nk = tmp.nk(eit, surface);
            if (stop) {
                if (type == CGAL::Nk::N) {
                    if (nk.n() == -2) {
                        stop = false;
                    }
                } else if (type == CGAL::Nk::K) {
                    if (nk.k() == -2) {
                        stop = false;
                    }
                }
            }
        }

        // ... and vertices
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 tmp.vertices_begin();
             vit != tmp.vertices_end(); vit++) { 
            CGAL::Nk nk = tmp.nk(vit, surface);
            //std::cout << "V1: " << nk << std::endl;
            if (nk.feature2() == CGAL::FACE) {
                CGAL_assertion(nk.mult() != -1);
                if (type == CGAL::Nk::N) {
                    if (!whole_plane && 
                        nk.n() == -2) {
                        //std::cout << "VN: " << value1 << std::endl;
                        nkd._set_nk_value(vit, surface, CGAL::Nk::N, value1);
                        set_value = true;
                    }
                } else if (type == CGAL::Nk::K) {
                    if (!whole_plane && 
                        nk.n() == value2 &&
                        nk.k() == -2) {
                        //std::cout << "VK: " << value1 << std::endl;
                        nkd._set_nk_value(vit, surface, CGAL::Nk::K, value1);
                        set_value = true;
                    }
                }
            }
            nk = tmp.nk(vit, surface);
            if (stop) {
                if (type == CGAL::Nk::N) {
                    if (nk.n() == -2) {
                        stop = false;
                    }
                } else if (type == CGAL::Nk::K) {
                    if (nk.k() == -2) {
                        stop = false;
                    }
                }
            }
        }
        
        nkd = tmp;
        
        return stop;
    }
};



} // namespace CGALi

/*!\brief
 * Functor class that provides creation of restriced_cads.
 */
template < class SurfaceZAtXyIsolatorTraits >
class Create_restricted_cad_3 { 
public:    
    
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of instantiated class
    typedef Create_restricted_cad_3< Surface_z_at_xy_isolator_traits > Self;
    
#if DOXYGEN_RUNNING
    //! type of projected curve
    typedef typename Surface_z_at_xy_isolator_traits::Curve_analysis_2 
    Curve_analysis_2;
    
    //! type of projected point
    typedef typename Surface_z_at_xy_isolator_traits::Point_2 Point_2;

    //! type of projected segment
    typedef typename Surface_z_at_xy_isolator_traits::X_monotone_curve_2
    X_monotone_curve_2;
#endif
    
    //! type of Restricted_cad_3
    typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits >
    Restricted_cad_3;

private:
    //! type of data
    typedef typename Restricted_cad_3::Data Data;
    

public:
    /*!\brief
     * creates cad for silhouette-curve of \c surface
     */
    Restricted_cad_3 operator()(const Surface_3& surface) const {
        return Restricted_cad_3::cad_cache()(surface);
    }
    
private:
    //! used to indicate special constructors
    class Non_cached_construction_tag {};

    /*!\brief
     * creates cad for silhouette-curve of \c surface (called in cache only)
     */
    Restricted_cad_3 operator()(const Surface_3& surface, 
                                Non_cached_construction_tag t) const {
        
        typedef CGAL::CGALi::Construct_nk_decomposition_2< 
            Surface_z_at_xy_isolator_traits >
            Construct_nk_decomposition_2;
        
        Construct_nk_decomposition_2 construct_nk_decomposition
            // TODO (traits.construct_nk_decomposition_2_object())
            ;
        
        Restricted_cad_3 nkd = construct_nk_decomposition(surface);
        
        nkd._finalize(surface);
        
        // and return it
        return nkd;
    }

public:    
    /*!\brief
     * creates cad for cut-curve of \c surface1 and \c surface2
     */
    Restricted_cad_3 operator()(const Surface_3& surface1, 
                                const Surface_3& surface2) const {
        return Restricted_cad_3::cad_cache()(surface1, surface2);
    }
    
private:
    /*!\brief
     * creates cad for cut-curve of \c surface1 and \c surface2
     * (called in cache only)
     */
    Restricted_cad_3 operator()(const Surface_3& surface1, 
                                const Surface_3& surface2,
                                Non_cached_construction_tag t) const {
        
        typename 
            Surface_z_at_xy_isolator_traits::
            Construct_projected_cuts_2
            construct_projected_cuts
            // TODO (traits.construct_projected_cuts_2_object())
            ;
        
        // store components of curve
        std::list< std::pair< Curve_analysis_2, int > > curves;
        
        construct_projected_cuts(
                surface1, surface2,
                std::back_inserter(curves)
        );
        

#if !NDEBUG
        std::cout << "#cuts-curves: " << curves.size() << std::endl;
        
        for (typename std::list< std::pair< Curve_analysis_2, int > >::
                 const_iterator it = curves.begin(); it != curves.end();
             it++) {
            std::cout << "cut-curve: " 
                      << it->first.polynomial_2() << std::endl;
        }
#endif     

        if (curves.empty()) {
            // empty arrangement 
            Restricted_cad_3 cad;
            cad._renew();

            // for all faces
            for (typename Restricted_cad_3::Face_const_iterator fit =
                     cad.faces_begin();
                 fit != cad.faces_end(); fit++) {
                if (!fit->data()) {
                    CGAL::Object obj = CGAL::make_object(fit);
                    cad.non_const_handle(fit)->data() = Data(CGAL::FACE, obj);
                    fit->data()->_set_rs_id(cad.id());
                    fit->data()->_add_surface(surface1);
                    fit->data()->_add_surface(surface2);
                    fit->data()->_set_dcel(surface1, surface2, CGAL::FACE, obj);
                }
            }
            // and return it
            return cad;
        }
        
        // for each part create RS_CAD_3
        std::vector< Restricted_cad_3 > 
            cads(std::distance(curves.begin(), curves.end()), 
                 Restricted_cad_3()
            );
        int i = 0;
        
        typedef CGAL::Restricted_cad_3_accessor< Restricted_cad_3 > 
            Accessor;

        for (typename std::list< std::pair< Curve_analysis_2, int > >::
                 const_iterator
                 cit = curves.begin(); cit != curves.end(); cit++, i++) {
            
            // make rep unique
            cads[i]._renew();
            
            // split curve and insert points/segments into the arrangement
            cads[i] = Accessor::construct_for_curve(cit->first);
            
            // assign values
            _set_cut(surface1, surface2, cit->second, cads[i]);
        }
        
        // merge the RS_CAD_3s
        Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > 
            overlay; // TODO single instance

        CGAL_assertion(cads.begin() != cads.end());
        Restricted_cad_3 tmp = overlay(cads.begin(), cads.end(), true);
        
        // it remains to set handles correctly
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 tmp.vertices_begin();
             vit != tmp.vertices_end(); vit++) {
            vit->data()->_set_rs_id(tmp.id());
            CGAL::Object obj = CGAL::make_object(vit);
            vit->data()->_set_dcel(CGAL::VERTEX, obj);
            vit->data()->_set_dcel(surface1, surface2, CGAL::VERTEX, obj);
        }
        for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                 tmp.halfedges_begin();
             hit != tmp.halfedges_end(); hit++) {
            bool correct_dir = (hit->direction() == CGAL::ARR_LEFT_TO_RIGHT);
            if (correct_dir) {
                hit->data()->_set_rs_id(tmp.id());
                CGAL::Object obj = CGAL::make_object(hit);
                hit->data()->_set_dcel(CGAL::EDGE, obj);
                hit->data()->_set_dcel(surface1, surface2, CGAL::EDGE, obj);
            }
        }
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 tmp.faces_begin();
             fit != tmp.faces_end(); fit++) {
            fit->data()->_set_rs_id(tmp.id());
            CGAL::Object obj = CGAL::make_object(fit);
            fit->data()->_set_dcel(CGAL::FACE, obj);
            fit->data()->_set_dcel(surface1, surface2, CGAL::FACE, obj);
#if 0
            // TASK benchmark: check whether this saves time
            // it may be even problematic, as we have to isolate
            // over non-face-features of the "virtual" silhouette-curve
            cad.slice(fit);
#endif
        }
        
        // and return it
        return tmp;
    }

public:
    /*!\brief
     * Compute cad for a range of surfaces defined by \c [begin,end)
     */
    // TASK move to Surface_cad??
    template < class InputIterator >
    Restricted_cad_3 operator()(InputIterator begin, InputIterator end) const {
        
        CGAL_precondition((boost::is_same< typename InputIterator::value_type,
                     Surface_3 >::value
                    ));
        
        std::vector< Restricted_cad_3 > boundaries;
        boundaries.reserve(std::distance(begin,end));
        for (InputIterator it = begin; it != end; it++) {
            boundaries.push_back(this->operator()(*it));
        }
        
        std::vector< Restricted_cad_3 > cuts;
        int d = std::distance(begin,end);
        cuts.reserve(d*(d-1));
        for (InputIterator oit = begin; oit != end; oit++) {
            for (InputIterator iit = oit; iit != end; iit++) {
                if (oit != iit) {
                    cuts.push_back(this->operator()(*oit,*iit));
                }
            }
        }
        
        Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > 
            overlay; // TODO single instance
        
        // TASK stream
#if NDEBUG
        std::cout << "Boundaries ... " << std::flush;
#endif
        Restricted_cad_3 bdry_overlay = overlay(boundaries.begin(),
                                                boundaries.end());
#if NDEBUG
        std::cout << "done." << std::endl;
#endif
        Restricted_cad_3 out;
        out._renew();
        
        if (static_cast< int >(boundaries.size()) > 1) {
#if NDEBUG
            std::cout << "Cuts ... " << std::flush;
#endif
            Restricted_cad_3 itxs_overlay = overlay(cuts.begin(),
                                                    cuts.end());
#if NDEBUG
            std::cout << "done." << std::endl;
#endif
            
#if NDEBUG
            std::cout << "Final ... " << std::flush;
#endif
            out = overlay(bdry_overlay, itxs_overlay);
#if NDEBUG
            std::cout << "done." << std::endl;
#endif
        } else {
            out = bdry_overlay;
#if NDEBUG
            std::cout << "Final consists of a single silhouette ... done." 
                      << std::endl;
#endif
        }
        
        return out;
    }

private:
    
    /*!\brief
     * sets features of dcel with silhouette of \c surface
     */
    void _set_silhouette(const Surface_3& surface, 
                       int multiplicity, bool vertical,
                       Restricted_cad_3& rscad) const {

        // for all vertices
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 rscad.vertices_begin();
             vit != rscad.vertices_end(); vit++) {
            CGAL_assertion(!vit->is_at_infinity());
            if (!vit->data()) {
                rscad.non_const_handle(vit)->data() = 
                    Data(CGAL::VERTEX, CGAL::make_object(vit));
                vit->data()->_add_surface(surface);
            }
            // store multiplicity of singular point
#if 0
            if (vit->point().curve().may_be_singular()) {
                int vmult = vit->point().curve().event_info_at_x(
                        vit->point().x().number().finite()
                ).multiplicity(
                        vit->point().arcno()
                );
                vit->data()->_add_silhouette(surface, multiplicity * vmult, 
                                             false, vertical);
            } else {
                vit->data()->_add_silhouette(surface, multiplicity, 
                                             false, vertical);
            }
#else
            vit->data()->_add_silhouette(surface, -1, 
                                         false, vertical);
#endif
        }

        // for all edges
        for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                 rscad.halfedges_begin();
             hit != rscad.halfedges_end(); hit++) {
            if (!hit->data()) {
                rscad.non_const_handle(hit)->data() = 
                    Data(CGAL::EDGE, CGAL::make_object(hit));
                hit->data()->_add_silhouette(surface, multiplicity, 
                                             false, vertical);
                hit->data()->_add_surface(surface);
            }
            if (!hit->twin()->data()) {
                rscad.non_const_handle(hit)->twin()->data() = *hit->data();
            }
        }

        // for all faces
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 rscad.faces_begin();
             fit != rscad.faces_end(); fit++) {
            if (!fit->data()) {
                rscad.non_const_handle(fit)->data() = 
                    Data(CGAL::FACE, CGAL::make_object(fit));
                fit->data()->_add_surface(surface);
            }
        }
    }
    
    /*!\brief
     * sets features of dcel with cut of \c surface1 and \c surface2
     */
    void _set_cut(const Surface_3& surface1, 
                  const Surface_3& surface2, 
                  int multiplicity,
                  Restricted_cad_3& rscad) const {
        
        // for all vertices
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 rscad.vertices_begin();
             vit != rscad.vertices_end(); vit++) {
            CGAL_assertion(!vit->is_at_infinity());
            if (!vit->data()) {
                rscad.non_const_handle(vit)->data() = 
                    Data(CGAL::VERTEX, CGAL::make_object(vit));
                vit->data()->_add_surface(surface1);
                vit->data()->_add_surface(surface2);
            }
            // store multiplicity of singular point
#if 0
            if (vit->point().curve().may_be_singular()) {
                int vmult = vit->point().curve().event_info_at_x(
                        vit->point().x().number().finite()
                ).multiplicity(
                        vit->point().arcno()
                );
                vit->data()->_add_cut(
                        surface1, surface2, multiplicity * vmult, false
                );
            } else {
                vit->data()->_add_cut(
                        surface1, surface2, multiplicity, false
                );
            }
#else
            vit->data()->_add_cut(
                    surface1, surface2, -1, false
            );
#endif
        }

        // for all edges
        for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                 rscad.halfedges_begin();
             hit != rscad.halfedges_end(); hit++) {
            if (!hit->data()) {
                rscad.non_const_handle(hit)->data() = 
                    Data(CGAL::EDGE, CGAL::make_object(hit));
                hit->data()->_add_cut(surface1, surface2, 
                                      multiplicity, false);
                hit->data()->_add_surface(surface1);
                hit->data()->_add_surface(surface2);
            }
            if (!hit->twin()->data()) {
                rscad.non_const_handle(hit)->twin()->data() = *hit->data();
            }
        }
        
        // for all faces
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 rscad.faces_begin();
             fit != rscad.faces_end(); fit++) {
            if (!fit->data()) {
                rscad.non_const_handle(fit)->data() = 
                    Data(CGAL::FACE, CGAL::make_object(fit));
                fit->data()->_add_surface(surface1);
                fit->data()->_add_surface(surface2);
            }
        }
    }

private:

    // friend
    // for non-cached construction
    friend class 
    CGAL::CGALi::Restricted_cad_3_cache< Surface_z_at_xy_isolator_traits >;
    
}; // Create


// TODO move Overlay_restricted_cad_3 to namespace CGALi?!?!
/*!\brief
 * Functor class that provides overlaying of restriced_cads.
 */
template < class SurfaceZAtXyIsolatorTraits >
class Overlay_restricted_cad_3 {
public:    

    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of instantiated class
    typedef Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > Self;
    
    //! type of Restricted_cad_3
    typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits >
    Restricted_cad_3;
    
    /*!\brief
     * merges two restricted cads
     */
    Restricted_cad_3 operator()(const Restricted_cad_3& rscad1, 
                                const Restricted_cad_3& rscad2) const {
        return this->operator()(rscad1, rscad2, false);
    }
    
private:
    /*!\brief
     * merges two restricted cads
     */
    Restricted_cad_3 operator()(const Restricted_cad_3& rscad1, 
                                const Restricted_cad_3& rscad2,
                                bool factors_of_same_curve) const {
        
        Restricted_cad_3 out; 
        out._renew();
        
        CGAL_precondition_code(
                for (typename Restricted_cad_3::Vertex_const_iterator vit =
                         rscad1.vertices_begin();
                     vit != rscad1.vertices_end(); vit++) {
                    CGAL_precondition(vit->data());
                }
                for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         rscad1.halfedges_begin();
                     hit != rscad1.halfedges_end(); hit++) {
                    CGAL_precondition(hit->data());
                }
                for (typename Restricted_cad_3::Face_const_iterator fit =
                         rscad1.faces_begin();
                     fit != rscad1.faces_end(); fit++) {
                    CGAL_precondition(fit->data());
                }
                for (typename Restricted_cad_3::Vertex_const_iterator vit =
                         rscad2.vertices_begin();
                     vit != rscad2.vertices_end(); vit++) {
                    CGAL_precondition(vit->data());
                }
                for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         rscad2.halfedges_begin();
                     hit != rscad2.halfedges_end(); hit++) {
                    CGAL_precondition(hit->data());
                }
                for (typename Restricted_cad_3::Face_const_iterator fit =
                         rscad2.faces_begin();
                     fit != rscad2.faces_end(); fit++) {
                    CGAL_precondition(fit->data());
                }
        );
        
        // overlay computation (unique place to do it)
        out._overlay(rscad1, rscad2, factors_of_same_curve);

        // TASK benchmark:
        // compute slices already here vs. delay their computation
        
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 out.vertices_begin();
             vit != out.vertices_end(); vit++) {
            vit->data()->_set_rs_id(out.id());
            CGAL::Object obj = CGAL::make_object(vit);
            vit->data()->_set_dcel(CGAL::VERTEX, obj);
            CGAL_postcondition(vit->data());
        }
        for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         out.halfedges_begin();
             hit != out.halfedges_end(); hit++) {
            bool correct_dir = (hit->direction() == CGAL::ARR_LEFT_TO_RIGHT);
            if (correct_dir) {
                hit->data()->_set_rs_id(out.id());
                CGAL::Object obj = CGAL::make_object(hit);
                hit->data()->_set_dcel(CGAL::EDGE, obj);
            }
            CGAL_postcondition(hit->data());
        }
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 out.faces_begin();
             fit != out.faces_end(); fit++) {
            fit->data()->_set_rs_id(out.id());
            CGAL::Object obj = CGAL::make_object(fit);
            fit->data()->_set_dcel(CGAL::FACE, obj);
            CGAL_postcondition(fit->data());
        }
        
        // cache id
        return Restricted_cad_3::cad_cache()(out);
    }


    Restricted_cad_3 operator()(const Restricted_cad_3& rscad1, 
                                const Restricted_cad_3& rscad2,
                                const Surface_3& surface,
                                CGAL::Nk::Value_type type) const {
        
        Restricted_cad_3 out; 
        out._renew();
        
        CGAL_precondition_code(
                for (typename Restricted_cad_3::Vertex_const_iterator vit =
                         rscad1.vertices_begin();
                     vit != rscad1.vertices_end(); vit++) {
                    CGAL_precondition(vit->data());
                }
                for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         rscad1.halfedges_begin();
                     hit != rscad1.halfedges_end(); hit++) {
                    CGAL_precondition(hit->data());
                }
                for (typename Restricted_cad_3::Face_const_iterator fit =
                         rscad1.faces_begin();
                     fit != rscad1.faces_end(); fit++) {
                    CGAL_precondition(fit->data());
                }
                for (typename Restricted_cad_3::Vertex_const_iterator vit =
                         rscad2.vertices_begin();
                     vit != rscad2.vertices_end(); vit++) {
                    CGAL_precondition(vit->data());
                }
                for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         rscad2.halfedges_begin();
                     hit != rscad2.halfedges_end(); hit++) {
                    CGAL_precondition(hit->data());
                }
                for (typename Restricted_cad_3::Face_const_iterator fit =
                         rscad2.faces_begin();
                     fit != rscad2.faces_end(); fit++) {
                    CGAL_precondition(fit->data());
                }
        );
        
        // overlay computation (unique place to do it)
        out._overlay(rscad1, rscad2, surface, type);

        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 out.vertices_begin();
             vit != out.vertices_end(); vit++) {
            CGAL_postcondition(vit->data());
        }
        for (typename Restricted_cad_3::Halfedge_const_iterator hit =
                         out.halfedges_begin();
             hit != out.halfedges_end(); hit++) {
            CGAL_postcondition(hit->data());
        }
        for (typename Restricted_cad_3::Face_const_iterator fit =
                 out.faces_begin();
             fit != out.faces_end(); fit++) {
            CGAL_postcondition(fit->data());
        }
        
        // cache id
        return Restricted_cad_3::cad_cache()(out);
    }
    
public: 
    /*!\brief
     * Computes overlay of a range of restricted cads
     */
    template < class InputIterator >
    Restricted_cad_3 operator()(InputIterator begin, InputIterator end) const {
        return this->operator()(begin, end, false);
    }
    
private:
    /*!\brief
     * Computes overlay of a range of restricted cads
     *
     * \precond range contains at least one element
     */
    template < class InputIterator >
    Restricted_cad_3 operator()(InputIterator begin, InputIterator end,
                                bool factors_of_same_curve) const {
        
        CGAL_precondition(begin != end);

        CGAL_precondition((boost::is_same< typename InputIterator::value_type,
                     Restricted_cad_3 >::value
                    ));
        
        std::vector< Restricted_cad_3 > cads;
        std::copy(begin, end, std::back_inserter(cads));
        
        int n = static_cast< int >(cads.size());
        CGAL_assertion(n > 0);
        while (n > 1) {
            // TASK better random generator
            std::random_shuffle(cads.begin(),cads.end());
            
            // TASK stream
#if NDEBUG
            std::cout << "Size " << n << std::flush;
#endif
            std::vector< Restricted_cad_3 > tmp;
            int i = 0;
            for (; i + 1 < n; i += 2) {
                tmp.push_back(this->operator()(cads[i],cads[i+1]));
            }
            // odd number of items
            if (i + 1 == n) {
                tmp.push_back(cads[i]);
            }
            cads.clear();
            cads = tmp;
            n = static_cast< int >(cads.size());
            if (n > 1) {
#if NDEBUG
                std::cout << "," << std::flush;
#endif
            }
#if NDEBUG
            std::cout << " " << std::flush;
#endif
        }
        CGAL_assertion(n == 1);
        return cads[0];
    }

    
    template < class InputIterator >
    Restricted_cad_3 operator()(InputIterator begin, InputIterator end,
                                const Surface_3& surface, 
                                CGAL::Nk::Value_type type) const {
        
        CGAL_precondition(begin != end);

        CGAL_precondition((boost::is_same< typename InputIterator::value_type,
                     Restricted_cad_3 >::value
                    ));
        
        std::vector< Restricted_cad_3 > cads;
        std::copy(begin, end, std::back_inserter(cads));
        
        int n = static_cast< int >(cads.size());
        CGAL_assertion(n > 0);
        while (n > 1) {
            // TASK better random generator
            std::random_shuffle(cads.begin(),cads.end());
            
            // TASK stream
#if NDEBUG
            std::cout << "Size " << n << std::flush;
#endif
            std::vector< Restricted_cad_3 > tmp;
            int i = 0;
            for (; i + 1 < n; i += 2) {
                tmp.push_back(this->operator()(cads[i], cads[i+1],
                                               surface, type));
            }
            // odd number of items
            if (i + 1 == n) {
                tmp.push_back(cads[i]);
            }
            cads.clear();
            cads = tmp;
            n = static_cast< int >(cads.size());
            if (n > 1) {
#if NDEBUG
                std::cout << "," << std::flush;
#endif
            }
#if NDEBUG
            std::cout << " " << std::flush;
#endif
        }
        CGAL_assertion(n == 1);
        return cads[0];
    }


private:
    
    // friends
    // for overlay with set factors_of_same_curve
    friend class Create_restricted_cad_3< Surface_z_at_xy_isolator_traits >;
    
    friend class 
    CGAL::CGALi::Construct_nk_decomposition_2< Surface_z_at_xy_isolator_traits >;
    
}; // Overlay

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_FUNCTORS_H
// EOF
