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

#ifndef CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_H
#define CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_H 1

/*!\file include/CGAL/Arrangement/Restricted_cad_3.h
 * \brief Contains class template Restricted_cad_3
 */

#include <CGAL/config.h>

#include <algorithm>
#include <vector>

#include <boost/optional/optional.hpp>
#include <boost/none.hpp>

#include <CGAL/Handle_with_policy.h>
#include <CGAL/Object.h>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_overlay_2.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/macros.h>
#include <CGAL/Arrangement_2l/P_dcel_info.h>
#include <CGAL/Arrangement_2l/Adjacencies_3.h>
#include <CGAL/Arrangement_2l/Arr_p_dcel_info_overlay_traits.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_accessor.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_functors.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief
 * Functor class that provides creation and overlaying of
 * restriced_cads.
 */
template < class SurfaceZAtXyIsolatorTraits >
class Restricted_cad_3_cache {
public:    

    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! type of instantiated class
    typedef Restricted_cad_3_cache< Surface_z_at_xy_isolator_traits > Self;

    //! type of Restricted_cad_3
    typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits >
    Restricted_cad_3;
    
    //! type of id
    typedef typename Restricted_cad_3::Id_type Id_type;

    //! type of Surface_3
    typedef typename Surface_z_at_xy_isolator_traits::Surface_3 Surface_3;

private:
    //! type of creator
    typedef CGAL::Create_restricted_cad_3< Surface_z_at_xy_isolator_traits >
    Creator;
    
    //! less for cad
    typedef CGAL::Handle_id_less_than< Restricted_cad_3 > Less;

    //! type of map for restricted cad 
    typedef std::map< Id_type, Restricted_cad_3, std::less< Id_type > > 
    Restricted_cad_map;
    
    //! type of Less for surfaces
    typedef typename Surface_3::Surface_less_than Surface_less_than;
  
    //! type of map
    typedef std::map< Surface_3, Id_type, Surface_less_than >
    Surface_cache_3;
    
    //! type of Surface_pair
    typedef std::pair< Surface_3, Surface_3 > Surface_pair_3;
    
    //! type of Less of surface pair
    typedef CGAL::Pair_lexicographical_less_than< Surface_3, Surface_3, 
                                Surface_less_than, Surface_less_than > 
    Surface_pair_less;
    
    //! type of pair map
    typedef std::map< Surface_pair_3, Id_type, Surface_pair_less >
    Surface_pair_cache_3;
    
    //! pair canonicalizer
    struct Canonicalize_pair {
        
        //! first in pair in always less than second
        Surface_pair_3 operator()(Surface_pair_3 pair) {
            Surface_less_than less;
            if (!less(pair.first, pair.second)) {
                std::swap(pair.first, pair.second);
            }
            return pair;
        }
    };

public:
    //!\name Caching for cads with ids
    //!@{
    
    //! returns cad of given \c id
    Restricted_cad_3 operator()(Id_type id) {
        typename Restricted_cad_map::iterator it = _m_cads.find(id);
        CGAL_assertion(it != _m_cads.end());
        return it->second;
    }

    //! cached \c cad
    Restricted_cad_3 operator()(Restricted_cad_3 cad) {
        typename Restricted_cad_map::iterator it = _m_cads.find(cad.id());
        CGAL_assertion(it == _m_cads.end());
        _m_cads.insert(it, std::make_pair(cad.id(), cad));
        return cad;
    }

    //!@}
    
    //!\name Caching for cads representing the silhouettes of single surface
    //!@{

public:
    //! returns \c cad of silhouette-curve for \c surface
    Restricted_cad_3 operator()(const Surface_3& surface) {
        typename Surface_cache_3::iterator it = _m_map.find(surface);
        if (it == _m_map.end()) {
            Creator creator; // TODO single instance
            typename Creator::Non_cached_construction_tag t;
            Restricted_cad_3 tmp = creator(surface, t);
            _m_cads.insert(std::make_pair(tmp.id(), tmp));
            it = _m_map.insert(it, std::make_pair(surface, tmp.id()));
        }
        return _m_cads.find(it->second)->second;
    }
    
private:
    //! return \c true if cad for \c surface's silhouette-curve is known
    bool is_cached(const Surface_3& surface) {
        return (_m_map.find(surface) != _m_map.end());
    }

    //!@}
    
    //!\name Caching for cads representing the cuts of surface pairs
    //!@{
    
public:
    //! return stored \c cad of cut-curve for \c surface1 and \c surface2
    Restricted_cad_3 operator()(const Surface_3& surface1, 
                                const Surface_3& surface2) {
        static Canonicalize_pair canonicalize;
        Surface_pair_3 pair = canonicalize(std::make_pair(surface1, surface2));
        typename Surface_pair_cache_3::iterator it = _m_pair_map.find(pair);
        if (it == _m_pair_map.end()) {
            Creator creator; // TODO single instance
            typename Creator::Non_cached_construction_tag t;
            Restricted_cad_3 tmp = creator(surface1, surface2, t);
            _m_cads.insert(std::make_pair(tmp.id(), tmp));
            it = _m_pair_map.insert(it, std::make_pair(pair, tmp.id()));
        }
        return _m_cads.find(it->second)->second;
    }
    
private:
    //! return \c true if cad for given surfaces cut-curve is known
    bool is_cached(const Surface_3& surface1, const Surface_3& surface2) {
        static Canonicalize_pair canonicalize;
        Surface_pair_3 pair = canonicalize(std::make_pair(surface1, surface2));
        return (_m_pair_map.find(pair) != _m_pair_map.end());
    }

    //!@}

    // Remark: We omit to cache sil+sil+cut - as this can be done in
    //         Surface_pair_3, and other classes can derive from it

private:
    
    // members
    //! cache for for cads
    Restricted_cad_map _m_cads;

    //! cache for silhouettes
    Surface_cache_3 _m_map;
    
    //! cache for cuts
    Surface_pair_cache_3 _m_pair_map;
};

} // namespace CGALi

// TODO document
template < class SurfaceZAtXyIsolatorTraits >
class Restricted_cad_3 : 
        public CGAL::Handle_with_policy < CGAL::Arrangement_2< 
typename SurfaceZAtXyIsolatorTraits::Arrangement_traits_2, 
        CGAL::Arr_extended_dcel< 
    typename SurfaceZAtXyIsolatorTraits::Arrangement_traits_2, 
        boost::optional< CGAL::P_dcel_info< SurfaceZAtXyIsolatorTraits > >, 
        boost::optional< CGAL::P_dcel_info< SurfaceZAtXyIsolatorTraits > >, 
        boost::optional< CGAL::P_dcel_info< SurfaceZAtXyIsolatorTraits > > 
> 
> 
> {
public:
    //! this instance' template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

#if DOXYGEN_RUNNING
    //! type of Surface_3
    typedef typename Surface_z_at_xy_isolator_traits::Surface_3 Surface_3;
#endif

    //! data type
    typedef CGAL::P_dcel_info< Surface_z_at_xy_isolator_traits > Data;

private:    
    //! optionalized data
    typedef boost::optional < Data > Opt_data;
    
#if DOXYGEN_RUNNING
    //! arrangment traits
    typedef typename Surface_z_at_xy_isolator_traits::Arrangement_traits_2 
    Arrangement_traits_2;
#endif

    //! dcel
    typedef CGAL::Arr_extended_dcel< Arrangement_traits_2, 
    Opt_data, Opt_data, Opt_data > Dcel;

public:
    //! rep class
    typedef CGAL::Arrangement_2< Arrangement_traits_2, Dcel > Rep;

private:
    //! base class
    typedef CGAL::Handle_with_policy< Rep > Base;
    
    typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits > Self;

    //! type of cache
    typedef CGALi::Restricted_cad_3_cache< Surface_z_at_xy_isolator_traits >
    Restricted_cad_3_cache;

public:

    //! Z_Stack type
    typedef typename Data::Z_stack Z_stack;
    
    // from arrangement
    typedef typename Dcel::Size Size;
    
#if DOXYGEN_RUNNING
    typedef typename Rep::Point_2 Point_2;
    typedef typename Rep::X_monotone_curve_2  X_monotone_curve_2;
#endif
    
private:
    typedef typename Rep::Vertex_iterator Vertex_iterator;
    typedef typename Rep::Edge_iterator Edge_iterator;
    typedef typename Rep::Halfedge_iterator Halfedge_iterator;
    typedef typename Rep::Face_iterator Face_iterator;

public:
    typedef typename Rep::Vertex_const_iterator Vertex_const_iterator;

    typedef typename Rep::Edge_const_iterator Edge_const_iterator;
private: 
    typedef typename Rep::Halfedge_const_iterator Halfedge_const_iterator;

public:
    typedef typename Rep::Face_const_iterator Face_const_iterator;

    typedef typename Rep::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
    
    typedef typename Rep::Outer_ccb_const_iterator Outer_ccb_const_iterator;
    
    typedef typename Rep::Inner_ccb_const_iterator Inner_ccb_const_iterator;
    
    typedef typename Rep::Isolated_vertex_const_iterator 
    Isolated_vertex_const_iterator;

public:
    typedef typename Rep::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
    
private:
    typedef Vertex_iterator Vertex_handle;
    typedef Halfedge_iterator Halfedge_handle;
    typedef Edge_iterator Edge_handle;
    typedef Face_iterator Face_handle;

public:
    typedef Vertex_const_iterator Vertex_const_handle;
private:
    typedef Halfedge_const_iterator Halfedge_const_handle;
public:
    typedef Edge_const_iterator Edge_const_handle;
    typedef Face_const_iterator Face_const_handle;

    // TODO choose better point location
    typedef CGAL::Arr_naive_point_location< Rep > Point_location;

public:
    //! type of X_coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;
    
    //! type of Boundary
    typedef typename Curve_kernel_2::Boundary Boundary;
    
    //!\name Static members
    //!@{
    
    //! returns static instance of cad_cache for silhouette- and cut-curves
    static Restricted_cad_3_cache& cad_cache() {
        static Restricted_cad_3_cache cache;
        return cache;
    }
    
    //!@}

private:

    //!\name Global
    //!@{

    //! make id unique and clear the rep
    void _renew() {
        this->copy_on_write();
        this->ptr()->clear();
    }

    //! construct cad for given curve
    static
    Restricted_cad_3 _construct_for_curve(const Curve_analysis_2& curve) {
        
        Restricted_cad_3 cad;
        cad._renew();
        
        // split curve into segments
        
        typename 
            Surface_z_at_xy_isolator_traits::Arrangement_traits_2
            geo_traits;
        ;
        
        std::list< Curve_analysis_2 > input;
        input.push_front(curve);
        
        std::list< X_monotone_curve_2 > xcurves;
        std::list< Point_2 > points;
        CGAL::make_x_monotone(
                input.begin(), input.end(),
                std::back_inserter(xcurves), 
                std::back_inserter(points),
                &geo_traits
        );
    
        // TODO use non-intersecting
        cad._insert_empty(xcurves.begin(), xcurves.end(),
                          points.begin(), points.end());
        
        return cad;
    }
    
    //! initializes a single feature
    template < class DcelConstHandle >
    void _init_feature(const Surface_3& surface, 
                       DcelConstHandle handle, 
                       CGAL::Dcel_feature feature) {
        
        CGAL::Object obj = CGAL::make_object(handle);
        this->non_const_handle(handle)->data() = Data(feature, obj);
        handle->data()->_set_dcel(feature, obj);
        handle->data()->_set_rs_id(this->id());
        
        handle->data()->_set_dcel(surface, feature, obj);
        
        handle->data()->_init_nk(surface);
    }
    
    //! initialize data objects
    void _init(const Surface_3& surface) {
        for (Vertex_const_iterator vit =
                 this->ptr()->vertices_begin();
             vit != this->ptr()->vertices_end(); vit++) {
            
            CGAL_precondition(!vit->data());
            
            _init_feature(surface, vit, CGAL::VERTEX);
            
            CGAL_postcondition(vit->data());
        }
        
        for (Halfedge_const_iterator hit =
                 this->ptr()->halfedges_begin();
             hit != this->ptr()->halfedges_end(); hit++) {
            
            if (!hit->data()) {
                
                CGAL_precondition(!hit->twin()->data());
                
                
                _init_feature(surface, hit, CGAL::EDGE);
                this->non_const_handle(hit->twin())->data() = *hit->data();
            }
            
            CGAL_postcondition(hit->data());
            CGAL_postcondition(hit->twin()->data());
        }
        
        for (Face_const_iterator fit =
                 this->ptr()->faces_begin();
             fit != this->ptr()->faces_end(); fit++) {
            
            CGAL_precondition(!fit->data());
            
            _init_feature(surface, fit, CGAL::FACE);
            
            CGAL_postcondition(fit->data());
        }
    }

    //! set final data
    void _finalize(const Surface_3& surface) {
        
        for (Vertex_const_iterator vit =
                 this->ptr()->vertices_begin();
             vit != this->ptr()->vertices_end(); vit++) {
            CGAL_precondition(nk(vit,surface).mult() != -1);
            vit->data()->_finalize_nk(surface);
        }
        
        for (Edge_const_iterator eit =
                 this->ptr()->edges_begin();
             eit != this->ptr()->edges_end(); eit++) {
            CGAL_precondition(nk(eit,surface).mult() != -1);
            eit->data()->_finalize_nk(surface);
        }
        
        for (Face_const_iterator fit =
                 this->ptr()->faces_begin();
             fit != this->ptr()->faces_end(); fit++) {
            CGAL_precondition(nk(fit,surface).mult() == -1);
            fit->data()->_finalize_nk(surface);
        }
    }
    
    //!@}
    
private:
    //!\name Arrangement methods 
    //!@{

    /*! insert a range of curves */
    template < class InputIterator >
    void _insert_non_intersecting_curves(
            InputIterator begin, InputIterator end
    ) {
        CGAL::insert_non_intersecting_curves(*(this->ptr()), begin, end);
    }

    /*! insert a range of curves */
    template < class CurveInputIterator, class PointInputIterator >
    void _insert_empty(
            CurveInputIterator cbegin, CurveInputIterator cend,
            PointInputIterator pbegin, PointInputIterator pend
            
    ) {
        CGAL::insert_empty(*(this->ptr()), 
                           cbegin, cend,
                           pbegin, pend);
    }
    
    /*! inserts an isolated point */
    void _insert_point(Point_2 pt) {
        CGAL::insert_point(*(this->ptr()),pt);
    }
    
    /*! remove edge from arrangement */
    void _remove_edge(Halfedge_handle he) {
        CGAL::remove_edge(*(this->ptr()), he);
    }

    /*! remove edge from arrangement */
    void _remove_vertex(Vertex_handle vh) {
        CGAL::remove_vertex(*(this->ptr()), vh);
    }

    /*! overlays two instances */
    void _overlay(Self rscad1, Self rscad2, bool factors_of_same_curve) {
        // use correct overlay instance
        CGAL::Arr_p_dcel_info_overlay_traits< Self > 
            ovltraits(*this, factors_of_same_curve);
        this->ptr()->clear();
        CGAL::overlay(*(rscad1.ptr()), *(rscad2.ptr()),
                      *(this->ptr()), ovltraits);
    }
    
    /*! overlays two instances */
    void _overlay(Self rscad1, Self rscad2, 
                  const Surface_3& surface, 
                  CGAL::Nk::Value_type type) {
        // use correct overlay instance
        CGAL::Arr_p_dcel_info_overlay_traits< Self > 
            ovltraits(surface, *this, type);
        this->ptr()->clear();
        CGAL::overlay(*(rscad1.ptr()), *(rscad2.ptr()),
                      *(this->ptr()), ovltraits);
    }
    

    //!@}

    //!\name Geometry and the Dcel
    //!@{
    
public:
    
    /*! perform point location */
    CGAL::Object locate(Point_2 pt) const {
        Point_location pl(*(this->ptr()));
        return pl.locate(pt);
    }
    
private:
    
    /*!\brief
     * use point location to check whether \c pt is part of given 
     * \c handle 
     */
    template < class DcelConstHandle >
    bool _point_on_dcel_handle(const Point_2& pt, 
                               DcelConstHandle handle) const {
        
        typedef DcelConstHandle Dcel_const_handle;
        // Perform the point-location query.
        CGAL::Object obj = this->locate(pt);
        
        // check if correct type
        Dcel_const_handle handle2;
        if (CGAL::assign (handle2, obj)) {
            if (handle->data()->id() == handle2->data()->id()) {
                return true;
            } else {
                return false;
            }
        } 
        // else
        return false;
    }

    /*!\brief 
     * constructs a rational point at coordinate \c x0, \c y0
     */
    static
    Point_2 _construct_point_with_rational_y(X_coordinate_1 x0, Boundary y0) {
        
        //! type of Polynomial 1
        typedef typename Polynomial_2::NT Polynomial_1;
        
        typedef CGAL::Fraction_traits< Boundary > FT;
        typedef typename FT::Numerator_type Integer;
        Integer yn;
        Integer yd;
        typename FT::Decompose decompose;
        decompose(y0, yn, yd);
        
        Polynomial_2 poly(Polynomial_1(-yn), Polynomial_1(yd));

        // construct curve analysis
        typename 
            Arrangement_traits_2::Curve_kernel_2::Construct_curve_2
            construct_curve = 
            Arrangement_traits_2::instance().kernel().
            construct_curve_2_object();
        typename 
            Arrangement_traits_2::Curve_kernel_2::Curve_analysis_2 hcurve = 
            construct_curve(poly);
        Point_2 pt(x0, hcurve, 0);
        return pt;
    }
    
    //! construct point, with rational x (or y) in interior of \c heh's curve
    static
    Point_2 _point_in_interior(const Halfedge_const_handle& heh) {
        Point_2 pt = _point_in_interior(heh->curve());
        CGAL_postcondition(
                Restricted_cad_3::cad_cache()(
                        heh->data()->ptr()->_m_rs_id
                )._point_on_dcel_handle(pt, heh)
        );
        return pt;
    }
    
    //! construct point, with rational x (or y) in interior of \c heh's curve
    static
    Point_2 _point_in_interior(const X_monotone_curve_2& cv) {

        Point_2 pt;
#if 1
        typename 
            Surface_z_at_xy_isolator_traits::Arrangement_traits_2
            geo_traits;

        typename Arrangement_traits_2::Construct_interior_vertex_2 
            construct_interior_vertex = 
            geo_traits.construct_interior_vertex_2_object();
        
        pt = construct_interior_vertex(cv);
#else 
        // TODO remove old code
        //! type of Event info
        typedef typename 
            Arrangement_traits_2::Curve_kernel_2::Curve_analysis_2::
            Status_line_1 
            Status_line_1;
        
        // find point in halfedge interior
        if (cv.is_vertical()) {
            // handle also vertical curves
            // find y0 between source.y and target.y
            
            Boundary y0(0);
            
            CGAL::Arr_parameter_space ps_y_min =
                cv.location(CGAL::ARR_MIN_END);
            CGAL::Arr_parameter_space ps_y_max =
                cv.location(CGAL::ARR_MAX_END);
            
            Point_2 min_p;
            if (ps_y_min == CGAL::ARR_INTERIOR) {
                min_p = cv.curve_end(CGAL::ARR_MIN_END);
            }
            Point_2 max_p;
            if (ps_y_max == CGAL::ARR_INTERIOR) {
                max_p = cv.curve_end(CGAL::ARR_MAX_END);
            }
            
            X_coordinate_1 x0 = cv.x();
            
            typename Arrangement_traits_2::Curve_kernel_2::Lower_boundary_y_2
                lower_boundary_y = 
                Arrangement_traits_2::instance().kernel().
                lower_boundary_y_2_object();
            typename Arrangement_traits_2::Curve_kernel_2::Upper_boundary_y_2
                upper_boundary_y = 
                Arrangement_traits_2::instance().kernel().
                upper_boundary_y_2_object();
            
            if (ps_y_min == CGAL::ARR_BOTTOM_BOUNDARY) { // -oo
                if (ps_y_max == CGAL::ARR_TOP_BOUNDARY) { // +oo
                    // nothing to do since y0 == 0 already
                } else {
                    Status_line_1 max_sl = 
                        max_p.curve().status_line_at_exact_x(x0);
                    int max_i = max_p.arcno();
                    y0 = lower_boundary_y(max_sl.algebraic_real_2(max_i)) - 
                        Boundary(1);
                }
            } else {
                if (ps_y_max == CGAL::ARR_TOP_BOUNDARY) { // +oo
                    Status_line_1 min_sl = 
                        min_p.curve().status_line_at_exact_x(x0);
                    int min_i = min_p.arcno();
                    y0 = upper_boundary_y(min_sl.algebraic_real_2(min_i)) + 
                        Boundary(1);
                } else {
                    // both endpoints are finite
                    Status_line_1 min_sl = 
                        min_p.curve().status_line_at_exact_x(x0);
                    int min_i = min_p.arcno();
                    
                    Status_line_1 max_sl = 
                        max_p.curve().status_line_at_exact_x(x0);
                    int max_i = max_p.arcno();
                    
                    typename 
                        Arrangement_traits_2::Curve_kernel_2::
                        Boundary_between_y_2 
                        boundary_between_y = 
                        Arrangement_traits_2::instance().kernel().
                        boundary_between_y_2_object();
                    
                    y0 = boundary_between_y(min_sl.algebraic_real_2(min_i),
                                            max_sl.algebraic_real_2(max_i));
                }
            }
            
            pt = _construct_point_with_rational_y(x0, y0);
        } else {
            // construct point on curve
            X_coordinate_1 x0(cv.boundary_in_x_range_interior());
            pt = Point_2(x0, cv.curve(), cv.arcno());
        }
#endif
        return pt;
    }

    //!@}

public:
    //!\name Counting
    //!@{
    /*! Check whether the arrangement is empty. */
    bool is_empty () const
    {
        return (this->ptr()->is_empty());
    }
    
    /*! Get the number of arrangement vertices. */
    Size number_of_vertices () const
    {
        return (this->ptr()->number_of_vertices());
    }
    
    /*! Get the number of isolated arrangement vertices. */
    Size number_of_isolated_vertices () const
    {
        return (this->ptr()->number_of_isolated_vertices());
    }
    
private:
    /*! Get the number of arrangement halfedges (the result is always even). */
    Size number_of_halfedges () const
    {
        return (this->ptr()->number_of_halfedges());
    }

public:
    /*! Get the number of arrangement edges. */
    Size number_of_edges () const
    {
        return (this->ptr()->number_of_halfedges() / 2);
    }

    /*! Get the number of arrangement faces. */
    Size number_of_faces () const
    {
        return (this->ptr()->number_of_faces());
    }
    
    /*! Get the number of arrangement faces. */
    Size number_of_unbounded_faces () const
    {
        return (this->ptr()->number_of_unbounded_faces());
    }
    
    //!@}

    //!\name Traversal functions for the arrangement vertices.
    //!@{

    /*! Get a const iterator for the first vertex in the arrangement. */
    Vertex_const_iterator vertices_begin() const
    { 
        return (this->ptr()->vertices_begin());
    }
    
    /*! Get a past-the-end const iterator for the arrangement vertices. */
    Vertex_const_iterator vertices_end() const
    {
        return (this->ptr()->vertices_end());
    }
    //!@}

private:
    //\name Traversal functions for the arrangement halfedges.
    //!@{

    /*! Get a const iterator for the first halfedge in the arrangement. */
    Halfedge_const_iterator halfedges_begin() const
    { 
        return (this->ptr()->halfedges_begin());
    }
  
    /*! Get a past-the-end const iterator for the arrangement halfedges. */
    Halfedge_const_iterator halfedges_end() const
    {
        return (this->ptr()->halfedges_end());
    }
    //!@}
    
public:
    //\name Traversal functions for the arrangement edges.
    //!@{

    /*! Get a const iterator for the first edge in the arrangement. */
    Edge_const_iterator edges_begin() const
    { 
        return (this->ptr()->edges_begin());
    }
  
    /*! Get a past-the-end const iterator for the arrangement halfedges. */
    Edge_const_iterator edges_end() const
    {
        return (this->ptr()->edges_end());
    }
    //!@}
    
    //!\name Traversal functions for the arrangement faces.
    //!@{

    /*! Get a const iterator for the first face in the arrangement. */
    Face_const_iterator faces_begin() const
    { 
        return (this->ptr()->faces_begin());
    }
    
    /*! Get a past-the-end const iterator for the arrangement faces. */
    Face_const_iterator faces_end() const
    {
        return (this->ptr()->faces_end());
    }
    //!@}

    //!\name Conversions of special halfedges to edges
    //!{
    
    //! convert to Edge_const_handle
    Edge_const_handle convert_to_edge_const_handle(
            const Halfedge_around_vertex_const_circulator& circ) const {
        Halfedge_const_handle heh(circ);
        Edge_const_handle eh = edges_begin();
        std::advance(eh, std::distance(halfedges_begin(),heh)/2);
        return eh;
    }

     //! convert to Edge_const_handle
    Edge_const_handle convert_to_edge_const_handle(
            const Ccb_halfedge_const_circulator& circ) const {
        Halfedge_const_handle heh(circ);
        Edge_const_handle eh = edges_begin();
        std::advance(eh, std::distance(halfedges_begin(),heh)/2);
        return eh;
    }
    
    //!@}

private:
    //!\name Casting away constness for handle types.
    //!@{
    //! converts const vertex handle to non-const vertex handle
    typename Rep::Vertex_handle non_const_handle (Vertex_const_handle vh)
    {
        return (this->ptr()->non_const_handle(vh));
    }
    
    //! converts const halfedge handle to non-const halfedge handle
    typename Rep::Halfedge_handle non_const_handle (Halfedge_const_handle hh)
    {
        return (this->ptr()->non_const_handle(hh));
    }
    
    //! converts const face handle to non-const face handle
    typename Rep::Face_handle non_const_handle (Face_const_handle fh)
    {
        return (this->ptr()->non_const_handle(fh));
    }
    //!@}
    
    //!\name Nk
    //!@{

private:
    //! set value for vertex handle
    void _set_nk_value(const Vertex_const_handle& vh, 
                      const Surface_3& surface,
                      CGAL::Nk::Value_type type, 
                      int value) const {
        switch (type) {
        case CGAL::Nk::MULT:
            vh->data()->_nk(surface).set_mult(value);
            break;
        case CGAL::Nk::N:
            vh->data()->_nk(surface).set_n(value);
            break;
        case CGAL::Nk::K:
            vh->data()->_nk(surface).set_k(value);
            break;
        }
    }

    //! set value for halfedge handle
    void _set_nk_value(const Halfedge_const_handle& heh, 
                      const Surface_3& surface,
                      CGAL::Nk::Value_type type, 
                      int value) const {
        switch (type) {
        case CGAL::Nk::MULT:
            heh->data()->_nk(surface).set_mult(value);
            break;
        case CGAL::Nk::N:
            heh->data()->_nk(surface).set_n(value);
            break;
        case CGAL::Nk::K:
            heh->data()->_nk(surface).set_k(value);
            break;
        }
    }

    //! set value for edge handle
    void _set_nk_value(const Edge_const_handle& eh, 
                      const Surface_3& surface,
                      CGAL::Nk::Value_type type, 
                      int value) const {
        switch (type) {
        case CGAL::Nk::MULT:
            eh->data()->_nk(surface).set_mult(value);
            break;
        case CGAL::Nk::N:
            eh->data()->_nk(surface).set_n(value);
            break;
        case CGAL::Nk::K:
            eh->data()->_nk(surface).set_k(value);
            break;
        }
    }

    //! set value for face handle
    void _set_nk_value(const Face_const_handle& fh, 
                      const Surface_3& surface,
                      CGAL::Nk::Value_type type, 
                       int value) const {
        switch (type) {
        case CGAL::Nk::MULT:
            fh->data()->_nk(surface).set_mult(value);
            break;
        case CGAL::Nk::N:
            fh->data()->_nk(surface).set_n(value);
            break;
        case CGAL::Nk::K:
            fh->data()->_nk(surface).set_k(value);
            break;
        }
    }
    
public:
    //! nk values for vertex handle
    const CGAL::Nk& nk(const Vertex_const_handle& vh, 
                      const Surface_3& surface) const {
        return vh->data()->_nk(surface);
    }

private:
    //! nk values for halfedge handle
    const CGAL::Nk& nk(const Halfedge_const_handle& heh, 
                      const Surface_3& surface) const {
        return heh->data()->_nk(surface);
    }

public:
    //! nk values for edge handle
    const CGAL::Nk& nk(const Edge_const_handle& eh, 
                      const Surface_3& surface) const {
        return eh->data()->_nk(surface);
    }

    //! nk values for face handle
    const CGAL::Nk& nk(const Face_const_handle& fh, 
                      const Surface_3& surface) const {
        return fh->data()->_nk(surface);
    }
    
    //!@}

public:
    //!\name Silhouette/Cut
    //!@{
    
    //! returns whether silhouette at given vertex handle exists
    bool has_silhouette(const Vertex_const_handle& vh) const {
        return vh->data()->_has_silhouette();
    }
private:
    //! returns whether silhouette at given halfedge handle exists
    bool has_silhouette(const Halfedge_const_handle& heh) const {
        return heh->data()->_has_silhouette();
    }
public:
    //! returns whether silhouette at given edge handle exists
    bool has_silhouette(const Edge_const_handle& eh) const {
        return eh->data()->_has_silhouette();
    }

    //! returns whether silhouette at given face handle exists
    bool has_silhouette(const Face_const_handle& fh) const {
        return fh->data()->_has_silhouette();
    }


    //! returns whether cut at given vertex handle exists
    bool has_cut(const Vertex_const_handle& vh) const {
        return vh->data()->_has_cut();
    }
private:
    //! returns whether cut at given halfedge handle exists
    bool has_cut(const Halfedge_const_handle& heh) const {
        return heh->data()->_has_cut();
    }
public:
    //! returns whether cut at given edge handle exists
    bool has_cut(const Edge_const_handle& eh) const {
        return eh->data()->_has_cut();
    }
    
    //! returns whether cut at given face handle exists
    bool has_cut(const Face_const_handle& fh) const {
        return fh->data()->_has_cut();
    }
    
    //! returns whether silhouette of \c surface at \c vh exisst
    bool has_silhouette(const Vertex_const_handle& vh, 
                        const Surface_3& surface) const {
        return vh->data()->_has_silhouette(surface);
    }
private:
    //! returns whether silhouette of \c surface at \c heh exists
    bool has_silhouette(const Halfedge_const_handle& heh, 
                        const Surface_3& surface) const {
        return heh->data()->_has_silhouette(surface);
    }
public:
    //! returns whether silhouette of \c surface at \c eh exisst
    bool has_silhouette(const Edge_const_handle& eh, 
                        const Surface_3& surface) const {
        return eh->data()->_has_silhouette(surface);
    }
    
    //! returns whether silhouette of \c surface at \c fh exisst
    bool has_silhouette(const Face_const_handle& fh, 
                        const Surface_3& surface) const {
        return fh->data()->_has_silhouette(surface);
    }
    
    //! returns whether \c vh supports vertical line of \c surface
    bool supports_vertical_line(const Vertex_const_handle& vh,
                                const Surface_3& surface) const {
        return vh->data()->_supports_vertical_line(surface);
    }
    
    //! returns whether \c surface supports a vertical line at \c point
    bool supports_vertical_line_at(const Point_2& point, 
                                   const Surface_3& surface) {
        CGAL::Object obj = this->locate(point);
        Face_const_handle fh;
        if (CGAL::assign(fh,obj)) {
            return false;
        }
        Halfedge_const_handle heh;
        if (CGAL::assign(heh, obj)) {
            return (this->nk(heh, surface).n() == -1);
        }
        Vertex_const_handle vh;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vh, obj);
        CGAL_assertion(check);
        return (this->nk(vh, surface).n() == -1);
    };
    
    //! returns whether  \c surface is vertical over \c vh
    bool is_vertical(const Vertex_const_handle& vh, 
                     const Surface_3& surface) const {
        return vh->data()->_is_surface_vertical(surface);
    }
private:
    //! returns whether  \c surface is vertical over \c heh
    bool is_vertical(const Halfedge_const_handle& heh, 
                        const Surface_3& surface) const {
        return heh->data()->_is_surface_vertical(surface);
    }
public:
    //! returns whether  \c surface is vertical over \c eh
    bool is_vertical(const Edge_const_handle& eh, 
                        const Surface_3& surface) const {
        return eh->data()->_is_surface_vertical(surface);
    }
    
    //! returns whether  \c surface is vertical over \c fh
    bool is_vertical(const Face_const_handle& fh, 
                        const Surface_3& surface) const {
        return false;
    }


    //! returns whether cut the two surface at \c vh exisst
    bool has_cut(const Vertex_const_handle& vh, 
                 const Surface_3& surface1,
                 const Surface_3& surface2) const {
        return vh->data()->_has_cut(surface1, surface2);
    }
private:
    //! returns whether cut the two surface at \c heh exisst
    bool has_cut(const Halfedge_const_handle& heh, 
                 const Surface_3& surface1,
                 const Surface_3& surface2) const {
        return heh->data()->_has_cut(surface1, surface2);
    }
public:
    //! returns whether cut the two surface at \c eh exisst
    bool has_cut(const Edge_const_handle& eh, 
                 const Surface_3& surface1,
                 const Surface_3& surface2) const {
        return eh->data()->_has_cut(surface1, surface2);
    }

    //! returns whether cut the two surface at \c fh exisst
    bool has_cut(const Face_const_handle& fh, 
                 const Surface_3& surface1,
                 const Surface_3& surface2) const {
        return fh->data()->_has_cut(surface1, surface2);
    }

    //!@}

    //!\name Sample Points
    //!@{
    
    //! returns sample point for vertex handle
    Point_2 sample_point(const Vertex_const_handle& vh) const {
        CGAL_precondition(vh->data()->_rs_id() == this->id());
        return vh->data()->_sample_point_for_vertex_handle(vh);
    }
private:
    //! returns sample point for given halfedge handle
    Point_2 sample_point(const Halfedge_const_handle& heh) const {
        CGAL_precondition(heh->data()->_rs_id() == this->id());
        return heh->data()->_sample_point_for_halfedge_handle(heh);
    }
public:
    //! returns sample point for given edge handle
    Point_2 sample_point(const Edge_const_handle& eh) const {
        CGAL_precondition(eh->data()->_rs_id() == this->id());
        Halfedge_const_handle heh = eh;
        return heh->data()->_sample_point_for_halfedge_handle(heh);
    }

    //! returns sample point for given edge handle
    Point_2 sample_point(const Face_const_handle& fh) const {
        CGAL_precondition(fh->data()->_rs_id() == this->id());
        return fh->data()->_sample_point_for_face_handle(fh);
    }

    //!@}
    
    //!\name Isolators
    //!@{

    //! returns isolator (if existing) for given surface at given handle
    boost::optional< Z_at_xy_isolator> 
    isolator(const Vertex_const_handle& vh,
             const Surface_3& surface) const {
        Z_stack z_stackp = z_stack(vh);
        bool empty;
        if (z_stackp._knows_isolator(surface, empty)) {
            return z_stackp._isolator(surface);
        }
        // else
        return boost::none;
    }

private:
    //! returns isolator (if existing) for given surface at given handle
    boost::optional< Z_at_xy_isolator> 
    isolator(const Halfedge_const_handle& heh,
             const Surface_3& surface) const {
        Z_stack z_stackp = z_stack(heh);
        bool empty;
        if (z_stackp._knows_isolator(surface, empty)) {
            return z_stackp._isolator(surface);
        }
        // else
        return boost::none;
    }

public:
    //! returns isolator (if existing) for given surface at given handle
    boost::optional< Z_at_xy_isolator> 
    isolator(const Edge_const_handle& eh,
             const Surface_3& surface) const {
        Z_stack z_stackp = z_stack(eh);
        bool empty;
        if (z_stackp._knows_isolator(surface, empty)) {
            return z_stackp._isolator(surface);
        }
        // else
        return boost::none;
    }
    
    //! returns isolator (if existing) for given surface at given handle
    boost::optional< Z_at_xy_isolator> 
    isolator(const Face_const_handle& fh,
             const Surface_3& surface) const {
        Z_stack z_stackp = z_stack(fh);
        bool empty;
        if (z_stackp._knows_isolator(surface, empty)) {
            return z_stackp._isolator(surface);
        }
        // else
        return boost::none;
    }
    
    //! returns isolator (if existing) for \c surface for given point \c pt
    boost::optional< Z_at_xy_isolator >
    isolator_for(const Point_2& pt,
                const Surface_3& surface) const {
        Z_stack z_stackp = z_stack_for(pt);
        bool empty;
        if (z_stackp._knows_isolator(surface, empty)) {
            return z_stackp._isolator(surface);
        }
        // else
        return boost::none;
    }
    
#if 1
    // TODO reimplement avoiding cache!
    //! returns isolator (if existing) for \c surface for given point \c pt
    boost::optional< Z_at_xy_isolator >
    isolator_at(const Point_2& pt,
                const Surface_3& surface) const {

        
        typename Surface_z_at_xy_isolator_traits::
            Construct_isolator
            construct_isolator; 
        // TODO (construct_isolator_object())
        // construct isolator
        
        CGAL::Object obj = (cad_cache()(surface)).locate(pt);

        CGAL::Nk nk;
        CGAL::Dcel_feature feature;

        Face_const_handle fh;
        Halfedge_const_handle heh;
        Vertex_const_handle vh;
        if (CGAL::assign(fh,obj)) {
            nk = fh->data()->_nk(surface);
            feature = CGAL::FACE;
        } else if (CGAL::assign(heh, obj)) {
            nk = heh->data()->_nk(surface);
            feature = CGAL::EDGE;
        } else {
            CGAL_assertion_code(bool check =)
                CGAL::assign (vh, obj);
            CGAL_assertion(check);
            nk = vh->data()->_nk(surface);
            feature = CGAL::VERTEX;
        }
        
        Z_at_xy_isolator isolator = 
            construct_isolator(surface, pt, nk, feature);
        
        return isolator;
    }
#endif
    
    //!@}
    
    //!\name Z_Stacks
    //!@{
    
    //! returns z_stack for given vertex handle
    Z_stack z_stack(const Vertex_const_handle& vh) const {
        CGAL_precondition(vh->data()->_rs_id() == this->id());
        return vh->data()->_z_stack_for_vertex_handle(vh);
    }
private:
    //! returns z_stack for given halfedge handle
    Z_stack z_stack(const Halfedge_const_handle& heh) const {
        CGAL_precondition(heh->data()->_rs_id() == this->id());
        return heh->data()->_z_stack_for_halfedge_handle(heh);
    }
public:
    //! returns z_stack for given edge handle
    Z_stack z_stack(const Edge_const_handle& eh) const {
        CGAL_precondition(eh->data()->_rs_id() == this->id());
        Halfedge_const_handle heh = eh;
        return heh->data()->_z_stack_for_halfedge_handle(heh);
    }

    //! returns z_stack for given edge handle
    Z_stack z_stack(const Face_const_handle& fh) const {
        CGAL_precondition(fh->data()->_rs_id() == this->id());
        return fh->data()->_z_stack_for_face_handle(fh);
    }

    //! returns z_stack for given point
    Z_stack z_stack_for(const Point_2& point) const {
        CGAL::Object obj = this->locate(point);
        Face_const_handle fh;
        if (CGAL::assign(fh,obj)) {
            return this->z_stack(fh);
        }
        Halfedge_const_handle heh;
        if (CGAL::assign(heh, obj)) {
            return this->z_stack(heh);
        }
        Vertex_const_handle vh;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vh, obj);
        CGAL_assertion(check);
        return this->z_stack(vh);
    };

#if 0
    //! returns z_stack for given point
    std::pair< Z_stack, CGAL::Dcel_feature >
    z_stack_at(const Point_2& point) const {
        CGAL::Object obj = this->locate(point);
        Face_const_handle fh;
        if (CGAL::assign(fh,obj)) {
            return std::make_pair(fh->data()->_z_stack(point, CGAL::FACE),
                                  CGAL::FACE);
        }
        Halfedge_const_handle heh;
        if (CGAL::assign(heh, obj)) {
            return std::make_pair(heh->data()->_z_stack(point, CGAL::EDGE),
                                  CGAL::EDGE);
        }
        Vertex_const_handle vh;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vh, obj);
        CGAL_assertion(check);
        return std::make_pair(this->z_stack(vh), CGAL::VERTEX);
    };
#endif
    
    //! returns z_stack of silhouette-curve of \c surface at given vertex handle
    std::pair< Z_stack, CGAL::Dcel_feature > 
    z_stack(const Vertex_const_handle& vh, 
            const Surface_3& surface) const {
        return vh->data()->_z_stack_of_surface(surface);
    }
private:
    //! returns z_stack of silhouette-curve of \c surface at given edge handle
    std::pair< Z_stack, CGAL::Dcel_feature > 
    z_stack(const Halfedge_const_handle& heh, 
            const Surface_3& surface) const {
        return heh->data()->_z_stack_of_surface(surface);
    }
public:
    //! returns z_stack of silhouette-curve of \c surface at given edge handle
    std::pair < Z_stack, CGAL::Dcel_feature > 
    z_stack(const Edge_const_handle& eh, 
            const Surface_3& surface) const {
        return eh->data()->_z_stack_of_surface(surface);
    }
    
    //! returns z_stack of silhouette-curve of \c surface at given face handle
    std::pair < Z_stack, CGAL::Dcel_feature > 
    z_stack(const Face_const_handle& fh, 
            const Surface_3& surface) const {
        return fh->data()->_z_stack_of_surface(surface);
    }

    //! returns z_stack for given point
    std::pair< Z_stack, CGAL::Dcel_feature >
    z_stack_for(const Point_2& point, const Surface_3& surface) const {
        CGAL::Object obj = this->locate(point);
        Face_const_handle fh;
        if (CGAL::assign(fh,obj)) {
            return this->z_stack(fh, surface);
        }
        Halfedge_const_handle heh;
        if (CGAL::assign(heh, obj)) {
            return this->z_stack(heh, surface);
        }
        Vertex_const_handle vh;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vh, obj);
        CGAL_assertion(check);
        return this->z_stack(vh, surface);
    };

#if 0
    //! returns z_stack for given point
    std::pair< Z_stack, CGAL::Dcel_feature >
    z_stack_at(const Point_2& point, const Surface_3& surface) const {
        CGAL::Object obj = this->locate(point);
        Face_const_handle fh;
        if (CGAL::assign(fh,obj)) {
            return fh->data()->_z_stack_of_surface_at(point, surface);
        }
        Halfedge_const_handle heh;
        if (CGAL::assign(heh, obj)) {
            return heh->data()->_z_stack_of_surface_at(point, surface);
        }
        Vertex_const_handle vh;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vh, obj);
        CGAL_assertion(check);
        return this->z_stack(vh, surface);
    };
#endif
    //!@}
    
    //!\name Adjacencies
    //!{@

    /*\brief Computes adjacency when going from \c surface sheet \c sheet in 
     * \c from to \c to, and return the corresponding z-cell as \c .first
     * (std::pair< Z_stack, std::pair<int,int> >) and the sheet numbers 
     * of \c surface at \c to as \c .second .
     */
    template < class DcelConstHandle1, class DcelConstHandle2 >
    std::pair< std::pair< Z_stack, std::pair< int, int > >, 
               std::pair< int, int > > adjacency(
                       DcelConstHandle1 from, 
                       const Surface_3& surface, int sheet, 
                       DcelConstHandle2 to) const {
        
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_from = 
            z_stack(from, surface);
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_to = 
            z_stack(to, surface);
        
        CGAL_assertion(from->data()->_rs_id() == to->data()->_rs_id());

        boost::optional< 
            std::pair< CGAL::Dcel_feature,CGAL::Object >
            > sffrom = from->data()->_dcel(surface);        
        CGAL_assertion(sffrom);
        
        boost::optional< 
            std::pair< CGAL::Dcel_feature,CGAL::Object >
            > sfto = to->data()->_dcel(surface);        
        CGAL_assertion(sfto);

        CGAL::Object fromobj = sffrom->second;
        CGAL::Object toobj = sfto->second;
        
        std::pair< int, int > sheet_to =
            z_stack_from.first.adjacency(
                    surface, sheet, 
                    z_stack_from.second, 
                    fromobj,
                    z_stack_to.first, 
                    z_stack_to.second, 
                    toobj
            );
        
        Z_stack goal = z_stack(to);

        if (sheet_to.first > sheet_to.second) {
            // not involved
            return std::make_pair(std::make_pair(goal, sheet_to), sheet_to);
        }
        if (sheet_to.first == -1) {
            if (sheet_to.second == -1) {
                // TODO handling at minus inf
                return std::make_pair(std::make_pair(goal, sheet_to), 
                                      sheet_to);
            } else if (sheet_to.second ==
                       z_stack_to.first.number_of_z_cells()) {
                int z = goal.number_of_z_cells();
                return std::make_pair(
                        std::make_pair(goal, 
                                       std::make_pair(-1, z)), 
                        sheet_to
                );
            } else {
                int z_to_high = 
                    goal.z_level_of_sheet(surface, sheet_to.second);
                return std::make_pair(
                        std::make_pair(goal, 
                                       std::make_pair(-1, z_to_high)), 
                        sheet_to
                );
            }
        } else if (sheet_to.second == z_stack_to.first.number_of_z_cells()) {
            int z = goal.number_of_z_cells();
            if (sheet_to.first == z_stack_to.first.number_of_z_cells()) {
                // TODO handling at plus inf
                return std::make_pair(
                        std::make_pair(goal, std::make_pair(z,z)), 
                        sheet_to
                );
            } else {
                int z_to_low = goal.z_level_of_sheet(surface, sheet_to.first);
                
                return std::make_pair(
                        std::make_pair(goal, std::make_pair(z_to_low, z)), 
                        sheet_to
                );
            }
        }
        
        // else
        int z_to_low = goal.z_level_of_sheet(surface, sheet_to.first);
        int z_to_high = goal.z_level_of_sheet(surface, sheet_to.second);
        
        return std::make_pair(
                std::make_pair(goal, std::make_pair(z_to_low, z_to_high)), 
                sheet_to
        );
    }
    
    //!@}

    //!\name Surface Adjacency
    //!@{

    /*!\brief Computes adjacency when going from \c surface sheet \c sheet in 
     * \c from to \c to, and return the corresponding in sil-arr as .first
     * and the sheet number of \c surface 
     * at \c to as \c .second .
     */
    template < class DcelConstHandle1, class DcelConstHandle2 >
    std::pair< Z_stack, std::pair< int, int > > surface_adjacency(
            DcelConstHandle1 from, 
            const Surface_3& surface, int sheet, 
            DcelConstHandle2 to) const {
        
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_from = 
            z_stack(from, surface);
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_to = 
            z_stack(to, surface);
        
        CGAL_assertion(from->data()->_rs_id() == to->data()->_rs_id());

        boost::optional< 
            std::pair< CGAL::Dcel_feature,CGAL::Object >
            > sffrom = from->data()->_dcel(surface);        
        CGAL_assertion(sffrom);
        
        boost::optional< 
            std::pair< CGAL::Dcel_feature,CGAL::Object >
            > sfto = to->data()->_dcel(surface);        
        CGAL_assertion(sfto);

        CGAL::Object fromobj = sffrom->second;
        CGAL::Object toobj = sfto->second;
        
        std::pair< int, int > sheet_to =
            z_stack_from.first.adjacency(
                    surface, sheet, 
                    z_stack_from.second, 
                    fromobj,
                    z_stack_to.first, 
                    z_stack_to.second, 
                    toobj
            );
        
        return std::make_pair(z_stack_to.first, sheet_to);
    }
    
    template < class DcelConstHandle1, class DcelConstHandle2 >
    CGAL::Adjacencies_3 adjacency(DcelConstHandle1 from, 
                                 const Surface_3& surface,
                                 DcelConstHandle2 to) const {
        
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_from = 
            z_stack(from, surface);
        std::pair< Z_stack, CGAL::Dcel_feature > z_stack_to = 
            z_stack(to, surface);
        
        CGAL_assertion(from->data()->_rs_id() == to->data()->_rs_id());
        
        boost::optional< 
        std::pair< CGAL::Dcel_feature,CGAL::Object > > sffrom = 
            from->data()->_dcel(surface);        
        CGAL_assertion(sffrom);
        
        boost::optional< 
        std::pair< CGAL::Dcel_feature,CGAL::Object > > sfto = 
            to->data()->_dcel(surface);        
        CGAL_assertion(sfto);
        
        CGAL::Object fromobj = sffrom->second;
        CGAL::Object toobj = sfto->second;
        
        typename Surface_z_at_xy_isolator_traits::Adjacency 
            compute_adj; // TODO single instance
        
        CGAL::Adjacencies_3 adj = 
            compute_adj(surface,
                        // is valid as we just constructed z-stack
                        *isolator(from, surface),
                        z_stack_from.second, fromobj, 
                        from->data()->_supports_vertical_line(surface),
                        // is valid as we just constructed z-stack
                        *isolator(to, surface),
                        z_stack_to.second, toobj,
                        to->data()->_supports_vertical_line(surface)
            );
        return adj;
    }

    //!@}
    
    // friends
    // for non_const_handle, _insert*, _overlay, halfedge-stuff
    friend class 
    CGAL::Create_restricted_cad_3< Surface_z_at_xy_isolator_traits >;
    friend class 
    CGAL::Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits >;
    friend class CGAL::P_dcel_info< Surface_z_at_xy_isolator_traits >;
    friend class CGAL::Arr_p_dcel_info_overlay_traits< Self >;
    // for types
    friend class CGAL::Z_stack< Surface_z_at_xy_isolator_traits, Data >;
    // for all private members
    friend class CGAL::Restricted_cad_3_accessor< Self >;
    friend class CGAL::CGALi::Construct_nk_decomposition_2< 
    Surface_z_at_xy_isolator_traits >;
    
};

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_H
// EOF
