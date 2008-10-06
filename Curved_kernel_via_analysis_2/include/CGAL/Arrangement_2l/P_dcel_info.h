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
// File          : include/SoX/GAPS/P_dcel_info.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file SoX/GAPS/P_dcel_info.h
 * \brief definition of P_dcel_info class template
 */

#ifndef SoX_GAPS_P_DCEL_INFO_H
#define SoX_GAPS_P_DCEL_INFO_H 1

#include <CGAL/config.h>

#include <algorithm>
#include <set>
#include <map>

#include <boost/optional.hpp>
#include <boost/none.hpp>

#include <CGAL/function_objects.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/algorithm.h>

#include <CGAL/ipower.h>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/macros.h>
#include <CGAL/Arrangement_2l/Z_stack.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_enums.h>
#include <CGAL/Arrangement_2l/Arr_p_dcel_info_overlay_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief
 * Representation class for P_dcel_info
 */
template < class SurfaceZAtXyIsolatorTraits >
class P_dcel_info_rep {
    
public:
    // types
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

    //! the class itself
    typedef P_dcel_info_rep< Surface_z_at_xy_isolator_traits > Self;
    
    //! type of Less for surfaces
    typedef typename Surface_3::Surface_less_than Surface_less_than;

    //! type of surface container
    typedef std::set< Surface_3, Surface_less_than > Surface_container;

    //! type of nk map for silhouettes
    typedef std::map< Surface_3, CGAL::Nk, Surface_less_than > 
    Silhouettes_nk_map;
    
    //! type of Silhouette map
    typedef std::map< Surface_3, int, Surface_less_than > Silhouettes_map;
    
    //! type of Surface_pair
    typedef std::pair< Surface_3, Surface_3 > Surface_pair;
    
    //! type of Less of surface pair
    typedef CGAL::Pair_lexicographical_less_than< Surface_3, Surface_3, 
                                Surface_less_than, Surface_less_than > 
    Surface_pair_less;
    
    //! type of Cuts map
    typedef std::map< Surface_pair, int, Surface_pair_less > Cuts_map;

    //! type of Dcel_data
    typedef CGAL::P_dcel_info< Surface_z_at_xy_isolator_traits > Dcel_data;
    
    //! type of Z_stack
    typedef CGAL::Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data > 
    Z_stack;

    //! type of id
    typedef std::ptrdiff_t Rs_id;
    
private:
    //! originating dcel-handles for single surfaces
    typedef std::map< Surface_3, std::pair< CGAL::Dcel_feature, CGAL::Object >, 
    Surface_less_than > Surface_dcel_handle_map;

    //! originating dcel-handles for single surfaces
    typedef std::map< Surface_pair, 
    std::pair< CGAL::Dcel_feature,CGAL::Object >, Surface_pair_less >
    Surface_pair_dcel_handle_map;


    //! Standard constructor
    P_dcel_info_rep (const CGAL::Dcel_feature feature, 
                     const CGAL::Object dcel_handle) :
        _m_rs_id(0),
        _m_dcel_feature(feature),
        _m_dcel_handle(dcel_handle) {
    }
    
private:
    // data members
    //! stored id of handle
    mutable Rs_id _m_rs_id;

    //! stored feature
    mutable CGAL::Dcel_feature _m_dcel_feature;

    //! stored for this dcel-handle
    mutable CGAL::Object _m_dcel_handle;

    //! stored point
    mutable boost::optional< Point_2 > _m_point;
    
    //! stored z-stack
    mutable boost::optional < Z_stack > _m_z_stack;
    
    //! list of surfaces
    mutable Surface_container _m_surfaces;
    
    //! map for Nk objects for silhouettes
    mutable Silhouettes_nk_map _m_silhouettes_nk;

    //! indication whether silhouette is stored
    mutable bool _m_has_silhouette_curve;

    //! indication whether vertical surface exists
    mutable bool _m_has_vertical_surface;
    
    //! indication whether vertical line exists
    mutable bool _m_has_vertical_line;
    
    //! list of involved cut curves
    mutable Cuts_map _m_cuts;

    //! stored boundary class in x
    mutable boost::optional < CGAL::Arr_boundary_type > _m_boundary_type_in_x;

    //! stored boundary class in y
    mutable boost::optional < CGAL::Arr_boundary_type > _m_boundary_type_in_y;
    
    //! stored data
    mutable std::list< Z_stack > _m_z_stacks;

    // TODO _m_intermediate_stack_cells
    
    //! stores originating dcel handles for surface
    mutable Surface_dcel_handle_map _m_surface_dcel_handle_map;

    //! stores originating dcel handles for surface pairs
    mutable Surface_pair_dcel_handle_map _m_surface_pair_dcel_handle_map;

    //! friends
    friend class P_dcel_info< Surface_z_at_xy_isolator_traits >;
    friend class Restricted_cad_3< Surface_z_at_xy_isolator_traits>;
};

} // namespace CGALi


/*!\brief
 * Attach information about surface to planar dcel.
 */
template < class SurfaceZAtXyIsolatorTraits >
class P_dcel_info : public 
::CGAL::Handle_with_policy< CGALi::P_dcel_info_rep< SurfaceZAtXyIsolatorTraits > > {
public:

     //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

    //! type of rep;
    typedef CGALi::P_dcel_info_rep< Surface_z_at_xy_isolator_traits > Rep;
    
    //! type of Base
    typedef CGAL::Handle_with_policy< Rep > Base;

    //! the class itself
    typedef P_dcel_info< Surface_z_at_xy_isolator_traits > Self;

    //! type of Z_stack
    typedef typename Rep::Z_stack Z_stack;

    //! type of Rs_id
    typedef typename Rep::Rs_id Rs_id;
    
protected:

    //! type of Silhouettes map
    typedef typename Rep::Silhouettes_map Silhouettes_map;

    //! type of Silhouettes nk map
    typedef typename Rep::Silhouettes_nk_map Silhouettes_nk_map;

    //! type of Cuts map
    typedef typename Rep::Cuts_map Cuts_map;
    
    //! type of Surface_pair
    typedef typename Rep::Surface_pair Surface_pair;

    //! type of Surface_dcel_handle_map
    typedef typename Rep::Surface_dcel_handle_map Surface_dcel_handle_map;

    //! type of Surface_pair_dcel_handle_map
    typedef typename Rep::Surface_pair_dcel_handle_map
    Surface_pair_dcel_handle_map;
    
    //! pair canonicalizer
    struct Canonicalize_pair {

        //! first in pair in always less than second
        Surface_pair operator()(Surface_pair pair) {
            typename Rep::Surface_less_than less;
            if (!less(pair.first, pair.second)) {
                std::swap(pair.first, pair.second);
            }
            return pair;
        }
    };

private:
    //! type of X_coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;

    //! type of Boundary
    typedef typename Curve_kernel_2::Boundary Boundary;

    //! type of Event info
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    //! type of Polynomial 1
    typedef typename Polynomial_2::NT Polynomial_1;

private:
    //! Standard constructor
    P_dcel_info(const CGAL::Dcel_feature feature, 
                const CGAL::Object& dcel_handle) :
        Base(Rep(feature, dcel_handle)) {
    }

    //!\name Id of cads
    //!@{

    //! returns id of underlying cad
    Rs_id _rs_id() const {
        return this->ptr()->_m_rs_id;
    }
    
    //! sets id of underlying cad
    void _set_rs_id(Rs_id rs_id) const {
        this->ptr()->_m_rs_id = rs_id;
    }

    //!@}

    //!\name Dcel Handles
    //!@{

    //! the dcel-handle 
    CGAL::Object _dcel_handle() const {
        return this->ptr()->_m_dcel_handle;
    }
    
    //! sets the dcel-handle 
    void _set_dcel(const CGAL::Dcel_feature feature, 
                   const CGAL::Object& dcel_handle) const {
        this->ptr()->_m_dcel_feature = feature;
        this->ptr()->_m_dcel_handle = dcel_handle;
    }
    
    //! return (if exiting) a pair of dcel-feature + handle for \c surface
    boost::optional< std::pair< CGAL::Dcel_feature, CGAL::Object > >
    _dcel(const Surface_3& surface) const {
        typedef
            boost::optional< std::pair< CGAL::Dcel_feature, CGAL::Object > >
            Return_type;
                typename Surface_dcel_handle_map::iterator it = 
            this->ptr()->_m_surface_dcel_handle_map.find(surface);
        if (it != this->ptr()->_m_surface_dcel_handle_map.end()) {
            Return_type tmp = it->second;
            return tmp;
        }
        // else
        return boost::none;
    }


    /*!\brief
     * stores for given \c surface the dcel in originating rscad
     */
    void _set_dcel(const Surface_3& surface,
                   const CGAL::Dcel_feature feature, 
                   CGAL::Object obj) const {

        CGAL_assertion_code((
        {
            typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
                Restricted_cad_3;
            switch (feature) {
            case CGAL::FACE: {
                typename Restricted_cad_3::Face_const_iterator fit;
                CGAL_assertion_code(bool check = )
                    CGAL::assign (fit, obj);
                CGAL_assertion(check);
                break;
            }
            case CGAL::EDGE: {
                typename Restricted_cad_3::Halfedge_const_iterator hit;
                CGAL_assertion_code(bool check = )
                    CGAL::assign (hit, obj);
                CGAL_assertion(check);
                break;
            }
            case CGAL::VERTEX: {
                typename Restricted_cad_3::Vertex_const_iterator vit;
                CGAL_assertion_code(bool check = )
                    CGAL::assign (vit, obj);
                CGAL_assertion(check);
                break;
            }
            }
        })
        );
        
        typename Surface_dcel_handle_map::iterator it = 
            this->ptr()->_m_surface_dcel_handle_map.find(surface);
        CGAL_precondition(it == 
                          this->ptr()->_m_surface_dcel_handle_map.end());
        this->ptr()->_m_surface_dcel_handle_map.insert(
                std::make_pair(surface,std::make_pair(feature,obj))
        );
    }
    
    //! return (if exiting) a pair of dcel-feature + handle for \c surfaces
    boost::optional< std::pair< CGAL::Dcel_feature, CGAL::Object > >
    _dcel(const Surface_3& surface1, const Surface_3& surface2) const {
        typedef 
            boost::optional< std::pair< CGAL::Dcel_feature, CGAL::Object > >
            Return_type;
        static Canonicalize_pair canonicalize;
        Surface_pair pair = canonicalize(std::make_pair(surface1, surface2));
        typename Surface_pair_dcel_handle_map::iterator it = 
            this->ptr()->_m_surface_pair_dcel_handle_map.find(pair);
        if (it != this->ptr()->_m_surface_pair_dcel_handle_map.end()) {
            Return_type tmp = it->second;
            return tmp;
        }
        // else
        return boost::none;
    }

    /*!\brief
     * stores for given \c surface the dcel in originating rscad
     */
    void _set_dcel(const Surface_3& surface1,
                   const Surface_3& surface2,
                   const CGAL::Dcel_feature feature, 
                   CGAL::Object obj) const {
        static Canonicalize_pair canonicalize;
        Surface_pair pair = canonicalize(std::make_pair(surface1, surface2));
        typename Surface_pair_dcel_handle_map::iterator it = 
            this->ptr()->_m_surface_pair_dcel_handle_map.find(pair);
        CGAL_precondition(
                it == this->ptr()->_m_surface_pair_dcel_handle_map.end()
        );
        this->ptr()->_m_surface_pair_dcel_handle_map.insert(
                std::make_pair(pair,std::make_pair(feature,obj))
        );
    }
    
    //!@}

    //!\name Adding 
    //!@{

private:
    /*! adds surface 
     */
    // TODO remove with adding init_nk(s1,s2,nk)
    void _add_surface(const Surface_3& surface) const {
        this->ptr()->_m_surfaces.insert(surface);
    }

    
    //! initialized nk-object for surface
    void _init_nk(const Surface_3& surface) const {
        this->ptr()->_m_surfaces.insert(surface);
        
        CGAL::Nk nk;
        
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        if (it == this->ptr()->_m_silhouettes_nk.end()) {
            this->ptr()->_m_silhouettes_nk.insert(
                    it, std::make_pair(surface, nk)
            );
        } 
    }
    
    /*! finalizes internal structures
     */
    void _finalize_nk(const Surface_3& surface) const {
        
        CGAL::Nk nk = _nk(surface);
        
#if !NDEBUG
        std::cout << "Finalizing Nk for " << this->ptr()->_m_dcel_feature 
                  << " with " << nk << " ... " << std::flush;
#endif
        
        if (nk.mult() != -1) {
            if (!this->ptr()->_m_has_silhouette_curve) {
                this->ptr()->_m_has_silhouette_curve = true;
            }
            if (nk.n() == -1) {
                if (this->ptr()->_m_dcel_feature == CGAL::VERTEX) {
                    // TODO be carefull if point is at infinity
                    if (!this->ptr()->_m_has_vertical_line) {
                        this->ptr()->_m_has_vertical_line = true;
                    }
                } else if (this->ptr()->_m_dcel_feature == CGAL::EDGE) {
                    if (!this->ptr()->_m_has_vertical_surface) {
                        this->ptr()->_m_has_vertical_surface = true;
                    }
                }
            } 
        }
#if !NDEBUG
        std::cout << "done." << std::endl;
#endif
    }
    
    /*!\brief
     * Add cut of \c surface1 and \c surface2 with \c multiplicity
     */
    void _add_cut(const Surface_3& surface1, 
                  const Surface_3& surface2, 
                  int multiplicity,
                  bool factors_of_same_curve) const {
        static Canonicalize_pair canonicalize;
        Surface_pair pair = canonicalize(std::make_pair(surface1, surface2));
        typename Cuts_map::iterator it = 
            this->ptr()->_m_cuts.find(pair);
        if (it == this->ptr()->_m_cuts.end()) {
            this->ptr()->_m_cuts.insert(
                    it, std::make_pair(pair, multiplicity)
            );
        } else {
            if (factors_of_same_curve) {
                it->second += multiplicity;
            }
        }
    }
    
    //!@}
    
private:
    //!\name Query members
    //!@{

    /*!\brief
     * returns \c true if carries positive number of stored silhouette curves
     */
    bool _has_silhouette() const {
        return (this->ptr()->_m_has_silhouette_curve);
    }

    /*!\brief
      * returns \c true if carries vertical silhouettes
     */
    bool _has_vertical_surface() const {
        return (this->ptr()->_m_has_vertical_surface);
    }
    
    /*!\brief
     * returns \c true if carries positive number of stored cut curves
     */
    bool _has_cut() const {
        return (!this->ptr()->_m_cuts.empty());
    }
    
    /*!\brief
     * returns \c true if handle is silhouette of \c surface
     */
    bool _has_silhouette(const Surface_3& surface) const {
        if (!this->ptr()->_m_has_silhouette_curve) {
            return false;
        }
        // else
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        if (it == this->ptr()->_m_silhouettes_nk.end()) {
            return false;
        }
        // else
        return it->second.mult() != -1;
    }

    /*! return nk-instance of given surface 
     */
    const CGAL::Nk& _nk(const Surface_3& surface) const {
        
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        CGAL_assertion(it != this->ptr()->_m_silhouettes_nk.end());
        return (it->second);
    } 

    /*!\brief
     * returns \c true of projected curve belongs to vertical part of 
     * \c surface
     */
    // TODO new name with n
    bool _is_surface_vertical(
            const Surface_3& surface
    ) const {
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        CGAL_assertion(it != this->ptr()->_m_silhouettes_nk.end());
        return (it->second.n() == -1);
    }

    /*!\brief
     * returns whether supports vertical line
     */
    // TODO new name with n
    bool _supports_vertical_line(const Surface_3& surface) const {
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        CGAL_assertion(it != this->ptr()->_m_silhouettes_nk.end());
        return (it->second.n() == -1);
    }
    
private:
    /*!\brief
     * returns \c true if carries positive number of stored cut curves
     */
    bool _has_cut(const Surface_3& surface1, 
                  const Surface_3& surface2) const {
        static Canonicalize_pair canonicalize;
        Surface_pair pair = canonicalize(std::make_pair(surface1, surface2));
        typename Cuts_map::iterator it = 
            this->ptr()->_m_cuts.find(pair);
        return (it != this->ptr()->_m_cuts.end());
    }

public:
    
    /*!\brief
     * returns multiplicity of silhouette curve of \c surface
     *
     * \pre Is info for an edge.
     */
    boost::optional< int > multiplicity_of_silhouette(
            const Surface_3& surface
    ) const {
        CGAL_precondition(this->ptr()->_m_dcel_feature == CGAL::EDGE);
        typename Silhouettes_nk_map::iterator it = 
            this->ptr()->_m_silhouettes_nk.find(surface);
        if (it == this->ptr()->_m_silhouettes_nk.end()) {
            return boost::none;
        }
        // else
        return boost::optional< int >(it->second.mult());
    }

    /*!\brief
     * returns multiplicity of cut curve of \c surface1 and 
     * \c surface2
     *
     * \pre Is info for an edg.e
     */
    boost::optional< int > multiplicity_of_cut(
            const Surface_3& surface1,
            const Surface_3& surface2
    ) const {
        CGAL_precondition(this->ptr()->_m_dcel_feature == CGAL::EDGE);
        static Canonicalize_pair canonicalize;
        Surface_pair pair = canonicalize(std::make_pair(surface1, surface2));
        typename Cuts_map::iterator it = 
            this->ptr()->_m_cuts.find(pair);
        if (it == this->ptr()->_m_cuts.end()) {
            return boost::none;
        }
        // else
        return boost::optional< int >(it->second);
    }

    //!@}

    // TOOD make private and accessible only through RS_3
    //!\name Boundary classes
    //!@{

    /*!\brief
     * returns known boundary class for object
     */
    CGAL::Arr_boundary_type boundary_type_in_x() const {
        if (this->ptr()->_m_boundary_type_in_x) {
            return *this->ptr()->_m_boundary_type_in_x;
        }
        // else
        return CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES;
    }

    /*!\brief
     * returns known boundary class for object
     */
    CGAL::Arr_boundary_type boundary_type_in_y() const {
        if (this->ptr()->_m_boundary_type_in_y) {
            return *this->ptr()->_m_boundary_type_in_y;
        }
        // else
        return CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES;
    }

    /*!\brief
     * Sets boundary class for data to given value \c cl
     */
    void set_boundary_type_in_x(CGAL::Arr_boundary_type cl) const {
        this->ptr()->_m_boundary_type_in_x = cl;
    }

    /*!\brief
     * Sets boundary class for data to given value \c cl
     */
    void set_boundary_type_in_y(CGAL::Arr_boundary_type cl) const {
        this->ptr()->_m_boundary_type_in_y = cl;
    }
    
    //!@}
    

private:
    // TASK move the following members to RS_3?; make static? make "global"

    //!\name Stored point
    //!@{
    
    //! returns the sample point for a vertex handle
    template < class VertexHandle >
    Point_2 _sample_point_for_vertex_handle(VertexHandle vh) const {
        CGAL_precondition(!vh->is_at_infinity());

        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;
        
#if 0 && SoX_Z_STACK_Z_TIMERS

        bool running = z_stack_time[0].is_running();
        if (!running) {
            z_stack_time[0].start();
        }
#endif

        Restricted_cad_3 cad = 
            Restricted_cad_3::cad_cache()(
                    vh->data()->ptr()->_m_rs_id
            );
        
        // check existing stacks
        boost::optional< Z_stack > helper;
        
        for (typename std::list< Z_stack >::iterator sit =
                 this->ptr()->_m_z_stacks.begin(); 
             sit != this->ptr()->_m_z_stacks.end();
             sit++) {
            
            if (cad._point_on_dcel_handle(sit->point(), vh)) {
                // if helper is not set or 
                // number of z_cell in current z_stack is greater
                // than number of z_stack cells in current helper
                if (!helper || 
                    (sit->number_of_z_cells() > 
                     helper->number_of_z_cells())) {
                    // store new values.
                    helper = *sit;
                }
            } // TASK implement assert
        }
        
        // to obtain a final candidate
        this->ptr()->_m_z_stacks.clear();
        if (helper) {
            this->ptr()->_m_z_stacks.push_front(*helper);
        }
        
        if (!this->ptr()->_m_point) {
            CGAL_assertion(!this->ptr()->_m_z_stack);
            
            Point_2 pt;
            // either helper gives point
            if (helper) {
                pt = helper->point();
                // _m_z_stacks store also z-stack candidate
            } else {
                // or we compute it 
                pt = vh->point();
            }

            _set_point(pt);
        } else {
            // second case: point is already set
            if (this->ptr()->_m_z_stack) {
                // and also z-stack
                if (helper) {
                    // -> delete candidate
                    this->ptr()->_m_z_stacks.clear();
                }
            } else {
                if (helper) {
                    // reset point as helper is set and z_stack is unset
                    Point_2 pt = helper->point();
                    _set_point(pt);
                }
            }
        }

#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[0].stop();
        }
#endif
        CGAL_postcondition(this->ptr()->_m_point);
        CGAL_postcondition(
                cad._point_on_dcel_handle(*this->ptr()->_m_point, vh)
        );
                       
        // return stored z_stack
        return *(this->ptr()->_m_point);
    }

    //! returns the sample point for a halfedge handle
    template < class HalfedgeHandle >
    Point_2 _sample_point_for_halfedge_handle(HalfedgeHandle heh) const {
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;

#if 0 && SoX_Z_STACK_Z_TIMERS
        bool running = z_stack_time[1].is_running();
        if (!running) {
            z_stack_time[1].start();
        }
#endif
        Restricted_cad_3 cad = 
            Restricted_cad_3::cad_cache()(
                    heh->data()->ptr()->_m_rs_id
            );
        
        // check existing stacks
        boost::optional< Z_stack > helper;
        
        for (typename std::list< Z_stack >::iterator sit =
                 this->ptr()->_m_z_stacks.begin(); 
             sit != this->ptr()->_m_z_stacks.end();
             sit++) {
            if (cad._point_on_dcel_handle(sit->point(), heh)) {
                // if helper is not set or 
                // number of z_cell in current z_stack is greater
                // than number of z_stack cells in current helper
                if (!helper || 
                    (sit->number_of_z_cells() > 
                     helper->number_of_z_cells())) {
                    // store new values.
                    helper = *sit;
                }
            } // TASK implement assert
        }
        
        // to obtain a final candidate
        this->ptr()->_m_z_stacks.clear();
        if (helper) {
            this->ptr()->_m_z_stacks.push_front(*helper);
        }

        
        if (!this->ptr()->_m_point) {
            CGAL_assertion(!this->ptr()->_m_z_stack);
            
            Point_2 pt;
            // either helper gives point
            if (helper) {
                pt = helper->point();
                // _m_z_stacks store also z-stack candidate
            } else {
                // or we compute it
                pt = Restricted_cad_3::_point_in_interior(heh);
            }

            _set_point(pt);
        } else {
            // second case: point is already set
            if (this->ptr()->_m_z_stack) {
                // and also z-stack
                if (helper) {
                    // -> delete candidate
                    this->ptr()->_m_z_stacks.clear();
                }
            } else {
                if (helper) {
                    // reset point as helper is set and z_stack is unset
                    Point_2 pt = helper->point();
                    _set_point(pt);
                }
            }
        }

#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[1].stop();
        }
#endif
        CGAL_postcondition(this->ptr()->_m_point);
        CGAL_postcondition(
                cad._point_on_dcel_handle(*this->ptr()->_m_point, heh)
        );
        // return stored z_stack
        return *(this->ptr()->_m_point);
    }

    //! returns sample point for face handle
    template < class FaceHandle >
    Point_2 _sample_point_for_face_handle(FaceHandle fh) const {
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;
        
#if 0 && SoX_Z_STACK_Z_TIMERS
        bool running = z_stack_time[2].is_running();
        if (!running) {
            z_stack_time[2].start();
        }
#endif

        Restricted_cad_3 cad = 
            Restricted_cad_3::cad_cache()(
                    fh->data()->ptr()->_m_rs_id
            );
        
        // check existing stacks
        boost::optional< Z_stack > helper;
        
        for (typename std::list< Z_stack >::iterator sit =
                 this->ptr()->_m_z_stacks.begin(); 
             sit != this->ptr()->_m_z_stacks.end();
             sit++) {
            
            if (cad._point_on_dcel_handle(sit->point(), fh)) {
                // if helper is not set or 
                // number of z_cell in current z_stack is greater
                // than number of z_stack cells in current helper
                if (!helper || 
                    (sit->number_of_z_cells() > 
                     helper->number_of_z_cells())) {
                    // store new values.
                    helper = *sit;
                }
            } // TASK implement assert
        }
        
        // to obtain a final candidate
        this->ptr()->_m_z_stacks.clear();
        if (helper) {
            this->ptr()->_m_z_stacks.push_front(*helper);
        }

        
        if (!this->ptr()->_m_point) {
            CGAL_assertion(!this->ptr()->_m_z_stack);
            
            Point_2 pt;
            // either helper gives point
            if (helper) {
                pt = helper->point();
                // _m_z_stacks store also z-stack candidate
            } else {
                // or we compute it 
                pt = _rational_point_in_face(fh, cad);
            }
            
            _set_point(pt);
        } else {
            // second case: point is already set
            if (this->ptr()->_m_z_stack) {
                // and also z-stack
                if (helper) {
                    // -> delete candidate
                    this->ptr()->_m_z_stacks.clear();
                }
            } else {
                if (helper) {
                    // reset point as helper is set and z_stack is unset
                    Point_2 pt = helper->point();
                    _set_point(pt);
                }
            }
        }
#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[2].stop();
        }
#endif
        CGAL_postcondition(this->ptr()->_m_point); 
        CGAL_postcondition(
                cad._point_on_dcel_handle(*this->ptr()->_m_point, fh)
        );
        
        // return stored z_stack
        return *(this->ptr()->_m_point);
    }
    
    //!@}
    
private:
    //!\name Stored z-stack
    //!@{

    /*!\brief
     * returns \c true if z-stack is known
     */
    bool _knows_z_stack() const {
        return this->ptr()->_m_z_stack;
    }
    
    /*!\brief
     * compute Z_stack for a given Vertex_handle
     */
    template < class VertexHandle >
    Z_stack _z_stack_for_vertex_handle(VertexHandle vh) const {
        CGAL_precondition(!vh->is_at_infinity());
#if 0 && SoX_Z_STACK_Z_TIMERS
        bool running = z_stack_time[0].is_running();
        if (!running) {
            z_stack_time[0].start();
        }
#endif
        if (!this->ptr()->_m_z_stack) {
            if (!this->ptr()->_m_point) {
                _sample_point_for_vertex_handle(vh);
            }
            
            // store z_stack
            _set_z_stack(CGAL::VERTEX);
        }

        // delete all old z_stacks
        this->ptr()->_m_z_stacks.clear();
      
#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[0].stop();
        }
#endif
  
        // return stored z_stack
        return *(this->ptr()->_m_z_stack);
    }

    /*!\brief
     * compute Z_stack for a given Halfedge_handle
     */
    template < class HalfedgeHandle >
    Z_stack _z_stack_for_halfedge_handle(HalfedgeHandle heh) const {
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;

#if 0 && SoX_Z_STACK_Z_TIMERS
        bool running = z_stack_time[1].is_running();
        if (!running) {
            z_stack_time[1].start();
        }
#endif
        if (!this->ptr()->_m_z_stack) {
            if (!this->ptr()->_m_point) {
                _sample_point_for_halfedge_handle(heh);
            }
            
            // store z_stack
            _set_z_stack(CGAL::EDGE);
        }

        // delete all old z_stacks
        this->ptr()->_m_z_stacks.clear();
        
#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[1].stop();
        }
#endif
        
        // return stored z_stack
        return *(this->ptr()->_m_z_stack);
    }

    /*!\brief
     * compute Z_stack for a given Face_handle
     */
    template < class FaceHandle >
    Z_stack _z_stack_for_face_handle(FaceHandle fh) const {
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;
        
#if 0 && SoX_Z_STACK_Z_TIMERS
        bool running = z_stack_time[2].is_running();
        if (!running) {
            z_stack_time[2].start();
        }
#endif
        if (!this->ptr()->_m_z_stack) {
            if (!this->ptr()->_m_point) {
                _sample_point_for_face_handle(fh);
            }

            // store z_stack
            _set_z_stack(CGAL::FACE);
        }

        // delete all old z_stacks
        this->ptr()->_m_z_stacks.clear();
        
#if 0 && SoX_Z_STACK_Z_TIMERS
        if (!running) {
            z_stack_time[2].stop();
        }
#endif        
        
        // return stored z_stack
        return *(this->ptr()->_m_z_stack);
    }

    //! returns for given \c surface and \cad the store z_stack in the sil-cad
    std::pair< Z_stack, CGAL::Dcel_feature > _z_stack_of_surface(
            const Surface_3& surface
    ) const {
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
            Restricted_cad_3;
        
        typename Surface_dcel_handle_map::iterator it = 
            this->ptr()->_m_surface_dcel_handle_map.find(surface);
        CGAL_precondition(it != this->ptr()->_m_surface_dcel_handle_map.end());
        
        if (it->second.first == CGAL::FACE) {
            typename Restricted_cad_3::Face_const_iterator fit;
            CGAL_assertion_code(bool check = )
                CGAL::assign (fit, it->second.second);
            CGAL_assertion(check);
            return std::make_pair(
                    fit->data()->_z_stack_for_face_handle(fit),
                    CGAL::FACE
            );
        } else if (it->second.first == CGAL::EDGE) {
            typename Restricted_cad_3::Halfedge_const_iterator hit;
            CGAL_assertion_code(bool check = )
                CGAL::assign (hit, it->second.second);
            CGAL_assertion(check);
            return std::make_pair(
                    hit->data()->_z_stack_for_halfedge_handle(hit),
                    CGAL::EDGE
            );
        } 
        // else
        CGAL_assertion(it->second.first == CGAL::VERTEX);
        typename Restricted_cad_3::Vertex_const_iterator vit;
        CGAL_assertion_code(bool check =)
            CGAL::assign (vit, it->second.second);
        CGAL_assertion(check);
        return std::make_pair(vit->data()->_z_stack_for_vertex_handle(vit),
                              CGAL::VERTEX);
    }
    
    /*!\brief
     * scans a connected component of the given ccb indicated by \c start
     */
    template < class HalfedgeHandle > 
    static
    std::pair< HalfedgeHandle, HalfedgeHandle >
    _scan_ccb(HalfedgeHandle start) {
        
        HalfedgeHandle invalid_he;
        HalfedgeHandle he;
        HalfedgeHandle he_v;
        
        HalfedgeHandle curr = start;
        do {
            ++curr;
            if (!curr->is_fictitious()) {
                //std::cout << "curve: " << curr->curve() << std::endl;
                //std::cout << "magickey" << std::endl;
                if (curr->curve().is_vertical()) {
                    if (he_v == invalid_he) {
                        he_v = curr;
                    }
                } else {
                    if (he == invalid_he) {
                        he = curr;
                    }
                }
            }
            if (he != invalid_he) {
                break;
            }
        } while (curr != start);
        
        // finally return 
        return std::make_pair(he, he_v);
    }
    
    /*!\brief
     * computes rational point in given face
     */
    template < class FaceHandle, class RestrictedCad_3>
    Point_2 _rational_point_in_face(FaceHandle fh, 
                                    RestrictedCad_3& cad) const {
        CGAL_precondition(fh != FaceHandle());
        
        typedef RestrictedCad_3 Restricted_cad_3;

        typedef typename Restricted_cad_3::Ccb_halfedge_const_circulator 
            Ccb_halfedge_const_circulator; 
        typedef typename Restricted_cad_3::Outer_ccb_const_iterator 
            Outer_ccb_const_iterator; 
        typedef typename Restricted_cad_3::Inner_ccb_const_iterator 
            Inner_ccb_const_iterator; 
        
        Ccb_halfedge_const_circulator invalid_he;
        
        std::pair< Ccb_halfedge_const_circulator, 
            Ccb_halfedge_const_circulator > ocbpair, icbpair;
        
        Ccb_halfedge_const_circulator he, he_v;
        
        // select halfedge on boundary of face 
        // avoid vertical halfedges
        for (Outer_ccb_const_iterator ocb = fh->outer_ccbs_begin();
             ocb != fh->outer_ccbs_end(); 
             ocb++) {
            
            ocbpair = _scan_ccb(*ocb);
            
            if (ocbpair.first != invalid_he && he == invalid_he) {
                he = ocbpair.first;
                break;
            }
            if (ocbpair.second != invalid_he && he_v == invalid_he) {
                he_v = ocbpair.second;
            }
        }
        
        if (he == invalid_he) {
            // search on inner ccbs
            for (Inner_ccb_const_iterator icb = fh->inner_ccbs_begin();
                 icb != fh->inner_ccbs_end(); 
                 icb++) {
                
                icbpair = _scan_ccb(*icb);
                
                if (icbpair.first != invalid_he && he == invalid_he) {
                    he = icbpair.first;
                    break;
                }
                if (icbpair.second != invalid_he && he_v == invalid_he) {
                    he_v = icbpair.second;
                }
            }
            
        }
        
        if (he == invalid_he) {
            // try whether vertical one has been found
            if (he_v != invalid_he) {
                he = he_v;
            }
        }
        // he is either invalid, or we've found a normal edge or a vertical one

        // helping variables
        int type = -1;
        Boundary xr(0), y0(0);
        X_coordinate_1 x0(xr);
        int count = 0;
        int arc = -1;
        Status_line_1 sl;

        // type init phase
        if (he == invalid_he) {
            type = 0;
        } else if (!he->curve().is_vertical()) {
            type = 1;
            Boundary xr1 = he->curve().boundary_in_x_range_interior();
            x0 = X_coordinate_1(xr1);
            
            arc = he->curve().arcno();
            sl = he->curve().curve().status_line_at_exact_x(x0);
        } else {
            CGAL_assertion(he->curve().is_vertical());
            type = 2;
        }
        

        typename Curve_kernel_2::Lower_boundary_y_2 lower_boundary_y = 
            Arrangement_traits_2::instance().kernel().
            lower_boundary_y_2_object();
        typename Curve_kernel_2::Upper_boundary_y_2 upper_boundary_y = 
            Arrangement_traits_2::instance().kernel().
            upper_boundary_y_2_object();
        typename Curve_kernel_2::Refine_y_2 refine_y = 
            Arrangement_traits_2::instance().kernel().
            refine_y_2_object();
        
        while (true) {
            // desired object
            Point_2 pt;
            
            Point_2 ptv;

            // compute value phase
            switch (type) {
            case 0: {
                x0 = X_coordinate_1(xr);
                y0 = xr * xr;
                break;
            }
            case 1:
                // we found a normal halfedge
                // (needs direction)
                y0 = (he->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                      upper_boundary_y(sl.algebraic_real_2(arc)) :
                      lower_boundary_y(sl.algebraic_real_2(arc))
                );
                break;
            case 2: {
                ++count;
                // in this case we could only find a 
                // vertical halfedge for the face 
                // that has to go from y=-oo to y =+00
                CGAL_assertion(he->curve().is_vertical());
                X_coordinate_1 x1 = he->curve().x();
                Boundary xt2 = (he->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                      x1.low() : x1.high());
                if (x1.is_rational()) {
                    xr = xt2 + 
                        (he->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                         Boundary(-1) : Boundary(1))/
                        CGAL::ipower(Boundary(2),count);
                } else {
                    // TODO use ak
                    xr = x1.rational_between(X_coordinate_1(xt2));
                }
                x0 = X_coordinate_1(xr);

                ptv = Restricted_cad_3::_point_in_interior(he);
                
                break;
            }
            default:
                CGAL_error_msg("Only 0,1,2 are allowed as types");
            }

            if (type == 2) {
                pt = Point_2(X_coordinate_1(x0),ptv.curve(),0);
            } else {
                // construct point phase
                pt = Restricted_cad_3::_construct_point_with_rational_y(
                        x0, y0
                );
            }
            
            // check phase
            // test whether point lies in the face with point location
            if (cad._point_on_dcel_handle(pt, fh)) {
                // found face matches given one, 
                // i.e., constructed point lies within correct face
                return pt;
            }

            // refine phase
            switch (type) {
            case 0:
                xr += Boundary(1);
                break;
            case 1:
                CGAL_assertion(arc != -1);
                // refine y-coordinate
                refine_y(sl.algebraic_real_2(arc));
                break;
            case 2:
                CGAL_assertion(arc == -1);
                // refine x-coordinate
                x0.strong_refine(xr);
                break;
            }
        }
    }
    
private:

    /*!\brief
     * sets point for z-stack construction
     */
    void _set_point(const Point_2& pt) const {
        CGAL_precondition(!this->ptr()->_m_point);
        this->ptr()->_m_point = pt;
    }

    /*!\brief
     * sets z_stack for object
     */
    void _set_z_stack(CGAL::Dcel_feature feature) const {
        CGAL_precondition(this->ptr()->_m_point);        
        CGAL_precondition(!this->ptr()->_m_z_stack);        

        Point_2 pt = *this->ptr()->_m_point;
        
#if !NDEBUG
        std::cout << "feature: " << feature << std::endl;
        std::cout << "Computing z-stack for " << pt << " ... " << std::flush;
#endif

        // tasks: 
        // - check whether all surface in stored z_stack exist in _m_surfaces
        if (static_cast< int >(this->ptr()->_m_z_stacks.size()) > 1) {
            CGAL_assertion(
                    static_cast< int >(this->ptr()->_m_z_stacks.size()) == 2
            );
            CGAL_assertion(pt == this->ptr()->_m_z_stacks.begin()->point());
            CGAL_assertion(
                    pt == (++this->ptr()->_m_z_stacks.begin())->point()
            );
            Z_stack z_stack(
                    this->ptr()->_m_z_stacks.begin()->_merge(
                            *(this->ptr()->_m_z_stacks.begin()++),
                            this
                    ),
                    feature, this
            );
            
            this->ptr()->_m_z_stack = z_stack;
        } else {

            CGAL_assertion_code(
                    if (static_cast< int >(this->ptr()->_m_z_stacks.size()) 
                        == 1) {
                        CGAL_assertion(
                                pt == this->ptr()->_m_z_stacks.begin()->point()
                        );
                    }
            );
            const Z_stack& z_stack = 
                (static_cast< int >(this->ptr()->_m_z_stacks.size()) == 1 ?
                 Z_stack(*this->ptr()->_m_z_stacks.begin(), feature, this) :
                 Z_stack(feature,this,pt));
            
            typename Rep::Surface_container erase;

            for (typename Rep::Surface_container::iterator it =
                     this->ptr()->_m_surfaces.begin();
                 it != this->ptr()->_m_surfaces.end(); it++) {
                
                if (this->ptr()->_m_dcel_feature == CGAL::EDGE && 
                    _is_surface_vertical(*it)) {
                    continue;
                }
                
                bool empty_isolator;
                if (!z_stack._knows_isolator(*it, empty_isolator)) {
                    // construct traits
                    typename Surface_z_at_xy_isolator_traits::
                        Construct_isolator
                        construct_isolator; 
                    // TODO (construct_isolator_object())
                    // construct isolator
                    boost::optional< 
                        std::pair< CGAL::Dcel_feature, CGAL::Object > >
                        dcelinfo = _dcel(*it);
                    CGAL_assertion(dcelinfo);
                    Z_at_xy_isolator isolator = 
                        construct_isolator(
                                *it, pt, _nk(*it), dcelinfo->first
                        );
                    // add roots in z_stack sequence
                    z_stack._add_surface(*it, isolator);
                    if (isolator.number_of_real_roots() == 0) {
                        // if isolator has no root remove 
                        // surface from _m_surfaces
                        erase.insert(*it);
                    }
                } else {
                    CGAL_assertion_code(
                            Z_at_xy_isolator isolator = z_stack._isolator(*it);
                    );
                    CGAL_assertion(isolator.traits().point() == pt);
                }
            }
            
            for (typename Rep::Surface_container::iterator it =
                     erase.begin();
                 it != erase.end(); it++) {
                this->ptr()->_m_surfaces.erase(*it);
            }
            
            this->ptr()->_m_z_stack = z_stack;
        }
        
        this->ptr()->_m_z_stacks.clear();
        
        CGAL_postcondition(this->ptr()->_m_z_stack);
        
        // check that each involved surface shows up in z_stack
        CGAL_postcondition_code((
        {
            std::list< int > levels;
            for (typename Rep::Surface_container::iterator it =
                     this->ptr()->_m_surfaces.begin();
                 it != this->ptr()->_m_surfaces.end(); it++){ 
                if (!_is_surface_vertical(*it)) {
                    levels.clear();
                    typename 
                        Surface_z_at_xy_isolator_traits::Construct_isolator
                        construct_isolator; 
                    // TODO (construct_isolator_object())
                    boost::optional< 
                        std::pair< CGAL::Dcel_feature, CGAL::Object > >
                        dcelinfo = _dcel(*it);
                    CGAL_assertion(dcelinfo);
                    Z_at_xy_isolator isolator = 
                        construct_isolator(
                                *it, pt, _nk(*it), dcelinfo->first
                        );
                    if (isolator.number_of_real_roots() > 0) {
                        this->ptr()->_m_z_stack->levels_of(
                                *it, std::back_inserter(levels)
                        );
                        CGAL_postcondition(!levels.empty());
                    }
                }
            }
        }
        ));

#if !NDEBUG
        std::cout << "done." << std::endl;
#endif
    }

    //!@}
    
    //!\name Overlay
    //!@{
    
    Self _overlay_with(const Self& info, bool factors_of_same_curve) const {
        Self tmp(info);
        tmp.copy_on_write();
        
        // set rs_id to zero
        tmp.ptr()->_m_rs_id = 0;

        // empty object
        CGAL::Object obj;
        tmp.ptr()->_m_dcel_handle = obj;
        
        // involved surfaces
        tmp.ptr()->_m_surfaces.insert(this->ptr()->_m_surfaces.begin(),
                                      this->ptr()->_m_surfaces.end());
        
        // silhouettes and ...
        for (typename Silhouettes_nk_map::iterator it = 
                 this->ptr()->_m_silhouettes_nk.begin(); 
             it != this->ptr()->_m_silhouettes_nk.end(); it++) {
            tmp.ptr()->_m_silhouettes_nk.insert(*it);
        }
        
        if (this->ptr()->_m_has_silhouette_curve) {
            tmp.ptr()->_m_has_silhouette_curve = true;
        }

        // cuts .. note that add_cut may also be additive
        for (typename Cuts_map::iterator it = 
                 this->ptr()->_m_cuts.begin(); 
             it != this->ptr()->_m_cuts.end(); it++) {
            tmp._add_cut(it->first.first, it->first.second, 
                                  it->second, factors_of_same_curve);
        }
        
        // overlay of boundary classes in x
        if (tmp.ptr()->_m_boundary_type_in_x) {
            if (this->ptr()->_m_boundary_type_in_x) {
                CGAL_assertion(tmp.boundary_type_in_x() == 
                               this->boundary_type_in_x());
            } 
        } else {
            if (this->ptr()->_m_boundary_type_in_x) {
                tmp.ptr()->_m_boundary_type_in_x = 
                    this->boundary_type_in_x();
            }
        }
        // overlay of boundary classes in y
        if (tmp.ptr()->_m_boundary_type_in_y) {
            if (this->ptr()->_m_boundary_type_in_y) {
                CGAL_assertion(tmp.boundary_type_in_y() == 
                               this->boundary_type_in_y());
            } 
        } else {
            if (this->ptr()->_m_boundary_type_in_y) {
                tmp.ptr()->_m_boundary_type_in_y = 
                    this->boundary_type_in_y();
            }
        }
        
        // new construction may benefit from original z_stacks 
        // store to _m_z_stacks 
        if (tmp.ptr()->_m_z_stack) {
            tmp.ptr()->_m_z_stacks.push_front(*(tmp.ptr()->_m_z_stack));
        }
        if (this->ptr()->_m_z_stack) {
            tmp.ptr()->_m_z_stacks.push_front(*(this->ptr()->_m_z_stack));
        }
        // this list will be deleted when asking for z_stack of *this
        
        // but first, we invalidate current z_stack
        tmp.ptr()->_m_point = boost::none;
        tmp.ptr()->_m_z_stack = boost::none;
        
        if (!factors_of_same_curve) {
            // surface dcel handles
            for (typename Surface_dcel_handle_map::iterator it = 
                     this->ptr()->_m_surface_dcel_handle_map.begin(); 
                 it != this->ptr()->_m_surface_dcel_handle_map.end(); it++) {
                tmp._set_dcel(it->first, it->second.first, it->second.second);
            }
            
            // surface pair dcel handles
            for (typename Surface_pair_dcel_handle_map::iterator it = 
                     this->ptr()->_m_surface_pair_dcel_handle_map.begin(); 
                 it != this->ptr()->_m_surface_pair_dcel_handle_map.end(); 
                 it++) {
                tmp._set_dcel(
                        it->first.first, it->first.second, 
                        it->second.first, it->second.second
                );
            }
        }
        
        return tmp;
    }
    
    //!@}

public:
    //!\name IO members
    //!@{

    /*!\brief
     * prints pretty-formated version of P_dcel_info
     */
    void pretty_print(std::ostream& os) const {
        os << "PDcelInfo(id: " << this->id() 
           << ", rsid: " <<  this->_rs_id() 
           << ", feat=" << this->ptr()->_m_dcel_feature << "): " 
           << std::endl;
        if (this->ptr()->_m_dcel_feature == CGAL::VERTEX) {
            typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
                Restricted_cad_3;
            typename Restricted_cad_3::Vertex_const_iterator vit;
            CGAL_assertion_code(bool check = )
                CGAL::assign(vit, this->ptr()->_m_dcel_handle);
            CGAL_assertion(check);
            os << "Point: " << vit->point() << std::endl;
        }
        if (this->ptr()->_m_dcel_feature == CGAL::EDGE) {
            typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits> 
                Restricted_cad_3;
            typename Restricted_cad_3::Halfedge_const_iterator hit;
            CGAL_assertion_code(bool check = )
                CGAL::assign(hit, this->ptr()->_m_dcel_handle);
            CGAL_assertion(check);
            os << "Curve: " << hit->curve() << std::endl;
        }
        os << "Surfaces ";
        for (typename Rep::Surface_container::iterator it =
                 this->ptr()->_m_surfaces.begin();
             it != this->ptr()->_m_surfaces.end(); it++){ 
            os << "<" << it->id() << ">" << std::flush;
        }
        if (this->ptr()->_m_dcel_feature != CGAL::FACE) {
            os << std::endl;
        }
        if (this->_has_silhouette()) {
            os << "Silhouettes ";
            for (typename Silhouettes_nk_map::iterator it = 
                     this->ptr()->_m_silhouettes_nk.begin(); 
                 it != this->ptr()->_m_silhouettes_nk.end(); it++) {
                if (it->second.mult() != -1) {
                    os << "<" << it->first.id() << "," << it->second;
                    if (it->second.n() == 0) {
                        if (this->ptr()->_m_dcel_feature == CGAL::EDGE) {
                            os << ",V";
                        } else if (this->ptr()->_m_dcel_feature == 
                                   CGAL::VERTEX) {
                            os << ",VL";
                        }
                    }
                    os << "> ";
                }
            }
            os << " ";
        }
        if (this->_has_cut()) {
            if (this->_has_silhouette()) {
                os << " ";
            }
            os << "Cuts ";
            for (typename Cuts_map::iterator it = 
                     this->ptr()->_m_cuts.begin(); 
                 it != this->ptr()->_m_cuts.end(); it++) {
                os << "<(" << it->first.first.id() << "," 
                   << it->first.second.id() << ")," << it->second 
                    //<< " (" << it->first.first.f() << "," 
                    //<< it->first.second.f() << ")" 
                   << "> ";
            }
        }
        if (this->ptr()->_m_boundary_type_in_x ||
            this->ptr()->_m_boundary_type_in_y) {
            if (this->ptr()->_m_boundary_type_in_x) {
                os << "BdryX " << *this->ptr()->_m_boundary_type_in_x;
            }
            if (this->ptr()->_m_boundary_type_in_y) {
                if (this->ptr()->_m_boundary_type_in_x) {
                    os << " ";
                }
                os << "BdryY " << *this->ptr()->_m_boundary_type_in_y;
            }
            if (this->ptr()->_m_dcel_feature != CGAL::FACE) {
                os << std::endl;
            }
        }
        if (this->ptr()->_m_z_stack) {
            os << *this->ptr()->_m_z_stack;
        }
    }
    
    //!@}
    
private:
    
    // friends
    // for z_stacks
    friend class CGAL::Restricted_cad_3< Surface_z_at_xy_isolator_traits>;
    friend class CGAL::Restricted_cad_3_accessor< Restricted_cad_3< 
    Surface_z_at_xy_isolator_traits> >;
    
    // for add members
    friend class CGAL::Create_restricted_cad_3< Surface_z_at_xy_isolator_traits >;
    // for _set_rs_id
    friend class CGAL::Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits >;

    // for _overlay_with
    friend class CGAL::Arr_p_dcel_info_overlay_traits< 
    Restricted_cad_3< Surface_z_at_xy_isolator_traits> >;
    
    // for _z_stack_of_surface
    friend class CGAL::Z_stack< Surface_z_at_xy_isolator_traits, Self>;
    
    // data members
};


/*!\relates P_dcel_info
 * \brief
 * outputs P_dcel_instance object to stream \c os
 */
template < class SurfaceZAtXyIsolatorTraits >
std::ostream& operator<<(
        std::ostream& os, 
        const P_dcel_info< SurfaceZAtXyIsolatorTraits >& data
) {
    data.pretty_print(os);
    return os;
}
    

CGAL_END_NAMESPACE

#endif // SoX_GAPS_P_DCEL_INFO_H
// EOF
