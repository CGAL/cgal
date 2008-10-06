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
// File          : include/SoX/GAPS/Surface_pair_3.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef SoX_GAPS_SURFACE_PAIR_3_H
#define SoX_GAPS_SURFACE_PAIR_3_H 1

/*!\file SoX/GAPS/Surface_pair_3.h
 * \brief 
 * definition of \c Surface_pair_3<>
 */

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arrangement_2l/macros.h>
#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_functors.h>

namespace SoX {

namespace Intern {

template < class SurfaceZAtXyIsolatorTraits >
class Surface_pair_3_rep {

public:
    //! this instance's template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! this instance itself
    typedef Surface_pair_3_rep< Surface_z_at_xy_isolator_traits > Self;

    //! type of restricted cad 
    typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits > 
    Restricted_cad_3;
    
    //! type of creator
    typedef Create_restricted_cad_3< Surface_z_at_xy_isolator_traits > Creator;

    //! type of creator
    typedef Overlay_restricted_cad_3< Surface_z_at_xy_isolator_traits > 
    Overlayer;
    
public:
    //!\name Constructors
    //!@{
    
    Surface_pair_3_rep(const Surface_3& surface1, const Surface_3& surface2) :
        _m_surface1(surface1),
        _m_surface2(surface2) {
    }
    
    //!@}
    
// TODO make again private:
    
    //! surface1
    mutable Surface_3 _m_surface1;

    //! surface2
    mutable Surface_3 _m_surface2;

    //! restricted cad of surface1
    mutable boost::optional< Restricted_cad_3 > _m_silhouette1;

    //! restricted cad of surface2
    mutable boost::optional< Restricted_cad_3 > _m_silhouette2;

    //! restricted cad of surface2
    mutable boost::optional< Restricted_cad_3 > _m_silhouettes;
    
    //! restricted cad of intersection
    mutable boost::optional< Restricted_cad_3 > _m_cut;

    //! restricted cad of first silhouette with x
    mutable boost::optional< Restricted_cad_3 > _m_silhouette1cut;

    //! restricted cad of second silhouette with x
    mutable boost::optional< Restricted_cad_3 > _m_silhouette2cut;
    
    //! restricted cad
    mutable boost::optional< Restricted_cad_3 > _m_silhouettescut;
    
    //! friends
    friend class Surface_pair_3< Surface_z_at_xy_isolator_traits, Self > ;
};

} // namespace Intern

template < 
class SurfaceZAtXyIsolatorTraits, 
class Rep_ = Intern::Surface_pair_3_rep < SurfaceZAtXyIsolatorTraits >
>
class Surface_pair_3 : 
        public ::CGAL::Handle_with_policy< Rep_ > {
    
public:
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;
    
    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of Base
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //! type of restricted cad
    typedef typename Rep::Restricted_cad_3 Restricted_cad_3;

    //! type of creator
    typedef typename Rep::Creator Creator;

    //! type of overlayer
    typedef typename Rep::Overlayer Overlayer;

    //! type of Z_stack
    typedef typename Restricted_cad_3::Z_stack Z_stack;

    //!\name Constructors
    //!@{
    
    Surface_pair_3(const Surface_3& surface1, const Surface_3& surface2) :
        Base(Rep(surface1, surface2)) {
    }
    
    //!@}
    
    //!\name Access members
    //!@{

    //! the first surface
    Surface_3 surface1() const {
        return this->ptr()->_m_surface1;
    }

    //! the second surface
    Surface_3 surface2() const {
        return this->ptr()->_m_surface2;
    }
    
    //!@}
    
     //!\name Cads
    //!@{

    //! rscad of first silhouette
    Restricted_cad_3 silhouette1() const {
        if (!this->ptr()->_m_silhouette1) {  
            Creator creator; // TODO static??
            // create rs_cad_3 for surface1
            this->ptr()->_m_silhouette1 = creator(this->ptr()->_m_surface1);
        }
        return *this->ptr()->_m_silhouette1;
    }

    //! rscad of second silhouette
    Restricted_cad_3 silhouette2() const {
        if (!this->ptr()->_m_silhouette2) {  
            Creator creator; // TODO static??
            // create rs_cad_3 for surface2
            this->ptr()->_m_silhouette2 = creator(this->ptr()->_m_surface2);
        }
        return *this->ptr()->_m_silhouette2;
    }

protected:
    //! rscad of both silhouettes
    Restricted_cad_3 _silhouettes() const {
        if (!this->ptr()->_m_silhouettes) {  
            Overlayer overlay; // TODO static?
            this->ptr()->_m_silhouettes = 
                overlay(silhouette1(), silhouette2());
        }
        return *this->ptr()->_m_silhouettes;
    }
    
public:
    //! rscad of cut
    Restricted_cad_3 cut() const {
        if (!this->ptr()->_m_cut) {  
            // create rs_cad_3 for intersection of surface1 and surface2
            Creator creator; // TODO static??
            this->ptr()->_m_cut = 
                creator(this->ptr()->_m_surface1, this->ptr()->_m_surface2);
        }
        return *this->ptr()->_m_cut;
    }

protected:    
    //! rscad of first silhouette with cut
    Restricted_cad_3 _silhouette1_cut() const {
        if (!this->ptr()->_m_silhouette1cut) {  
            Overlayer overlay; // TODO static?
            this->ptr()->_m_silhouette1cut = overlay(silhouette1(), cut());
        }
        return *this->ptr()->_m_silhouette1cut;
    }
    
    //! rscad of first silhouette with cut
    Restricted_cad_3 _silhouette2_cut() const {
        if (!this->ptr()->_m_silhouette2cut) {  
            Overlayer overlay; // TODO static?
            this->ptr()->_m_silhouette2cut = overlay(silhouette2(), cut());
        }
        return *this->ptr()->_m_silhouette2cut;
    }

public:
    //! rscad of first silhouette with cut
    Restricted_cad_3 silhouettes_cut() const {
        if (!this->ptr()->_m_silhouettescut) {  
            Overlayer overlay; // TODO static?
            this->ptr()->_m_silhouettescut = overlay(_silhouettes(), cut());
        }
        return *this->ptr()->_m_silhouettescut;
    }

    //!@}

    //!\name Curves
    //!@{

    /*!\brief
     * returns boundary curves and points for a given \c surface in pair
     */
    template < class CurveOutputIterator, class PointOutputIterator >
    void silhouette_objects(
            const Surface_3& surface, 
            CurveOutputIterator coi, PointOutputIterator poi) const {
        
        const Restricted_cad_3& rsc = (surface == this->surface1() ?
                                       this->silhouette1() :
                                       this->silhouette2());
        
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 rsc.vertices_begin();
             vit != rsc.vertices_end(); vit++) {
            if (vit->is_isolated()) {
                *poi++ = vit->point();
            }
        }
        for (typename Restricted_cad_3::Edge_const_iterator eit =
                 rsc.edges_begin();
             eit != rsc.edges_end(); eit++) {
            *coi++ = eit->curve();
        }
    }

    /*!\brief
     * returns projected intersection curves and points of the pair
     */
    template < class CurveOutputIterator, class PointOutputIterator >
    void cut_objects(
            CurveOutputIterator coi, PointOutputIterator poi) const {
        
        const Restricted_cad_3& rsc = this->cut();
        
        for (typename Restricted_cad_3::Vertex_const_iterator vit =
                 rsc.vertices_begin();
             vit != rsc.vertices_end(); vit++) {
            if (vit->is_isolated()) {
                *poi++ = vit->point();
            }
        }
        for (typename Restricted_cad_3::Edge_const_iterator eit =
                 rsc.edges_begin();
             eit != rsc.edges_end(); eit++) {
            *coi++ = eit->curve();
        }
    }
    
    //!@}

    //!\name Z_stack for compare_xyz
    //!@{
    
    //! returns z-stack of silhouttes_cut() for \c point
    Z_stack z_stack_for(const Point_2& point) const {
        return silhouettes_cut().z_stack_for(point);
    }
    
    //!@}
};

} // namespace SoX

#endif // SoX_GAPS_SURFACE_PAIR_3_H
// EOF
