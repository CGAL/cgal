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

#ifndef CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ACCESSOR_H
#define CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ACCESSOR_H 1

/*!\file include/CGAL/Arrangement_2l/Restricted_cad_3_accessor.h
 * \brief Contains class template Restricted_cad_3_accessor
 */

#include <CGAL/config.h>

#include <CGAL/Arrangement_2l/Restricted_cad_3_enums.h>

CGAL_BEGIN_NAMESPACE

/*!\brief
 * Accessor class for private types and members of RestrictedCad_3
 */
template < class RestrictedCad_3 >
class Restricted_cad_3_accessor {
public:
    
    //! this instance's first template parameter
    typedef RestrictedCad_3 Restricted_cad_3;
    
    //! type of instantiated class
    typedef Restricted_cad_3_accessor< Restricted_cad_3 > Self;

    //! standard constructor
    Restricted_cad_3_accessor(const Restricted_cad_3& rscad) :
        _m_rscad(rscad) {
    }
    
private:
    typedef Restricted_cad_3 RSC_3;
    
public:
    //! rep class
    typedef typename RSC_3::Rep Rep;
    
    //! base class
    typedef typename RSC_3::Base Base;

    //! type of Surface
    typedef typename RSC_3::Surface_3 Surface_3;
    
    //! Z_Stack type
    typedef typename RSC_3::Z_stack Z_stack;

    //! Z_at_xy_isolator type
    typedef typename RSC_3::Z_at_xy_isolator Z_at_xy_isolator;
    
    // from arrangement
    typedef typename RSC_3::Size Size;
   
    typedef typename RSC_3::Curve_analysis_2 Curve_analysis_2;
    
    typedef typename RSC_3::Point_2 Point_2;
    typedef typename RSC_3::X_monotone_curve_2  X_monotone_curve_2;
    
    typedef typename RSC_3::Vertex_const_iterator Vertex_const_iterator;

    typedef typename RSC_3::Edge_const_iterator Edge_const_iterator;

    typedef typename RSC_3::Halfedge_const_iterator Halfedge_const_iterator;
    
    typedef typename RSC_3::Face_const_iterator Face_const_iterator;

    typedef typename RSC_3::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
    
    typedef typename RSC_3::Outer_ccb_const_iterator Outer_ccb_const_iterator;
    
    typedef typename RSC_3::Inner_ccb_const_iterator Inner_ccb_const_iterator;
    
    typedef typename RSC_3::Vertex_const_handle Vertex_const_handle;
    typedef typename RSC_3::Halfedge_const_handle Halfedge_const_handle;
    typedef typename RSC_3::Edge_const_handle Edge_const_handle;
    typedef typename RSC_3::Face_const_handle Face_const_handle;

    typedef typename RSC_3::Vertex_handle Vertex_handle;
    typedef typename RSC_3::Halfedge_handle Halfedge_handle;
    typedef typename RSC_3::Edge_handle Edge_handle;
    typedef typename RSC_3::Face_handle Face_handle;

    typedef typename RSC_3::Point_location Point_location;

    //! type of X_coordinate
    typedef typename RSC_3::X_coordinate_1 X_coordinate_1;
    
    //! type of Boundary
    typedef typename RSC_3::Boundary Boundary;
    
    //!\name Global
    //!@{

    //! Returns a pointer to the representation of the cad
    inline const Rep* rep() const {
        return _m_rscad.ptr();
    }

    //! make id unique and clear the rep
    inline void renew() {
        _m_rscad._renew();
    }
    
    //! construct cad for curve
    static
    Restricted_cad_3 construct_for_curve(const Curve_analysis_2& curve) {
        return Restricted_cad_3::_construct_for_curve(curve);
    }

    //! initialize data objects
    inline void init(const Surface_3& surface) {
        _m_rscad._init(surface);
    }
    
    //! set final data
    inline void finalize(const Surface_3& surface) {
        _m_rscad._finalize(surface);
    }
    
    //!@}
    
    //!\name Arrangement methods 
    //!@{

    /*! insert a range of curves */
    template < class InputIterator >
    inline
    void insert_non_intersecting_curves(
            InputIterator begin, InputIterator end
    ) {
        _m_rscad._insert_non_intersecting_curves(begin, end);
    }
    
    /*! insert a range of curves */
    template < class CurveInputIterator, class PointInputIterator >
    inline
    void insert_empty(
            CurveInputIterator cbegin, CurveInputIterator cend,
            PointInputIterator pbegin, PointInputIterator pend
            
    ) {
        _m_rscad._insert_empty(cbegin, cend, pbegin, pend);
    }
    

    /*! inserts an isolated point */
    inline
    void insert_point(Point_2 pt) {
        _m_rscad._insert_point(pt);
    }

    /*! remove edge from arrangement */
    void remove_edge(Halfedge_handle he) {
        _m_rscad._remove_edge(he);
    }

    /*! remove edge from arrangement */
    void remove_vertex(Vertex_handle vh) {
        _m_rscad._remove_vertex(vh);
    }
    
    /*! overlays two instances */
    inline
    void overlay(Restricted_cad_3 rscad1, 
                 Restricted_cad_3 rscad2, bool factors_of_same_curve) {
        _m_rscad._overlay(rscad1, rscad2, factors_of_same_curve);
    }

    /*! overlays two instances */
    inline
    void overlay(Restricted_cad_3 rscad1, 
                 Restricted_cad_3 rscad2, 
                 Surface_3 surface, CGAL::Nk::Value_type type) {
        _m_rscad._overlay(rscad1, rscad2, surface, type);
    }
    
    //!@}

    //!\name Geometry and the Dcel
    //!@{
    
    /*!\brief
     * use point location to check whether \c pt is part of given 
     * \c handle 
     */
    template < class DcelConstHandle >
    inline
    bool point_on_dcel_handle(const Point_2& pt, 
                              DcelConstHandle handle) const {
        return _m_rscad._point_on_dcel_handle(pt, handle);
    }

    /*!\brief 
     * constructs a rational point at coordinate \c x0, \c y0
     */
    static
    Point_2 construct_point_with_rational_y(X_coordinate_1 x0, Boundary y0) {
        return RSC_3::_construct_point_with_rational_y(x0, y0);
    }
    
    //! construct point, with rational x (or y) in interior of \c heh's curve
    static
    Point_2 point_in_interior(const Halfedge_const_handle& heh) {
        return RSC_3::_point_in_interior(heh);
    }

    //! construct point, with rational x (or y) in interior of \c cv
    static
    Point_2 point_in_interior(const X_monotone_curve_2& cv) {
        return RSC_3::_point_in_interior(cv);
    }
    
    //!@}

    //!\name Counting
    //!@{
    /*! Check whether the arrangement is empty. */

    /*! Get the number of arrangement halfedges (the result is always even). */
    inline
    Size number_of_halfedges () const
    {
        return (_m_rscad.number_of_halfedges());
    }
    
    //!@}

    //\name Traversal functions for the arrangement halfedges.
    //!@{

    /*! Get a const iterator for the first halfedge in the arrangement. */
    inline
    Halfedge_const_iterator halfedges_begin() const
    { 
        return (_m_rscad.halfedges_begin());
    }
  
    /*! Get a past-the-end const iterator for the arrangement halfedges. */
    inline
    Halfedge_const_iterator halfedges_end() const
    {
        return (_m_rscad.halfedges_end());
    }
    //!@}
    
    //!\name Casting away constness for handle types.
    //!@{
    
    //! converts const vertex handle to non-const vertex handle
    inline
    typename Rep::Vertex_handle non_const_handle (Vertex_const_handle vh)
    {
        return (_m_rscad.non_const_handle(vh));
    }
    
    //! converts const halfedge handle to non-const halfedge handle
    inline
    typename Rep::Halfedge_handle non_const_handle (Halfedge_const_handle hh)
    {
        return (_m_rscad.non_const_handle(hh));
    }

    //! converts const face handle to non-const face handle
    inline
    typename Rep::Face_handle non_const_handle (Face_const_handle fh)
    {
        return (_m_rscad.non_const_handle(fh));
    }
    //!@}

    //!\name Nk
    //!@{

    
    //! set value for vertex handle
    inline void set_nk_value(const Vertex_const_handle& vh, 
                             const Surface_3& surface,
                             CGAL::Nk::Value_type type, 
                             int value) const {
        _m_rscad._set_nk_value(vh, surface, type, value);
    }

    //! set value for halfedge handle
    inline void set_nk_value(const Halfedge_const_handle& heh, 
                             const Surface_3& surface,
                             CGAL::Nk::Value_type type, 
                             int value) const {
        _m_rscad._set_nk_value(heh, surface, type, value);
    }

    //! set value for edge handle
    inline void set_nk_value(const Edge_const_handle& eh, 
                             const Surface_3& surface,
                             CGAL::Nk::Value_type type, 
                             int value) const {
        _m_rscad._set_nk_value(eh, surface, type, value);
    }

    //! set value for face handle
    inline void set_nk_value(const Face_const_handle& fh, 
                             const Surface_3& surface,
                             CGAL::Nk::Value_type type, 
                             int value) const {
        _m_rscad._set_nk_value(fh, surface, type, value);
    }
    
    //! nk for halfedge handle
    const CGAL::Nk& nk(const Halfedge_const_handle& heh, 
                      const Surface_3& surface) const {
        return _m_rscad.nk(heh, surface);
    }

    //!@}


    //!\name Silhouette/Cut
    //!@{
    
    //! returns whether silhouette at given halfedge handle exists
    bool has_silhouette(const Halfedge_const_handle& heh) const {
        return (_m_rscad.has_silhouette(heh));
    }
    
    //! returns whether cut at given halfedge handle exists
    bool has_cut(const Halfedge_const_handle& heh) const {
        return (_m_rscad.has_cut(heh));
    }
    
    //! returns whether silhouette of \c surface at \c heh exists
    bool has_silhouette(const Halfedge_const_handle& heh, 
                        const Surface_3& surface) const {
        return (_m_rscad.has_silhouette(heh,surface));
    }
    
    //! returns whether cut the two surface at \c heh exisst
    bool has_cut(const Halfedge_const_handle& heh, 
                 const Surface_3& surface1,
                 const Surface_3& surface2) const {
        return (_m_rscad.has_cut(heh,surface1, surface2));
    }
    
    //!@}

    //!@{
    
    //! returns sample point for halfedge handle
    inline
    Point_2 sample_point(const Halfedge_const_handle& heh) const {
        return _m_rscad.sample_point(heh);
    }

    //!@}

    //!\name Isolators
    //!@{
    
     //! returns isolator (if existing) for given surface at given handle
    boost::optional< Z_at_xy_isolator> 
    isolator(const Halfedge_const_handle& heh,
             const Surface_3& surface) const {
        return _m_rscad_isolator(heh, surface);
    }
    
    //!@}

    //!\name Z_Stacks
    //!@{
    
    //! returns whether z_stack is known
    bool knows_z_stack(const Vertex_const_handle& vh) const {
        return vh->data()->_knows_z_stack();
    }

    //! returns whether z_stack is known
    bool knows_z_stack(const Halfedge_const_handle& heh) const {
        return heh->data()->_knows_z_stack();
    }

    //! returns whether z_stack is known
    bool knows_z_stack(const Edge_const_handle& eh) const {
        return eh->data()->_knows_z_stack();
    }

      //! returns whether z_stack is known
    bool knows_z_stack(const Face_const_handle& fh) const {
        return fh->data()->_knows_z_stack();
    }

    //! returns z_stack for given edge handle
    inline
    Z_stack z_stack(const Halfedge_const_handle& heh) const {
        return _m_rscad.z_stack(heh);
    }

    //! returns z_stack of silhouette-curve of \c surface at halfedge handle
    inline
    std::pair< Z_stack, CGAL::Dcel_feature > 
    z_stack(const Halfedge_const_handle& heh, 
          const Surface_3& surface) const {
        return _m_rscad.z_stack(heh, surface);
    }

    //!@}
  
private:
    Restricted_cad_3 _m_rscad;
    
};

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ACCESSOR_H
// EOF
