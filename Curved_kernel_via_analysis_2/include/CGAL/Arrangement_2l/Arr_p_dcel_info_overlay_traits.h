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

/*!\file include/CGAL/Arrangement_2l/Arr_p_dcel_info_overlay_traits.h
 * \brief definition of Arr_p_dcel_info_overlay_traits class template
 */

#ifndef CGAL_ARRANGEMENT_2l_ARR_P_DCEL_INFO_OVERLAY_TRAITS_H
#define CGAL_ARRANGEMENT_2l_ARR_P_DCEL_INFO_OVERLAY_TRAITS_H 1

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Arr_default_overlay_traits.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_accessor.h>

CGAL_BEGIN_NAMESPACE

/*!\brief
 * Model of CGAL::ArrangementOverlayTraits_2 for Arrangements
 * attached with CGAL::P_dcel_info<> data.
 */
template < class RestrictedCad_3 >
class Arr_p_dcel_info_overlay_traits : public 
CGAL::Arr_default_overlay_traits< typename RestrictedCad_3::Rep > {
    
public:
    //! this instance's first template parameter
    typedef RestrictedCad_3 Restricted_cad_3;
    
    //! type of arrangement
    typedef typename Restricted_cad_3::Rep Arrangement_2;

    //! type of surface
    typedef typename Restricted_cad_3::Surface_3 Surface_3;
    
    //! base type
    //typedef CGAL::Arr_default_overlay_traits< Arrangement_2 > Base;

    typedef Restricted_cad_3 RS_3;

    //! type of accesor
    typedef CGAL::Restricted_cad_3_accessor< Restricted_cad_3 > Accessor;

    //! the class itself
    typedef Arr_p_dcel_info_overlay_traits< Arrangement_2 > Self;
    
    typedef typename RS_3::Vertex_const_handle Vertex_const_handle_A;
    typedef typename RS_3::Vertex_const_handle Vertex_const_handle_B;
    typedef typename RS_3::Vertex_const_handle Vertex_const_handle_R;

    typedef typename RS_3::Halfedge_const_handle Halfedge_const_handle_A;
    typedef typename RS_3::Halfedge_const_handle Halfedge_const_handle_B;
    typedef typename RS_3::Halfedge_const_handle Halfedge_const_handle_R;

    typedef typename RS_3::Face_const_handle Face_const_handle_A;
    typedef typename RS_3::Face_const_handle Face_const_handle_B;
    typedef typename RS_3::Face_const_handle Face_const_handle_R;

    
    /*!\brief
     * Constructs a new overlay traits
     */
    Arr_p_dcel_info_overlay_traits(const Surface_3& surface, 
                                   const Restricted_cad_3& cad,
                                   CGAL::Nk::Value_type type) :
        _m_surface_mode(true),
        _m_surface(surface),
        _m_cad(cad),
        _m_acc(cad),
        _m_type(type),
        _m_factors_of_same_curve(false)
    {
        
    }


    /*!\brief
     * Constructs a new overlay traits where \c factors_of_same_curve defines
     * whether the two parts originate from two factors of the same curve
     */
    Arr_p_dcel_info_overlay_traits(const Restricted_cad_3& cad, 
                                   bool factors_of_same_curve) :
        _m_surface_mode(false),
        _m_cad(cad),
        _m_acc(cad),
        _m_factors_of_same_curve(factors_of_same_curve) {
    }

private:

    // for vertex/face
    inline void merge(CGAL::Nk nk1, CGAL::Nk nk2, const CGAL::Nk& nk, 
                      const CGAL::Dcel_feature& feature) const {

        if (this->_m_type == CGAL::Nk::MULT) {
            // mults!
            if (feature == CGAL::VERTEX) {
                nk.set_mult(0);
            } else if (feature == CGAL::EDGE) {
                nk.set_mult(
                        (nk1.mult() > 0 ? nk1.mult() : 0) +
                        (nk2.mult() > 0 ? nk2.mult() : 0)
                );
                CGAL_assertion(nk.mult() > 0);
            }
        } else { 
            // copy already set values!
            if (nk1.mult() >= 0) {
                if (feature == CGAL::VERTEX) {
                    nk.set_mult(0);
                } else {
                    CGAL_assertion(feature == CGAL::EDGE);
                    nk.set_mult(nk1.mult());
                }
            }
            if (nk1.n() >= -1) {
                nk.set_n(nk1.n());
            }
            if (nk1.k() >= -1) {
                nk.set_k(nk1.k());
            }
        }
    }

public:    


    /*!
     * Create a vertex v that corresponds to the coinciding vertices v1 and v2.
     */
    virtual void create_vertex (Vertex_const_handle_A v1,
                                Vertex_const_handle_B v2,
                                Vertex_const_handle_R v)
    {
        CGAL_precondition(v1->point() == v->point());
        CGAL_precondition(v2->point() == v->point());
        CGAL_assertion(v1->data());
        CGAL_assertion(v2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!v->data());
            _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
            
            v1->data()->_nk(_m_surface);
            v2->data()->_nk(_m_surface);
            v->data()->_nk(_m_surface);

            merge(v1->data()->_nk(_m_surface), 
                  v2->data()->_nk(_m_surface), 
                  v->data()->_nk(_m_surface), CGAL::VERTEX);
            CGAL_postcondition(v->data());
            
            v->data()->_nk(_m_surface).set_feature1(CGAL::VERTEX);
            v->data()->_nk(_m_surface).set_feature2(CGAL::VERTEX);
        } else {
            _m_acc.non_const_handle(v)->data() = 
                v1->data()->_overlay_with(*v2->data(), 
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create a vertex v that mathces v1, which lies of the edge e2.
     */
    virtual void create_vertex (Vertex_const_handle_A v1,
                                Halfedge_const_handle_B e2,
                                Vertex_const_handle_R v)
    {
        CGAL_precondition(v1->point() == v->point());
        CGAL_assertion(v1->data());
        CGAL_assertion(e2->data());
        
        if (_m_surface_mode) {
            CGAL_precondition(!v->data());
            _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
            merge(v1->data()->_nk(_m_surface), 
                  e2->data()->_nk(_m_surface), 
                  v->data()->_nk(_m_surface), CGAL::VERTEX);
            CGAL_postcondition(v->data());

            v->data()->_nk(_m_surface).set_feature1(CGAL::VERTEX);
            v->data()->_nk(_m_surface).set_feature2(CGAL::EDGE);
        } else {
            _m_acc.non_const_handle(v)->data() = 
                v1->data()->_overlay_with(*e2->data(), 
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create a vertex v that mathces v1, contained in the face f2.
     */
    virtual void create_vertex (Vertex_const_handle_A v1,
                                Face_const_handle_B f2,
                                Vertex_const_handle_R v)
    {
        CGAL_precondition(v1->point() == v->point());
        CGAL_assertion(v1->data());
        CGAL_assertion(f2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!v->data());
            _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
            merge(v1->data()->_nk(_m_surface), 
                  f2->data()->_nk(_m_surface), 
                  v->data()->_nk(_m_surface), CGAL::VERTEX);
            CGAL_postcondition(v->data());

            v->data()->_nk(_m_surface).set_feature1(CGAL::VERTEX);
            v->data()->_nk(_m_surface).set_feature2(CGAL::FACE);

        } else {
            _m_acc.non_const_handle(v)->data() = 
                v1->data()->_overlay_with(*f2->data(),
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create a vertex v that mathces v2, which lies of the edge e1.
     */
    virtual void create_vertex (Halfedge_const_handle_A e1,
                                Vertex_const_handle_B v2,
                                Vertex_const_handle_R v)
    {
        CGAL_precondition(v2->point() == v->point());
        CGAL_assertion(e1->data());
        CGAL_assertion(v2->data());
        
        if (_m_surface_mode) {
              CGAL_precondition(!v->data());
              _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
              merge(e1->data()->_nk(_m_surface), 
                    v2->data()->_nk(_m_surface), 
                    v->data()->_nk(_m_surface), CGAL::VERTEX);
              CGAL_postcondition(v->data());

              v->data()->_nk(_m_surface).set_feature1(CGAL::EDGE);
              v->data()->_nk(_m_surface).set_feature2(CGAL::VERTEX);
        } else {
            _m_acc.non_const_handle(v)->data() = 
                v2->data()->_overlay_with(*e1->data(),
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create a vertex v that mathces v2, contained in the face f1.
     */
    virtual void create_vertex (Face_const_handle_A f1,
                                Vertex_const_handle_B v2,
                                Vertex_const_handle_R v)
    {
        CGAL_precondition(v2->point() == v->point());
        CGAL_assertion(f1->data());
        CGAL_assertion(v2->data());

        if (_m_surface_mode) {
              CGAL_precondition(!v->data());
              _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
              merge(f1->data()->_nk(_m_surface), 
                    v2->data()->_nk(_m_surface), 
                    v->data()->_nk(_m_surface), CGAL::VERTEX);
              CGAL_postcondition(v->data());

              v->data()->_nk(_m_surface).set_feature1(CGAL::FACE);
              v->data()->_nk(_m_surface).set_feature2(CGAL::VERTEX);
        } else {
            _m_acc.non_const_handle(v)->data() = 
                v2->data()->_overlay_with(*f1->data(),
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create a vertex v that mathces the intersection of the edges e1 and e2.
     */
    virtual void create_vertex (Halfedge_const_handle_A e1,
                                Halfedge_const_handle_B e2,
                                Vertex_const_handle_R v)
    {
        CGAL_assertion(e1->data());
        CGAL_assertion(e2->data());
        
        if (_m_surface_mode) {
              CGAL_precondition(!v->data());
              _m_cad._init_feature(_m_surface, v, CGAL::VERTEX);
              merge(e1->data()->_nk(_m_surface), 
                    e2->data()->_nk(_m_surface), 
                    v->data()->_nk(_m_surface), CGAL::VERTEX);
              CGAL_postcondition(v->data());

              v->data()->_nk(_m_surface).set_feature1(CGAL::EDGE);
              v->data()->_nk(_m_surface).set_feature2(CGAL::EDGE);
        } else {
            _m_acc.non_const_handle(v)->data() = 
                e1->data()->_overlay_with(*e2->data(),
                                          _m_factors_of_same_curve);
        }
    }
    
    /*!
     * Create an edge e that matches the overlap between e1 and e2.
     */
    virtual void create_edge (Halfedge_const_handle_A e1,
                              Halfedge_const_handle_B e2,
                              Halfedge_const_handle_R e)
    {
        CGAL_assertion(e1->data());
        CGAL_assertion(e2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!e->data());
            _m_cad._init_feature(_m_surface, e, CGAL::EDGE);
            merge(e1->data()->_nk(_m_surface), 
                  e2->data()->_nk(_m_surface), 
                  e->data()->_nk(_m_surface), CGAL::EDGE);
            CGAL_postcondition(e->data());

            e->data()->_nk(_m_surface).set_feature1(CGAL::EDGE);
            e->data()->_nk(_m_surface).set_feature2(CGAL::EDGE);
        } else {
            _m_acc.non_const_handle(e)->data() = 
                e1->data()->_overlay_with(*e2->data(),
                                          _m_factors_of_same_curve);
        }
        CGAL_assertion(!e->twin()->data());
        _m_acc.non_const_handle(e->twin())->data() = *e->data();
    }
    
    /*!
     * Create an edge e that matches the edge e1, contained in the face f2.
     */
    virtual void create_edge (Halfedge_const_handle_A e1,
                              Face_const_handle_B f2,
                              Halfedge_const_handle_R e)
    {
        CGAL_assertion(e1->data());
        CGAL_assertion(f2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!e->data());
            _m_cad._init_feature(_m_surface, e, CGAL::EDGE);
            merge(e1->data()->_nk(_m_surface), 
                  f2->data()->_nk(_m_surface), 
                  e->data()->_nk(_m_surface), CGAL::EDGE);
            CGAL_postcondition(e->data());

            e->data()->_nk(_m_surface).set_feature1(CGAL::EDGE);
            e->data()->_nk(_m_surface).set_feature2(CGAL::FACE);
        } else {
            _m_acc.non_const_handle(e)->data() = 
                e1->data()->_overlay_with(*f2->data(),
                                          _m_factors_of_same_curve);
        }
        CGAL_assertion(!e->twin()->data());
        _m_acc.non_const_handle(e->twin())->data() = *e->data();
    }
    
    /*!
     * Create an edge e that matches the edge e2, contained in the face f1.
     */
    virtual void create_edge (Face_const_handle_A f1,
                              Halfedge_const_handle_B e2,
                              Halfedge_const_handle_R e)
    {
        CGAL_assertion(f1->data());
        CGAL_assertion(e2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!e->data());
            _m_cad._init_feature(_m_surface, e, CGAL::EDGE);
            merge(f1->data()->_nk(_m_surface), 
                  e2->data()->_nk(_m_surface), 
                  e->data()->_nk(_m_surface), CGAL::EDGE);
            CGAL_postcondition(e->data());

            e->data()->_nk(_m_surface).set_feature1(CGAL::FACE);
            e->data()->_nk(_m_surface).set_feature2(CGAL::EDGE);
        } else {
            _m_acc.non_const_handle(e)->data() = 
                e2->data()->_overlay_with(*f1->data(),
                                            _m_factors_of_same_curve);
        }
        CGAL_assertion(!e->twin()->data());
        _m_acc.non_const_handle(e->twin())->data() = *e->data();
    }
    
    /*!
     * Create a face f that matches the overlapping region between f1 and f2.
     */
    virtual void create_face (Face_const_handle_A f1,
                              Face_const_handle_B f2,
                              Face_const_handle_R f)
    {
        CGAL_assertion(f1->data());
        CGAL_assertion(f2->data());

        if (_m_surface_mode) {
            CGAL_precondition(!f->data());
            _m_cad._init_feature(_m_surface, f, CGAL::FACE);
            CGAL_assertion(f1->data()->_nk(_m_surface).mult() == -1);
            CGAL_assertion(f2->data()->_nk(_m_surface).mult() == -1);
            merge(f1->data()->_nk(_m_surface), 
                  f2->data()->_nk(_m_surface), 
                  f->data()->_nk(_m_surface), CGAL::FACE);
            CGAL_postcondition(f->data());

            f->data()->_nk(_m_surface).set_feature1(CGAL::FACE);
            f->data()->_nk(_m_surface).set_feature2(CGAL::FACE);
        } else {
            _m_acc.non_const_handle(f)->data() = 
                f1->data()->_overlay_with(*f2->data(),
                                          _m_factors_of_same_curve);
        }
    }

private:
    // members
    //! surface mode?
    bool _m_surface_mode;
    
    //! surface
    Surface_3 _m_surface;

    //! cad
    Restricted_cad_3 _m_cad;

    //! accessor
    Accessor _m_acc;

    //! overlay type
    CGAL::Nk::Value_type _m_type;
    
    //! factors of same curve
    bool _m_factors_of_same_curve;
};

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_ARR_P_DCEL_INFO_OVERLAY_TRAITS_H
// EOF
