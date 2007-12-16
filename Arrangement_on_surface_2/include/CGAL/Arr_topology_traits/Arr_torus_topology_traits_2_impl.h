// Copyright (c) 2006,2007 Max-Planck-Institute Saarbruecken (Germany).
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
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_TORUS_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_TORUS_TOPOLOGY_TRAITS_2_IMPL_H

/*! \file
 * Member-function definitions for the
 * Arr_torus_topology_traits_2<GeomTraits> class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Default constructor.
//
template <class GeomTraits, class Dcel_>
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
Arr_torus_topology_traits_2() :
    _m_own_traits (true), 
    _m_f_top(NULL)
{
    // status: correct
    _m_traits = new Traits_adaptor_2;

    _m_identification_WE = Identification_WE(Point_2_less_WE(_m_traits));
    _m_identification_NS = Identification_NS(Point_2_less_NS(_m_traits));
}

//-----------------------------------------------------------------------------
// Constructor with a geometry-traits class.
//
template <class GeomTraits, class Dcel_>
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
Arr_torus_topology_traits_2 (Geometry_traits_2 *tr) : 
    _m_own_traits(false),  
    _m_f_top(NULL)
{
    // status: correct
    _m_traits = static_cast<Traits_adaptor_2*>(tr);

    _m_identification_WE = Identification_WE(Point_2_less_WE(_m_traits));
    _m_identification_NS = Identification_NS(Point_2_less_NS(_m_traits));
}

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <class GeomTraits, class Dcel_>
void Arr_torus_topology_traits_2<GeomTraits, Dcel_>::assign (const Self& other)
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 assign"  << std::endl;

    // Assign the class.
    // Clear the current DCEL and duplicate the other DCEL.
    _m_dcel.delete_all();
    _m_dcel.assign(other._m_dcel);
    
    // Take care of the traits object.
    if (_m_own_traits && _m_traits != NULL) {
        delete _m_traits;
    }
 
    if (other._m_own_traits) {
        _m_traits = new Traits_adaptor_2;
    } else {
        _m_traits = other._m_traits;
    }
    _m_own_traits = other._m_own_traits;
 
    // Update the special properties of the topology traits.
    dcel_updated();
}

//-----------------------------------------------------------------------------
// Make the necessary updates after the DCEL structure have been updated.
//
template <class GeomTraits, class Dcel_>
void Arr_torus_topology_traits_2<GeomTraits, Dcel_>::dcel_updated ()
{
    // status: missing location of f_top

    // Go over the DCEL vertices and locate all points with boundary condition
    typename Dcel::Vertex_iterator       vit;
    Arr_parameter_space                        bx, by;

    for (vit = this->_m_dcel.vertices_begin();
         vit != this->_m_dcel.vertices_end(); ++vit) {
        // First check whether the vertex has a boundary condition in x.
        bx = vit->parameter_space_in_x();
        if (bx != CGAL::ARR_INTERIOR) {
            
            std::pair< typename Identification_WE::iterator, bool > res =
                _m_identification_WE.insert (std::make_pair(vit->point(),
                                                            &(*vit)));
            CGAL_assertion(! res.second);
            
            _m_vertices_on_identification_WE[&(*vit)] = res.first;
        }
        

        // First check whether the vertex has a boundary condition in y.
        by = vit->parameter_space_in_y();
        if (by != CGAL::ARR_INTERIOR) {
            
            std::pair< typename Identification_NS::iterator, bool > res =
                _m_identification_NS.insert (std::make_pair(vit->point(),
                                                            &(*vit)));
            CGAL_assertion(! res.second);
            
            _m_vertices_on_identification_NS[&(*vit)] = res.first;
        }
    }
    
    // Go over the DCEL faces and locate the top face, which is the only
    // face with no outer CCB.
    typename Dcel::Face_iterator         fit;
    
#if 0 // TODO
    _m_f_top = NULL;
    if (this->_m_dcel.number_of_faces() == 1) {
        _m_f_top = this->_m_dcel->faces_begin();
    } else {
        for (fit = this->_m_dcel.faces_begin();
             fit != this->_m_dcel.faces_end(); ++fit) {
            
            if (fit->number_of_outer_ccbs() == 0) {
                // includes the case that it touches the pole!
                CGAL_assertion (_m_f_top == NULL);
                
                _m_f_top = &(*fit);
                break;
            } else {
                if (fit->number_of_outer_ccbs() == 2) {
                    // search lowest/leftmost ident points
                    
                    // find halfedges
                    Halfedge *e1 = *(fit->outer_ccbs_begin());
                    Halfedge *e2 = *(++fit->outer_ccbs_begin());
                    
                    // collect data of perimetric paths
                    CGAL::Sign sign1 = _sign_of_paths(e1, e1);
                    
                    CGAL::Sign sign2 = _sign_of_paths(e2, e2);
                    
                    bool check_lowest = false;
                    
                    if (counters1.first % 2 != 0) {
                        CGAL_assertion(counters2.first % 2 != 0);
                        CGAL_assertion(counters1.second % 2 == 0);
                        CGAL_assertion(counters2.second % 2 == 0);
                        
                        if (less_we(bottommost1->point(), 
                                    bottommost2->point())) {
                            if (bottommost_crossing1 == BEFORE_TO_AFTER) {
                                _m_f_top = &(*fit);
                                break;
                            } 
                        } else {
                            CGAL_assertion(
                                    less_we(bottommost2->point(), 
                                            bottommost1->point())
                            );
                            if (crossing2 == BEFORE_TO_AFTER) {
                                _m_f_top = &(*fit);
                                break;
                            }
                        }
                        
                    } else if (counters1.second % 2 != 0) {
                        CGAL_assertion(counters2.second % 2 != 0);
                        CGAL_assertion(counters1.first % 2 == 0);
                        CGAL_assertion(counters2.first % 2 == 0);
                        
                        if (less_ns(leftmost1->point(), 
                                    leftmost2->point())) {
                            if (crossing1 == AFTER_TO_BEFORE) {
                                _m_f_top = &(*fit);
                                break;
                            } 
                        } else {
                            CGAL_assertion(
                                    less_ns(leftmost2->point(), 
                                            leftmost1->point())
                            );
                            if (crossing2 == AFTER_TO_BEFORE) {
                                _m_f_top = &(*fit);
                                break;
                            }
                        }
                        
                    } else {
                        CGAL_assertion(counters1.first % 2 == 0 && 
                                       counters1.second % 2 == 0 &&
                                       counters2.first % 2 == 0 && 
                                       counters2.second % 2 == 0);
                        
                        if (less_we(bottommost1->point(), 
                                    bottommost2->point())) {
                            if (crossing1 == AFTER_TO_BEFORE) {
                                _m_f_top = &(*fit);
                                break;
                            } 
                        } else {
                            CGAL_assertion(
                                    less_we(bottommost2->point(), 
                                            bottommost1->point())
                            );
                            if (crossing2 == AFTER_TO_BEFORE) {
                                _m_f_top = &(*fit);
                                break;
                            }
                        }
                    }
                }
                // cannot have a single outer ccb
            }
        }
    }
    CGAL_assertion (_m_f_top != NULL);
#endif
    
    return;
}

//-----------------------------------------------------------------------------
// Initialize an empty DCEL structure.
//
template <class GeomTraits, class Dcel_>
void Arr_torus_topology_traits_2<GeomTraits, Dcel_>::init_dcel ()
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 init_dcel"  
    //          << std::endl;

    // Clear the current DCEL.
    this->_m_dcel.delete_all();

    // create the face
    this->_m_f_top = this->_m_dcel.new_face();
    
    // bounded
    this->_m_f_top->set_unbounded (false);   

    // set not fictious
    this->_m_f_top->set_fictitious (false);
    
    // identifications
    this->_m_identification_WE.clear();
    this->_m_identification_NS.clear();
    
    this->_m_vertices_on_identification_WE.clear();
    this->_m_vertices_on_identification_NS.clear();

}

//-----------------------------------------------------------------------------
// Compare the relative y-position of the given point and the given edge
//
template <class GeomTraits, class Dcel_>
Comparison_result
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
compare_y_at_x (const Point_2& p, const Halfedge* he) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 compare_y_at_x"
    //          << std::endl;
    
    // all edges are valid, therefore just compare p to its associated curve.
    return (this->_m_traits->compare_y_at_x_2_object() (p, he->curve()));
}

//-----------------------------------------------------------------------------
// Check if the given vertex is associated with the given curve end.
//
template <class GeomTraits, class Dcel_>
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
are_equal (const Vertex *v,
           const X_monotone_curve_2& cv, Arr_curve_end ind,
           Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 are_equal"  
    //          << std::endl;
    
    CGAL_precondition(_valid(ps_x, ps_y));
    
    // check wether the two concrete points are equal
    if (v->parameter_space_in_x() != CGAL::ARR_INTERIOR) {
        if (ps_x == CGAL::ARR_INTERIOR) {
            //std::cout << "false2" << std::endl;
            return false;
        }
        // else curve-end lies on NS
        bool res = 
            (this->_m_traits->compare_y_on_identification_2_object() (
                    v->point(),
                    (ind == CGAL::ARR_MIN_END ?
                     this->_m_traits->construct_min_vertex_2_object()(cv) :
                     this->_m_traits->construct_max_vertex_2_object()(cv))) 
             == CGAL::EQUAL
            );
        //std::cout << "res1: " << res << std::endl;
        return res;
    }
    CGAL_assertion(v->parameter_space_in_y() != CGAL::ARR_INTERIOR);
    if (ps_y == CGAL::ARR_INTERIOR) {
        return false;
        //std::cout << "false3" << std::endl;
    }
    // else curve-end lies on WE
    bool res = 
        (this->_m_traits->compare_x_on_identification_2_object() (
                v->point(),
                (ind == CGAL::ARR_MIN_END ?
                 this->_m_traits->construct_min_vertex_2_object()(cv) :
                 this->_m_traits->construct_max_vertex_2_object()(cv))) 
         == EQUAL
        );
    //std::cout << "res2: " << res << std::endl;
    return res;
}

//-----------------------------------------------------------------------------
// Given a curve end with boundary conditions and a face that contains the
// interior of the curve, find a place for a boundary vertex that will
// represent the curve end along the face boundary.
//
template <class GeomTraits, class Dcel_>
CGAL::Object
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
place_boundary_vertex (Face *f,
                       const X_monotone_curve_2& cv, Arr_curve_end ind,
                       Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 place_boundary_vertex"  
    //          << std::endl;

    CGAL_precondition(_valid(ps_x, ps_y));

    // this topology return either an empty object or a DCEL vertex,
    // but never a fictious edge!!!
    
    Vertex *v = NULL;
    
    const Point_2& key = 
        (ind == CGAL::ARR_MIN_END ?
         this->_m_traits->construct_min_vertex_2_object()(cv) :
         this->_m_traits->construct_max_vertex_2_object()(cv));
    
    if (ps_x != CGAL::ARR_INTERIOR) {
        // locate curve-end (here a concrete point) 
        // in local structure of for points on identification_WE
        v = vertex_WE(key);
    } else {
        CGAL_assertion(ps_y != CGAL::ARR_INTERIOR);
        // locate curve-end (here a concrete point) 
        // in local structure of for points on 
        // identification_NS
        v = vertex_NS(key);
    }
    // if there is no vertex found, return empty object
    if (v == NULL) {
        //std::cout << "no vertex found" << std::endl;
        return CGAL::Object();
    }
    
    // else we return the vertex we have located.
    CGAL_assertion(v->parameter_space_in_x() == ps_x && 
                   v->parameter_space_in_y() == ps_y);

    CGAL_assertion(!v->has_null_point());
    return (CGAL::make_object (v));
}

//-----------------------------------------------------------------------------
// Locate the predecessor halfedge for the given curve around a given
// vertex with boundary conditions.
//
template <class GeomTraits, class Dcel_>
typename Arr_torus_topology_traits_2<GeomTraits,Dcel_>::Halfedge* 
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
locate_around_boundary_vertex (Vertex *v,
                               const X_monotone_curve_2& cv, Arr_curve_end ind,
                               Arr_parameter_space ps_x,
                               Arr_parameter_space ps_y) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 locate_around_boundary_vertex"  
    //          << std::endl;
    
    CGAL_precondition(_valid(ps_x, ps_y));
    
    // If the vertex is isolated, there is no predecssor halfedge.
    if (v->is_isolated()) {
        //std::cout << "no edge" << std::endl;
        return NULL;
    }
    
    // Get the first incident halfedge around v and the next halfedge.
    Halfedge * first = v->halfedge();
    Halfedge * curr = first;
    CGAL_assertion(curr != NULL);
    Halfedge * next = curr->next()->opposite();
    
    // If is only one halfedge incident to v, return this halfedge as xc's
    // predecessor:
    if (curr == next) {
        //std::cout << "single edge" << std::endl;
        return curr;
    }
    
    // Otherwise, we traverse the halfedges around v until we find the pair
    // of adjacent halfedges between which we should insert xc.
    typename Traits_adaptor_2::Is_between_cw_2 is_between_cw =
        _m_traits->is_between_cw_2_object();
    bool eq_curr, eq_next;
    
#if 0
    std::cout << "??????????????????????????????????????????" << std::endl;
    std::cout << "search: " << std::endl;

    std::cout << "curr: " << curr->curve() << std::endl;
    std::cout << "dir: " << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << (next->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    
    std::cout << "******************************************" << std::endl;
#endif
    
    while (!is_between_cw(cv, (ind == CGAL::ARR_MIN_END),
                          curr->curve(), 
                          (curr->direction() == CGAL::ARR_RIGHT_TO_LEFT),
                          next->curve(), 
                          (next->direction() == CGAL::ARR_RIGHT_TO_LEFT),
                          v->point(), eq_curr, eq_next))
    {
        // The curve must not be equal to one of the curves 
        // already incident to v.
        CGAL_assertion(!eq_curr && !eq_next);
        
#if 0
        std::cout << "circle: " << std::endl;
        std::cout << "??????????????????????????????????????????" << std::endl;
        std::cout << "search: " << std::endl;
        
        std::cout << "curr: " << curr->curve() << std::endl;
        std::cout << "dir: " << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                                 "L2R" : "R2L") << std::endl;
        std::cout << "next: " << next->curve() << std::endl;
        std::cout << "dir: " << (next->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                                 "L2R" : "R2L") << std::endl;
        
        std::cout << "******************************************" << std::endl;
#endif
        
        // Move to the next pair of incident halfedges.
        curr = next;
        next = curr->next()->opposite();
        
        // Make sure we have not completed a full traversal around v without
        // locating a place for the new curve xc.
        CGAL_assertion (curr != first);
    }

#if 0
    std::cout << "after circle: " << std::endl;
    
    std::cout << "??????????????????????????????????????????" << std::endl;
    std::cout << "search: " << std::endl;

    std::cout << "curr: " << curr->curve() << std::endl;
    std::cout << "dir: " << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << (next->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    
    std::cout << "******************************************" << std::endl;
#endif

    return curr;
}

//-----------------------------------------------------------------------------
// Notifies on the creation of a boundary vertex
//
template <class GeomTraits, class Dcel_>
void 
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
notify_on_boundary_vertex_creation (Vertex *v,
                                    const X_monotone_curve_2& cv,
                                    Arr_curve_end ind,
                                    Arr_parameter_space ps_x,
                                    Arr_parameter_space ps_y) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2::" 
    //          << "notify_on_boundary_vertex_creation"
    //           << std::endl;       

    CGAL_precondition(_valid(ps_x, ps_y));
    
    CGAL_assertion(v->parameter_space_in_x() == ps_x);
    CGAL_assertion(v->parameter_space_in_y() == ps_y);

    CGAL_assertion(!v->has_null_point());

    const Point_2& key = 
        (ind == CGAL::ARR_MIN_END ?
         this->_m_traits->construct_min_vertex_2_object()(cv) :
         this->_m_traits->construct_max_vertex_2_object()(cv));
    
    if (ps_x != CGAL::ARR_INTERIOR) {
        
        CGAL_assertion_code(
                int lod_size = 
                static_cast< int >(this->_m_identification_WE.size())
        );
        CGAL_assertion(
                static_cast< int >(_m_vertices_on_identification_WE.size()) ==
                lod_size
        );
        
        // update the local structure for points on the curve of identification
        typename Identification_WE::iterator it = 
            this->_m_identification_WE.find(key);
        // not existing so far
        CGAL_assertion(it == this->_m_identification_WE.end());
        // therefore insert it
        this->_m_identification_WE.insert(it, std::make_pair(key, v));
        CGAL_assertion(
                static_cast< int >(_m_identification_WE.size()) ==
                lod_size + 1
        );
        
        // store iterator for vertex 
        // -> needed to delete vertex if becoming redundant
        typename Vertices_on_identification_WE::iterator vit = 
            _m_vertices_on_identification_WE.find(v);
        if (vit == _m_vertices_on_identification_WE.end()) {
            _m_vertices_on_identification_WE.insert(
                    vit, std::make_pair(v,it)
            );
            CGAL_assertion(
                    static_cast< int >(_m_vertices_on_identification_WE.size())
                    ==
                    lod_size + 1
            );
        }
        return;
    } 

    // else
    CGAL_assertion(ps_y != CGAL::ARR_INTERIOR);
    CGAL_assertion_code(
            int lod_size = 
            static_cast< int >(this->_m_identification_NS.size())
    );
    CGAL_assertion(
            static_cast< int >(_m_vertices_on_identification_NS.size()) ==
            lod_size
    );
    
    // update the local structure for points on the curve of identification
    typename Identification_NS::iterator it = 
        this->_m_identification_NS.find(key);
    // not existing so far
    CGAL_assertion(it == this->_m_identification_NS.end());
    // therefore insert it
    this->_m_identification_NS.insert(it, std::make_pair(key, v));
    CGAL_assertion(
            static_cast< int >(_m_identification_NS.size()) ==
            lod_size + 1
    );
    
    // store iterator for vertex 
    // -> needed to delete vertex if becoming redundant
    //CGAL_assertion(_m_vertices_on_identification_NS.find(*v) ==
    //               _m_vertices_on_identification_NS.end());
    typename Vertices_on_identification_NS::iterator vit = 
        _m_vertices_on_identification_NS.find(v);
    if (vit == _m_vertices_on_identification_NS.end()) {
        _m_vertices_on_identification_NS.insert(
                vit, std::make_pair(v,it)
        );
        CGAL_assertion(
                static_cast< int >(_m_vertices_on_identification_NS.size())
                ==
                lod_size + 1
        );
    }
    return;
}

//-----------------------------------------------------------------------------
// Locate curve end with respect to dcel
template <class GeomTraits, class Dcel_>
CGAL::Object 
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
locate_curve_end (const X_monotone_curve_2& cv, Arr_curve_end ind,
                  Arr_parameter_space ps_x, Arr_parameter_space ps_y) const 
{
    // status: to implement
    //std::cout << "Arr_torus_topology_traits_2 locate_curve_end"  
    //          << std::endl;
    CGAL_precondition(_valid(ps_x, ps_y));

    // torus does not contain unbounded curves
    
    Vertex* v = NULL;
    typename Identification_NS::iterator ns_it;
    typename Identification_WE::iterator we_it;

    Point_2 key = (ind == CGAL::ARR_MIN_END ?
                   this->m_traits->construct_min_vertex_2_object()(cv) :
                   this->m_traits->construct_max_vertex_2_object()(cv)); 

    if (ps_x != CGAL::ARR_INTERIOR) {
        
        we_it = _m_identification_WE.find (key);
        if (we_it != _m_identification_WE.end()) {
            v = we_it->second;
            return CGAL::make_object(v);
        }
        // else
        we_it = _m_identification_WE.lower_bound(key);
        
        if (we_it != _m_identification_WE.end()) {
            v = we_it->second;
        }
        // TODO return top face?
        
    } else if (v == NULL || ps_y != CGAL::ARR_INTERIOR) {
        
        if (ps_x == CGAL::ARR_INTERIOR) {
            ns_it = _m_identification_NS.find (key);
            if (ns_it != _m_identification_NS.end()) {
                v = ns_it->second;
                return CGAL::make_object(v);
            }
            // else
            ns_it = _m_identification_NS.lower_bound(key);
        } else {
            ns_it = _m_identification_NS.begin();
        }
        
        if (ns_it == _m_identification_NS.end()) {
            return CGAL::make_object(_m_f_top);
        }
        // else 
        v = ns_it->second;
    }
    CGAL_assertion(v != NULL);
    
    return CGAL::make_object(_face_before_vertex_on_discontinuity(v));
}

//-----------------------------------------------------------------------------
// Given two predecessor halfedges that belong to the same inner CCB of
// a face, determine what happens when we insert an edge connecting the
// target vertices of the two edges.
//
template <class GeomTraits, class Dcel_>
std::pair<bool, bool>
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
face_split_after_edge_insertion (const Halfedge *prev1,
                                 const Halfedge *prev2,
                                 const X_monotone_curve_2& cv) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 face_split"  << std::endl;
    
    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());
    
    //std::cout << "Path1: " << std::endl;
    CGAL::Sign sign_12 = _sign_of_path(prev1, prev2, cv);
    //std::cout << "sign1: " << sign_12 << std::endl;
    //std::cout << "Path2: " << std::endl;
    CGAL::Sign sign_21 = _sign_of_path(prev2, prev1, cv);
    //std::cout << "sign2: " << sign_21 << std::endl;
    
    // TODO use arr function for to check perimetry
    bool perimetric = (sign_12 != CGAL::ZERO && sign_21 != CGAL::ZERO);
    
    // on a torus except for one case, there is a face split
    if (perimetric) {
        // must be topface
        if (prev1->inner_ccb()->face() == top_face()) {
            if (prev1->inner_ccb()->face()->number_of_outer_ccbs() == 0) {
                // the special case is when the initial perimetric path is
                // found, this juste creates two outer ccbs for the face
                // that contained a "perimetric" hole before
                //std::cout << "face_split false, false" << std::endl;
                return std::make_pair(false, false);
            }
        }
        
        // else
        // there is a face split, but no hole is created
        //std::cout << "face_split true, false" << std::endl;
        return std::make_pair(true, false);
    }
    // else
    // face is splitted and it forms a new hole in the old
    //std::cout << "face_split true, true";
    return std::make_pair(true, true);
}

//-----------------------------------------------------------------------------
// Determine whether the removal of the given edge will cause the creation
// of a hole.
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
hole_creation_after_edge_removal (const Halfedge *he) const
{
    // status: to implement
    std::cout << "Arr_torus_topology_traits_2 hole_creation"  << std::endl;

    CGAL_error(); // hole_creation not finally implemented for torus

    CGAL_precondition (! he->is_on_inner_ccb());
    CGAL_precondition (! he->opposite()->is_on_inner_ccb());

    // TODO hole_creation for torus
    
    // TODO use arr function for to check perimetry

    // Check whether the halfedge and its twin belong to the same outer CCB
    // (and are therefore incident to the same face).
    if (he->outer_ccb() == he->opposite()->outer_ccb()) {
        // precondition is: does not form an antenna, or a simply to remove
        // halfedge
        
         // Check the two cycles that will be created once we remove he and its
        // twin (from he->next() to he's twin, not inclusive, and from the
        // successor of he's twin to he, not inclusive).
        if (_sign_of_path(he->next(), he->opposite()) != CGAL::ZERO
            &&
            _sign_of_path(he->opposite()->next(), he) != CGAL::ZERO
        ) {
            // Both paths are perimetric, so the two cycles become two separate
            // outer CCBs of the same face, and no hole is created.
            return (false);
        } else {
            // At least one cyclic path is non-perimetic. 
            // This cycle will become
            // an inner CCB representing a hole in the face.
            return (true);
        }
    } else {
        // The edge to be removed separates two faces.
        // Check the cyclic path from he and back, and from its twin and back.
        if (_sign_of_path(he, he) != CGAL::ZERO &&
            _sign_of_path(he->opposite(), he->opposite()) != CGAL::ZERO) {
            if (dcel().number_of_faces() == 1) {
                CGAL_assertion_code(
                        Face *f = dcel()->faces_begin();
                        CGAL_assertion(f->number_of_outer_ccbs() == 2);
                );
                // there is no face merge in this case ... but the remaining
                // path consists of a hole in the interior
                return (true);
            } 
            // else
            // In this case we disconnect a perimetric cycle around the torus,
            // causing two perimetric faces to merge. 
            // The remainder of the cycle
            // becomes an inner CCB (a hole) in the merged face.
            return (true);
        } else {
            // In this case we are about to merge to incident faces, so their
            // outer CCBs are merged and no new hole is created.
            return (false);
        }
    }
}

//-----------------------------------------------------------------------------
// checks whether halfedges are on a new perimetric face boundary
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
is_on_new_perimetric_face_boundary (const Halfedge *prev1,
                                    const Halfedge *prev2,
                                    const X_monotone_curve_2& cv) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2::" 
    //          << "is_on_new_perimetric_face_boundary"
    //          << std::endl;

    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());
    
    // maintain the invariant that the pole is always in the top_face,
    // i.e, it is the face that contains everything and has now outer ccb
    // If pole is part of a ccb itself, it incident face is the face that 
    // contains everything.
    
    CGAL::Sign sign = _sign_of_path(prev2, prev1, cv);
    CGAL_assertion(sign != CGAL::ZERO);
    
    return (sign == CGAL::POSITIVE);
}

//-----------------------------------------------------------------------------
// checks whether halfedges are boundaries of the same face
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
boundaries_of_same_face (const Halfedge *e1, const Halfedge *e2) const {
    
    // status: correct
    //std::cout << " Arr_torus_topology_traits_2::boundaries_of_same_face" 
    //          << std::endl;
    // This predicate is only used for case 3.3.2 of the insertion process

#if 0    
    std::cout << "e1: " << e1->curve() << std::endl;
    std::cout << "dir1: " 
            << (e1->direction() == CGAL::CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L") 
              << std::endl;
    std::cout << "e1->occbf: " << &(*e1->outer_ccb()->face()) << std::endl;
    std::cout << "e2: " << e2->curve() << std::endl;
    std::cout << "dir2: " 
            << (e2->direction() == CGAL::CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L") 
              << std::endl;
    std::cout << "e2->occbf: " << &(*e2->outer_ccb()->face()) << std::endl;
#endif

    // compute signs of path (both must be non-zero)
    CGAL::Sign sign1 = _sign_of_path(e1, e1);
    CGAL_assertion(sign1 != CGAL::ZERO);
    CGAL::Sign sign2 = _sign_of_path(e2, e2);
    CGAL_assertion(sign2 != CGAL::ZERO);
    
    return (sign1 != sign2);
}

//-----------------------------------------------------------------------------
// Determine whether the given vertex lies in the interior of the given face.
//
template <class GeomTraits, class Dcel_>
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
is_in_face (const Face *f, const Point_2& p, const Vertex *v) const
{
    // status: not implemented
    std::cout << "TODO: Arr_torus_topology_traits_2::is_in_face" 
              << std::endl;

    CGAL_error(); // is_in_face not implemented for torus
    // TODO is_in_face NEEDED for incremental insertion
    
    return false;
}

//-----------------------------------------------------------------------------
// Determine whether a boundary vertex is redundant
//
template <class GeomTraits, class Dcel_>
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::is_redundant
(const Vertex *v) const
{
    // status: correct
    //std::cout << "Arr_torus_topology_traits_2 is_redundant"  << std::endl;
    CGAL_precondition(_valid(v->parameter_space_in_x(),v->parameter_space_in_y()));
    
    // if there are not incident edges just remove it
    // TASK: check whether isolated or degree == 0 is needed!
    return (v->is_isolated());
}

//-----------------------------------------------------------------------------
// Determine whether a boundary vertex is redundant
//
template <class GeomTraits, class Dcel_>
typename Arr_torus_topology_traits_2<GeomTraits, Dcel_>::Halfedge* 
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
erase_redundant_vertex (Vertex *v) 
{
    // status: correct
    
    //std::cout << "Arr_torus_topology_traits_2 erase_redundant_vertex"  
    //          << std::endl;
    CGAL_precondition(_valid(v->parameter_space_in_x(),v->parameter_space_in_y()));
    
    // no incident curve-end can give us the key
    // -> but we stored something else useful: find iterator
    if (v->parameter_space_in_x() != CGAL::ARR_INTERIOR) {

        typename Vertices_on_identification_WE::iterator 
            vit = _m_vertices_on_identification_WE.find(v);
        
        // and delete this item
        _m_identification_WE.erase(vit->second);
        _m_vertices_on_identification_WE.erase(vit);
        
    } else {
        CGAL_assertion(v->parameter_space_in_y() != CGAL::ARR_INTERIOR);
        typename Vertices_on_identification_NS::iterator 
            vit = _m_vertices_on_identification_NS.find(v);
        
        // and delete this item
        _m_identification_NS.erase(vit->second);
        _m_vertices_on_identification_NS.erase(vit);
    }
    
    // a valid halfedge-pointer is only requested for if vertex
    // has been connecting fictiuous halfedges, this is not the case here,
    // so we 
    return NULL;
}

// protected:

//-----------------------------------------------------------------------------
// Number of crossing with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_path(const Halfedge* he1, const Halfedge* he2) const {
    
    // status: move to arr
    
    //std::cout << "Arr_torus_topology_traits: "
    //          << "_sign_of_paths" << std::endl;
    
    typedef typename Geometry_traits_2::Curve_kernel_2::Curve_2::Asymptote_y
        Asymptote_y;
    
    int x_counter = 0;
    int y_counter = 0;

    if (he1->next() == he2 && he2->next () == he1) {
        return CGAL::ZERO;
    }

    typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
        _m_traits->parameter_space_in_x_2_object();
    typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
        _m_traits->parameter_space_in_y_2_object();
    
    // Start with the next of prev1:
    const Halfedge * curr = he1->next();
    // Save its src condition
    Arr_curve_end curr_src_ind;
    Arr_curve_end curr_trg_ind;
    if (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
        curr_src_ind = CGAL::ARR_MIN_END;
        curr_trg_ind = CGAL::ARR_MAX_END;
    } else {
        curr_src_ind = CGAL::ARR_MAX_END;
        curr_trg_ind = CGAL::ARR_MIN_END;
    }
    CGAL_assertion(!curr->has_null_curve());
    CGAL_assertion_code(Arr_curve_end first_src_ind = curr_src_ind;);

    Arr_parameter_space first_src_bcx =
      parameter_space_in_x(curr->curve(), curr_src_ind);
    Arr_parameter_space curr_trg_bcx =
      parameter_space_in_x(curr->curve(), curr_trg_ind);  
    Arr_parameter_space first_src_bcy =
      parameter_space_in_y(curr->curve(), curr_src_ind);
    Arr_parameter_space curr_trg_bcy =
      parameter_space_in_y(curr->curve(), curr_trg_ind);  

    while (curr != he2) {

        //std::cout << "curr: " << curr->curve() << std::endl;
        //std::cout << "dir: " 
        //          << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L") 
        //          << std::endl;

        const Halfedge * next = curr->next();
        
        Arr_curve_end next_src_ind;
        Arr_curve_end next_trg_ind;
        if (next->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
            next_src_ind = CGAL::ARR_MIN_END;
            next_trg_ind = CGAL::ARR_MAX_END;
        } else {
            next_src_ind = CGAL::ARR_MAX_END;
            next_trg_ind = CGAL::ARR_MIN_END;
        }
        Arr_parameter_space next_src_bcx = 
            parameter_space_in_x(next->curve(), next_src_ind);
        Arr_parameter_space next_trg_bcx = 
            parameter_space_in_x(next->curve(), next_trg_ind);
        Arr_parameter_space next_src_bcy = 
            parameter_space_in_y(next->curve(), next_src_ind);
        Arr_parameter_space next_trg_bcy = 
            parameter_space_in_y(next->curve(), next_trg_ind);

        if (curr_trg_bcx != next_src_bcx) {
            CGAL_assertion(curr_trg_bcx != CGAL::ARR_INTERIOR);
            CGAL_assertion(next_src_bcx != CGAL::ARR_INTERIOR);
            if (curr_trg_bcx == CGAL::ARR_RIGHT_BOUNDARY) {
                ++x_counter;
                //std::cout << "+x1" << std::endl;
            } else {
                --x_counter;
                //std::cout << "-x1" << std::endl;
            }
        } else {
            // can influence pole
            if (curr_trg_bcx != CGAL::ARR_INTERIOR) {
                CGAL_assertion(next_src_bcx != CGAL::ARR_INTERIOR);
                
                Asymptote_y yc, yn;
                
                if (curr_trg_bcx == CGAL::ARR_LEFT_BOUNDARY) {
                    
                    CGAL_assertion(curr_trg_ind == CGAL::ARR_MIN_END);
                    CGAL_assertion(next_src_ind == CGAL::ARR_MIN_END);
                    
                    yc = curr->curve().curve().
                        horizontal_asymptote_for_arc_to_minus_infinity(
                                curr->curve().arcno()
                        );
                    yn = next->curve().curve().
                        horizontal_asymptote_for_arc_to_minus_infinity(
                                next->curve().arcno()
                        );
                    
                } else {

                    CGAL_assertion(curr_trg_ind == CGAL::ARR_MAX_END);
                    CGAL_assertion(next_src_ind == CGAL::ARR_MAX_END);
                    
                    yc = curr->curve().curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                curr->curve().arcno()
                        );
                    yn = next->curve().curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                next->curve().arcno()
                        );
                    
                }
                
                if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                    
                    if (yc.infty() == NiX::PLUS_INFTY) {
                        // bottom to top
                        //std::cout << "+y1" << std::endl;
                        ++y_counter;
                    } else {
                        // top to bottom
                        --y_counter;
                        //std::cout << "-y1" << std::endl;
                    }
                }
            }
            
        }
        if (curr_trg_bcy != next_src_bcy) {
            CGAL_assertion(curr_trg_bcy != CGAL::ARR_INTERIOR);
            CGAL_assertion(next_src_bcy != CGAL::ARR_INTERIOR);
            if (curr_trg_bcy == CGAL::ARR_TOP_BOUNDARY) {
                ++y_counter;
                //std::cout << "+y2" << std::endl;
            } else {
                --y_counter;
                //std::cout << "-y2" << std::endl;
            }
        }
        curr = next;
        curr_src_ind = next_src_ind;
        curr_trg_ind = next_trg_ind;
        curr_trg_bcx = next_trg_bcx;
        curr_trg_bcy = next_trg_bcy;
    }

    if (he1 == he2) {
        CGAL_assertion_code(Arr_curve_end last_trg_ind = curr_trg_ind;);
        Arr_parameter_space last_trg_bcx = curr_trg_bcx;
        if (last_trg_bcx != first_src_bcx) {
            if (last_trg_bcx == CGAL::ARR_RIGHT_BOUNDARY) {
                ++x_counter;
                //std::cout << "+x2" << std::endl;
            } else {
                --x_counter;
                //std::cout << "-x2" << std::endl;
            }
        } else {
            // can influence pole
            if (last_trg_bcx != CGAL::ARR_INTERIOR) {
                CGAL_assertion(first_src_bcx != CGAL::ARR_INTERIOR);
                
                Asymptote_y yc, yn;
                
                if (curr_trg_bcx == CGAL::ARR_LEFT_BOUNDARY) {
                    
                    CGAL_assertion(last_trg_ind == CGAL::ARR_MIN_END);
                    CGAL_assertion(first_src_ind == CGAL::ARR_MIN_END);
                    
                    yc = curr->curve().curve().
                        horizontal_asymptote_for_arc_to_minus_infinity(
                                curr->curve().arcno()
                        );
                    yn = he1->curve().curve().
                        horizontal_asymptote_for_arc_to_minus_infinity(
                                he1->curve().arcno()
                        );
                    
                } else {

                    CGAL_assertion(last_trg_ind == CGAL::ARR_MAX_END);
                    CGAL_assertion(first_src_ind == CGAL::ARR_MAX_END);
                    
                    yc = curr->curve().curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                curr->curve().arcno()
                        );
                    yn = he1->curve().curve().
                        horizontal_asymptote_for_arc_to_plus_infinity(
                                he1->curve().arcno()
                        );
                }
                
                if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                    
                    if (yc.infty() == NiX::PLUS_INFTY) {
                        // bottom to top
                        ++y_counter;
                        //std::cout << "+y3" << std::endl;
                    } else {
                        // top to bottom
                        --y_counter;
                        //std::cout << "-y3" << std::endl;
                    }
                }
            }
        }

        Arr_parameter_space last_trg_bcy = curr_trg_bcy;

        if (last_trg_bcy != first_src_bcy) {
            if (last_trg_bcy == CGAL::ARR_TOP_BOUNDARY) {
                ++y_counter;
                //std::cout << "+y4" << std::endl;
            } else {
                --y_counter;
                //std::cout << "-y4" << std::endl;
            }
        }
    }
    
    //std::cout << "x: " << x_counter << std::endl;
    //std::cout << "y: " << y_counter << std::endl;

    return (CGAL::sign((x_counter + y_counter) % 2));
}

//-----------------------------------------------------------------------------
// Number of crossing with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_path (const Halfedge* he1, const Halfedge* he2, 
               const X_monotone_curve_2& cv) const {
    
    // status: move to arr

    typedef typename Geometry_traits_2::Curve_kernel_2::Curve_2::Asymptote_y
        Asymptote_y;
    
    CGAL::Sign sign = _sign_of_path(
            he2, he1
    );
    
    int s = sign;

    //std::cout << "sinit: " << s << std::endl;
    
    const Halfedge* prev1 = he1;
    const Halfedge* prev2 = he2;
    
    //std::cout << "prev1: " << prev1->curve() << std::endl;
    //std::cout << "dir1: " 
    //          << (prev1->direction() == CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L")
    //          << std::endl;
    //std::cout << "prev2: " << prev2->curve() << std::endl;
    //std::cout << "dir2: " 
    //          << (prev2->direction() == CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L")
    //          << std::endl;
    //std::cout << "cv: " << cv << std::endl;
    typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
        _m_traits->parameter_space_in_x_2_object();
    typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
        _m_traits->parameter_space_in_y_2_object();
    
    // check whether cv can influence the counters

    Arr_parameter_space bcv1x = parameter_space_in_x(cv, CGAL::ARR_MIN_END);
    Arr_parameter_space bcv1y = parameter_space_in_y(cv, CGAL::ARR_MIN_END);

    Arr_parameter_space bcv2x = parameter_space_in_x(cv, CGAL::ARR_MAX_END);  
    Arr_parameter_space bcv2y = parameter_space_in_y(cv, CGAL::ARR_MAX_END);  
    
    if (bcv1x != CGAL::ARR_INTERIOR || bcv1y != CGAL::ARR_INTERIOR || 
        bcv2x != CGAL::ARR_INTERIOR || bcv2y != CGAL::ARR_INTERIOR) {
        // sign can change!

        bool equal = false;

        Point_2 minp = this->_m_traits->construct_min_vertex_2_object()(cv);
        Arr_parameter_space ps_x = 
            this->_m_traits->parameter_space_in_x_2_object()
            (cv, CGAL::ARR_MIN_END);
        Arr_parameter_space ps_y = 
            this->_m_traits->parameter_space_in_y_2_object()
            (cv, CGAL::ARR_MIN_END);

        bool v1_on_boundary = 
            (prev1->vertex()->parameter_space_in_x() != CGAL::ARR_INTERIOR ||
             prev1->vertex()->parameter_space_in_y() != CGAL::ARR_INTERIOR);
        bool min_on_boundary = 
            (ps_x != CGAL::ARR_INTERIOR || ps_y != CGAL::ARR_INTERIOR);
        
        if (v1_on_boundary == min_on_boundary) {
            if (v1_on_boundary) {
                // compare at boundary
                equal = this->are_equal(prev1->vertex(), cv, 
                                        CGAL::ARR_MIN_END, 
                                        ps_x, ps_y);
            } else {
                equal = (this->_m_traits->compare_xy_2_object()(
                                 prev1->vertex()->point(), minp
                         ) == CGAL::EQUAL);
            }
        }

        if (!equal) {
            std::swap(bcv1x, bcv2x);
            std::swap(bcv1y, bcv2y);
        }

        // orders are now with respect to prev1 and prev2
        
        if (bcv1x != CGAL::ARR_INTERIOR || bcv1y != CGAL::ARR_INTERIOR) {
            // the counter can change at prev1
            
            Arr_curve_end prev1_trg_ind;
            if (prev1->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                prev1_trg_ind = CGAL::ARR_MAX_END;
            } else {
                prev1_trg_ind = CGAL::ARR_MIN_END;
            }
            
            CGAL_assertion(!prev1->has_null_curve());
            Arr_parameter_space prev1_trg_bcx = 
                parameter_space_in_x(prev1->curve(), prev1_trg_ind);
            Arr_parameter_space prev1_trg_bcy = 
                parameter_space_in_y(prev1->curve(), prev1_trg_ind);
            
            if (prev1_trg_bcx != ARR_INTERIOR) {
                if (prev1_trg_bcx != bcv1x) {
                    bool modify = true;

                    if (prev1->vertex()->halfedge() !=
                        prev1->vertex()->halfedge()->opposite()->prev()) {
                        const Halfedge* next1 = prev1->next();
                        
                        Arr_curve_end next1_src_ind;
                        if (next1->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                            next1_src_ind = CGAL::ARR_MIN_END;
                        } else {
                            next1_src_ind = CGAL::ARR_MAX_END;
                        }
                        
                        CGAL_assertion(!next1->has_null_curve());
                        CGAL::Arr_parameter_space next1_src_bcx = 
                            parameter_space_in_x(next1->curve(), 
                                                 next1_src_ind);
                        
                        CGAL_assertion(next1_src_bcx != ARR_INTERIOR);
                        
                        modify = (next1_src_bcx != bcv1x);
                    }
                    
                    if (modify) {
                        if (prev1_trg_bcx == CGAL::ARR_RIGHT_BOUNDARY) {
                            s = (s + 1) % 2;
                            //std::cout << "+s1" << std::endl;
                        } else {
                            s = (s - 1) % 2;
                            //std::cout << "-s1" << std::endl;
                        }
                    }
                    
                } else {
                    
                    CGAL_assertion(bcv1x != CGAL::ARR_INTERIOR);
                    
                    Asymptote_y yc, yn;
                    
                    if (prev1_trg_bcx == CGAL::ARR_LEFT_BOUNDARY) {
                        
                        yc = prev1->curve().curve().
                            horizontal_asymptote_for_arc_to_minus_infinity(
                                    prev1->curve().arcno()
                            );
                        yn = cv.curve().
                            horizontal_asymptote_for_arc_to_minus_infinity(
                                    cv.arcno()
                            );
                        
                    } else {
                        
                        yc = prev1->curve().curve().
                            horizontal_asymptote_for_arc_to_plus_infinity(
                                    prev1->curve().arcno()
                            );
                        yn = cv.curve().
                            horizontal_asymptote_for_arc_to_plus_infinity(
                                    cv.arcno()
                            );
                    }
                    
                    if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                        
                        if (yc.infty() == NiX::PLUS_INFTY) {
                            // bottom to top
                            s = (s + 1) % 2;
                            //std::cout << "+s2" << std::endl;
                        } else {
                            // top to bottom
                            s = (s - 1) % 2;
                            //std::cout << "-s2" << std::endl;
                        }
                    }
                }
            }
            
            if (prev1_trg_bcy != CGAL::ARR_INTERIOR) {
                if (prev1_trg_bcy != bcv1y) {
                    
                    bool modify = true;
                    
                    if (prev1->vertex()->halfedge() !=
                        prev1->vertex()->halfedge()->opposite()->prev()) {
                        const Halfedge* next1 = prev1->next();

                        Arr_curve_end next1_src_ind;
                        if (next1->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                            next1_src_ind = CGAL::ARR_MIN_END;
                        } else {
                            next1_src_ind = CGAL::ARR_MAX_END;
                        }
                        
                        CGAL_assertion(!next1->has_null_curve());
                        CGAL::Arr_parameter_space next1_src_bcy = 
                            parameter_space_in_y(next1->curve(), 
                                                 next1_src_ind);
                        
                        CGAL_assertion(next1_src_bcy != CGAL::ARR_INTERIOR);
        
                        modify = (next1_src_bcy != bcv1y);
                    }
                    
                    if (modify) {
                        if (prev1_trg_bcy == ARR_TOP_BOUNDARY) {
                            s = (s + 1) % 2;
                            //std::cout << "+s3" << std::endl;
                        } else {
                            s = (s - 1) % 2;
                            //std::cout << "-s3" << std::endl;
                        }
                    }
                }
            }
        }
        // now for prev2
        if (bcv2x != CGAL::ARR_INTERIOR || bcv2y != CGAL::ARR_INTERIOR) {
            // the counter can change at prev2

            Arr_curve_end prev2_trg_ind;
            if (prev2->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                prev2_trg_ind = CGAL::ARR_MAX_END;
            } else {
                prev2_trg_ind = CGAL::ARR_MIN_END;
            }

            CGAL_assertion(!prev2->has_null_curve());
            Arr_parameter_space prev2_trg_bcx = 
                parameter_space_in_x(prev2->curve(), prev2_trg_ind);
            Arr_parameter_space prev2_trg_bcy = 
                parameter_space_in_y(prev2->curve(), prev2_trg_ind);
            
            if (prev2_trg_bcx != CGAL::ARR_INTERIOR) {
                if (prev2_trg_bcx != bcv2x) {

                    bool modify = true;

                    if (prev2->vertex()->halfedge() !=
                        prev2->vertex()->halfedge()->opposite()->prev()) {
                        const Halfedge* next2 = prev2->next();
                        
                        Arr_curve_end next2_src_ind;
                        if (next2->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                            next2_src_ind = CGAL::ARR_MIN_END;
                        } else {
                            next2_src_ind = CGAL::ARR_MAX_END;
                        }
                        
                        CGAL_assertion(!next2->has_null_curve());
                        CGAL::Arr_parameter_space next2_src_bcx = 
                            parameter_space_in_x(next2->curve(), 
                                                 next2_src_ind);
                        
                        CGAL_assertion(next2_src_bcx != CGAL::ARR_INTERIOR);
                        
                        modify = (next2_src_bcx != bcv2x);
                    }
                    
                    if (modify) {
                        if (prev2_trg_bcx == CGAL::ARR_RIGHT_BOUNDARY) {
                            s = (s + 1) % 2;
                            //std::cout << "+s4" << std::endl;
                        } else {
                            s = (s - 1) % 2;
                            //std::cout << "-s4" << std::endl;
                        }
                    }

                } else {
                    
                    CGAL_assertion(bcv2x != CGAL::ARR_INTERIOR);
                    
                    Asymptote_y yc, yn;
                    
                    if (prev2_trg_bcx == CGAL::ARR_LEFT_BOUNDARY) {
                        
                        yc = prev2->curve().curve().
                            horizontal_asymptote_for_arc_to_minus_infinity(
                                    prev2->curve().arcno()
                            );
                        yn = cv.curve().
                            horizontal_asymptote_for_arc_to_minus_infinity(
                                    cv.arcno()
                            );
                        
                    } else {
                        
                        yc = prev2->curve().curve().
                            horizontal_asymptote_for_arc_to_plus_infinity(
                                    prev2->curve().arcno()
                            );
                        yn = cv.curve().
                            horizontal_asymptote_for_arc_to_plus_infinity(
                                    cv.arcno()
                            );
                    }
                    
                    if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                        
                        if (yc.infty() == NiX::PLUS_INFTY) {
                            // bottom to top
                            s = (s + 1) % 2;
                            //std::cout << "+s5" << std::endl;
                        } else {
                            // top to bottom
                            s = (s - 1) % 2;
                            //std::cout << "-s5" << std::endl;
                        }
                    }
                }
            }

            if (prev2_trg_bcy != CGAL::ARR_INTERIOR) {
                if (prev2_trg_bcy != bcv2y) {
                    
                    bool modify = true;

                    if (prev2->vertex()->halfedge() !=
                        prev2->vertex()->halfedge()->opposite()->prev()) {

                        const Halfedge* next2 = prev2->next();
                        
                        Arr_curve_end next2_src_ind;
                        if (next2->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                            next2_src_ind = CGAL::ARR_MIN_END;
                        } else {
                            next2_src_ind = CGAL::ARR_MAX_END;
                        }
                        
                        CGAL_assertion(!next2->has_null_curve());
                        CGAL::Arr_parameter_space next2_src_bcy = 
                            parameter_space_in_y(next2->curve(), 
                                                 next2_src_ind);
                        
                        CGAL_assertion(next2_src_bcy != CGAL::ARR_INTERIOR);
                        
                        modify = (next2_src_bcy != bcv2y);
                    }
                    
                    if (modify) {
                        if (prev2_trg_bcy == CGAL::ARR_TOP_BOUNDARY) {
                            s = (s + 1) % 2;
                            //std::cout << "+s6" << std::endl;
                        } else {
                            s = (s - 1) % 2;
                            //std::cout << "-s6" << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    //std::cout << "s: " << s << std::endl;

    return (CGAL::sign(s % 2));
}

/*! \brief Return the face that lies before the given vertex, which lies
 * on the line of discontinuity.
 */
template <class GeomTraits, class Dcel>
typename Arr_torus_topology_traits_2<GeomTraits, Dcel>::Face *
Arr_torus_topology_traits_2<GeomTraits, Dcel>::
_face_before_vertex_on_identifications (Vertex * v) const {
    
    // status: correct

    // If the vertex is isolated, just return the face that contains it.
    if (v->is_isolated()) {
        return (v->isolated_vertex()->face());
    }
    
    // Get the first incident halfedge around v and the next halfedge.
    Halfedge  *first = v->halfedge();
    Halfedge  *curr = first;
    CGAL_assertion(curr != NULL);
    Halfedge  *next = curr->next()->opposite();
    
    // If there is only one halfedge incident to v, return its incident
    // face.
    if (curr == next) {
        if (curr->is_on_inner_ccb()) {
            return (curr->inner_ccb()->face());
        } else {
            return (curr->outer_ccb()->face());
        }
    }

    if (v->parameter_space_in_x() != CGAL::ARR_INTERIOR) {

        // Otherwise, we traverse the halfedges around v and locate the first
        // halfedge we encounter if we go from "6 o'clock" clockwise.
        // First locate the lower left and the top right halfedges around v.
        typename Traits_adaptor_2::Compare_y_at_x_right_2 
            compare_y_at_x_right =
            _m_traits->compare_y_at_x_right_2_object();
        typename Traits_adaptor_2::Compare_y_at_x_left_2 compare_y_at_x_left =
            _m_traits->compare_y_at_x_left_2_object();
        
        Halfedge  *lowest_left = NULL;
        Halfedge  *top_right = NULL;
        
        do {
            // Check whether the current halfedge is defined 
            // to the left or to the right of the given vertex.
            if (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
                // The curve associated with the current halfedge 
                // is defined to the left of v.
                if (lowest_left == NULL ||
                    compare_y_at_x_left (curr->curve(),
                                         lowest_left->curve(), 
                                         v->point()) == SMALLER) {
                    lowest_left = curr;
                }
            } else {
                // The curve associated with the current halfedge 
                // is defined to the right of v.
                if (top_right == NULL ||
                    compare_y_at_x_right (curr->curve(),
                                          top_right->curve(), 
                                          v->point()) == LARGER) {
                    top_right = curr;
                }
            }
            
            // Move to the next halfedge around the vertex.
            curr = curr->next()->opposite();
            
        } while (curr != first);
        
        // The first halfedge we encounter is the lowest to the left, 
        // but if there is no edge to the left, we first encounter the 
        // topmost halfedge to the right. Note that as the halfedge we 
        // located has v as its target, we now have to return its twin.
        if (lowest_left != NULL) {
            first = lowest_left->opposite();
        } else {
            first = top_right->opposite();
        }
        
    } else {
        CGAL_assertion(v->parameter_space_in_y() != CGAL::ARR_INTERIOR);
        
        // Otherwise, we traverse the halfedges around v and locate the first
        // halfedge we encounter if we go from "3 o'clock" clockwise.
        // First locate the lower left and the top right halfedges around v.
        typename Traits_adaptor_2::Compare_x_2 compare_x =
            _m_traits->compare_x_2_object();
        
        CGAL::Arr_curve_end leftmost_top_end = CGAL::ARR_MIN_END;
        Halfedge  *leftmost_top = NULL;
        CGAL::Arr_curve_end rightmost_bottom_end = CGAL::ARR_MIN_END;
        Halfedge  *rightmost_bottom = NULL;
        
        do {
            typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
                _m_traits->parameter_space_in_x_2_object();
            typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
                _m_traits->parameter_space_in_y_2_object();
            
            CGAL::Arr_curve_end ind = CGAL::ARR_MIN_END;
            
            Arr_parameter_space bd_x = 
                parameter_space_in_x(curr->curve(), CGAL::ARR_MAX_END);
            Arr_parameter_space bd_y = 
                parameter_space_in_y(curr->curve(), CGAL::ARR_MAX_END);
            if (are_equal(v, curr->curve(), CGAL::ARR_MAX_END, bd_x, bd_y)) {
                ind = CGAL::ARR_MAX_END;
            }
            
            if (parameter_space_in_y(curr->curve(),ind) == 
                CGAL::ARR_BOTTOM_BOUNDARY) {
                // TOP-side
                if ((leftmost_top == NULL) || 
                    (leftmost_top->direction() == CGAL::ARR_LEFT_TO_RIGHT &&
                     leftmost_top->direction() != curr->direction()) ||
                    (leftmost_top->direction() == curr->direction() &&
                     compare_x(curr->curve(), ind, 
                               leftmost_top->curve(), leftmost_top_end) == 
                     CGAL::SMALLER)) {
                    leftmost_top_end = ind;
                    leftmost_top = curr;
                } 
            } else {
                // same for BOTTOM-side
                
                if ((rightmost_bottom == NULL) || 
                    (rightmost_bottom->direction() == CGAL::ARR_RIGHT_TO_LEFT &&
                     rightmost_bottom->direction() != curr->direction()) ||
                    (rightmost_bottom->direction() == curr->direction() &&
                     compare_x(curr->curve(), ind, 
                               rightmost_bottom->curve(), 
                               rightmost_bottom_end) 
                     == CGAL::LARGER)) {
                    rightmost_bottom_end = ind;
                    rightmost_bottom = curr;
                } 
            }
            
            // Move to the next halfedge around the vertex.
            curr = curr->next()->opposite();
            
        } while (curr != first);
        
        // The first halfedge we encounter is the leftmost to top, but if there
        // is no edge to the left, we first encounter the righmost halfedge 
        // to the bottom. Note that as the halfedge we located has v as 
        // its target, we now have to return its twin.
        if (leftmost_top != NULL) {
            first = leftmost_top->opposite();
        } else {
            first = rightmost_bottom->opposite();
        }
    }

    // Return the incident face.
    if (first->is_on_inner_ccb()) {
        return (first->inner_ccb()->face());
    } else {
        return (first->outer_ccb()->face());
    }
}

CGAL_END_NAMESPACE

#endif
// EOF
