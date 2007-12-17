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
    std::cout << "dir: " << curr->direction() << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << next->direction() << std::endl;
    
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
        std::cout << "dir: " << curr->direction() << std::endl;
        std::cout << "next: " << next->curve() << std::endl;
        std::cout << "dir: " << next->direction() << std::endl;
        
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
    std::cout << "dir: " << curr->direction() << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << next->direction() << std::endl;
    
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
    std::cout << "dir1: " << e1->direction() << std::endl;
    std::cout << "e1->occbf: " << &(*e1->outer_ccb()->face()) << std::endl;
    std::cout << "e2: " << e2->curve() << std::endl;
    std::cout << "dir2: " << e2->direction() << std::endl;
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
// Returns sign of "crossing" with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_subpath(const Halfedge* he1, const Halfedge* he2) const {

    CGAL::Sign result = CGAL::ZERO;
    
    CGAL::Arr_curve_end end1 = 
        (he1->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
         CGAL::ARR_MAX_END : CGAL::ARR_MIN_END);
    
    CGAL::Arr_curve_end end2 = 
        (he2->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
         CGAL::ARR_MIN_END : CGAL::ARR_MAX_END);
    
    typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
        _m_traits->parameter_space_in_x_2_object();
    typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
        _m_traits->parameter_space_in_y_2_object();
    
    CGAL::Arr_parameter_space he1_psx =
        parameter_space_in_x(he1->curve(), end1);

    CGAL::Arr_parameter_space he1_psy =
        parameter_space_in_y(he1->curve(), end1);
    
    CGAL::Arr_parameter_space he2_psx =
        parameter_space_in_x(he2->curve(), end2);

    CGAL::Arr_parameter_space he2_psy =
        parameter_space_in_y(he2->curve(), end2);
    
    if (he1_psx == CGAL::ARR_INTERIOR &&
        he1_psy == CGAL::ARR_INTERIOR &&
        he2_psx == CGAL::ARR_INTERIOR &&
        he2_psy == CGAL::ARR_INTERIOR) {
        return result;
    }
    
    if (he1_psx != he2_psx) {
        
        if (he1_psx == CGAL::ARR_RIGHT_BOUNDARY) {
            result = CGAL::POSITIVE;
            //std::cout << "+xp1" << std::endl;
        } else {
            result = CGAL::NEGATIVE;
            //std::cout << "-xn1" << std::endl;
        }
        
    } else { 
        
        if (he1_psx != CGAL::ARR_INTERIOR) {
            // can influence pole
            
            CGAL_assertion(he2_psx != CGAL::ARR_INTERIOR);
            
            typename Geometry_traits_2::Curve_kernel_2::Curve_2::Asymptote_y
                yc, yn;
            
            if (he1_psx == CGAL::ARR_LEFT_BOUNDARY) {
                
                CGAL_assertion(end1 == CGAL::ARR_MIN_END);
                CGAL_assertion(end2 == CGAL::ARR_MIN_END);
                
                yc = he1->curve().curve().
                    horizontal_asymptote_for_arc_to_minus_infinity(
                            he1->curve().arcno()
                    );
                yn = he2->curve().curve().
                    horizontal_asymptote_for_arc_to_minus_infinity(
                            he2->curve().arcno()
                    );
                
            } else {
                
                CGAL_assertion(end1 == CGAL::ARR_MAX_END);
                CGAL_assertion(end2 == CGAL::ARR_MAX_END);
                
                yc = he1->curve().curve().
                    horizontal_asymptote_for_arc_to_plus_infinity(
                            he1->curve().arcno()
                    );
                yn = he2->curve().curve().
                    horizontal_asymptote_for_arc_to_plus_infinity(
                            he2->curve().arcno()
                    );
            }
            
            if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                
                if (yc.infty() == NiX::PLUS_INFTY) {
                    // bottom to top
                    result = CGAL::POSITIVE;
                    //std::cout << "+xp2" << std::endl;
                } else {
                    // top to bottom
                    result = CGAL::NEGATIVE;
                    //std::cout << "-xp2" << std::endl;
                }
            }
        } else {
            
            CGAL_assertion(he1_psy != CGAL::ARR_INTERIOR);
            CGAL_assertion(he2_psy != CGAL::ARR_INTERIOR);
            
            if (he1_psy != he2_psy) {
                if (he1_psy == CGAL::ARR_TOP_BOUNDARY) {
                    result = CGAL::POSITIVE;
                    //std::cout << "+yp1" << std::endl;
                } else {
                    result = CGAL::NEGATIVE;
                    //std::cout << "-yn1" << std::endl;
                }
            }
        }
    }
        
    return result;
}

//-----------------------------------------------------------------------------
// Returns sign of "crossing" with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_subpath (const Halfedge* he1,
                  const X_monotone_curve_2& cv2,
                  const CGAL::Arr_curve_end& end2) const {
    
    //std::cout << "he1: " << he1->curve() << std::endl;
    //std::cout << "dir1: " << he1->direction() << std::endl;
    //std::cout << "cv2: " << cv2 << std::endl;
    //std::cout << "end: " << end2 << std::endl;

    CGAL_assertion(!he1->has_null_curve());
    
    CGAL::Sign result = CGAL::EQUAL;
    
    typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
        _m_traits->parameter_space_in_x_2_object();
    typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
        _m_traits->parameter_space_in_y_2_object();
    
    // check whether cv can influence the counters
    
    CGAL::Arr_parameter_space ps_x = 
        parameter_space_in_x(cv2, end2);
    CGAL::Arr_parameter_space ps_y = 
        parameter_space_in_y(cv2, end2);
    
    CGAL::Arr_curve_end he1_trg_end =
        (he1->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
         CGAL::ARR_MAX_END : CGAL::ARR_MIN_END
        );
    
    CGAL::Arr_parameter_space he1_trg_ps_x = 
        parameter_space_in_x(he1->curve(), he1_trg_end);
    CGAL::Arr_parameter_space he1_trg_ps_y = 
        parameter_space_in_y(he1->curve(), he1_trg_end);
        
    if (he1_trg_ps_x != CGAL::ARR_INTERIOR) {
        
        if (he1_trg_ps_x != ps_x) {
            // possible jump over x

            bool modify = true;
            
            // check next
            if (he1->vertex()->halfedge() !=
                he1->vertex()->halfedge()->opposite()->prev()) {
                
                const Halfedge* next1 = he1->next();
                CGAL_assertion(!next1->has_null_curve());
                
                CGAL::Arr_curve_end next1_src_end =
                    (next1->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END);
                
                CGAL::Arr_parameter_space next1_src_ps_x = 
                    parameter_space_in_x(next1->curve(), 
                                         next1_src_end);
                
                CGAL_assertion(next1_src_ps_x != CGAL::ARR_INTERIOR);
                
                modify = (next1_src_ps_x != ps_x);
            }
            
            if (modify) {
                if (he1_trg_ps_x == CGAL::ARR_RIGHT_BOUNDARY) {
                    result = CGAL::POSITIVE;
                    //std::cout << "+xp3" << std::endl;
                } else {
                    result = CGAL::NEGATIVE;
                    //std::cout << "-xn3" << std::endl;
                }
            }
            
        } else {
            
            CGAL_assertion(ps_x != CGAL::ARR_INTERIOR);
            
            typename Geometry_traits_2::Curve_kernel_2::Curve_2::Asymptote_y
                yc, yn;
            
            if (he1_trg_ps_x == CGAL::ARR_LEFT_BOUNDARY) {
                
                yc = he1->curve().curve().
                    horizontal_asymptote_for_arc_to_minus_infinity(
                            he1->curve().arcno()
                    );
                yn = cv2.curve().
                    horizontal_asymptote_for_arc_to_minus_infinity(
                            cv2.arcno()
                    );
                
            } else {
                
                yc = he1->curve().curve().
                    horizontal_asymptote_for_arc_to_plus_infinity(
                            he1->curve().arcno()
                    );
                yn = cv2.curve().
                    horizontal_asymptote_for_arc_to_plus_infinity(
                            cv2.arcno()
                    );
            }
            
            if (!yc.is_finite() && !yn.is_finite() && yc != yn) {
                
                if (yc.infty() == NiX::PLUS_INFTY) {
                    // bottom to top
                    result = CGAL::POSITIVE;
                    //std::cout << "+xp4" << std::endl;
                } else {
                    // top to bottom
                    result = CGAL::NEGATIVE;
                    //std::cout << "-xn4" << std::endl;
                }
            }
        }
    } else if (he1_trg_ps_y != CGAL::ARR_INTERIOR) {
        
        if (he1_trg_ps_y != ps_y) {
            // possible jump over y

            bool modify = true;
            
            // check next
            if (he1->vertex()->halfedge() !=
                he1->vertex()->halfedge()->opposite()->prev()) {
                
                const Halfedge* next1 = he1->next();
                CGAL_assertion(!next1->has_null_curve());
                
                CGAL::Arr_curve_end next1_src_end =
                    (next1->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                     CGAL::ARR_MIN_END : CGAL::ARR_MAX_END);
                
                CGAL::Arr_parameter_space next1_src_ps_y = 
                    parameter_space_in_y(next1->curve(), 
                                         next1_src_end);
                
                CGAL_assertion(next1_src_ps_y != CGAL::ARR_INTERIOR);
                
                modify = (next1_src_ps_y != ps_y);
            }
            
            if (modify) {
                if (he1_trg_ps_y == CGAL::ARR_TOP_BOUNDARY) {
                    result = CGAL::POSITIVE;
                    //std::cout << "+yp2" << std::endl;
                } else {
                    result = CGAL::NEGATIVE;
                    //std::cout << "-yn2" << std::endl;
                }
            }
        }
    }
    
    //std::cout << "result: " << result << std::endl;
    return result;
}

//-----------------------------------------------------------------------------
// Returns sign of crossings with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_path(const Halfedge* he1, const Halfedge* he2) const {
    
    // status: move to arr
    
    //std::cout << "Arr_torus_topology_traits: "
    //          << "_sign_of_subpath" << std::endl;
    
    CGAL::Sign result = CGAL::ZERO;
    
    if (he1->next() == he2 && he2->next () == he1) {
        return result;
    }
    
    // Start with the next of prev1:
    const Halfedge* curr = he1->next();
    
    while (curr != he2) {

        CGAL_assertion(!curr->has_null_curve());
        
        const Halfedge* next = curr->next();

        CGAL_assertion(!next->has_null_curve());
        
        CGAL::Sign tmp = _sign_of_subpath(curr, next);
        
        if (tmp != CGAL::ZERO) {
            switch (result) {
            case ZERO:
                result = tmp;
                break;
            default:
                CGAL_assertion(result == -tmp || result == tmp);
                result = CGAL::ZERO;
            }
        }
        
        curr = next;
    }

    if (he1 == he2) {
        CGAL::Sign tmp = _sign_of_subpath(he1, he1->next());
        
        if (tmp != CGAL::ZERO) {
            switch (result) {
            case ZERO:
                result = tmp;
                break;
            default:
                CGAL_assertion(result == -tmp || result == tmp);
                result = CGAL::ZERO;
            }
        }

    }
    
    return result;
}

//-----------------------------------------------------------------------------
// Returns sign of crossings with the curve of identification
//
template <class GeomTraits, class Dcel_>
CGAL::Sign
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_sign_of_path (const Halfedge* he1, const Halfedge* he2, 
               const X_monotone_curve_2& cv) const {
    
    // status: move to arr
    CGAL::Sign result = _sign_of_path(he2, he1);


    typename Traits_adaptor_2::Parameter_space_in_x_2 parameter_space_in_x =
        _m_traits->parameter_space_in_x_2_object();
    typename Traits_adaptor_2::Parameter_space_in_y_2 parameter_space_in_y =
        _m_traits->parameter_space_in_y_2_object();
    
    // check whether cv can influence the counters
    
    CGAL::Arr_parameter_space ps_min_x = 
        parameter_space_in_x(cv, CGAL::ARR_MIN_END);
    CGAL::Arr_parameter_space ps_min_y = 
        parameter_space_in_y(cv, CGAL::ARR_MIN_END);
    
    CGAL::Arr_parameter_space ps_max_x = 
        parameter_space_in_x(cv, CGAL::ARR_MAX_END);  
    CGAL::Arr_parameter_space ps_max_y = 
        parameter_space_in_y(cv, CGAL::ARR_MAX_END);  
    
    if (ps_min_x != CGAL::ARR_INTERIOR || ps_min_y != CGAL::ARR_INTERIOR || 
        ps_max_x != CGAL::ARR_INTERIOR || ps_max_y != CGAL::ARR_INTERIOR) {
        
        // sign can change!
        
        bool equalmin = false;
        
        Point_2 minp = this->_m_traits->construct_min_vertex_2_object()(cv);
        
        CGAL::Arr_curve_end he1_trg_ind =
            (he1->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
             CGAL::ARR_MAX_END : CGAL::ARR_MIN_END
            );
        
        CGAL::Arr_parameter_space he1_trg_psx = 
            parameter_space_in_x(he1->curve(), he1_trg_ind);
        CGAL::Arr_parameter_space he1_trg_psy = 
            parameter_space_in_y(he1->curve(), he1_trg_ind);
        
        bool v1_on_boundary = 
            (he1_trg_psx != CGAL::ARR_INTERIOR ||
             he1_trg_psy != CGAL::ARR_INTERIOR);
        
        bool min_on_boundary = 
            (ps_min_x != CGAL::ARR_INTERIOR || ps_min_y != CGAL::ARR_INTERIOR);
        
        if (v1_on_boundary == min_on_boundary) {
            if (v1_on_boundary) {
                // compare at boundary
                equalmin = this->are_equal(he1->vertex(), cv, 
                                           CGAL::ARR_MIN_END, 
                                           ps_min_x, ps_min_y);
            } else {
                equalmin = (this->_m_traits->compare_xy_2_object()(
                                    he1->vertex()->point(), minp
                            ) == CGAL::EQUAL);
            }
        }

        if (ps_min_x != CGAL::ARR_INTERIOR || ps_min_y != CGAL::ARR_INTERIOR) {
            
            CGAL::Sign tmp1 = 
                _sign_of_subpath(
                        he1, cv, 
                        (equalmin ? CGAL::ARR_MIN_END : CGAL::ARR_MAX_END)
                );
            
            if (tmp1 != CGAL::ZERO) {
                switch (result) {
                case ZERO:
                    result = tmp1;
                    break;
                default:
                    CGAL_assertion(result == -tmp1 || result == tmp1);
                    result = CGAL::ZERO;
                }
            }
        }
        
        if (ps_min_x != CGAL::ARR_INTERIOR || ps_min_y != CGAL::ARR_INTERIOR) {
            
            CGAL::Sign tmp2 = 
                _sign_of_subpath(
                        he2, cv, 
                        (equalmin ? CGAL::ARR_MAX_END : CGAL::ARR_MIN_END)
                );
            
            if (tmp2 != CGAL::ZERO) {
                switch (result) {
                case ZERO:
                    result = tmp2;
                    break;
                default:
                    CGAL_assertion(result == -tmp2 || result == tmp2);
                    result = CGAL::ZERO;
                }
            }
        }
    }
    
    return result;
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
