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
void Arr_torus_topology_traits_2<GeomTraits, Dcel_>::assign
    (const Self& other)
{
    // status: missing dcel-assign
    std::cout << "Arr_torus_topology_traits_2 assign"  << std::endl;
    // Assign the class.
    // Clear the current DCEL and duplicate the other DCEL.
    _m_dcel.delete_all();
    _m_dcel.assign (other._m_dcel);
    
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
    // Go over the DCEL vertices and locate all points with boundary condition
    typename Dcel::Vertex_iterator       vit;
    Boundary_type                        bx, by;

    for (vit = this->m_dcel.vertices_begin();
         vit != this->m_dcel.vertices_end(); ++vit) {
        // First check whether the vertex has a boundary condition in x.
        bx = vit->boundary_in_x();
        if (bx != CGAL::NO_BOUNDARY) {
            
            std::pair< typename Identification_WE::iterator, bool > res =
                _m_identification_WE.insert (std::make_pair(vit->point(),
                                                            &(*vit)));
            CGAL_assertion(! res.second);
            
            _m_vertices_on_identification_WE[&(*vit)] = res.first;
        }
        

        // First check whether the vertex has a boundary condition in y.
        by = vit->boundary_in_y();
        if (by != CGAL::NO_BOUNDARY) {
            
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
                    Identification_crossing leftmost_crossing1;
                    Identification_crossing bottommost_crossing1;
                    
                    bool touching1 = false;
                    bool crossing1 = false;
                    
                    Vertex *leftmost1;
                    Vertex *bottommost1;
                    
                    std::pair< int, int > counters1 = 
                        _crossings_with_identifications(
                                e1, e1, touching1, crossing1,
                                leftmost1, leftmost_crossing1, 
                                bottommost1, bottommost_crossing1);
                    
                    CGAL_assertion_code(bool tpole1;);
                    CGAL_assertion(
                            _is_perimetric_data(counters1, touching1, 
                                                crossing1, tpole1)
                    );
                    
                    
                    Identification_crossing leftmost_crossing2;
                    Identification_crossing bottommost_crossing2;
                    
                    bool touching2 = false;
                    bool crossing2 = false;
                    
                    Vertex *leftmost2;
                    Vertex *bottommost2;
                    
                    std::pair< int, int > counters2 = 
                        _crossings_with_identifications(
                                e2, e2, touching2, crossing2,
                                leftmost2, leftmost_crossing2, 
                                bottommost2, bottommost_crossing2);
                    
                    CGAL_assertion_code(bool tpole2;);
                    CGAL_assertion(
                            _is_perimetric_data(counters2, touching2, 
                                                crossing2, tpole2)
                    );
                    
                    bool check_lowest = false;

                    Point_2_less_NS less_ns(_m_traits);
                    Point_2_less_WE less_we(_m_traits);
                    
                    if (counters1.first % 2 != 0) {
                        CGAL_assertion(counters2.first % 2 != 0);
                        CGAL_assertion(counters1.second % 2 == 0);
                        CGAL_assertion(counters2.second % 2 == 0);
                        
                        if (less_we(bottommost1->point(), 
                                    bottommost2->point())) {
                            if (crossing1 == BEFORE_TO_AFTER) {
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
    
    return;
}

//-----------------------------------------------------------------------------
// Initialize an empty DCEL structure.
//
template <class GeomTraits, class Dcel_>
void Arr_torus_topology_traits_2<GeomTraits, Dcel_>::init_dcel ()
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2 init_dcel"  
              << std::endl;

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
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::compare_y_at_x
(const Point_2& p, const Halfedge* he) const
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
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::are_equal
(const Vertex *v,
 const X_monotone_curve_2& cv, Curve_end ind,
 CGAL::Boundary_type bound_x, CGAL::Boundary_type bound_y) const
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2 are_equal"  
              << std::endl;
    
    CGAL_precondition(_valid(bound_x, bound_y));
    
    // In case the given boundary conditions do not match those of the given
    // vertex, v cannot represent the curve end.
    if (bound_x != v->boundary_in_x() || bound_y != v->boundary_in_y()) {
        return false;
    }
    
    // check wether the two concrete points are equal
    return (this->_m_traits->compare_xy_2_object() (
                    v->point(),
                    (ind == CGAL::MIN_END ?
                     this->_m_traits->construct_min_vertex_2_object()(cv) :
                     this->_m_traits->construct_max_vertex_2_object()(cv))) 
            == CGAL::EQUAL
    );
}

//-----------------------------------------------------------------------------
// Given a curve end with boundary conditions and a face that contains the
// interior of the curve, find a place for a boundary vertex that will
// represent the curve end along the face boundary.
//
template <class GeomTraits, class Dcel_>
CGAL::Object
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::place_boundary_vertex
    (Face *f,
     const X_monotone_curve_2& cv, CGAL::Curve_end ind,
     Boundary_type bound_x, Boundary_type bound_y)
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2 place_boundary_vertex"  
              << std::endl;

    CGAL_precondition(_valid(bound_x, bound_y));

    // this topology return either an empty object or a DCEL vertex,
    // but never a fictious edge!!!
    
    Vertex *v = NULL;
    
    const Point_2& key = 
        (ind == CGAL::MIN_END ?
         this->_m_traits->construct_min_vertex_2_object()(cv) :
         this->_m_traits->construct_max_vertex_2_object()(cv));
    
    if (bound_x != CGAL::NO_BOUNDARY) {
        // locate curve-end (here a concrete point) 
        // in local structure of for points on identification_WE
        v = vertex_WE(key);
    } else {
        CGAL_assertion(bound_y != CGAL::NO_BOUNDARY);
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
    CGAL_assertion(v->boundary_in_x() == bound_x && 
                   v->boundary_in_y() == bound_y);

    CGAL_assertion(!v->has_null_point());
    return (CGAL::make_object (v));
}

//-----------------------------------------------------------------------------
// Locate the predecessor halfedge for the given curve around a given
// vertex with boundary conditions.
//
template <class GeomTraits, class Dcel_>
typename Arr_torus_topology_traits_2<GeomTraits,Dcel_>::Halfedge* 
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::locate_around_boundary_vertex
    (Vertex *v,
     const X_monotone_curve_2& cv, Curve_end ind,
     Boundary_type bound_x, Boundary_type bound_y) const
{
    // status: correct
    CGAL_precondition(_valid(bound_x, bound_y));

    std::cout << "Arr_torus_topology_traits_2 locate_around_boundary_vertex"  
              << std::endl;
    
     // If the vertex is isolated, there is no predecssor halfedge.
    if (v->is_isolated()) {
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
    std::cout << "dir: " << (curr->direction() == LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << (next->direction() == LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    
    std::cout << "******************************************" << std::endl;
#endif
    
    while (!is_between_cw(cv, (ind == MIN_END),
                          curr->curve(), 
                          (curr->direction() == RIGHT_TO_LEFT),
                          next->curve(), 
                          (next->direction() == RIGHT_TO_LEFT),
                          v->point(), eq_curr, eq_next))
    {
        // The curve must not be equal to one of the curves 
        // already incident to v.
        CGAL_assertion(!eq_curr && !eq_next);
        
        // Move to the next pair of incident halfedges.
        curr = next;
        next = curr->next()->opposite();
        
        // Make sure we have not completed a full traversal around v without
        // locating a place for the new curve xc.
        CGAL_assertion (curr != first);
    }
    return curr;
}

//-----------------------------------------------------------------------------
// Notifies on the creation of a boundary vertex
//
template <class GeomTraits, class Dcel_>
void 
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::
notify_on_boundary_vertex_creation
(Vertex *v,
 const X_monotone_curve_2& cv,
 Curve_end ind,
 Boundary_type bound_x,
 Boundary_type bound_y) const
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2::" 
              << "notify_on_boundary_vertex_creation"
              << std::endl;       

    CGAL_precondition(_valid(bound_x, bound_y));
    
    CGAL_assertion(v->boundary_in_x() == bound_x);
    CGAL_assertion(v->boundary_in_y() == bound_y);

    CGAL_assertion(!v->has_null_point());

    const Point_2& key = 
        (ind == CGAL::MIN_END ?
         this->_m_traits->construct_min_vertex_2_object()(cv) :
         this->_m_traits->construct_max_vertex_2_object()(cv));
    
    if (bound_x != CGAL::NO_BOUNDARY) {
        
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
    CGAL_assertion(bound_y != CGAL::NO_BOUNDARY);
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
// checks whether two halfedges form a perimetric path
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::_is_perimetric_path
(const Halfedge *e1,
 const Halfedge *e2,
 bool& touching_pole) const
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2::is_perimetric_path" 
              << std::endl;
    
    Vertex *leftmost;
    Vertex *bottommost;
    
    Identification_crossing leftmost_crossing;
    Identification_crossing bottommost_crossing;
    
    bool touching = false;
    bool crossing = false;
    
    std::pair< int, int > counters = 
        _crossings_with_identifications(e1, e2, touching, crossing,
                                        leftmost, leftmost_crossing, 
                                        bottommost, bottommost_crossing);
    
    return _is_perimetric_data(counters, touching, crossing, touching_pole);
}

//-----------------------------------------------------------------------------
// checks whether given data indicates a perimetric path
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::_is_perimetric_data
(const std::pair< int, int >& counters, bool touching, bool crossing,
 bool& touching_pole) const
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2::is_perimetric_data" 
              << std::endl;
    
    if (!crossing) {
        return false;
    }
    
    // path crosses identification, which includes "crossing" at pole

    int x_counter = counters.first;
    int y_counter = counters.second;
    
    touching_pole = false;

    // it is perimetric if it crosses NS or WE an odd number of times
    // or both identification an even number of times ("% 2" includes 
    // "degenerate" crossing at pole)
    if (x_counter % 2 != 0 && y_counter == 0 ||
        x_counter == 0 && y_counter % 2 != 0) {
        return true;
    } else if (touching && x_counter == 0 && y_counter == 0) {
        CGAL_assertion(crossing);
        // touching path are holes touching_pole = (touching && !crossing);
        return true;
    }
    // else 
    return false;
}

//-----------------------------------------------------------------------------
// Given two predecessor halfedges that belong to the same inner CCB of
// a face, determine what happens when we insert an edge connecting the
// target vertices of the two edges.
//
template <class GeomTraits, class Dcel_>
std::pair<bool, bool>
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::face_split_after_edge_insertion
    (const Halfedge *prev1,
     const Halfedge *prev2) const
{
    // status: correct
    std::cout << "Arr_torus_topology_traits_2 face_split"  << std::endl;
    
    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());
    
    bool tpole1, tpole2;

    bool perimetric = 
        (_is_perimetric_path (prev1->next(), prev2->next(), tpole1) ||
         _is_perimetric_path (prev2->next(), prev1->next(), tpole2));
    
    // on a torus except for one case, there is a face split
    if (perimetric) {
#if 0
        if (tpole1 && tpole2) {
            std::cout << "face_split D true, true" << std::endl;
            return std::make_pair(true, true);
        } else 
#endif
            if (this->number_of_valid_faces() == 1) {
                // must be topface
                CGAL_assertion(prev1->inner_ccb()->face() == top_face());
                if (top_face()->number_of_outer_ccbs() == 0) {
                    // the special case is when the initial perimetric path is
                    // found, this juste creates two outer ccbs for the face
                    // that contained a "perimetric" hole before
                    std::cout << "face_split A false, false" << std::endl;
                    return std::make_pair(false, false);
                } 
            }
        
        // else
        // there is a face split, but no hole is created
        std::cout << "face_split A true, false" << std::endl;
        return std::make_pair(true, false);
    }
    // else
    // face is splitted and it forms a new hole in the old
    std::cout << "face_split A true, true" << std::endl;
    return std::make_pair(true, true);
}

//-----------------------------------------------------------------------------
// Determine whether the removal of the given edge will cause the creation
// of a hole.
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::hole_creation_after_edge_removal
(const Halfedge *he) const
{
    // status: check implementation
    std::cout << "Arr_torus_topology_traits_2 hole_creation"  << std::endl;

    CGAL_assertion(false); // hole_creation not finally implemented for torus

    CGAL_precondition (! he->is_on_inner_ccb());
    CGAL_precondition (! he->opposite()->is_on_inner_ccb());

    // TODO hole_creation for torus
    
    // Check whether the halfedge and its twin belong to the same outer CCB
    // (and are therefore incident to the same face).
    if (he->outer_ccb() == he->opposite()->outer_ccb()) {
        // precondition is: does not form an antenna, or a simply to remove
        // halfedge
        
        bool tpole1, tpole2;
        
        // Check the two cycles that will be created once we remove he and its
        // twin (from he->next() to he's twin, not inclusive, and from the
        // successor of he's twin to he, not inclusive).
        if (_is_perimetric_path (he->next(), he->opposite(), tpole1) &&
            _is_perimetric_path (he->opposite()->next(), he, tpole2)) {
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
        bool tpole1, tpole2;
        // The edge to be removed separates two faces.
        // Check the cyclic path from he and back, and from its twin and back.
        if (_is_perimetric_path (he, he, tpole1) &&
            _is_perimetric_path (he->opposite(), he->opposite(), tpole2)) {
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
is_on_new_perimetric_face_boundary
(const Halfedge *prev1,
 const Halfedge *prev2,
 const X_monotone_curve_2& cv) const
{
    // status: check correctness of implementation
    std::cout << "Arr_torus_topology_traits_2::" 
              << "is_on_new_perimetric_face_boundary"
              << std::endl;

    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());
    
    Identification_crossing leftmost_crossing;
    Identification_crossing bottommost_crossing;

    CGAL_assertion_code(bool tpole);
    CGAL_assertion(_is_perimetric_path(prev2, prev1, tpole));
    
    
    // maintain the invariant that the pole is always in the top_face,
    // i.e, it is the face that contains everything and has now outer ccb
    // If pole is part of a ccb itself, it incident face is the face that 
    // contains everything.
    
    if (prev1->inner_ccb()->face() != top_face()) {
        // actual answer does not matter, as a perimetric face not containing
        // the pole of the torus is split into two perimetric faces
        // it is not important whether prev1 or prev2 belongs to it
        return false;
    }
    
    bool touching = false;
    bool crossing = false;
    
    Vertex *leftmost;
    Vertex *bottommost;

    std::pair< int, int > counters = 
        _crossings_with_identifications(prev2, prev1, touching, crossing,
                                        leftmost, leftmost_crossing, 
                                        bottommost, bottommost_crossing);
    
    CGAL_assertion(crossing);
    
    int x_counter = counters.first;
    int y_counter = counters.second;
    
    if (x_counter % 2 != 0) {

        CGAL_assertion(bottommost != NULL);
        
        bool smallest_perimetric = true;
        
        typename Identification_WE::iterator it = 
            this->_m_identification_WE.find(bottommost->point());
        
        while (smallest_perimetric && 
               it != this->_m_identification_WE.begin()) {
            it--;
            Vertex *v = it->second;
            smallest_perimetric = !_has_crossing_perimetric_path(v);
        }
        
        if (smallest_perimetric) {
            return (bottommost_crossing == BEFORE_TO_AFTER);
        } else {
            return (bottommost_crossing == AFTER_TO_BEFORE);
        }
    } else if (y_counter % 2 != 0) {
      
        CGAL_assertion(leftmost != NULL);

        bool smallest_perimetric = true;
        
        typename Identification_NS::iterator it = 
            this->_m_identification_NS.find(leftmost->point());
        
        while (smallest_perimetric && 
               it != this->_m_identification_NS.begin()) {
            it--;
            Vertex *v = it->second;
            smallest_perimetric = !_has_crossing_perimetric_path(v);
        }
        
        if (smallest_perimetric) {
            return (leftmost_crossing == AFTER_TO_BEFORE);
        } else {
            return (leftmost_crossing == BEFORE_TO_AFTER);
        }
    } 
    // else 
    CGAL_assertion(touching && x_counter == 0 && y_counter == 0);

    return (bottommost_crossing == BEFORE_TO_AFTER);
}

//-----------------------------------------------------------------------------
// checks whether halfedges are boundaries of the same face
//
template <class GeomTraits, class Dcel_>
bool
Arr_torus_topology_traits_2<GeomTraits,Dcel_>::boundaries_of_same_face
(const Halfedge *e1,
 const Halfedge *e2) const
{
    // status: check correctness of implementation
    std::cout << " Arr_torus_topology_traits_2::boundaries_of_same_face" 
              << std::endl;
    // This predicate is only used for case 3.3.2 of the insertion process
    
    // first check it for e2 ...
    Identification_crossing leftmost_crossing2;
    Identification_crossing bottommost_crossing2;
    
    bool touching2 = false;
    bool crossing2 = false;
    
    Vertex *leftmost2;
    Vertex *bottommost2;

    std::pair< int, int > counters2 = 
        _crossings_with_identifications(e2, e2, touching2, crossing2,
                                        leftmost2, leftmost_crossing2, 
                                        bottommost2, bottommost_crossing2);
    
    CGAL_assertion_code(bool tpole2;);
    if (!_is_perimetric_data(counters2, touching2, crossing2, tpole2)) {
        // ... as there is a quick exit
        return false;
    }
    
    // if e2 is perimetric e1 must be as well
    
    int x_counter2 = counters2.first;
    int y_counter2 = counters2.second;
    
    Identification_crossing leftmost_crossing1;
    Identification_crossing bottommost_crossing1;
    
    bool touching1 = false;
    bool crossing1 = false;

    Vertex *leftmost1;
    Vertex *bottommost1;

    std::pair< int, int > counters1 = 
        _crossings_with_identifications(e1, e1, touching1, crossing1,
                                        leftmost1, leftmost_crossing1, 
                                        bottommost1, bottommost_crossing1);
    
    CGAL_assertion_code(bool tpole1;);
    CGAL_assertion(
            _is_perimetric_data(counters1, touching1, crossing1, tpole1)
    );
    
    int x_counter1 = counters1.first;
    int y_counter1 = counters1.second;
    
    if (x_counter1 % 2 != 0) {
        CGAL_assertion(crossing1);
        CGAL_assertion(crossing2);
        CGAL_assertion(x_counter2 % 2 != 0);
        return (bottommost_crossing1 != bottommost_crossing2);
    } else if (y_counter1 % 2 != 0) {
        CGAL_assertion(crossing1);
        CGAL_assertion(crossing2);
        CGAL_assertion(y_counter2 % 2 != 0);
        return (leftmost_crossing1 != leftmost_crossing2);
    } 
    // else 
    
    CGAL_assertion(touching1 && x_counter1 == 0 && y_counter1 == 0 &&
                   touching2 && x_counter2 == 0 && y_counter2 == 0);
    
    CGAL_assertion(crossing1);
    CGAL_assertion(crossing2);
    CGAL_assertion((leftmost_crossing1 == AFTER_TO_BEFORE && 
                    bottommost_crossing1 == BEFORE_TO_AFTER) ||
                   (leftmost_crossing1 == BEFORE_TO_AFTER && 
                    bottommost_crossing1 == AFTER_TO_BEFORE));
    CGAL_assertion((leftmost_crossing2 == AFTER_TO_BEFORE && 
                    bottommost_crossing2 == BEFORE_TO_AFTER) ||
                   (leftmost_crossing2 == BEFORE_TO_AFTER && 
                    bottommost_crossing2 == AFTER_TO_BEFORE));
    
    bool res = (bottommost_crossing1 != bottommost_crossing2);
    CGAL_assertion(res == (leftmost_crossing1 != leftmost_crossing2));
    return res;
}

//-----------------------------------------------------------------------------
// Determine whether the given vertex lies in the interior of the given face.
//
template <class GeomTraits, class Dcel_>
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::is_in_face
(const Face *f, const Point_2& p, const Vertex *v) const
{
    // status: not implemented
    std::cout << "TODO: Arr_torus_topology_traits_2::is_in_face" 
              << std::endl;

    CGAL_assertion(false); // is_in_face not implemented for torus
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
    
    std::cout << "Arr_torus_topology_traits_2 is_redundant"  << std::endl;
    CGAL_precondition(_valid(v->boundary_in_x(),v->boundary_in_y()));
    
    // if there are not incident edges just remove it
    // TASK: check whether isolated or degree == 0 is needed!
    return (v->is_isolated());
}

//-----------------------------------------------------------------------------
// Determine whether a boundary vertex is redundant
//
template <class GeomTraits, class Dcel_>
typename Arr_torus_topology_traits_2<GeomTraits, Dcel_>::Halfedge* 
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::erase_redundant_vertex
(Vertex *v) 
{
    // status: correct
    
    std::cout << "Arr_torus_topology_traits_2 erase_redundant_vertex"  
              << std::endl;
    CGAL_precondition(_valid(v->boundary_in_x(),v->boundary_in_y()));
    
    // no incident curve-end can give us the key
    // -> but we stored something else useful: find iterator
    if (v->boundary_in_x() != CGAL::NO_BOUNDARY) {

        typename Vertices_on_identification_WE::iterator 
            vit = _m_vertices_on_identification_WE.find(v);
        
        // and delete this item
        _m_identification_WE.erase(vit->second);
        _m_vertices_on_identification_WE.erase(vit);
        
    } else {
        CGAL_assertion(v->boundary_in_y() != CGAL::NO_BOUNDARY);
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

//-----------------------------------------------------------------------------
// Number of crossing with the curve of identification
//
template <class GeomTraits, class Dcel_>
std::pair< int, int >
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_crossings_with_identifications(
        const Halfedge* he1, const Halfedge* he2, 
        bool& touching, bool& crossing, 
        const Vertex *leftmost_vertex, 
        Identification_crossing& leftmost_crossing,
        const Vertex *bottommost_vertex, 
        Identification_crossing& bottommost_crossing) const {
    
    // status: check implementation
    
    std::cout << "Arr_torus_topology_traits: "
              << "_crossings_with_identifications" << std::endl;

    touching = false;
    crossing = false;

    leftmost_vertex = NULL;
    bottommost_vertex = NULL;
    
    const Vertex *leftmost_tvertex = NULL;
    const Vertex *bottommost_tvertex = NULL;
    
    Point_2_less_NS less_ns(_m_traits);
    Point_2_less_WE less_we(_m_traits);
    
    int x_counter = 0;
    int y_counter = 0;

    if (he1->next() == he2 && he2->next () == he1) {
        return std::make_pair(x_counter, y_counter);
    }

    typename Traits_adaptor_2::Boundary_in_x_2 boundary_in_x =
        _m_traits->boundary_in_x_2_object();
    typename Traits_adaptor_2::Boundary_in_y_2 boundary_in_y =
        _m_traits->boundary_in_y_2_object();
    
    // Start with the next of prev1:
    const Halfedge * curr = he1->next();
    // Save its src condition
    Curve_end curr_src_ind;
    Curve_end curr_trg_ind;
    if (curr->direction() == LEFT_TO_RIGHT) {
        curr_src_ind = MIN_END;
        curr_trg_ind = MAX_END;
    } else {
        curr_src_ind = MAX_END;
        curr_trg_ind = MIN_END;
    }
    CGAL_assertion(!curr->has_null_curve());
    Boundary_type first_src_bcx = boundary_in_x(curr->curve(), curr_src_ind);
    Boundary_type curr_trg_bcx = boundary_in_x(curr->curve(), curr_trg_ind);  
    Boundary_type first_src_bcy = boundary_in_y(curr->curve(), curr_src_ind);
    Boundary_type curr_trg_bcy = boundary_in_y(curr->curve(), curr_trg_ind);  
    while (curr != he2) {
        const Halfedge * next = curr->next();
        
        Curve_end next_src_ind;
        Curve_end next_trg_ind;
        if (next->direction() == LEFT_TO_RIGHT) {
            next_src_ind = MIN_END;
            next_trg_ind = MAX_END;
        } else {
            next_src_ind = MAX_END;
            next_trg_ind = MIN_END;
        }
        Boundary_type next_src_bcx = 
            boundary_in_x(next->curve(), next_src_ind);
        Boundary_type next_trg_bcx = 
            boundary_in_x(next->curve(), next_trg_ind);
        Boundary_type next_src_bcy = 
            boundary_in_y(next->curve(), next_src_ind);
        Boundary_type next_trg_bcy = 
            boundary_in_y(next->curve(), next_trg_ind);
        if (curr_trg_bcx != CGAL::NO_BOUNDARY) {
            if (bottommost_tvertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_ident
                less_we(curr->vertex()->point(),
                        bottommost_tvertex->point())) {
                bottommost_tvertex = curr->vertex();
            }
        }
        if (curr_trg_bcy != CGAL::NO_BOUNDARY) {
            if (leftmost_tvertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_ident
                less_ns(curr->vertex()->point(),
                        leftmost_tvertex->point())) {
                leftmost_tvertex = curr->vertex();
            }
        }
        if (curr_trg_bcx != next_src_bcx) {
            CGAL_assertion(curr_trg_bcx != CGAL::NO_BOUNDARY);
            CGAL_assertion(next_src_bcx != CGAL::NO_BOUNDARY);
            if (curr_trg_bcx == BEFORE_DISCONTINUITY) {
                if (bottommost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_we(curr->vertex()->point(),
                            bottommost_vertex->point())) {
                    bottommost_vertex = curr->vertex();
                    bottommost_crossing = BEFORE_TO_AFTER;
                }
                ++x_counter;
            } else {
                if (bottommost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_we(curr->vertex()->point(),
                            bottommost_vertex->point())) {
                    bottommost_vertex = curr->vertex();
                    bottommost_crossing = AFTER_TO_BEFORE;
                }
                --x_counter;
            }
        }
        if (curr_trg_bcy != next_src_bcy) {
            CGAL_assertion(curr_trg_bcy != CGAL::NO_BOUNDARY);
            CGAL_assertion(next_src_bcy != CGAL::NO_BOUNDARY);
            if (curr_trg_bcy == BEFORE_DISCONTINUITY) {
                if (leftmost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_ns(curr->vertex()->point(),
                            leftmost_vertex->point())) {
                    leftmost_vertex = curr->vertex();
                    leftmost_crossing = BEFORE_TO_AFTER;
                }
                ++y_counter;
            } else {
                if (leftmost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_ns(curr->vertex()->point(),
                            leftmost_vertex->point())) {
                    leftmost_vertex = curr->vertex();
                    leftmost_crossing = AFTER_TO_BEFORE;
                }
                --y_counter;
            }
        }
        curr = next;
        curr_trg_bcx = next_trg_bcx;
        curr_trg_bcy = next_trg_bcy;
    }
    if (he1 == he2) {
        Boundary_type last_trg_bcx = curr_trg_bcx;
        Boundary_type last_trg_bcy = curr_trg_bcy;
        if (last_trg_bcx != CGAL::NO_BOUNDARY) {
            if (bottommost_tvertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_ident
                less_we(curr->vertex()->point(),
                        bottommost_tvertex->point())) {
                bottommost_tvertex = curr->vertex();
            }
        }
        if (last_trg_bcy != CGAL::NO_BOUNDARY) {
            if (leftmost_tvertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_ident
                less_ns(curr->vertex()->point(),
                        leftmost_tvertex->point())) {
                leftmost_tvertex = curr->vertex();
            }
        }
        if (last_trg_bcx != first_src_bcx) {
            //CGAL_assertion(last_trg_bcx != CGAL::NO_BOUNDARY);
            //CGAL_assertion(first_src_bcx != CGAL::NO_BOUNDARY);
            if (last_trg_bcx == BEFORE_DISCONTINUITY) {
                if (bottommost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_we(curr->vertex()->point(),
                            bottommost_vertex->point())) {
                    bottommost_vertex = curr->vertex();
                    bottommost_crossing = BEFORE_TO_AFTER;
                }
                ++x_counter;
            } else {
                if (bottommost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_we(curr->vertex()->point(),
                            bottommost_vertex->point())) {
                    bottommost_vertex = curr->vertex();
                    bottommost_crossing = AFTER_TO_BEFORE;
                }
                --x_counter;
            }
        }
        if (last_trg_bcy != first_src_bcy) {
            //CGAL_assertion(last_trg_bcy != CGAL::NO_BOUNDARY);
            //CGAL_assertion(first_src_bcy != CGAL::NO_BOUNDARY);
            if (last_trg_bcy == BEFORE_DISCONTINUITY) {
                if (leftmost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_ns(curr->vertex()->point(),
                            leftmost_vertex->point())) {
                    leftmost_vertex = curr->vertex();
                    leftmost_crossing = BEFORE_TO_AFTER;
                }
                ++y_counter;
            } else {
                if (leftmost_vertex == NULL || 
                    // TASK avoid real comparisons, ask _m_vertices_on_ident
                    less_ns(curr->vertex()->point(),
                            leftmost_vertex->point())) {
                    leftmost_vertex = curr->vertex();
                    leftmost_crossing = AFTER_TO_BEFORE;
                }
                --y_counter;
            }
        }
    }
    
    if (leftmost_vertex != NULL || bottommost_vertex != NULL) {
        crossing = true;
    } else {
        if (leftmost_tvertex != NULL || bottommost_tvertex != NULL) {
            touching = true;
            if (leftmost_vertex == NULL) {
                leftmost_vertex = leftmost_tvertex;
            }
            if (bottommost_vertex == NULL) {
                bottommost_vertex = bottommost_tvertex;
            }
        }
    }

    return (std::make_pair(x_counter, y_counter));
}

CGAL_END_NAMESPACE

#endif
// EOF
