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
 
    // TODO assign the vertices and the face correctly
    _m_identification_WE = other._m_identification_WE;
    _m_identification_NS = other._m_identification_NS;
    
    _m_f_top = other._m_f_top;
    
    CGAL_assertion(_m_f_top != NULL);

#if 0
    // Go over the DCEL vertices and locate all points with boundary condition
    typename Dcel::Vertex_iterator       vit;
    //Boundary_type                        bx, by;
    //Halfedge                            *first_he, *next_he;

    for (vit = this->_m_dcel.vertices_begin(); vit != 
             this->_m_dcel.vertices_end(); ++vit)
    {
        if (vit->boundary_in_x() != CGAL::NO_BOUNDARY) {
            
            // fill _m_identification_WE
            std::pair< typename Identification_WE::iterator, bool > p =
                _m_identification_WE.insert(
                        // each point on loc is concrete
                        std::make_pair(vit->f, &(*vit))
                );
            CGAL_assertion(!p.second);
            
            // vertices on idendification
            _m_vertices_on_identification_WE[(*vit)] = p.first;
            
        } else {
            if (vit->boundary_in_y() != CGAL::NO_BOUNDARY) {
                // fill _m_identification_NS
                std::pair< typename Identification_NS::iterator, bool > p =
                    _m_identification_NS.insert(
                            // each point on loc is concrete
                            std::make_pair(vit->f, &(*vit))
                    );
                CGAL_assertion(!p.second);
                
                // vertices on idendification
                _m_vertices_on_identification_NS[(*vit)] = p.first;
            }
        }
    }
#endif
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
        
        // update the local structure for points on the line of discontinuity
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
    
    // update the local structure for points on the line of discontinuity
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
 const Halfedge *e2) const
{
    // status: check implementation
    std::cout << "TODO: Arr_torus_topology_traits_2::is_perimetric_path" 
              << std::endl;
    
    int x_counter = 0;
    int y_counter = 0;
    typename Traits_adaptor_2::Boundary_in_x_2 boundary_in_x =
        _m_traits->boundary_in_x_2_object();
    typename Traits_adaptor_2::Boundary_in_y_2 boundary_in_y =
        _m_traits->boundary_in_y_2_object();
    
    // Start with the next of prev1:
    const Halfedge * curr = e1->next();
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
    Boundary_type first_src_bcx = boundary_in_x(curr->curve(), curr_src_ind);
    Boundary_type curr_trg_bcx = boundary_in_x(curr->curve(), curr_trg_ind);  
    Boundary_type first_src_bcy = boundary_in_y(curr->curve(), curr_src_ind);
    Boundary_type curr_trg_bcy = boundary_in_y(curr->curve(), curr_trg_ind);  
    while (curr != e2) {
        std::cout << "curve: " << curr->curve() << std::endl;
        const Halfedge * next = curr->next();
        std::cout << "next: " << next->curve() << std::endl;
        
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
        if (curr_trg_bcx != next_src_bcx) {
            CGAL_assertion(curr_trg_bcx != CGAL::NO_BOUNDARY);
            CGAL_assertion(next_src_bcx != CGAL::NO_BOUNDARY);
            if (curr_trg_bcx == BEFORE_DISCONTINUITY) {
                ++x_counter;
            } else {
                --x_counter;
            }
        }
        if (curr_trg_bcy != next_src_bcy) {
            CGAL_assertion(curr_trg_bcy != CGAL::NO_BOUNDARY);
            CGAL_assertion(next_src_bcy != CGAL::NO_BOUNDARY);
            if (curr_trg_bcy == BEFORE_DISCONTINUITY) {
                ++y_counter;
            } else {
                --y_counter;
            }
        }
        curr = next;
        curr_trg_bcx = next_trg_bcx;
        curr_trg_bcy = next_trg_bcy;
    }
    std::cout << "curve: " << curr->curve() << std::endl;
    Boundary_type last_trg_bcx = curr_trg_bcx;
    Boundary_type last_trg_bcy = curr_trg_bcy;
    if (last_trg_bcx != first_src_bcx) {
        CGAL_assertion(last_trg_bcx != CGAL::NO_BOUNDARY);
        CGAL_assertion(first_src_bcx != CGAL::NO_BOUNDARY);
        if (last_trg_bcx == BEFORE_DISCONTINUITY) {
            ++x_counter;
        } else {
            --x_counter;
        }
    }
    if (last_trg_bcy != first_src_bcy) {
        std::cout << "last_trg_bcy: " << last_trg_bcy << std::endl;
        std::cout << "first_src_bcy: " << first_src_bcy << std::endl;
        //CGAL_assertion(last_trg_bcy != CGAL::NO_BOUNDARY);
        //CGAL_assertion(first_src_bcy != CGAL::NO_BOUNDARY);
        if (last_trg_bcy == BEFORE_DISCONTINUITY) {
            ++y_counter;
        } else {
            --y_counter;
        }
    }
    // Path must be perimetric:
    return (x_counter % 2 == 1 && y_counter == 0 ||
            x_counter == 0 && y_counter % 2 == 0 ||
            x_counter == 0 && y_counter == 0);
    
#if 0
    Identification_crossing tmp;
    std::pair< unsigned int, unsigned int > crossings;
    // TODO check order
    if (e1 == e2) {
        crossings = _crossings_with_identification_NS(e1, NULL, tmp);
    } else {
        crossings = _crossings_with_identification_NS(e1, e2, tmp);
    }
    // return whether there has been an odd number of intersections
    bool res = ((crossings.first + crossings.second) % 2 == 1);
    //std::cout << "res: " << res << std::endl;
    return res;
#endif
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
    // status: check implementation
    std::cout << "Arr_torus_topology_traits_2 face_split"  << std::endl;


    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());
    
    // TODO face_split for torus

    // In case of a quadric topology, connecting two vertices on the same
    // inner CCB always causes a face split. We need to check if at least one
    // of the paths from prev1's target to prev2's target, or from prev2's to
    // prev1's target is perimetric. If so, we split two adjacent faces.
    // If both are not perimetric, then the split face becomes a hole in the
    // original face.
    bool    face_split = true;
    bool    is_hole = (! _is_perimetric_path (prev1->next(), prev2->next()) &&
                       ! _is_perimetric_path (prev2->next(), prev1->next()));

    return (std::make_pair (face_split, is_hole));
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

    CGAL_precondition (! he->is_on_inner_ccb());
    CGAL_precondition (! he->opposite()->is_on_inner_ccb());

    // TODO hole_creation for torus

    // Check whether the halfedge and its twin belong to the same outer CCB
    // (and are therefore incident to the same face).
    if (he->outer_ccb() == he->opposite()->outer_ccb())
    {
      // Check the two cycles that will be created once we remove he and its
      // twin (from he->next() to he's twin, not inclusive, and from the
      // successor of he's twin to he, not inclusive).
      if (_is_perimetric_path (he->next(), he->opposite()) &&
          _is_perimetric_path (he->opposite()->next(), he))
      {
        // Both paths are perimetric, so the two cycles become to separate
        // outer CCBs of the same face, and no hole is created.
        return (false);
      }
      else
      {
        // At least one cyclic path is non-perimetic. This cycle will become
        // an inner CCB representing a hole in the face.
        return (true);
      }
    }
    else
    {
      // The edge to be removed separates two faces.
      // Check the cyclic path from he and back, and from its twin and back.
      if (_is_perimetric_path (he, he) &&
          _is_perimetric_path (he->opposite(), he->opposite()))
      {
        // In this case we disconnect a perimetric cycle around the quadric,
        // causing two perimetric faces to merge. The remainder of the cycle
        // becomes an inner CCB (a hole) in the merged face.
        return (true);
      }
      else
      {
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
    // status: check implementation
    std::cout << "TODO: Arr_torus_topology_traits_2::" 
              << "is_on_new_perimetric_face_boundary"
              << std::endl;
#if 0
    
    // can only be possible for a paraboloid.
    if ((_m_left == CGAL::MINUS_INFINITY && _m_right == PLUS_INFINITY) ||
        (_m_left == CGAL::AFTER_SINGULARITY && _m_right == BEFORE_SINGULARITY)) {
        return false;
    }
    CGAL_assertion(_m_quadric.is_elliptic_paraboloid());
    
    // so we use leftmost disconti-crossing of 
    // path ending in prev1 (and starting in prev2)
    CGAL_assertion(is_perimetric_path(prev2, prev1));
    Identification_crossing leftmost;
    std::pair< unsigned int, unsigned int > crossings =
        _crossings_with_line_of_discontinuity(prev2, prev1, leftmost);
    // is perimetric test
    CGAL_assertion((crossings.first + crossings.second) % 2 == 1);

    // to check whether new face contains singular point
    // and therefore prev1 belongs to the outer_ccb of this new face
    if (leftmost == BEFORE_TO_AFTER) { 
        // the face is "on the left" side of the surfacs
        return (_m_left == CGAL::AFTER_SINGULARITY);
    } else {
        CGAL_assertion(leftmost == AFTER_TO_BEFORE);
        // the face is "on the rightt" side of the surfacs
        return (_m_right == CGAL::BEFORE_SINGULARITY);
    }
#endif
    CGAL_assertion(false);
    return false;
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
    // status: check implementation
    std::cout << "TODO: Arr_torus_topology_traits_2::boundaries_of_same_face" 
              << std::endl;
    // This predicate is only used for case 3.3.2 of the insertion process
    
    // TODO check whether asympotote of identification crossings are equal

#if 0
    _m_boundaries_of_same_face++;
    std::cout << "mpWE: " << _m_boundaries_of_same_face << std::endl;
    if (_m_boundaries_of_same_face < 2) {
        std::cout << "return true" << std::endl;
        return true;
    }
    // else
    // otherwise return true
    std::cout << "return false" << std::endl;
#endif
    return false;
    

#if 0
    Identification_crossing leftmost1;
    std::pair< unsigned int, unsigned int > crossings1 =
        _crossings_with_identification_NS(e1, NULL, leftmost1);
    
    Identification_crossing leftmost2;
    std::pair< unsigned int, unsigned int > crossings2 =
        _crossings_with_identification_NS(e2, NULL, leftmost2);
    
    return (leftmost1 != leftmost2);
#endif
}

//-----------------------------------------------------------------------------
// Determine whether the given vertex lies in the interior of the given face.
//
template <class GeomTraits, class Dcel_>
bool Arr_torus_topology_traits_2<GeomTraits, Dcel_>::is_in_face
(const Face *f, const Point_2& p, const Vertex *v) const
{
    // status: check implementation
    // TODO is_in_face NEEDED for incremental insertion
    std::cout << "TODO: Arr_torus_topology_traits_2::is_in_face" 
              << std::endl;
#if 0
    CGAL_precondition (v == NULL || ! v->has_null_point());
    CGAL_precondition (v == NULL || 
                       _m_traits->equal_2_object()(p, v->point()));
    
    // In case the face is unbounded and has no outer ccbs, this is the single
    // (un)bounded face of an arrangement of bounded curves. 
    // This face obviously contains any point in its interior.
    int nu_m_occb = f->number_of_outer_ccbs();
    CGAL_assertion(num_occb <= 2);
    if (num_occb == 0) {
        std::cout << "IS_IN_FACE: number_of_outer_ccbs == 0" << std::endl;
        return (true);
    }
    
    // if f would be a perimetric face than it has at least one outer boundary
    
    Outer_ccb_const_iterator occb_it = f->outer_ccbs_begin();

    bool perimetic_face = false;

    if (num_occb == 1) {
        // then face contains v_left/v_right if its outer_ccb is perimetric
        
        if (is_perimetric_path(*occb_it, *occb_it)) {
            std::cout << "IS_IN_FACE: perimetric boundary face" << std::endl;
            perimetic_face = true;
        } else {
            std::cout << "IS_IN_FACE: normal face" << std::endl;
        }
        // otherwise, it is a normal face (maybe going over the curve of disc)
        
    }

    if (num_occb == 2) {
        // both outer_ccbs should be perimetric
        
        CGAL_assertion(is_perimetric_path(*occb_it, *occb_it));
        occb_it++;
        CGAL_assertion(is_perimetric_path(*occb_it, *occb_it));
        std::cout << "IS_IN_FACE: perimetric middle face" << std::endl;
        perimetic_face = true;
    }
    
    
    if (perimetic_face) {
        // iterate over all halfedges of outer_ccbs and locate the
        // one that is closest to p
        
        Halfedge *closest;
        closest = NULL;
        
        const Halfedge    *first = *(f->outer_ccbs_begin());
        const Halfedge    *curr = first;
        for(;curr != first; curr++) {

            // jump over antennas
            if (! curr->opposite()->is_on_inner_ccb() &&
                curr->outer_ccb()->face() == 
                curr->opposite()->outer_ccb()->face()) {
                if (curr == first || curr->next() == first) {
                    break;
                }
                // we skip, this and the next halfedge
                curr = curr->next()->next();
                continue;
            }
            


            if (this->_m_traits->is_in_x_range_2_object() (
                        curr->curve(), p
                ) &&
                this->_m_traits->compare_y_at_x_2_object() (
                        p, curr->curve()) == CGAL::SMALLER) {
                if (closest == NULL) {
                    *closest = *curr;
                } else {
                    // compare with closest
                    if (this->_m_traits->compare_y_position_2_object() (
                                curr->curve(), closest->curve()
                        ) == CGAL::SMALLER) {
                        *closest = *curr;
                    }
                }
            }
        }
        
        // if no halfedge found
        if (closest == NULL) {
            std::cout << "No halfedge found" << std::endl;
            // search the line of discontinuity for right most point
            // whose x is smaller than p's x
            typename Line_of_discontinuity::const_iterator
                lower_bound,
                it = curve_ends_and_vertices_on_line_of_discontinuity_begin();
            
            // TODO do not store points, do store x-coordinates
            bool lower_bound_valid = false;
            while (it != 
                   curve_ends_and_vertices_on_line_of_discontinuity_end()) {
                
                if (it->first.x() < p.x()) {
                    lower_bound = it;
                    lower_bound_valid = true;
                }
                if (it->first.x() > p.x()) {
                    break;
                }
                it++;
            }
            // TODO check whether valid!!!
            // TODO select adjacent edge
            if (lower_bound_valid) {
                *closest = *lower_bound->second->halfedge();
                return true; // TODO closest->face() == f;
            }
            // TODO check
            return true;

        } 
        
        CGAL_assertion(closest != NULL);

        return (closest->direction() == CGAL::RIGHT_TO_LEFT);
        
    }



    // Keep two counter that store how many segments are smaller/larger than
    // p, i.e. intersecting the rays shootet in vertical direction from p
    // emanating from p (except for some degenerate cases that are
    // explained below).
    unsigned int seg_smaller = 0;
    unsigned int seg_larger = 0;
    
    // Get a halfedge along the outer CCB of the given face, go over all curves
    // of the boundary component, and count those which are above p.
    // We begin by comparing p to the source vertex of the first halfedge.
    // Note that if p coincides with this vertex, p is obviously not in the
    // interior of the component.
    const Halfedge    *first = *(f->outer_ccbs_begin());
    const Halfedge    *curr = first;
    Comparison_result  res_source;
    Comparison_result  res_target;
    Comparison_result  res_y_at_x;
    
    if (curr->opposite()->vertex() == v) {
        return (false);
    }
    
    CGAL::Boundary_type bound_x = curr->opposite()->vertex()->boundary_in_x();
    if (bound_x < 0) {
        res_source = CGAL::LARGER;
    } else if (bound_x > 0) {
        res_source = CGAL::SMALLER;
    } else {
        res_source = this->_m_traits->compare_xy_2_object() (
                p, curr->opposite()->vertex()->point()
        );
    }
    CGAL_assertion(res_source != CGAL::EQUAL);
    
    CGAL::Boundary_type last_by = CGAL::NO_BOUNDARY;
    CGAL::Boundary_type curr_by = CGAL::NO_BOUNDARY;
    

    do {
#if 0
        // first jump over antennas
        // In case the current halfedge belongs to an "antenna", namely its
        // incident face is the same as its twin's, we can simply skip it
        // (in order not to count it twice).
        if (! curr->opposite()->is_on_inner_ccb() &&
            curr->outer_ccb()->face() == 
            curr->opposite()->outer_ccb()->face()) {
            if (curr == first || curr->next() == first) {
                break;
            }
            // we skip, this and the next halfedge
            curr = curr->next()->next();
            continue;
        }
#endif
        
        // Compare p to the target vertex of the current halfedge.
        // If the vertex v associated with p (if v is given and is not NULL)
        // on the boundary of the component, p is obviously not in the interior
        // the component.
        if (curr->vertex() == v) {
            std::cout << "IS_IN_FACE: false" << std::endl;
            return (false);
        }
        
        bound_x = curr->vertex()->boundary_in_x();
        if (bound_x < 0) {
            res_target = CGAL::LARGER;
        } else if (bound_x > 0) {
            res_target = CGAL::SMALLER;
        } else {
            res_target = this->_m_traits->compare_xy_2_object()
                (p, curr->vertex()->point());
        }
        CGAL_assertion(res_target != CGAL::EQUAL);  
        
        // read the boundary_type at src of curr
        curr_by = this->_m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MIN_END : CGAL::MAX_END)
        );
        
#if 1
        if ((last_by == CGAL::AFTER_DISCONTINUITY &&
             curr_by == CGAL::BEFORE_DISCONTINUITY) ||
            (last_by == CGAL::BEFORE_DISCONTINUITY &&
             curr_by == CGAL::AFTER_DISCONTINUITY)) {
            // "jumped over the line of discontinuity"
            std::swap(seg_smaller, seg_larger);
        }
#endif
#if 0
        if (curr->curve().get_level() != curr->next()->curve().get_level()) {
            // "jumped over the line of discontinuity"
            std::swap(seg_smaller, seg_larger);
            std::cout << "swap: less=" << seg_smaller << ", greater=" << seg_larger << std::endl;
        }
#endif        

        // Check that if we shoot a "tilted" vertical ray from p upward
        // (by "tilted" we mean the angle it forms with the x-axis is
        //  PI/2 + epsilon, where epsilon is arbitrarily small), then we hit
        // the x-monotone curve associated with curr once.
        if (res_source != res_target) {
            // if p is in x-range of the curr halfedge
            
            res_y_at_x = compare_y_at_x (p, curr);
            
            std::cout << "pt: " << p << std::endl;
            std::cout << "curr: " << curr->curve() << std::endl;
            std::cout << "dir: " << (curr->direction() == CGAL::LEFT_TO_RIGHT ?
                                     "LEFT_TO_RIGHT" : "RIGHT_TO_LEFT") 
                      << std::endl;
            if (res_y_at_x == CGAL::SMALLER) {
                seg_smaller++;
            } else if (res_y_at_x == CGAL::LARGER) {
                seg_larger++;
            } else {
                CGAL_assertion(res_y_at_x == CGAL::EQUAL);
                // In this case p lies on the current edge, 
                // so it is obviously not
                // contained in the interior of the component.
                std::cout << "IS_IN_FACE: false" << std::endl;
                return (false);
            }
        }
            
        // store as last_by the boundary at target of curr
        last_by = this->_m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MAX_END : CGAL::MIN_END)
        );
        
        // Proceed to the next halfedge along the component boundary.
        // Note that the source vertex of this halfedge is the current target.
        curr = curr->next();
        res_source = res_target;
        
    } while (curr != first);
    
    // The query point lies inside the connected components if and only if the
    // ray we shoot from it intersects the boundary an odd number of time.
    std::cout << "smaller: " << seg_smaller << std::endl;
    std::cout << "larger: " << seg_larger << std::endl;
    std::cout << "is_in_face: " << ((seg_smaller % 2) != 0) << std::endl;
    std::cout << "IS_IN_FACE: " << ((seg_smaller % 2) != 0) << std::endl;
    return ((seg_smaller % 2) != 0);
#endif
    CGAL_assertion(false);
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

    // std::cout << "Arr_torus_topology_traits_2 erase_redundant_vertex"  
    // << std::endl;
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


// TODO required??!
//-----------------------------------------------------------------------------
// Number of crossing with the line of discontiniuty
//
template <class GeomTraits, class Dcel_>
std::pair< unsigned int, unsigned int >
Arr_torus_topology_traits_2<GeomTraits, Dcel_>::
_crossings_with_identification_NS(
        const Halfedge* he1, const Halfedge* he2, 
        Identification_crossing& leftmost) const {

    std::cout << "TODO: Arr_torus_topology_traits: "
              << "_crossings_with_identification_NS" << std::endl;
    
    CGAL::Boundary_type thistgt_by = CGAL::NO_BOUNDARY;
    CGAL::Boundary_type nextsrc_by = CGAL::NO_BOUNDARY;
    
    const Halfedge *curr = he1;
    const Vertex *leftmost_vertex = NULL;
    
    // we count the number of crossings with the line of disc
    unsigned int n_crossings_before_to_after = 0;
    unsigned int n_crossings_after_to_before = 0;

    if (he2 == NULL) {
        // also check prev()->tgt with curr->src()
        // read the boundary_type at tgt of curr
        thistgt_by = this->_m_traits->boundary_in_y_2_object()(
                curr->prev()->curve(), 
                (curr->prev()->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MAX_END : CGAL::MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->_m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MIN_END : CGAL::MAX_END)
        );
        
        if (thistgt_by == CGAL::AFTER_DISCONTINUITY &&
            nextsrc_by == CGAL::BEFORE_DISCONTINUITY) {
#if 0
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_lod
                Vertex_key_less_NS(_m_traits)(curr->vertex()->point(),
                                              leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
#endif
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::BEFORE_DISCONTINUITY &&
            nextsrc_by == CGAL::AFTER_DISCONTINUITY) {
#if 0
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_lod
                Vertex_key_less_NS(_m_traits)(curr->vertex()->point(),
                                              leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = BEFORE_TO_AFTER;
            }
#endif
            n_crossings_before_to_after++;
        }
    } 
    
    const Halfedge *last = (he2 == NULL ? he1 : he2);

    do {
        // note that boundary conditions at beginning vertex of path
        // and at end vertex of task do not count, since they offer
        // always possibilities to connect the path on both sides of the
        // discontinuity. 

#if 0 
        // TASK antennas are not so critical ... 
        // since counted an even number of times, therefore not changing the
        // result
        // does only make sense to avoid checking edges, but makes handling
        // of loop more complicated

        
        // first jump over antennas (the following CODE is NOT tested)
        // In case the current halfedge belongs to an "antenna", namely its
        // incident face is the same as its twin's, we can simply skip it
        // (in order not to count it twice).
        if ((curr->in_on_inner_ccb() && 
             curr->inner_ccb()->face() == 
             curr->opposite()->inner_ccb()->face())
            ||
            (curr->in_on_outer_ccb() && 
             curr->outer_ccb()->face() == 
             curr->opposite()->outer_ccb()->face())) {
            if (curr == first || curr->next() == first) {
                break;
            }
            // we skip, this and the next halfedge
            curr = curr->next()->next();
            continue;
        }
#endif
        
        // read the boundary_type at tgt of curr
        thistgt_by = this->_m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MAX_END : CGAL::MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->_m_traits->boundary_in_y_2_object()(
                curr->next()->curve(), 
                (curr->next()->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MIN_END : CGAL::MAX_END)
        );
        
        if (thistgt_by == CGAL::AFTER_DISCONTINUITY &&
            nextsrc_by == CGAL::BEFORE_DISCONTINUITY) {
#if 0
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_lod
                Point_2_NS_less(_m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
#endif
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::BEFORE_DISCONTINUITY &&
            nextsrc_by == CGAL::AFTER_DISCONTINUITY) {
#if 0
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask _m_vertices_on_lod
                Point_2_NS_less(_m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = BEFORE_TO_AFTER;
            }
#endif
            n_crossings_before_to_after++;
        }
        
        // iterate
        curr = curr->next();
        
    } while (curr != last);
    
    return std::make_pair(n_crossings_after_to_before,
                          n_crossings_before_to_after);
}

CGAL_END_NAMESPACE

#endif
// EOF
