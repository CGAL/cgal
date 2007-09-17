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

#ifndef CGAL_ARR_QDX_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_QDX_TOPOLOGY_TRAITS_2_IMPL_H

/*! \file
 * Member-function definitions for the
 * Arr_qdx_topology_traits_2<GeomTraits> class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Default constructor.
//
template <class GeomTraits, class Dcel_>
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::
Arr_qdx_topology_traits_2() :
    m_own_traits (true), 
    m_left(CGAL::NO_BOUNDARY),
    m_right(CGAL::NO_BOUNDARY),
    v_left(NULL),
    v_right(NULL),
    f_top(NULL)
{
    //std::cout << "Arr_qdx_topology_traits_2 ctor"  << std::endl;
    m_traits = new Traits_adaptor_2;
    m_line_of_discontinuity = Line_of_discontinuity(Point_2_less(m_traits));
}

//-----------------------------------------------------------------------------
// Constructor with a geometry-traits class.
//
template <class GeomTraits, class Dcel_>
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::
Arr_qdx_topology_traits_2 (Geometry_traits_2 *tr) : 
    m_own_traits(false),  
    m_left(CGAL::NO_BOUNDARY),
    m_right(CGAL::NO_BOUNDARY),
    v_left(NULL),
    v_right(NULL),
    f_top(NULL)
{
    //std::cout << "Arr_qdx_topology_traits_2 ctor(*tr)"  << std::endl;
    m_traits = static_cast<Traits_adaptor_2*>(tr);
    m_line_of_discontinuity = Line_of_discontinuity(Point_2_less(m_traits));

    _initialize_with_quadric(tr->base_quadric());
}

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <class GeomTraits, class Dcel_>
void Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::assign // open
    (const Self& other)
{
    // Assign the class.
    // Clear the current DCEL and duplicate the other DCEL.
    m_dcel.delete_all();
    m_dcel.assign (other.m_dcel);
    
    // Take care of the traits object.
    if (m_own_traits && m_traits != NULL) {
        delete m_traits;
    }
    
    if (other.m_own_traits) {
        m_traits = new Traits_adaptor_2;
    } else {
        m_traits = other.m_traits;
    }
    m_own_traits = other.m_own_traits;
    
    m_left = other.m_left;
    m_right = other.m_right;

    // Update the special properties of the topology traits.
    dcel_updated();

    return;
}

//-----------------------------------------------------------------------------
// Make the necessary updates after the DCEL structure have been updated.
//
template <class GeomTraits, class Dcel_>
void Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::dcel_updated ()
{
    // Go over the DCEL vertices and locate all points with boundary condition
    typename Dcel::Vertex_iterator       vit;
    Boundary_type                        bx, by;

    // TODO what about m_left/m_right if called from accessor?

    v_left = v_right = NULL;
    for (vit = this->m_dcel.vertices_begin();
         vit != this->m_dcel.vertices_end(); ++vit) {
        // First check whether the vertex has a boundary condition in x.
        // If so, then a negative boundary condition indicates it is the left
        // vertex, and a positive boundary condition indicates it is the right
        // vertex.
        bx = vit->boundary_in_x();
        if (bx != CGAL::NO_BOUNDARY) {
            if (CGAL::sign (bx) == CGAL::NEGATIVE) {
                CGAL_assertion (bx == MINUS_INFINITY ||
                                bx == AFTER_SINGULARITY);
                v_left = &(*vit);
            } else {
                CGAL_assertion (bx == PLUS_INFINITY ||
                                bx == BEFORE_SINGULARITY);
                v_right = &(*vit);
            }
        }
        
        // In case the vertex lies on the line of dicontinuity, it is
        // associated with a concrete point. Map this point to the vertex.
        by = vit->boundary_in_y();
        if (by != CGAL::NO_BOUNDARY) {
            
            std::pair< typename Line_of_discontinuity::iterator, bool > res =
                m_line_of_discontinuity.insert (std::make_pair(vit->point(),
                                                               &(*vit)));
            CGAL_assertion(! res.second);
            
            m_vertices_on_line_of_discontinuity[&(*vit)] = res.first;
        }
    }
    //CGAL_assertion (v_left != NULL);
    //CGAL_assertion (v_right != NULL);

    // Go over the DCEL faces and locate the top face, which is the only
    // face with no outer CCB.
    typename Dcel::Face_iterator         fit;
    
    f_top = NULL;
    for (fit = this->m_dcel.faces_begin();
         fit != this->m_dcel.faces_end(); ++fit) {
        if (fit->number_of_outer_ccbs() == 0) {
            CGAL_assertion (f_top == NULL);
            
            f_top = &(*fit);
            break;
        }
    }
    CGAL_assertion (f_top != NULL);
    
    return;
}

//-----------------------------------------------------------------------------
// Initialize an empty DCEL structure.
//
template <class GeomTraits, class Dcel_>
void Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::init_dcel () // open
{
    //std::cout << "Arr_qdx_topology_traits_2 init_dcel"  << std::endl;
    // Clear the current DCEL.
    this->m_dcel.delete_all();

#if 1
    // TODO workaround for missing toptraits constructor (default ctor of Arr)
    // assume elliposid
    if (this->m_left == CGAL::NO_BOUNDARY &&
        this->m_right == CGAL::NO_BOUNDARY) {
        //this->m_left = CGAL::AFTER_SINGULARITY;
        //this->m_right = CGAL::BEFORE_SINGULARITY;
        this->m_left = CGAL::MINUS_INFINITY;
        this->m_right = CGAL::PLUS_INFINITY;
    } 
#endif
    CGAL_precondition(this->m_left != CGAL::NO_BOUNDARY &&
                      this->m_right != CGAL::NO_BOUNDARY);
    CGAL_precondition(this->m_left < this->m_right);

    // create the face
    this->f_top = this->m_dcel.new_face();
    
    // set not fictious
    this->f_top->set_fictitious (false);

    // bounded or unbounded
    if (this->m_left == CGAL::MINUS_INFINITY || 
        this->m_right == CGAL::PLUS_INFINITY) {
        this->f_top->set_unbounded (true);
    } else {
        this->f_top->set_unbounded (false);   
    }

    v_left = NULL;
    v_right = NULL;
    m_line_of_discontinuity.clear();
    m_vertices_on_line_of_discontinuity.clear();
}

// TASK: check wther two following two methods are still needed
//-----------------------------------------------------------------------------
// Compare the relative y-position of the given point and the given edge
//
template <class GeomTraits, class Dcel_>
Comparison_result
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::compare_y_at_x
(const Point_2& p, const Halfedge* he) const
{
    //std::cout << "Arr_qdx_topology_traits_2 compare_y_at_x"  << std::endl;
    // all edges are valid, therefore just compare p to its associated curve.
    return (this->m_traits->compare_y_at_x_2_object() (p, he->curve()));
}

//-----------------------------------------------------------------------------
// Check if the given vertex is associated with the given curve end.
//
template <class GeomTraits, class Dcel_>
bool Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::are_equal
(const Vertex *v,
 const X_monotone_curve_2& cv, Curve_end ind,
 CGAL::Boundary_type bound_x, CGAL::Boundary_type bound_y) const
{
    CGAL_precondition (bound_x == m_left || 
                       bound_x == m_right ||
                       bound_y == CGAL::AFTER_DISCONTINUITY ||
                       bound_y == CGAL::BEFORE_DISCONTINUITY
    );
    
    // In case the given boundary conditions do not match those of the given
    // vertex, v cannot represent the curve end.
    if (bound_x != v->boundary_in_x()) {
        return false;
    }
    
    if (bound_x != CGAL::NO_BOUNDARY) {
        return (bound_x == v->boundary_in_x());
    } else  {
        CGAL_assertion (bound_y != CGAL::NO_BOUNDARY);
        // check wether the two concrete points are equal
        return (this->m_traits->compare_xy_2_object() (
                // since bound_x == CGAL::NO_BOUNDARY we immediately have
                // a valid vertex v (with point)
                        v->point(),
                        (ind == CGAL::MIN_END ?
                         this->m_traits->construct_min_vertex_2_object()(cv) :
                         this->m_traits->construct_max_vertex_2_object()(cv))) 
                == CGAL::EQUAL
        );
    }
}

//-----------------------------------------------------------------------------
// Given a curve end with boundary conditions and a face that contains the
// interior of the curve, find a place for a boundary vertex that will
// represent the curve end along the face boundary.
//
template <class GeomTraits, class Dcel_>
CGAL::Object
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::place_boundary_vertex // done
    (Face *f,
     const X_monotone_curve_2& cv, CGAL::Curve_end ind,
     Boundary_type bound_x, Boundary_type bound_y)
{
    //std::cout << "Arr_qdx_topology_traits_2 place_boundary_vertex"  
    //          << std::endl;

    CGAL_precondition(
            (bound_x == CGAL::NO_BOUNDARY && 
             (bound_y == CGAL::NO_BOUNDARY ||
              bound_y == CGAL::AFTER_DISCONTINUITY ||
              bound_y == CGAL::BEFORE_DISCONTINUITY))
            ||
            ((bound_x == CGAL::MINUS_INFINITY ||
              bound_x == CGAL::PLUS_INFINITY ||
              bound_x == CGAL::AFTER_SINGULARITY ||
              bound_x == CGAL::BEFORE_SINGULARITY) &&
             bound_y == CGAL::NO_BOUNDARY)
    );
    
    // this topology return either an empty object or a DCEL vertex,
    // but never a fictious edge!!!
    
    Vertex *v = NULL;
    
    if (bound_x != CGAL::NO_BOUNDARY) {
        // for points at infinity/singularity
        // curve-end goes to v_left or to v_right
        v = (bound_x < CGAL::NO_BOUNDARY ? this->v_left : this->v_right);
    } else {
        CGAL_assertion(bound_y != CGAL::NO_BOUNDARY);
        // locate curve-end (here a concrete point) 
        // in local structure of for points on 
        // line of discontiuity!
        typename Line_of_discontinuity::iterator it = 
            this->m_line_of_discontinuity.find(
                    (ind == CGAL::MIN_END ?
                     this->m_traits->construct_min_vertex_2_object()(cv) :
                     this->m_traits->construct_max_vertex_2_object()(cv)
                    )
            );

        if (it != this->m_line_of_discontinuity.end()) {
            // found vertex for curve-end 
            v = it->second;
            CGAL_assertion(!v->has_null_point());
        }
    }
    // if there is no vertex found, return empty object
    if (v == NULL) {
        //std::cout << "no vertex found" << std::endl;
        return CGAL::Object();
    }

    // else we return the vertex we have located.
    CGAL_assertion(v->boundary_in_x() == bound_x && 
                   v->boundary_in_y() == bound_y);

    CGAL_assertion(!v->has_null_point() || 
                   v->boundary_in_x() == CGAL::AFTER_SINGULARITY ||
                   v->boundary_in_x() == CGAL::BEFORE_SINGULARITY);
    
    return (CGAL::make_object (v));
}

//-----------------------------------------------------------------------------
// Locate the predecessor halfedge for the given curve around a given
// vertex with boundary conditions.
//
template <class GeomTraits, class Dcel_>
typename Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::Halfedge* 
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::locate_around_boundary_vertex
    (Vertex *v,
     const X_monotone_curve_2& cv, Curve_end ind,
     Boundary_type bound_x, Boundary_type bound_y) const
{
    CGAL_precondition(
            (bound_x == CGAL::NO_BOUNDARY || bound_y == CGAL::NO_BOUNDARY) &&
            !(bound_x != CGAL::NO_BOUNDARY && bound_y != CGAL::NO_BOUNDARY)
    );

    // std::cout << "locate_around_boundary_vertex()" << std::endl;
    if (bound_x == m_left && ind == CGAL::MIN_END) {
        CGAL_assertion(bound_y == CGAL::NO_BOUNDARY);
        CGAL_assertion(v == v_left);
        bool dummy;
        return (_locate_around_vertex_with_boundary_at_x(
                        v_left, cv, ind, dummy
                )
        );
    }
    if (bound_x == m_right && ind == CGAL::MAX_END) {
        CGAL_assertion(bound_y == CGAL::NO_BOUNDARY);
        CGAL_assertion(v == v_right);
        bool dummy;
        return (_locate_around_vertex_with_boundary_at_x(
                        v_right, cv, ind, dummy
                )
        );
    }
    
    CGAL_assertion(bound_x == CGAL::NO_BOUNDARY && 
                   (bound_y == CGAL::AFTER_DISCONTINUITY) ||
                   (bound_y == CGAL::BEFORE_DISCONTINUITY));
    
    return (_locate_around_vertex_on_discontinuity (v, cv, ind));
}


template <class GeomTraits, class Dcel_>
typename Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::Halfedge* 
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::
_locate_around_vertex_with_boundary_at_x
(Vertex *v,
 const X_monotone_curve_2& cv, Curve_end ind, bool& equal,
 bool allow_equal = false) const {

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

    const Curve_end curve_end = 
        (v == v_left) ? CGAL::MIN_END : CGAL::MAX_END;

    // If we compare a curve and its successor around the left/right
    // event, the result LARGER (left) [SMALLER (right)] 
    // indicates that the line of
    // discontinuity is located in between the two curves.
    const Comparison_result cross_res = 
        (curve_end == CGAL::MAX_END ? CGAL::LARGER : CGAL::SMALLER);
    
    // Traverse all other halfedges, and compare their y-positions next to the
    // pole with the query curve xc.
    typename Traits_adaptor_2::Compare_y_at_x_2 cmp_y_at_x = 
        m_traits->compare_y_at_x_2_object();
    
    Comparison_result curr_res, next_res;
    Comparison_result curr_next_res;
    
    curr_res = cmp_y_at_x (cv, curr->curve(), curve_end);
    do {
        if (allow_equal && curr_res == CGAL::EQUAL) {
            return (curr);
        }
        next_res = cmp_y_at_x (cv, next->curve(), curve_end);
        if (allow_equal && next_res == CGAL::EQUAL) {
            return (next);
        }
        curr_next_res = 
            cmp_y_at_x(curr->curve(), next->curve(), curve_end);
        if (curr_next_res == cross_res) {
            // The line of discontinuity must lie between curr and next, so the
            // comparison result of cv with the two curves should be equal:
            if (curr_res == next_res) {
                return curr;
            }
        } else {
            // The line of discontinuity does not lie between curr and next, 
            // so the comparison results must be different if cv 
            // lies in between.
            if (curr_res != next_res) {
                return curr;
            }
        }
        
        // Move to the next halfedge around the pole.
        curr = next;
        curr_res = next_res;
        next = curr->next()->opposite();
    } while (curr != first);

    // We sould never reach here:
    CGAL_assertion(0);
    return NULL;
}


/*! \brief returns the halfedge, the target vertex of which is given, that is
 * the predecessor of a halfedge, the curve of which is given, that is about
 * to be inserted into the dcel.
 */
template <class GeomTraits, class Dcel>
typename Arr_qdx_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_qdx_topology_traits_2<GeomTraits, Dcel>::
_locate_around_vertex_on_discontinuity(Vertex* v,
                                       const X_monotone_curve_2 & xc,
                                       Curve_end ind) const
{
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
        m_traits->is_between_cw_2_object();
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
    
    while (!is_between_cw(xc, (ind == MIN_END),
                          curr->curve(), 
                          (curr->direction() == RIGHT_TO_LEFT),
                          next->curve(), 
                          (next->direction() == RIGHT_TO_LEFT),
                          v->point(), eq_curr, eq_next))
    {
        // The curve must not be equal to one of the curves 
        // already incident to v.
        CGAL_assertion(!eq_curr && !eq_next);
        

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
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::notify_on_boundary_vertex_creation
(Vertex *v,
 const X_monotone_curve_2& cv,
 Curve_end ind,
 Boundary_type bound_x,
 Boundary_type bound_y) const
{
    // In the planar-topology traits this function should never be invoked:
    //std::cout << "Arr_qdx_topology_traits_2::" 
    //          << "notify_on_boundary_vertex_creation"
    //          << std::endl;       

    CGAL_assertion(v->boundary_in_x() == bound_x);
    CGAL_assertion(v->boundary_in_y() == bound_y);
    
    // update locate structures
    if (bound_x < CGAL::NO_BOUNDARY) {
        //std::cout << "LEFT vertex created" << std::endl;
        this->v_left = v; 
        CGAL_assertion(this->v_left->boundary_in_x() == m_left);
        return;
    }
    if (bound_x > CGAL::NO_BOUNDARY) {
        //std::cout << "RIGHT vertex created" << std::endl;
        this->v_right = v;
        CGAL_assertion(this->v_right->boundary_in_x() == m_right);
        return;
    }

    // else: boundary in y
    CGAL_assertion(bound_y != NO_BOUNDARY);
    CGAL_assertion(!v->has_null_point());

    CGAL_assertion_code(
            int lod_size = 
            static_cast< int >(this->m_line_of_discontinuity.size())
    );
    CGAL_assertion(
            static_cast< int >(m_vertices_on_line_of_discontinuity.size()) ==
            lod_size
    );
    
    // update the local structure for points on the line of discontinuity
    typename Line_of_discontinuity::iterator it = 
        this->m_line_of_discontinuity.find(v->point());
    // not existing so far
    CGAL_assertion(it == this->m_line_of_discontinuity.end());
    // therefore insert it
    this->m_line_of_discontinuity.insert(it, std::make_pair(v->point(), v));
    CGAL_assertion(
            static_cast< int >(m_line_of_discontinuity.size()) ==
            lod_size + 1
    );

    //std::cout << "points on LoD: " 
    //          << m_line_of_discontinuity.size() << std::endl;
    
    // store iterator for vertex 
    // -> needed to delete vertex if becoming redundant
    m_vertices_on_line_of_discontinuity[v] = it; 
    CGAL_assertion(
            static_cast< int >(m_vertices_on_line_of_discontinuity.size()) ==
            lod_size + 1
    );
}

//-----------------------------------------------------------------------------
// Locate a DCEL feature that contains the given curve end.
//
template <class GeomTraits, class Dcel_>
CGAL::Object Arr_qdx_topology_traits_2<GeomTraits, Dcel_>:: // open
locate_curve_end (const X_monotone_curve_2& cv, Curve_end ind,
                  Boundary_type bound_x, Boundary_type bound_y)
{
    // \todo RWRW: Add support for all boundary conditions, not just
    //             for unbounded curve-ends!

    // NEEDED for incremental insertion
    std::cout << "Arr_qdx_topology_traits_2 locate_unb_curve_end"  
              << std::endl;
    CGAL_precondition(bound_x == CGAL::MINUS_INFINITY ||
                      bound_x == CGAL::PLUS_INFINITY);

    Vertex *v = (bound_x == CGAL::MINUS_INFINITY ? v_left : v_right);
    if (v == NULL) {
        // return the "initial" unbounded face
        // TODO locate correct face when perimetric paths exists
        CGAL_assertion(this->f_top->is_unbounded());
        return CGAL::make_object(this->f_top);
    } 

    // otherwise 
    // search for face that contains the curve-end
    bool overlaps;
    Halfedge *pred = 
        _locate_around_vertex_with_boundary_at_x(v, cv, ind, overlaps, true);
    
    if (overlaps) {
        // or the half-edge overlapping with this curve-end
        return CGAL::make_object(pred);
    } 

    // else
    // retusrn correct incident face
    // TODO check wether outer_ccb or inner_ccb must be returned
    return CGAL::make_object(pred->outer_ccb()->face());
}

//-----------------------------------------------------------------------------
// Given two predecessor halfedges that belong to the same inner CCB of
// a face, determine what happens when we insert an edge connecting the
// target vertices of the two edges.
//
template <class GeomTraits, class Dcel_>
std::pair<bool, bool>
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::face_split_after_edge_insertion
    (const Halfedge *prev1,
     const Halfedge *prev2) const
{
    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());

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
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::hole_creation_after_edge_removal
    (const Halfedge *he) const
{
    CGAL_precondition (! he->is_on_inner_ccb());
    CGAL_precondition (! he->opposite()->is_on_inner_ccb());

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
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::is_on_new_perimetric_face_boundary
(const Halfedge *prev1,
 const Halfedge *prev2,
 const X_monotone_curve_2& cv) const
{
    //std::cout << "Arr_qdx_topology_traits_2::" 
    //          << "is_on_new_perimetric_face_boundary"
    //          << std::endl;
    
    // can only be possible for a paraboloid.
    if ((m_left == CGAL::MINUS_INFINITY && m_right == PLUS_INFINITY) ||
        (m_left == CGAL::AFTER_SINGULARITY && m_right == BEFORE_SINGULARITY)) {
        return false;
    }
    CGAL_assertion(m_quadric.is_elliptic_paraboloid());
    
    // so we use leftmost disconti-crossing of 
    // path ending in prev1 (and starting in prev2)
    CGAL_assertion(_is_perimetric_path(prev2, prev1));
    Discontinuity_crossing leftmost;
    std::pair< unsigned int, unsigned int > crossings =
        _crossings_with_line_of_discontinuity(prev2, prev1, leftmost);
    // is perimetric test
    CGAL_assertion((crossings.first + crossings.second) % 2 == 1);

    // to check whether new face contains singular point
    // and therefore prev1 belongs to the outer_ccb of this new face
    if (leftmost == BEFORE_TO_AFTER) { 
        // the face is "on the left" side of the surfacs
        return (m_left == CGAL::AFTER_SINGULARITY);
    } else {
        CGAL_assertion(leftmost == AFTER_TO_BEFORE);
        // the face is "on the rightt" side of the surfacs
        return (m_right == CGAL::BEFORE_SINGULARITY);
    }
}

//-----------------------------------------------------------------------------
// checks whether halfedges are boundaries of the same face
//
template <class GeomTraits, class Dcel_>
bool
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::boundaries_of_same_face
(const Halfedge *e1,
 const Halfedge *e2) const
{
    //std::cout << "Arr_qdx_topology_traits_2::boundaries_of_same_face" 
    //          << std::endl;
    // This predicate is only used for case 3.3.2 of the insertion process
    
    Discontinuity_crossing leftmost1;
    std::pair< unsigned int, unsigned int > crossings1 =
        _crossings_with_line_of_discontinuity(e1, NULL, leftmost1);

    Discontinuity_crossing leftmost2;
    std::pair< unsigned int, unsigned int > crossings2 =
        _crossings_with_line_of_discontinuity(e2, NULL, leftmost2);

    return (leftmost1 != leftmost2);
}

//-----------------------------------------------------------------------------
// Determine whether the given vertex lies in the interior of the given face.
//
template <class GeomTraits, class Dcel_>
bool Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::is_in_face // open
(const Face *f, const Point_2& p, const Vertex *v) const
{
    // NEEDED for incremental insertion
    // TODO open is_in_face
    std::cout << "Arr_qdx_topology_traits_2::is_in_face" 
              << std::endl;

    CGAL_precondition (v == NULL || ! v->has_null_point());
    CGAL_precondition (v == NULL || 
                       m_traits->equal_2_object()(p, v->point()));
    
    // In case the face is unbounded and has no outer ccbs, this is the single
    // (un)bounded face of an arrangement of bounded curves. 
    // This face obviously contains any point in its interior.
    int num_occb = f->number_of_outer_ccbs();
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
        
        if (_is_perimetric_path(*occb_it, *occb_it)) {
            std::cout << "IS_IN_FACE: perimetric boundary face" << std::endl;
            perimetic_face = true;
        } else {
            std::cout << "IS_IN_FACE: normal face" << std::endl;
        }
        // otherwise, it is a normal face (maybe going over the curve of disc)
        
    }

    if (num_occb == 2) {
        // both outer_ccbs should be perimetric
        
        CGAL_assertion(_is_perimetric_path(*occb_it, *occb_it));
        occb_it++;
        CGAL_assertion(_is_perimetric_path(*occb_it, *occb_it));
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
            


            if (this->m_traits->is_in_x_range_2_object() (
                        curr->curve(), p
                ) &&
                this->m_traits->compare_y_at_x_2_object() (
                        p, curr->curve()) == CGAL::SMALLER) {
                if (closest == NULL) {
                    *closest = *curr;
                } else {
                    // compare with closest
                    if (this->m_traits->compare_y_position_2_object() (
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
        res_source = this->m_traits->compare_xy_2_object() (
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
            res_target = this->m_traits->compare_xy_2_object()
                (p, curr->vertex()->point());
        }
        CGAL_assertion(res_target != CGAL::EQUAL);  
        
        // read the boundary_type at src of curr
        curr_by = this->m_traits->boundary_in_y_2_object()(
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
        last_by = this->m_traits->boundary_in_y_2_object()(
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
}

//-----------------------------------------------------------------------------
// Determine whether the given face is unbounded.
//
template <class GeomTraits, class Dcel_>
bool Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::is_unbounded
(const Face *f) const
{
    //std::cout << "Arr_qdx_topology_traits_2 is_unbounded"  << std::endl;
    
    // if ellipsoid then naturally all faces are bounded
    if (this->m_left != CGAL::MINUS_INFINITY &&
        this->m_right != CGAL::PLUS_INFINITY) {
        return (false);
    }
    
    // Go over the outer CBB of the given face and 
    // look for a non-concrete point
    
    // we are on a cylinder or paraboloid

    switch (f->number_of_outer_ccbs()) {
    case 0:
        // if there is no edge then it must be unbounded
        return true;
    case 1: {
        // check whether it contains a vertex at inf
        const Halfedge *first = *(f->outer_ccbs_begin());
        const Halfedge *curr = first;
        do {
            if (curr->vertex()->has_null_point()) {
                // Found a non-concrete vertex along the boundary: 
                // f is unbounded.
                return (true);
            }
            curr = curr->next();
        } while (curr != first);
        
        // if not
        if (this->m_quadric.is_elliptic_cylinder()) {
            // each face consiting of a single outer_ccb that does not
            // touch a point at inf contains a "inf" 
            // therefore:
            return true;
        }
        //CGAL_assertion(this->m_quadric.is_elliptic_paraboloid());
        
        // check wether the face contains the singularity -> return false;
        // otherwise it contains a point at infinity -> return true;

        // we use the perimetric path defined by the outer ccb
        CGAL_assertion(_is_perimetric_path(first, first->next()));
        Discontinuity_crossing leftmost;
        std::pair< unsigned int, unsigned int > crossings =
            _crossings_with_line_of_discontinuity(first, NULL, leftmost);
        // is perimetric test
        CGAL_assertion((crossings.first + crossings.second) % 2 == 1);
        
        // and check how the leftmost crossing of the path crosses the line of
        /// discontiuity
        if (leftmost == BEFORE_TO_AFTER) { 
            // the face is "on the left" side of the surfacs
            return (m_left == CGAL::MINUS_INFINITY);
        } else {
            CGAL_assertion(leftmost == AFTER_TO_BEFORE);
            // the face is "on the rightt" side of the surfacs
            return (m_right == CGAL::PLUS_INFINITY);
        }
        
        /* NOT REACHED */
        CGAL_assertion(false);
        return false;
    }
    case 2:
        // these two outer_ccbs are perimetric and therefore the face
        // must be bounded
        return false;
    default:
        //std::cout << "More than two outer_ccbs! Not nice!" << std::endl;
        CGAL_assertion(false);
    }
    
    /* should not be reached */
    CGAL_assertion(false);
    return (false);
}


//-----------------------------------------------------------------------------
// Determine whether a boundary vertex is redundant
//
template <class GeomTraits, class Dcel_>
bool Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::is_redundant
(const Vertex *v) const
{
    //std::cout << "Arr_qdx_topology_traits_2 is_redundant"  << std::endl;
    CGAL_precondition(v->boundary_in_x() != CGAL::NO_BOUNDARY || 
                      v->boundary_in_y() != CGAL::NO_BOUNDARY);
    
    // if there are not incident edges just remove it
    // TASK: check whether isolated or degree == 0 is needed!
    return (v->is_isolated());
}

//-----------------------------------------------------------------------------
// Determine whether a boundary vertex is redundant
//
template <class GeomTraits, class Dcel_>
typename Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::Halfedge* 
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::erase_redundant_vertex // open
(Vertex *v) 
{
    std::cout << "Arr_qdx_topology_traits_2 erase_redundant_vertex"  
              << std::endl;
    if (v == v_left) {
        v_left = NULL;
    }
    if (v == v_right) {
        v_right = NULL;
    }
    // if there is no unbounded face connected to a special vertex
    // TASK check whether this is needed
    if (v_left == NULL && v_right == 0 && (m_left == CGAL::MINUS_INFINITY ||
                                           m_right == PLUS_INFINITY)) {
        // ensure that f_top is set to the correct one;
        if (!f_top->is_unbounded()) {
            // go over the DCEL faces and locate the unbounded face.
            for (typename Dcel::Face_iterator fit = this->m_dcel.faces_begin();
                 fit != this->m_dcel.faces_end(); ++fit) {
                if (fit->is_unbounded()) {
                    this->f_top = &(*fit);
                    break;
                }
            }
        }
    }

    // else v is a vertex on the line of discontinuity
    // unfortnunately no incident curve-end can give us the key
    // -> but we stored something else useful: find iterator
    typename Line_of_discontinuity::iterator it =
        m_vertices_on_line_of_discontinuity[v];
    // and delete this item
    m_line_of_discontinuity.erase(it);

    // valid a valid halfedge-pointer is only requested for if vertex
    // has been connecting fictiuous halfedges
    return NULL;
}

//-----------------------------------------------------------------------------
// Number of crossing with the line of discontiniuty
//
template <class GeomTraits, class Dcel_>
    std::pair< unsigned int, unsigned int >
Arr_qdx_topology_traits_2<GeomTraits, Dcel_>::
_crossings_with_line_of_discontinuity(
        const Halfedge* he1, const Halfedge* he2, 
        Discontinuity_crossing& leftmost) const {
    
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
        thistgt_by = this->m_traits->boundary_in_y_2_object()(
                curr->prev()->curve(), 
                (curr->prev()->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MAX_END : CGAL::MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MIN_END : CGAL::MAX_END)
        );
        
        if (thistgt_by == CGAL::AFTER_DISCONTINUITY &&
            nextsrc_by == CGAL::BEFORE_DISCONTINUITY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::BEFORE_DISCONTINUITY &&
            nextsrc_by == CGAL::AFTER_DISCONTINUITY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = BEFORE_TO_AFTER;
            }
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
        thistgt_by = this->m_traits->boundary_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MAX_END : CGAL::MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->m_traits->boundary_in_y_2_object()(
                curr->next()->curve(), 
                (curr->next()->direction() == CGAL::LEFT_TO_RIGHT ? 
                 CGAL::MIN_END : CGAL::MAX_END)
        );
        
        if (thistgt_by == CGAL::AFTER_DISCONTINUITY &&
            nextsrc_by == CGAL::BEFORE_DISCONTINUITY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::BEFORE_DISCONTINUITY &&
            nextsrc_by == CGAL::AFTER_DISCONTINUITY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = BEFORE_TO_AFTER;
            }
            n_crossings_before_to_after++;
        }
        
        // iterate
        curr = curr->next();
        
    } while (curr != last);
    
    return std::make_pair(n_crossings_after_to_before,
                          n_crossings_before_to_after);
}

//-----------------------------------------------------------------------------
// Check whether the path between two halfedges is perimetric.
//
template <class GeomTraits, class Dcel_>
bool
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::_is_perimetric_path
(const Halfedge *e1,
 const Halfedge *e2) const
{
    Discontinuity_crossing tmp;
    std::pair< unsigned int, unsigned int > crossings;
    // TODO check order
    if (e1 == e2) {
        crossings = _crossings_with_line_of_discontinuity(e1, NULL, tmp);
    } else {
        crossings = _crossings_with_line_of_discontinuity(e1, e2, tmp);
    }
    // return whether there has been an odd number of intersections
    bool res = ((crossings.first + crossings.second) % 2 == 1);
    return res;
}

CGAL_END_NAMESPACE

#endif
// EOF
