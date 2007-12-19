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
    m_left(CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES),
    m_right(CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES),
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
    m_left(CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES),
    m_right(CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES),
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
    typename Dcel::Vertex_iterator vit;
    CGAL::Arr_parameter_space ps_x, ps_y;

    // TODO what about m_left/m_right if called from accessor?

    v_left = v_right = NULL;
    for (vit = this->m_dcel.vertices_begin();
         vit != this->m_dcel.vertices_end(); ++vit) {
        // First check whether the vertex has a boundary condition in x.
        // If so, then a negative boundary condition indicates it is the left
        // vertex, and a positive boundary condition indicates it is the right
        // vertex.
        ps_x = vit->parameter_space_in_x();
        if (ps_x != CGAL::ARR_INTERIOR) {
            if (ps_x == CGAL::ARR_LEFT_BOUNDARY) {
                v_left = &(*vit);
            } else {
                v_right = &(*vit);
            }
        }
        
        // In case the vertex lies on the line of dicontinuity, it is
        // associated with a concrete point. Map this point to the vertex.
        ps_y = vit->parameter_space_in_y();
        if (ps_y != CGAL::ARR_INTERIOR) {
            
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
    if (this->m_left == CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES &&
        this->m_right == CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES) {
        //this->m_left = CGAL::ARR_CONTRACTION;
        //this->m_right = CGAL::ARR_CONTRACTION;
        this->m_left = CGAL::ARR_UNBOUNDED;
        this->m_right = CGAL::ARR_UNBOUNDED;
    } 
#endif
    CGAL_precondition(this->m_left != CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES &&
                      this->m_right != CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES);
    
    // create the face
    this->f_top = this->m_dcel.new_face();
    
    // set not fictious
    this->f_top->set_fictitious (false);

    // bounded or unbounded
    if (this->m_left == CGAL::ARR_UNBOUNDED || 
        this->m_right == CGAL::ARR_UNBOUNDED) {
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
 const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
 CGAL::Arr_parameter_space ps_x, CGAL::Arr_parameter_space ps_y) const
{
    CGAL_precondition (ps_x != CGAL::ARR_INTERIOR || 
                       ps_y != CGAL::ARR_INTERIOR
    );
    
    // In case the given boundary conditions do not match those of the given
    // vertex, v cannot represent the curve end.
    if (ps_x != v->parameter_space_in_x()) {
        return false;
    }
    
    if (ps_x != CGAL::ARR_INTERIOR) {
        return (ps_x == v->parameter_space_in_x());
    } else  {
        CGAL_assertion (ps_y != CGAL::ARR_INTERIOR);
        // check wether the two concrete points are equal
        return (this->m_traits->compare_xy_2_object() (
                // since ps_x == CGAL::ARR_INTERIOR we immediately have
                // a valid vertex v (with point)
                        v->point(),
                        (ind == CGAL::ARR_MIN_END ?
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
     const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
     CGAL::Arr_parameter_space ps_x, CGAL::Arr_parameter_space ps_y)
{
    //std::cout << "Arr_qdx_topology_traits_2 place_boundary_vertex"  
    //          << std::endl;

    CGAL_precondition(
            (ps_x != CGAL::ARR_INTERIOR ||
             ps_y != CGAL::ARR_INTERIOR)
    );
    
    // this topology return either an empty object or a DCEL vertex,
    // but never a fictious edge!!!
    
    Vertex *v = NULL;
    
    if (ps_x != CGAL::ARR_INTERIOR) {
        // for points at infinity/singularity
        // curve-end goes to v_left or to v_right
        v = (ps_x == CGAL::ARR_LEFT_BOUNDARY ? this->v_left : this->v_right);
    } else {
        CGAL_assertion(ps_y != CGAL::ARR_INTERIOR);
        // locate curve-end (here a concrete point) 
        // in local structure of for points on 
        // line of discontiuity!
        typename Line_of_discontinuity::iterator it = 
            this->m_line_of_discontinuity.find(
                    (ind == CGAL::ARR_MIN_END ?
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
    CGAL_assertion(v->parameter_space_in_x() == ps_x && 
                   v->parameter_space_in_y() == ps_y);

    CGAL_assertion(!v->has_null_point() || 
                   (ps_x == CGAL::ARR_LEFT_BOUNDARY && 
                    m_left == CGAL::ARR_CONTRACTION) || 
                   (ps_x == CGAL::ARR_RIGHT_BOUNDARY &&
                    m_right == CGAL::ARR_CONTRACTION)
    );
    
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
     const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
     CGAL::Arr_parameter_space ps_x, CGAL::Arr_parameter_space ps_y) const
{
    CGAL_precondition(
            (ps_x == CGAL::ARR_INTERIOR || ps_y == CGAL::ARR_INTERIOR) &&
            !(ps_x != CGAL::ARR_INTERIOR && ps_y != CGAL::ARR_INTERIOR)
    );

    // std::cout << "locate_around_boundary_vertex()" << std::endl;
    if (ps_x == CGAL::ARR_LEFT_BOUNDARY && ind == CGAL::ARR_MIN_END) {
        CGAL_assertion(ps_y == CGAL::ARR_INTERIOR);
        CGAL_assertion(v == v_left);
        bool dummy;
        return (_locate_around_vertex_with_boundary_at_x(
                        v_left, cv, ind, dummy
                )
        );
    }
    if (ps_x == CGAL::ARR_RIGHT_BOUNDARY && ind == CGAL::ARR_MAX_END) {
        CGAL_assertion(ps_y == CGAL::ARR_INTERIOR);
        CGAL_assertion(v == v_right);
        bool dummy;
        return (_locate_around_vertex_with_boundary_at_x(
                        v_right, cv, ind, dummy
                )
        );
    }
    
    CGAL_assertion(ps_x == CGAL::ARR_INTERIOR && 
                   ps_y != CGAL::ARR_INTERIOR);
    
    return (_locate_around_vertex_on_discontinuity (v, cv, ind));
}


template <class GeomTraits, class Dcel_>
typename Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::Halfedge* 
Arr_qdx_topology_traits_2<GeomTraits,Dcel_>::
_locate_around_vertex_with_boundary_at_x
(Vertex *v,
 const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind, bool& equal,
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

    const CGAL::Arr_curve_end curve_end = 
        (v == v_left) ? CGAL::ARR_MIN_END : CGAL::ARR_MAX_END;

    // If we compare a curve and its successor around the left/right
    // event, the result LARGER (left) [SMALLER (right)] 
    // indicates that the line of
    // discontinuity is located in between the two curves.
    const Comparison_result cross_res = 
        (curve_end == CGAL::ARR_MAX_END ? CGAL::LARGER : CGAL::SMALLER);
    
    // Traverse all other halfedges, and compare their y-positions next to the
    // pole with the query curve xc.
    typename Traits_adaptor_2::Compare_y_near_boundary_2 
        compare_y_near_boundary = 
        m_traits->compare_y_near_boundary_2_object();
    
    Comparison_result curr_res, next_res;
    Comparison_result curr_next_res;
    
    curr_res = compare_y_near_boundary (cv, curr->curve(), curve_end);
    do {
        if (allow_equal && curr_res == CGAL::EQUAL) {
            return (curr);
        }
        next_res = compare_y_near_boundary(cv, next->curve(), curve_end);
        if (allow_equal && next_res == CGAL::EQUAL) {
            return (next);
        }
        curr_next_res = 
            compare_y_near_boundary(curr->curve(), next->curve(), curve_end);
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
    CGAL_error();
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
                                       CGAL::Arr_curve_end ind) const
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
    std::cout << "dir: " << (curr->direction() == ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << (next->direction() == ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    
    std::cout << "******************************************" << std::endl;
#endif
    
    while (!is_between_cw(xc, (ind == ARR_MIN_END),
                          curr->curve(), 
                          (curr->direction() == ARR_RIGHT_TO_LEFT),
                          next->curve(), 
                          (next->direction() == ARR_RIGHT_TO_LEFT),
                          v->point(), eq_curr, eq_next))
    {
        // The curve must not be equal to one of the curves 
        // already incident to v.
        CGAL_assertion(!eq_curr && !eq_next);
        

#if 0
    std::cout << "??????????????????????????????????????????" << std::endl;
    std::cout << "search: " << std::endl;

    std::cout << "curr: " << curr->curve() << std::endl;
    std::cout << "dir: " << (curr->direction() == ARR_LEFT_TO_RIGHT ?
                             "L2R" : "R2L") << std::endl;
    std::cout << "next: " << next->curve() << std::endl;
    std::cout << "dir: " << (next->direction() == ARR_LEFT_TO_RIGHT ?
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
 CGAL::Arr_curve_end ind,
 CGAL::Arr_parameter_space ps_x,
 CGAL::Arr_parameter_space ps_y) const
{
    // In the planar-topology traits this function should never be invoked:
    //std::cout << "Arr_qdx_topology_traits_2::" 
    //          << "notify_on_boundary_vertex_creation"
    //          << std::endl;       

    CGAL_assertion(v->parameter_space_in_x() == ps_x);
    CGAL_assertion(v->parameter_space_in_y() == ps_y);
    
    // update locate structures
    if (ps_x == CGAL::ARR_LEFT_BOUNDARY) {
        //std::cout << "LEFT vertex created" << std::endl;
        this->v_left = v; 
        return;
    }
    if (ps_x == CGAL::ARR_RIGHT_BOUNDARY) {
        //std::cout << "RIGHT vertex created" << std::endl;
        this->v_right = v;
        return;
    }

    // else: boundary in y
    CGAL_assertion(ps_y != ARR_INTERIOR);
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
locate_curve_end (const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
                  CGAL::Arr_parameter_space ps_x, 
                  CGAL::Arr_parameter_space ps_y)
{
    // NEEDED for incremental insertion
    std::cout << "Arr_qdx_topology_traits_2 locate_curve_end"  
              << std::endl;
    CGAL_precondition(ps_x != CGAL::ARR_INTERIOR || 
                      ps_y != CGAL::ARR_INTERIOR);
    
    Vertex* v = NULL;
    typename Line_of_discontinuity::iterator  it;
    bool locate = false;

    CGAL_assertion_code(
            bool search_face = false;
    );

    if (ps_x != CGAL::ARR_INTERIOR) {
        bool contraction = false;
        switch (ps_x) {
        case CGAL::ARR_LEFT_BOUNDARY: {
            contraction = (m_left == CGAL::ARR_CONTRACTION);
            v = v_left;
            if (v != NULL) {
                if (contraction) {
                    return CGAL::make_object(v);
                } else {
                    locate = true;
                }
            } else {
                // search for face
                CGAL_assertion_code(search_face = true);
                it = m_line_of_discontinuity.begin();
            }
            break;
        }
        case CGAL::ARR_RIGHT_BOUNDARY: {
            contraction = (m_right == CGAL::ARR_CONTRACTION);
            v = v_right;
            if (v != NULL) {
                if (contraction) {
                    return CGAL::make_object(v);
                } else {
                    locate = true;
                }
            } else {
                return CGAL::make_object(f_top);
            }
            break;
        }
        default:
            CGAL_error(); // cannot happen
            break;
        }
    }

    if (ps_y != CGAL::ARR_INTERIOR) {
        Point_2 key = (ind == CGAL::ARR_MIN_END ?
                       this->m_traits->construct_min_vertex_2_object()(cv) :
                       this->m_traits->construct_max_vertex_2_object()(cv)); 
        
        it = m_line_of_discontinuity.find (key);
        if (it != m_line_of_discontinuity.end()) {
            v = it->second;
            return CGAL::make_object(v);
        }
        // else
        CGAL_assertion_code(search_face = true;);
        it = m_line_of_discontinuity.lower_bound (key);
    }
    
    if (locate) {
        // search for face that contains the curve-end
        bool overlaps;
        Halfedge *pred = 
            _locate_around_vertex_with_boundary_at_x(
                    v, cv, ind, overlaps, true
            );
        
        if (overlaps) {
            // or the half-edge overlapping with this curve-end
            return CGAL::make_object(pred);
        } 
        
        // else
        // return correct incident face
        if (pred->is_on_inner_ccb()) {
            return CGAL::make_object(pred->inner_ccb()->face());
        } else {
            return CGAL::make_object(pred->outer_ccb()->face());
        }
    }

    // else
    CGAL_assertion(search_face);
    
    // At this point, the iterator it points to a vertex on the line of
    // discontinuity that is strictly after the curve end. If there is none,
    // we know the curve end is contained in the top face. Otherwise,
    // we return the face that lies below the vertex v.
    if (it == m_line_of_discontinuity.end()) {
        return CGAL::make_object(f_top);
    }
    v = it->second;
    return CGAL::make_object(_face_before_vertex_on_discontinuity(v));
    
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
     const Halfedge *prev2,
     const X_monotone_curve_2& cv) const
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
    if ((m_left == CGAL::ARR_UNBOUNDED && m_right == CGAL::ARR_UNBOUNDED) ||
        (m_left == CGAL::ARR_CONTRACTION && m_right == CGAL::ARR_CONTRACTION)
    ) {
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
        return (m_left == CGAL::ARR_CONTRACTION);
    } else {
        CGAL_assertion(leftmost == AFTER_TO_BEFORE);
        // the face is "on the rightt" side of the surfacs
        return (m_right == CGAL::ARR_CONTRACTION);
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

        return (closest->direction() == CGAL::ARR_RIGHT_TO_LEFT);
        
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
    
    CGAL::Arr_parameter_space ps_x = 
        curr->opposite()->vertex()->parameter_space_in_x();
    if (ps_x == CGAL::ARR_LEFT_BOUNDARY) {
        res_source = CGAL::LARGER;
    } else if (ps_x == CGAL::ARR_RIGHT_BOUNDARY) {
        res_source = CGAL::SMALLER;
    } else {
        res_source = this->m_traits->compare_xy_2_object() (
                p, curr->opposite()->vertex()->point()
        );
    }
    CGAL_assertion(res_source != CGAL::EQUAL);
    
    CGAL::Arr_parameter_space last_by = CGAL::ARR_INTERIOR;
    CGAL::Arr_parameter_space curr_by = CGAL::ARR_INTERIOR;
    

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
        
        ps_x = curr->vertex()->parameter_space_in_x();
        if (ps_x == CGAL::ARR_LEFT_BOUNDARY) {
            res_target = CGAL::LARGER;
        } else if (ps_x == CGAL::ARR_RIGHT_BOUNDARY) {
            res_target = CGAL::SMALLER;
        } else {
            res_target = this->m_traits->compare_xy_2_object()
                (p, curr->vertex()->point());
        }
        CGAL_assertion(res_target != CGAL::EQUAL);  
        
        // read the boundary_type at src of curr
        curr_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MIN_END : CGAL::ARR_MAX_END)
        );
        
#if 1
        if ((last_by == CGAL::ARR_BOTTOM_BOUNDARY &&
             curr_by == CGAL::ARR_TOP_BOUNDARY) ||
            (last_by == CGAL::ARR_TOP_BOUNDARY &&
             curr_by == CGAL::ARR_BOTTOM_BOUNDARY)) {
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
            std::cout << "dir: " << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ?
                                     "ARR_LEFT_TO_RIGHT" : "ARR_RIGHT_TO_LEFT") 
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
        last_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MAX_END : CGAL::ARR_MIN_END)
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
    if (this->m_left != CGAL::ARR_UNBOUNDED &&
        this->m_right != CGAL::ARR_UNBOUNDED) {
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
            return (m_left == CGAL::ARR_UNBOUNDED);
        } else {
            CGAL_assertion(leftmost == AFTER_TO_BEFORE);
            // the face is "on the rightt" side of the surfacs
            return (m_right == CGAL::ARR_UNBOUNDED);
        }
        
        /* NOT REACHED */
        CGAL_error();
        return false;
    }
    case 2:
        // these two outer_ccbs are perimetric and therefore the face
        // must be bounded
        return false;
    default:
        //std::cout << "More than two outer_ccbs! Not nice!" << std::endl;
        CGAL_error();
    }
    
    /* should not be reached */
    CGAL_error();
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
    CGAL_precondition(v->parameter_space_in_x() != CGAL::ARR_INTERIOR || 
                      v->parameter_space_in_y() != CGAL::ARR_INTERIOR);
    
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
    if (v_left == NULL && v_right == 0 && (m_left == CGAL::ARR_UNBOUNDED ||
                                           m_right == CGAL::ARR_UNBOUNDED)) {
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

/*! \brief Return the face that lies before the given vertex, which lies
 * on the line of discontinuity.
 */
template <class GeomTraits, class Dcel>
typename Arr_qdx_topology_traits_2<GeomTraits, Dcel>::Face *
Arr_qdx_topology_traits_2<GeomTraits, Dcel>::
_face_before_vertex_on_discontinuity (Vertex * v) const {
    
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
    
    // Otherwise, we traverse the halfedges around v and locate the first
    // halfedge we encounter if we go from "3 o'clock" clockwise.
    // First locate the lower left and the top right halfedges around v.
    typename Traits_adaptor_2::Compare_x_on_identification_2 
        compare_x_on_identification = 
        m_traits->compare_x_on_identification_2_object();
    
    CGAL::Arr_curve_end leftmost_top_end = CGAL::ARR_MIN_END;
    Halfedge  *leftmost_top = NULL;
    CGAL::Arr_curve_end rightmost_bottom_end = CGAL::ARR_MIN_END;
    Halfedge  *rightmost_bottom = NULL;
    
    do {
        typename Traits_adaptor_2::Parameter_space_in_x_2 
            parameter_space_in_x =
            m_traits->parameter_space_in_x_2_object();
        typename Traits_adaptor_2::Parameter_space_in_y_2 
            parameter_space_in_y =
            m_traits->parameter_space_in_y_2_object();

        CGAL::Arr_curve_end ind = CGAL::ARR_MIN_END;
        
        CGAL::Arr_parameter_space ps_x = 
            parameter_space_in_x(curr->curve(), CGAL::ARR_MAX_END);
        CGAL::Arr_parameter_space ps_y = 
            parameter_space_in_y(curr->curve(), CGAL::ARR_MAX_END);
        if (are_equal(v, curr->curve(), CGAL::ARR_MAX_END, ps_x, ps_y)) {
            ind = CGAL::ARR_MAX_END;
        }
        
        if (parameter_space_in_y(curr->curve(), ind) == 
            // TODO check whether not bottom?
            CGAL::ARR_TOP_BOUNDARY) {
            if ((leftmost_top == NULL) || 
                (leftmost_top->direction() == CGAL::ARR_LEFT_TO_RIGHT &&
                 leftmost_top->direction() != curr->direction()) ||
                (leftmost_top->direction() == curr->direction() &&
                 compare_x_on_identification(
                         (ind == CGAL::ARR_MIN_END ?
                          this->m_traits->construct_min_vertex_2_object()(
                                  curr->curve()
                          ) :
                          this->m_traits->construct_max_vertex_2_object()(
                                  curr->curve()
                          )),
                         (leftmost_top_end == CGAL::ARR_MIN_END ?
                          this->m_traits->construct_min_vertex_2_object()(
                                  leftmost_top->curve()
                          ) :
                          this->m_traits->construct_max_vertex_2_object()(
                                  leftmost_top->curve()
                          ))
                 ) == CGAL::SMALLER
                )
            ) {
                leftmost_top_end = ind;
                leftmost_top = curr;
            } 
        } else {
            if ((rightmost_bottom == NULL) || 
                (rightmost_bottom->direction() == CGAL::ARR_RIGHT_TO_LEFT &&
                 rightmost_bottom->direction() != curr->direction()) ||
                (rightmost_bottom->direction() == curr->direction() &&
                 compare_x_on_identification(
                         (ind == CGAL::ARR_MIN_END ?
                          this->m_traits->construct_min_vertex_2_object()(
                                  curr->curve()
                          ) :
                          this->m_traits->construct_max_vertex_2_object()(
                                  curr->curve()
                          )),
                         (rightmost_bottom_end == CGAL::ARR_MIN_END ?
                          this->m_traits->construct_min_vertex_2_object()(
                                  rightmost_bottom->curve()
                          ) :
                          this->m_traits->construct_max_vertex_2_object()(
                                  rightmost_bottom->curve()
                          ))
                 ) == CGAL::LARGER
                )
            ) {
                rightmost_bottom_end = ind;
                rightmost_bottom = curr;
            } 
        }
        
        // Move to the next halfedge around the vertex.
        curr = curr->next()->opposite();
        
    } while (curr != first);
    
    // The first halfedge we encounter is the leftmost to top, but if there
    // is no edge to the left, we first encounter the righmost halfedge to the 
    // bottom. Note that as the halfedge we located has v as its target, we now
    // have to return its twin.
    if (leftmost_top != NULL) {
        first = leftmost_top->opposite();
    } else {
        first = rightmost_bottom->opposite();
    }
    
    // Return the incident face.
    if (first->is_on_inner_ccb()) {
        return (first->inner_ccb()->face());
    } else {
        return (first->outer_ccb()->face());
    }
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
    
    CGAL::Arr_parameter_space thistgt_by = CGAL::ARR_INTERIOR;
    CGAL::Arr_parameter_space nextsrc_by = CGAL::ARR_INTERIOR;
    
    const Halfedge *curr = he1;
    const Vertex *leftmost_vertex = NULL;
    
    // we count the number of crossings with the line of disc
    unsigned int n_crossings_before_to_after = 0;
    unsigned int n_crossings_after_to_before = 0;

    if (he2 == NULL) {
        // also check prev()->tgt with curr->src()
        // read the boundary_type at tgt of curr
        thistgt_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->prev()->curve(), 
                (curr->prev()->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MAX_END : CGAL::ARR_MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MIN_END : CGAL::ARR_MAX_END)
        );
        
        if (thistgt_by == CGAL::ARR_BOTTOM_BOUNDARY &&
            nextsrc_by == CGAL::ARR_TOP_BOUNDARY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::ARR_TOP_BOUNDARY &&
            nextsrc_by == CGAL::ARR_BOTTOM_BOUNDARY) {
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
        thistgt_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->curve(), 
                (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MAX_END : CGAL::ARR_MIN_END)
        );
        
        // read the boundary_type at src of next
        nextsrc_by = this->m_traits->parameter_space_in_y_2_object()(
                curr->next()->curve(), 
                (curr->next()->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 CGAL::ARR_MIN_END : CGAL::ARR_MAX_END)
        );
        
        if (thistgt_by == CGAL::ARR_BOTTOM_BOUNDARY &&
            nextsrc_by == CGAL::ARR_TOP_BOUNDARY) {
            if (leftmost_vertex == NULL || 
                // TASK avoid real comparisons, ask m_vertices_on_lod
                Point_2_less(m_traits)(curr->vertex()->point(),
                                       leftmost_vertex->point())) {
                leftmost_vertex = curr->vertex();
                leftmost = AFTER_TO_BEFORE;
            }
            n_crossings_after_to_before++;
        }
        if (thistgt_by == CGAL::ARR_TOP_BOUNDARY &&
            nextsrc_by == CGAL::ARR_BOTTOM_BOUNDARY) {
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
