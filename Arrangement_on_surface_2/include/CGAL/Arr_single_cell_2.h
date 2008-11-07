// Copyright (c) 2008  Tel-Aviv University (Israel).
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
// Author(s)     : Eric Berberich     <ericb@post.tau.ac.il>


// flags in this file:

#ifndef CGAL_SINGLE_CELL_RI_NAIVE
#define CGAL_SINGLE_CELL_RI_NAIVE 0
#endif


#ifndef CGAL_ARR_SINGLE_CELL_2_H
#define CGAL_ARR_SINGLE_CELL_2_H

#include <boost/optional.hpp>
#include <boost/none.hpp>

#include <CGAL/Timer.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

// TASK select best point location strategy
#include <CGAL/Arr_naive_point_location.h>

// DOES NOT WORK FOR UNBOUNDED
//#include <CGAL/Arr_simple_point_location.h>
//#include <CGAL/Arr_walk_along_line_point_location.h> 

// required for Constructor of RI_observer
#include <CGAL/Arr_simple_point_location.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief
 * class that observes a city-arrangement
 * 
 * Keeps eye on face- and edge-splits
 * and updates internal structures wrt to given query point
 */
template < class Arrangement_2_ >
struct RI_observer : CGAL::Arr_observer< Arrangement_2_ > {
    
    //! this class template parameter
    typedef Arrangement_2_ Arrangement_2;

    //! base class
    typedef CGAL::Arr_observer< Arrangement_2 > Base;
    
    //! geometric traits class
    typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
    
    //! traits adaptor
    typedef Arr_traits_basic_adaptor_2< Geometry_traits_2 > Traits_adaptor_2;
    
    //! type of point
    typedef typename Geometry_traits_2::Point_2 Point_2;
    
    //! type of x-monotone curve
    typedef typename Geometry_traits_2::X_monotone_curve_2 
    X_monotone_curve_2;
    
    //! type of curve
    typedef typename Geometry_traits_2::Curve_2 Curve_2;
    
    typedef typename Arrangement_2::Vertex_const_handle 
    Vertex_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle 
    Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle Face_const_handle;
    
    typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

    typedef typename Arrangement_2::Vertex_handle            Vertex_handle;
    typedef typename Arrangement_2::Halfedge_handle          Halfedge_handle;
    typedef typename Arrangement_2::Face_handle              Face_handle;
    
    //!\name Constructors
    //!@{
    
    // FUTURE TODO allow multiple points

    //! default constructor from arr and point
    RI_observer(Arrangement_2& arr, const Point_2& point) :
        Base(arr),
        _m_point(point),
        _m_check_split_edges(false) {
        
        // determine initial _m_edge_handle
        // TODO avoud ray_shoot up
        CGAL::Arr_simple_point_location< Arrangement_2 > 
            ray(*(Base::arrangement()));
        CGAL::Object obj = ray.ray_shoot_up(_m_point);
        Vertex_const_handle vh;
        
        if (CGAL::assign(_m_halfedge_handle, obj)) {
            
            if (_m_halfedge_handle->direction() != CGAL::ARR_RIGHT_TO_LEFT) {
                _m_halfedge_handle = _m_halfedge_handle->twin();
            }
            CGAL_assertion(
                    _m_halfedge_handle->direction() == CGAL::ARR_RIGHT_TO_LEFT
            );
            
            _m_cell_handle = CGAL::make_object(_m_halfedge_handle->face());
            
        } else if (CGAL::assign(vh, obj)) {
            
            // TODO correct assign
            std::cerr << "has been assigned point" << std::endl;
            
        }  else {
            
            // TODO does not work with unbounded faces!!!
            Face_const_handle fh;
            CGAL_assertion_code(bool check =)
                CGAL::assign(fh, obj);
            CGAL_assertion(check);
            
            Vertex_const_handle 
                vtr(Base::arrangement()->topology_traits()->
                    top_right_vertex());
            
            Halfedge_around_vertex_const_circulator heh  = 
                vtr->incident_halfedges();
            heh++;
            
            CGAL_assertion_code(
                    Vertex_const_handle 
                    vtl(Base::arrangement()->topology_traits()->
                        top_left_vertex())
            );
            CGAL_assertion(heh->source() == vtl);
            
            _m_halfedge_handle = heh->twin();
            CGAL_assertion(
                    _m_halfedge_handle->direction() == CGAL::ARR_RIGHT_TO_LEFT
            );
            _m_cell_handle = CGAL::make_object(_m_halfedge_handle->face());

        }
            
#if !NDEBUG
        if (_m_halfedge_handle->is_fictitious()) {
            std::cout << "Initial curve: fict" << std::endl;
        } else {
            std::cout << "Initial curve: " << _m_halfedge_handle->curve()
                      << std::endl;
        }
#endif
    }

    //!@}
    
    //!\name Notifications
    //!@{

    /*!
     * Notification before the creation of a new edge.
     * \param c The x-monotone curve to be associated with the edge.
     * \param v1 A handle to the first end-vertex of the edge.
     * \param v2 A handle to the second end-vertex of the edge.
     */
    virtual void before_create_edge (const X_monotone_curve_2& /* c */,
                                     Vertex_handle /* v1 */,
                                     Vertex_handle /* v2 */)
    {}
    
    /*!
     * Notification after the creation of a new edge.
     * \param e A handle to one of the twin halfedges that were created.
     */
    virtual void after_create_edge (Halfedge_handle e)
    {
        typename Geometry_traits_2::Is_vertical_2 is_vertical =
            Base::arrangement()->geometry_traits()->is_vertical_2_object();
        typename Traits_adaptor_2::Is_in_x_range_2 is_in_x_range =
            static_cast< Traits_adaptor_2* >
            (Base::arrangement()->geometry_traits())->is_in_x_range_2_object();
        typename Geometry_traits_2::Compare_y_at_x_2 compare_y_at_x =
            Base::arrangement()->geometry_traits()->compare_y_at_x_2_object();
        
        CGAL_assertion(!e->is_fictitious());
        
        X_monotone_curve_2 new_cv = e->curve();
        
        std::cout << "after create edge" << std::endl;
        
        if (!is_vertical(new_cv)) {
            // TODO symbolic pertubation for endpoints!
            if (is_in_x_range(new_cv, _m_point)) {
                CGAL::Comparison_result res = 
                    compare_y_at_x(_m_point, new_cv);
                
                if (res == CGAL::SMALLER) {
                    
                    bool change = false;
                    
                    if (_m_halfedge_handle->is_fictitious()) {
                        change = true;
                        //std::cout << "CHANGE 1 (old fict)" << std::endl;
                    } else {
                        X_monotone_curve_2 old_cv = 
                            _m_halfedge_handle->curve();
                        
                        //std::cout << "old_cv: " << old_cv << std::endl;
                        //std::cout << "new_cv: " << new_cv << std::endl;

                        typename Traits_adaptor_2::Compare_y_position_2 
                            compare_y_pos =
                            static_cast< Traits_adaptor_2* >
                            (Base::arrangement()->geometry_traits())->
                            compare_y_position_2_object();
                        
                        if (is_in_x_range(new_cv, old_cv)) {
                            if (compare_y_pos(new_cv, old_cv) ==
                                CGAL::SMALLER) {
                                change = true;
                                std::cout << "CHANGE 2" << std::endl;
                            }
                        } else {
                            typename Geometry_traits_2::Compare_x_2 compare_x =
                                Base::arrangement()->geometry_traits()->
                                compare_x_2_object();
                            typename Geometry_traits_2::Construct_min_vertex_2 
                                construct_min =
                                Base::arrangement()->geometry_traits()->
                                construct_min_vertex_2_object();
                            typename Geometry_traits_2::Construct_max_vertex_2 
                                construct_max =
                                Base::arrangement()->geometry_traits()->
                                construct_max_vertex_2_object();
                            
                            Point_2 min = construct_min(new_cv);
                            if (compare_x(min, _m_point) == CGAL::EQUAL) {
                                change = (compare_y_at_x(min, old_cv) 
                                          == CGAL::SMALLER);
                                if (change) {
                                    std::cout << "CHANGE 3" << std::endl;
                                }
                            } else {
                                Point_2 max = construct_max(new_cv);
                                CGAL_assertion(
                                        compare_x(max, _m_point) == CGAL::EQUAL
                                );
                                change = (compare_y_at_x(max, old_cv)
                                          == CGAL::SMALLER);
                                if (change) {
                                    std::cout << "CHANGE 4" << std::endl;
                                }
                            }
                        }
                    }

                    // finally
                    if (change) {
                        // the new edge replaces the old one
                        // we only have to ensure the correct order
                        _m_halfedge_handle = 
                            ((e->direction() == CGAL::ARR_RIGHT_TO_LEFT) ?
                             e : e->twin());
                        
                        CGAL_assertion(
                                _m_halfedge_handle->direction() == 
                                CGAL::ARR_RIGHT_TO_LEFT
                        );
                        
                        _m_cell_handle = 
                            CGAL::make_object(_m_halfedge_handle->face());

                        std::cout << "HE1 set to cv: " 
                                  << _m_halfedge_handle->curve() 
                                  << std::endl;
                    }
                } else if (res == CGAL::EQUAL) {
                    // TODO on curve!!!!
                }
            }
        } else {
            // TODO what if e->curve() is vertical
        }
    }

#if 0
    // TODO remove? as after_create_edge subsumes this case
    /*!
     * Notification before the splitting of a face into two.
     * \param f A handle to the existing face.
     * \param e The new edge whose insertion causes the face to split.
     */
    virtual void before_split_face (Face_handle /* f */,
                                    Halfedge_handle e) 
    {
        typename Geometry_traits_2::Is_vertical_2 is_vertical =
            Base::arrangement()->geometry_traits()->is_vertical_2_object();
        typename Traits_adaptor_2::Is_in_x_range_2 is_in_x_range =
            static_cast< Traits_adaptor_2* >
            (Base::arrangement()->geometry_traits())->is_in_x_range_2_object();
        
        typename Geometry_traits_2::Compare_y_at_x_2 compare_y_at_x =
            Base::arrangement()->geometry_traits()->compare_y_at_x_2_object();
        
        X_monotone_curve_2 cv = e->curve();
        
        std::cout << "before split face" << std::endl;
        std::cout << "cv: " << cv << std::endl;
        std::cout << "he: " << _m_halfedge_handle->curve() << std::endl;

        if (!is_vertical(cv)) {
            // TODO symbolic pertubation for endpoints!
            if (is_in_x_range(cv, _m_point)) {
                CGAL::Comparison_result res = compare_y_at_x(_m_point, cv);
                
                if (res == CGAL::SMALLER) {
                    typename Traits_adaptor_2::Compare_y_position_2 
                        compare_y_pos =
                        static_cast< Traits_adaptor_2* >
                        (Base::arrangement()->geometry_traits())->
                        compare_y_position_2_object();
                    
                    bool change = false;

                    if (is_in_x_range(cv, _m_halfedge_handle->curve())) {
                        if (compare_y_pos(cv, _m_halfedge_handle->curve()) ==
                            CGAL::SMALLER) {
                            change = true;
                        }
                    } else {
                        typename Geometry_traits_2::Compare_x_2 compare_x =
                            Base::arrangement()->geometry_traits()->
                            compare_x_2_object();
                        typename Geometry_traits_2::Construct_min_vertex_2 
                            construct_min =
                            Base::arrangement()->geometry_traits()->
                            construct_min_vertex_2_object();
                        typename Geometry_traits_2::Construct_max_vertex_2 
                            construct_max =
                            Base::arrangement()->geometry_traits()->
                            construct_max_vertex_2_object();
                        
                        Point_2 min = construct_min(cv);
                        if (compare_x(min, _m_point) == CGAL::EQUAL) {
                            change = (compare_y_at_x(
                                              min, 
                                              _m_halfedge_handle->curve()) 
                                      == CGAL::SMALLER);
                        } else {
                            Point_2 max = construct_max(cv);
                            CGAL_assertion(
                                    compare_x(max, _m_point) == CGAL::EQUAL
                            );
                            change = (compare_y_at_x(
                                              max, 
                                              _m_halfedge_handle->curve()) 
                                      == CGAL::SMALLER);
                        }
                        
                        if (change) {
                            // new curve is below old curve
                            // thus, face that contains point is 
                            // restricted by this
                            // edge - we only have to ensure the correct order
                            _m_halfedge_handle = 
                                ((e->direction() == CGAL::ARR_RIGHT_TO_LEFT) ?
                                 e : e->twin());
                            
                            CGAL_assertion(
                                    _m_halfedge_handle->direction() == 
                                    CGAL::ARR_RIGHT_TO_LEFT
                            );
                            
                            _m_cell_handle = 
                                CGAL::make_object(_m_halfedge_handle->face());
                            

                            std::cout << "HE2 set to cv: " 
                                      << _m_halfedge_handle->curve() 
                                      << std::endl;

                        }
                    }
                } else if (res == CGAL::EQUAL) {
                    // TODO on curve!!!!
                }
            } 
        } else {
            // TODO what if e->curve() is vertical
        }
        
        // else no action is required
    }
#endif
   
    /*!
     * Notification after a face was split.
     * \param f A handle to the face we have just split.
     * \param new_f A handle to the new face that has been created.
     * \param is_hole Whether the new face forms a hole inside f.
     */
    virtual void after_split_face (Face_handle /* f */,
                                   Face_handle /* new_f */,
                                   bool /* is_hole */)
    {
        std::cout << "after_split_face" << std::endl;
        CGAL_assertion(
                _m_halfedge_handle->direction() == CGAL::ARR_RIGHT_TO_LEFT
        );
        //_m_cell_handle = CGAL::make_object(_m_halfedge_handle->face());
    }
    
    /*!
     * Notification before the splitting of an edge into two.
     * \param e A handle to one of the existing halfedges.
     * \param v A vertex representing the split point.
     * \param c1 The x-monotone curve to be associated with the first edge.
     * \param c2 The x-monotone curve to be associated with the second edge.
     */
    virtual void before_split_edge (Halfedge_handle e,
                                    Vertex_handle /* v */,
                                    const X_monotone_curve_2& /* c1 */,
                                    const X_monotone_curve_2& /* c2 */)
    {
        Halfedge_handle curr = 
            Base::arrangement()->non_const_handle(_m_halfedge_handle);
        if (curr->is_fictitious()) {
            std::cout << "curr fict";
        } else {
            std::cout << "curr: " << curr->curve() << std::endl;
        }
        std::cout << std::endl;
        if (e == curr || e->twin() == curr) {
            std::cout << "It's maybe required UPDATE _m_halfedge_handle"
                      << std::endl;
            _m_check_split_edges = true;
        }
        
    }
    
    /*!
     * Notification after an edge was split.
     * \param e1 A handle to one of the twin halfedges forming the first edge.
     * \param e2 A handle to one of the twin halfedges forming the second edge.
     */
    virtual void after_split_edge (Halfedge_handle e1,
                                   Halfedge_handle e2)
    {
        typename Geometry_traits_2::Is_vertical_2 is_vertical =
            Base::arrangement()->geometry_traits()->is_vertical_2_object();
        typename Traits_adaptor_2::Is_in_x_range_2 is_in_x_range =
            static_cast< Traits_adaptor_2* >
            (Base::arrangement()->geometry_traits())->is_in_x_range_2_object();
        
        typename Geometry_traits_2::Compare_y_at_x_2 compare_y_at_x =
            Base::arrangement()->geometry_traits()->compare_y_at_x_2_object();
        
        X_monotone_curve_2 cv1 = e1->curve();
        X_monotone_curve_2 cv2 = e2->curve();
        
        if (_m_check_split_edges) {
            std::cout << "pt : " << _m_point << std::endl;
            std::cout << "cv1: " << cv1 << std::endl;
            std::cout << "cv2: " << cv2 << std::endl;
            
            // Remark:
            // as e1 and e2 originate from _m_halfedge_handle
            // we only have to check x-ranges ...
            if (!is_vertical(cv1)) {
                CGAL_assertion(!is_vertical(cv2));
                // TODO symbolic pertubation for endpoints!
                if (is_in_x_range(cv1, _m_point)) {
                    _m_halfedge_handle = 
                        ((e1->direction() == CGAL::ARR_RIGHT_TO_LEFT) ?
                         e1 : e1->twin());

                    CGAL_assertion(
                            _m_halfedge_handle->direction() == 
                            CGAL::ARR_RIGHT_TO_LEFT
                    );
                    
                    _m_cell_handle = 
                        CGAL::make_object(_m_halfedge_handle->face());
                    
                    std::cout << "HE3 set to cv: " 
                              << _m_halfedge_handle->curve() << std::endl;

                } else {
                    CGAL_assertion(is_in_x_range(cv2, _m_point));
                    _m_halfedge_handle = 
                        ((e2->direction() == CGAL::ARR_RIGHT_TO_LEFT) ?
                         e2 : e2->twin());
                    
                    CGAL_assertion(
                            _m_halfedge_handle->direction() == 
                            CGAL::ARR_RIGHT_TO_LEFT
                    );
                    
                    _m_cell_handle = 
                        CGAL::make_object(_m_halfedge_handle->face());
                    
                    std::cout << "HE4 set to cv: " 
                              << _m_halfedge_handle->curve() << std::endl;
                    
                }
            } else {
                // or something in vertical case
                // TODO what if e->curve() is vertical
            }
        }
        
        // finally reset trigger
        _m_check_split_edges = false;
    }
    
    

    /*!
     * Notification before the splitting of a fictitious edge into two.
     * \param e A handle to one of the existing halfedges.
     * \param v A vertex representing the unbounded split point.
     */
    virtual void before_split_fictitious_edge (Halfedge_handle /* e */,
                                               Vertex_handle /* v */)
    {}
    
    /*!
     * Notification after a fictitious edge was split.
     * \param e1 A handle to one of the twin halfedges forming the first edge.
     * \param e2 A handle to one of the twin halfedges forming the second edge.
     */
    virtual void after_split_fictitious_edge (Halfedge_handle /* e1 */,
                                              Halfedge_handle /* e2 */)
    {}

    //!@}
    
   //! returns cell_handle to cell of given point
    CGAL::Object cell_handle() {
        return _m_cell_handle;
    }
    
private:
    //! query point
    Point_2 _m_point;
    
    //! halfedge "immediately above point"
    mutable Halfedge_const_handle _m_halfedge_handle;
    
    //! cell handle to compute
    CGAL::Object _m_cell_handle;

    // helper for split edge
    //! if true has to choose among two split edges
    bool _m_check_split_edges;
};



/*! \class Functor to compute single cell of point */
template < class Arrangement_2_ >
class Construct_single_cell_2 {
    
public:

    //! this instance's first template parameter
    typedef Arrangement_2_ Arrangement_2;

    //! the class itself
    typedef Construct_single_cell_2< Arrangement_2 > Self;

    //! geometric traits class
    typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;

    //! type of point
    typedef typename Geometry_traits_2::Point_2 Point_2;

    //! type of x-monotone curve
    typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_curve_2;

    //! type of curve
    typedef typename Geometry_traits_2::Curve_2 Curve_2;


    typedef typename Arrangement_2::Vertex_const_handle Vertex_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle 
    Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle Face_const_handle;

    // TASK select best point location strategy
    //! type of point location strategy
    typedef CGAL::Arr_naive_point_location< Arrangement_2 > Point_location;
    
    // DO NOT WORK FOR UNBOUNDED!
    //typedef CGAL::Arr_simple_point_location< Arrangement_2 > Point_location;
    //typedef Arr_walk_along_line_point_location< Arrangement_2 > 
    //Point_location;

    //!\name Constructors
    //!@{
    
    template < class InputIterator >
    Construct_single_cell_2(InputIterator begin, InputIterator end) :
        _m_cell_handle_pl(boost::none),
        _m_cell_handle_ri(boost::none),
        _m_cell_handle_rbo(boost::none) {
        
        // TODO check that iteratortype is CGAL::Object
        _m_objects.reserve(std::distance(begin, end));
        
#if !NDEBUG 
        std::cout << "Created Construct_single_cell_2 instance for "
                  << std::distance(begin, end) << " input objects."
                  << std::endl;
#endif

        for (InputIterator it = begin; it != end; it++) {
            _m_objects.push_back(*it);
            
            X_monotone_curve_2 curr_xcurve;
            Point_2 curr_point;
            Curve_2 curr_curve;
            
            if (CGAL::assign(curr_curve, *it)) {
                std::list< CGAL::Object > tmp;
#if CGAL_USE_ACK_2
                typename Geometry_traits_2::Make_x_monotone_2
                    make_x_monotone =
                    Geometry_traits_2::instance().make_x_monotone_2_object();
#else
                typename Geometry_traits_2::Make_x_monotone_2
                    make_x_monotone;
#endif
                make_x_monotone(curr_curve, std::back_inserter(tmp));
                
                for (typename std::list< CGAL::Object >::const_iterator 
                         iit = tmp.begin();
                     iit != tmp.end(); iit++) {
                    if (CGAL::assign(curr_xcurve,*iit)) {
                        _m_xcvs.push_back(curr_xcurve);
                    } else {
                        CGAL_assertion_code(bool check = )
                            CGAL::assign(curr_point,*iit);
                        CGAL_assertion(check);
                        _m_pts.push_back(curr_point);
                    }
                }
            } else if (CGAL::assign(curr_xcurve,*it)) {
                _m_xcvs.push_back(curr_xcurve);
            } else {
                CGAL_assertion_code(bool check = )
                    CGAL::assign(curr_point,*it);
                CGAL_assertion(check);
                _m_pts.push_back(curr_point);
            }
        }
    }
    
    //!@}
    
public:
    
    //!\name cell localizations
    //!@{

    //! returns the cell using point location
    CGAL::Object cell_pl(const Point_2& pt) const {
        if (!_m_cell_handle_pl) {
#if !NDEBUG
            std::cout << "Computing cell with POINT-LOCATION-method ... " 
                      << std::flush;
#endif
            // compute full arr
            this->_m_full_arr = Arrangement_2();
            _mt_full_arr.start();
            CGAL::insert_empty(*_m_full_arr,
                               _m_xcvs.begin(), _m_xcvs.end(),
                               _m_pts.begin(), _m_pts.end()
            );
            _mt_full_arr.stop();
            
            std::cout << "The full_arr sizes:" << std::endl
                      << "   V = " << _m_full_arr->number_of_vertices()
                      << ",  E = " << _m_full_arr->number_of_edges() 
                      << ",  F = " << _m_full_arr->number_of_faces() 
                      << std::endl;
            
            // locate point
            Point_location pl(*_m_full_arr);
            _mt_pl.start();
            _m_cell_handle_pl = pl.locate(pt);
            _mt_pl.stop();
            
#if !NDEBUG
            std::cout << "done."   << std::endl << std::endl;
#endif
            std::cout << "tFullArr: " << _mt_full_arr.time() 
                      << " sec" << std::endl;
            std::cout << "tPL     : " << _mt_pl.time() 
                      << " sec" << std::endl;
        }
        return *_m_cell_handle_pl;
    }

public:
    //! returns the cell using random incremental
    CGAL::Object cell_ri(const Point_2& pt) const {
        if (!_m_cell_handle_ri) {
#if !NDEBUG
            std::cout << "Computing cell with RANDOM_INCREMENTAL-method ... " 
                      << std::flush;
#endif
            _m_arr_city = Arrangement_2();
            
            CGAL::Object cell_handle;            
            
            // let's start with all points
            CGAL::non_intersecting_insert_empty(
                    *_m_arr_city,
                    _m_xcvs.begin(), _m_xcvs.begin(), // NO CURVE
                    _m_pts.begin(), _m_pts.end()
            );
            
            Point_location pl(*_m_arr_city);
            cell_handle = pl.locate(pt);
            
            Face_const_handle fh;
            if (!CGAL::assign(fh, cell_handle)) {
                // simple case
                _m_cell_handle_ri = cell_handle; // found point

            } else {
                
                // TODO reactivate random shuffle
                //std::random_shuffle(_m_xcvs.begin(), _m_xcvs.end());
                
                Face_const_handle fh = _m_arr_city->faces_begin();
                cell_handle = CGAL::make_object(fh);
                
                // real RI-case
                for (typename 
                         std::vector< X_monotone_curve_2 >::const_iterator 
                         cit = _m_xcvs.begin(); cit != _m_xcvs.end(); cit++) {

                    std::cout << "----------------------------------------"
                              << std::endl;
                    std::cout << "insert: " << *cit << std::endl; 
                    std::cout << std::endl;
                    
#if !CGAL_SINGLE_CELL_RI_NAIVE
                    CGALi::RI_observer< Arrangement_2 > 
                        obs((*_m_arr_city), pt);
#endif
                    // add *cit using zone
                    CGAL::insert(*_m_arr_city, *cit);

#if CGAL_SINGLE_CELL_RI_NAIVE
                    // make point location
                    cell_handle = pl.locate(pt);
#else
                    // TODO remove point location and replace by
                    //      observer! observer should also try to determine
                    //      whether point lies on new curve 
                    //      (using e.g., before_new_vertex!)
                    cell_handle = obs.cell_handle();
#endif
                    
                    // simplify!
                    Arrangement_2 new_city;
                    cell_arr(cell_handle, new_city);
                    _m_arr_city = new_city;

                    // TODO alternative: 
                    //      delete all feature not defining "cell_handle"
                }

#if CGAL_SINGLE_CELL_RI_NAIVE
#else
                // TODO remove this location (ie, improve arr-maintenance)
                //cell_handle = pl.locate(pt);
#endif
                _m_cell_handle_ri = cell_handle;
            }
            
#if !NDEBUG
            std::cout << "done."   << std::endl;
#endif
        }
        return *_m_cell_handle_ri;
    }

    CGAL::Object cell_rbo_naive(const Point_2& pt) const {
        if (!_m_cell_handle_rbo) {
#if !NDEBUG
            std::cout << "Computing cell for " 
                      << (_m_xcvs.size() + _m_pts.size())
                      << " input objects " 
                      << "with NAIVE-RED-BLUE-OVERLAY-method ... " 
                      << std::flush;
#endif
            // TODO permute INPUT randomly!!
            
            if (_m_xcvs.size() + _m_pts.size() <= 4) {
                
#if !NDEBUG
                std::cout << "Anchor" << std::endl;
#endif
                
                _mt_rec_anchor.start();
                CGAL::Object cell_handle = cell_pl(pt);
                _mt_rec_anchor.stop();
                
                _m_cell_handle_rbo = cell_handle;
                
            } else {
                
#if !NDEBUG
                std::cout << "Red-blue split" << std::endl;
#endif
                // permute input
                // TODO use better random values
                std::random_shuffle(_m_xcvs.begin(), _m_xcvs.end());
                std::random_shuffle(_m_pts.begin(), _m_pts.end());

                // split input into two sets
                std::vector< X_monotone_curve_2 > xcvs[2];
                typename 
                    std::vector< X_monotone_curve_2 >::iterator 
                    xcvs_mid =
                    _m_xcvs.begin();
                std::advance(xcvs_mid, (_m_xcvs.size() / 2));
                xcvs[0].reserve(std::distance(_m_xcvs.begin(), xcvs_mid));
                std::copy(_m_xcvs.begin(), xcvs_mid, 
                          std::back_inserter(xcvs[0]));
                xcvs[1].reserve(std::distance(xcvs_mid, _m_xcvs.end()));
                std::copy(xcvs_mid, _m_xcvs.end(), 
                          std::back_inserter(xcvs[1]));
                
                std::vector< Point_2 > pts[2];
                typename std::vector< Point_2 >::iterator pts_mid =
                    _m_pts.begin();
                std::advance(pts_mid, (_m_pts.size() / 2));
                pts[0].reserve(std::distance(_m_pts.begin(), pts_mid));
                std::copy(_m_pts.begin(), pts_mid, std::back_inserter(pts[0]));
                pts[1].reserve(std::distance(pts_mid, _m_pts.end()));
                std::copy(pts_mid, _m_pts.end(), std::back_inserter(pts[1]));
                
                CGAL::Object cell_handle[2];
                
                Arrangement_2 cell[2];
                
                for (int i = 0; i < 2; i++) {
                    
                    std::list< CGAL::Object > objects;
                    for (typename 
                             std::vector< X_monotone_curve_2 >::const_iterator
                             it = xcvs[i].begin(); it != xcvs[i].end(); it++) {
                        objects.push_back(CGAL::make_object(*it));
                    }
                    for (typename std::vector< Point_2 >::const_iterator
                             it = pts[i].begin(); it != pts[i].end(); it++) {
                        objects.push_back(CGAL::make_object(*it));
                    }
                    
                    Self recursive(objects.begin(), objects.end());
                    
                    cell_handle[i] = recursive.cell_rbo_naive(pt);

                    cell_arr(cell_handle[i], cell[i]);
                }
                
#if !NDEBUG
                std::cout << "Start overlay ... " << std::flush;
#endif

                _mt_rec_overlay.start();
                CGAL::Arr_default_overlay_traits< Arrangement_2 > ovltraits;

                Arrangement_2 arr_purple;

                CGAL::overlay(cell[0], cell[1], 
                              arr_purple, ovltraits);
                _mt_rec_overlay.stop();
                
                _m_arr_purple = arr_purple;
                
#if !NDEBUG
                std::cout << "done." << std::endl;
#endif
                
                _mt_pl.start();
                Point_location pl_purple(*_m_arr_purple);
                CGAL::Object cell_handle_pl_purple = pl_purple.locate(pt);
                _mt_pl.stop();            
                
                _m_cell_handle_rbo = cell_handle_pl_purple;
            }
#if !NDEBUG
            std::cout << "done." << std::endl;
            std::cout << std::endl;
#endif
            std::cout << "tRecAnchor : " << _mt_rec_anchor.time() 
                      << " sec" << std::endl;
            std::cout << "tRecOverlay: " << _mt_rec_overlay.time() 
                      << " sec" << std::endl;
            std::cout << "tPLs       : " << _mt_pl.time() 
                      << " sec" << std::endl;
        }
        return *_m_cell_handle_rbo;
    }


    //! returns the cell using red-blue overlay
    CGAL::Object cell_rbo(const Point_2& pt) const {
        if (!_m_cell_handle_rbo) {
#if !NDEBUG
            std::cout << "Computing cell with RED-BLUE-OVERLAY-method ... " 
                      << std::flush;
#endif
            _m_cell_handle_rbo = CGAL::Object();
            
#if !NDEBUG
            std::cout << "done."   << std::endl;
#endif
        }
        return *_m_cell_handle_rbo;
    }
    
    //!@}

    //!\name Helpers
    //!@{

    /*!\brief 
     * converts a cell-handle (face, edge, vertex) into its induced arrangement
     */
    // FUTURE TODO allow multiple cell_handles!
    void cell_arr(CGAL::Object cell_handle, Arrangement_2& cell) const {
        
        _mt_cell.start();
        
        std::list< Point_2 > cell_pts;
        std::list< X_monotone_curve_2 > cell_xcvs;
        
        Vertex_const_handle vh;
        Halfedge_const_handle heh;
        Face_const_handle fh;
        if (CGAL::assign(vh, cell_handle)) {
            cell_pts.push_back(vh->point());
        } else if (CGAL::assign(heh, cell_handle)) {
            cell_xcvs.push_back(heh->curve());
        } else {
            CGAL_assertion_code(bool check =)
                CGAL::assign(fh, cell_handle);
            CGAL_assertion(check);
            
            // copy curves of CCBs of face
            typedef typename Arrangement_2::Outer_ccb_const_iterator 
                Outer_ccb_const_iterator;
            
            typedef typename Arrangement_2::Inner_ccb_const_iterator 
                Inner_ccb_const_iterator;
            
            typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                Ccb_halfedge_const_circulator;
            
            typedef typename Arrangement_2::Isolated_vertex_const_iterator 
                Isolated_vertex_const_iterator;
            
            for (Outer_ccb_const_iterator ocb = fh->outer_ccbs_begin();
                 ocb != fh->outer_ccbs_end(); 
                 ocb++) {
                Ccb_halfedge_const_circulator he = *ocb;
                if (!he->is_fictitious()) {
                    cell_xcvs.push_back(he->curve());
                }
                he++;
                for (; he != *ocb; he++) {
                    if (!he->is_fictitious()) {
                        cell_xcvs.push_back(he->curve());
                    }
                }
            }
            
            for (Inner_ccb_const_iterator icb = fh->inner_ccbs_begin();
                 icb != fh->inner_ccbs_end(); 
                 icb++) {
                Ccb_halfedge_const_circulator he = *icb;
                if (!he->is_fictitious()) {
                    cell_xcvs.push_back(he->curve());
                }
                he++;
                for (; he != *icb; he++) {
                    if (!he->is_fictitious()) {
                        cell_xcvs.push_back(he->curve());
                    }
                }
            }
            
            // copy isolated points of face to cell
            for (Isolated_vertex_const_iterator vt = 
                     fh->isolated_vertices_begin(); 
                 vt != fh->isolated_vertices_end(); vt++) {
                cell_pts.push_back(vt->point());
            }
        }
        
        std::cout << "#cell-curves:" << cell_xcvs.size() << std::endl;
        std::cout << "#cell-points:" << cell_pts.size() << std::endl;

        for (typename 
                 std::list< X_monotone_curve_2 >::const_iterator
                 it = cell_xcvs.begin(); it != cell_xcvs.end(); it++) {
            //std::cout << "CURVE: " << *it << std::endl;
            //std::cout << "poly: " << it->curve().polynomial_2() << std::endl;
            
        }
        
        CGAL::set_pretty_mode(std::cerr);
        // TODO use CGAL::non_intersecting_insert_empty(
        CGAL::insert_empty(
                cell, 
                cell_xcvs.begin(), cell_xcvs.end(),
                cell_pts.begin(), cell_pts.end()
        );
        
        _mt_cell.stop();
        
        std::cout << "tCell   : " << _mt_cell.time() 
                  << " sec" << std::endl;
    }

    //!@}

    
private:

    //////////////////////////////////////////////////////////////////////////
    // members

    //! input objects
    std::vector< CGAL::Object > _m_objects;
    
    //! input curves
    mutable std::vector< X_monotone_curve_2 > _m_xcvs;

    //! input points
    mutable std::vector< Point_2 > _m_pts;
    
    //! the cell
    mutable boost::optional< CGAL::Object > _m_cell_handle_pl; 
    mutable boost::optional< CGAL::Object > _m_cell_handle_ri; 
    mutable boost::optional< CGAL::Object > _m_cell_handle_rbo; 

    // helper for pl
    mutable boost::optional< Arrangement_2 > _m_full_arr; 

    // helper for ri
    mutable boost::optional< Arrangement_2 > _m_arr_city;

    // helper for rbo_naive
    mutable boost::optional< Arrangement_2 > _m_arr_purple; 

    // timers
    mutable CGAL::Timer _mt_full_arr;
    mutable CGAL::Timer _mt_pl;
    mutable CGAL::Timer _mt_cell;

    mutable CGAL::Timer _mt_rec_anchor;
    mutable CGAL::Timer _mt_rec_overlay;
};

} // namespace CGALi

/*!
 * Construct the single cell containing a point using point location
 * \param point The reference point
 * \param begin An iterator for the first input object defining the full 
 *              arrangement
 * \param end A past-the-end iterator for the input objects defining the
 *            full arrangement
 * \param cell Output: The cell of Arr(begin,end) containing point as 
 *                     arrangement
 * \pre The value-type of InputIterator is CGAL::Object which can be passed to
 *      GeoTraits_2::Make_x_monotone_2()
 */
template < typename InputIterator, typename GeoTraits_2 >
CGAL::Object single_cell_pl_2(
        typename GeoTraits_2::Point_2 point,
        InputIterator begin, InputIterator end,
        CGAL::Arrangement_2< GeoTraits_2 >& cell) {
    
    typedef GeoTraits_2                                       Geo_traits_2;
    typedef CGAL::Arrangement_2< Geo_traits_2 >               Arrangement_2;
    
    typedef CGALi::Construct_single_cell_2< Arrangement_2 >   
        Construct_single_cell_2;
    
    Construct_single_cell_2 single_cell(
            begin, end
    );

     CGAL::Object cell_handle = single_cell.cell_pl(point);
     single_cell.cell_arr(cell_handle, cell);
    return cell_handle;
 }

/*!
 * Construct the single cell containing a point using random incremental 
 * \param point The reference point
 * \param begin An iterator for the first input object defining the full 
 *              arrangement
 * \param end A past-the-end iterator for the input objects defining the
 *            full arrangement
 * \param cell Output: The cell of Arr(begin,end) containing point as 
 *                     arrangement
 * \pre The value-type of InputIterator is CGAL::Object which can be passed to
 *      GeoTraits_2::Make_x_monotone_2()
 */
template < typename InputIterator, typename GeoTraits_2 >
CGAL::Object single_cell_ri_2(
        typename GeoTraits_2::Point_2 point,
        InputIterator begin, InputIterator end,
        CGAL::Arrangement_2< GeoTraits_2 >& cell) {
    
    typedef GeoTraits_2                                       Geo_traits_2;
    typedef CGAL::Arrangement_2< Geo_traits_2 >               Arrangement_2;
    
    typedef CGALi::Construct_single_cell_2< Arrangement_2 >   
        Construct_single_cell_2;
    
    Construct_single_cell_2 single_cell(
            begin, end
    );

    CGAL::Object cell_handle = single_cell.cell_ri(point);
    std::cout << "located" << std::endl;
    single_cell.cell_arr(cell_handle, cell);
    return cell_handle;
 }


/*!
 * Construct the single cell containing a point using naive red-blue overlay
 * \param point The reference point
 * \param begin An iterator for the first input object defining the full 
 *              arrangement
 * \param end A past-the-end iterator for the input objects defining the
 *            full arrangement
 * \param cell Output: The cell of Arr(begin,end) containing point as 
 *                     arrangement
 * \pre The value-type of InputIterator is CGAL::Object which can be passed to
 *      GeoTraits_2::Make_x_monotone_2()
 */
template < typename InputIterator, typename GeoTraits_2 >
CGAL::Object single_cell_rbo_naive_2(
        typename GeoTraits_2::Point_2 point,
        InputIterator begin, InputIterator end,
        CGAL::Arrangement_2< GeoTraits_2 >& cell) {
    
    typedef GeoTraits_2                                       Geo_traits_2;
    typedef CGAL::Arrangement_2< Geo_traits_2 >               Arrangement_2;
    
    typedef CGALi::Construct_single_cell_2< Arrangement_2 >   
        Construct_single_cell_2;
    
    Construct_single_cell_2 single_cell(
            begin, end
    );
    
    CGAL::Object cell_handle = single_cell.cell_rbo_naive(point);
    single_cell.cell_arr(cell_handle, cell);
    return cell_handle;
}


CGAL_END_NAMESPACE

#endif // CGAL_ARR_SINGLE_CELL_2_H
// EOF

