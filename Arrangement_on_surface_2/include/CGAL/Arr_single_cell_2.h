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


#ifndef CGAL_ARR_SINGLE_CELL_2_H
#define CGAL_ARR_SINGLE_CELL_2_H

#include <boost/optional.hpp>
#include <boost/none.hpp>

#include <CGAL/Arrangement_2.h>

// replace point location strategy
#include <CGAL/Arr_naive_point_location.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*! \class Functor to compute single cell of point */
template < class Arrangement_2_ >
class Construct_single_cell_2 {
    
public:

    //! this instance's first template parameter
    typedef Arrangement_2_ Arrangement_2;

    //! geometric traits class
    typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;

    //! type of point
    typedef typename Geometry_traits_2::Point_2 Point_2;

    //! type of x-monotone curve
    typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_curve_2;

    //! type of curve
    typedef typename Geometry_traits_2::Curve_2 Curve_2;

    // TODO remove tags?
    // tags
    struct Construction_method {};
    
    struct Point_location : public Construction_method {};

    struct Random_incremental : public Construction_method {};
    
    struct Red_blue_overlay : public Construction_method {};
    
    
    template < class InputIterator >
    Construct_single_cell_2(InputIterator begin, InputIterator end) :
        _m_cell_pl(boost::none),
        _m_cell_ri(boost::none),
        _m_cell_rbo(boost::none) {


        for (InputIterator it = begin; it != end; it++) {
            X_monotone_curve_2 curr_xcurve;
            Point_2 curr_point;
            Curve_2 curr_curve;
            
            if (CGAL::assign(curr_curve, *it)) {
                std::list< CGAL::Object > tmp;
                typename Geometry_traits_2::Make_x_monotone_2(
                        Geometry_traits_2::instance().
                        make_x_monotone_2_object()
                )(
                        curr_curve, std::back_inserter(tmp)
                );
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
            // TODO deal with Curve_2 using Make_x_monotone_2
        }
    }
    
    //! returns the cell using point location
    const Arrangement_2& cell_pl(const Point_2& pt) const {
        if (!_m_cell_pl) {
            Point_location pl;
            this->_compute_cell(pt, pl);
        }
        return *_m_cell_pl;
    }

    //! returns the cell using random incremental
    const Arrangement_2& cell_ri(const Point_2& pt) const {
        if (!_m_cell_ri) {
            Random_incremental ric;
            this->_compute_cell(pt, ric);
        }
        return *_m_cell_ri;
    }

    //! returns the cell using red-blue overlay
    const Arrangement_2& cell_rbo(const Point_2& pt) const {
        if (!_m_cell_rbo) {
            Red_blue_overlay rbo;
            this->_compute_cell(pt, rbo);
        }
        return *_m_cell_rbo;
    }
    
private:

    //////////////////////////////////////////////////////////////////////////
    // compute cell 
    void _compute_cell(const Point_2& pt, Point_location method) const {
        
#if !NDEBUG
        std::cout << "Computing cell with POINT-LOCATION-method ... " 
                  << std::flush;
#endif
        
        // TODO replace point location strategy
        typedef CGAL::Arr_naive_point_location< Arrangement_2 > PL;
        
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
        PL pl(*_m_full_arr);
        _mt_pl.start();
        CGAL::Object obj = pl.locate(pt);
        _mt_pl.stop();
        
        _mt_cell.start();
        
        std::list< Point_2 > cell_pts;
        std::list< X_monotone_curve_2 > cell_xcvs;

        typename Arrangement_2::Vertex_const_handle vh;
        typename Arrangement_2::Halfedge_const_handle heh;
        typename Arrangement_2::Face_const_handle fh;
        if (CGAL::assign(vh, obj)) {
            cell_pts.push_back(vh->point());
        } else if (CGAL::assign(heh, obj)) {
            cell_xcvs.push_back(heh->curve());
        } else {
            CGAL_assertion_code(bool check =)
                CGAL::assign(fh, obj);
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
                    cell_xcvs.push_back (he->curve());
                }
                he++;
                for (; he != *ocb; he++) {
                    if (!he->is_fictitious()) {
                        cell_xcvs.push_back (he->curve());
                    }
                }
            }
            for (Inner_ccb_const_iterator icb = fh->inner_ccbs_begin();
                 icb != fh->inner_ccbs_end(); 
                 icb++) {
                Ccb_halfedge_const_circulator he = *icb;
                if (!he->is_fictitious()) {
                    cell_xcvs.push_back (he->curve());
                }
                he++;
                for (; he != *icb; he++) {
                    if (!he->is_fictitious()) {
                        cell_xcvs.push_back (he->curve());
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
        
        Arrangement_2 cell;
        
        CGAL::non_intersecting_insert_empty(
                cell, 
                cell_xcvs.begin(), cell_xcvs.end(),
                cell_pts.begin(), cell_pts.end()
        );
        
        _mt_cell.stop();

        _m_cell_pl = cell;
        
#if !NDEBUG
        std::cout << "done."   << std::endl;
#endif
        std::cout << "tFullArr: " << _mt_full_arr.time() 
                  << " sec" << std::endl;
        std::cout << "tPl     : " << _mt_pl.time() 
                  << " sec" << std::endl;
        std::cout << "tCell   : " << _mt_cell.time() 
                  << " sec" << std::endl;
        
    }

    void _compute_cell(const Point_2& pt, Random_incremental method) const {
        
#if !NDEBUG
        std::cout << "Computing cell with RANDOM_INCREMENTAL-method ... " 
                  << std::flush;
#endif
        _m_cell_ri = Arrangement_2();

#if !NDEBUG
        std::cout << "done."   << std::endl;
#endif
        
    }

    void _compute_cell(const Point_2& pt, Red_blue_overlay method) const {
        
#if !NDEBUG
        std::cout << "Computing cell with RED-BLUE-OVERLAY-method ... " 
                  << std::flush;
#endif
        _m_cell_rbo = Arrangement_2();

#if !NDEBUG
        std::cout << "done."   << std::endl;
#endif
        
    }

    //////////////////////////////////////////////////////////////////////////
    // members
    //! construction
    Construction_method _m_method;
    
    //! input curves
    std::list< X_monotone_curve_2 > _m_xcvs;

    //! input points
    std::list< Point_2 > _m_pts;
    
    //! the cell
    mutable boost::optional< Arrangement_2 > _m_cell_pl; 
    mutable boost::optional< Arrangement_2 > _m_cell_ri; 
    mutable boost::optional< Arrangement_2 > _m_cell_rbo; 

    // helper for pl
    mutable boost::optional< Arrangement_2 > _m_full_arr; 

    
    // timers
    mutable CGAL::Timer _mt_full_arr;
    mutable CGAL::Timer _mt_pl;
    mutable CGAL::Timer _mt_cell;
};

} // namespace CGALi

/*!
 * Construct the single cell containing a point
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
void single_cell_2(typename GeoTraits_2::Point_2 point,
                   InputIterator begin, InputIterator end,
                   CGAL::Arrangement_2< GeoTraits_2 >& cell) {
    typedef GeoTraits_2                                       Geo_traits_2;
    typedef CGAL::Arrangement_2< Geo_traits_2 >               Arrangement_2;
    
    typedef CGALi::Construct_single_cell_2< Arrangement_2 >   
        Construct_single_cell_2;
    
    Construct_single_cell_2 single_cell(
            begin, end
    );

    // TODO write three functions!
    cell = single_cell.cell_pl(point);
    
    //cell = single_cell.cell_ri(point);
    
    //cell = single_cell.cell_rbo(point);
}


CGAL_END_NAMESPACE


#if 0


#include <CGAL/Envelope_3/Envelope_divide_and_conquer_3.h>
#include <CGAL/Envelope_3/Envelope_pm_dcel.h>
#include <CGAL/Envelope_3/Envelope_overlay_functor.h>
#include <CGAL/Arr_accessor.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of an envelope diagram (a minimization diagram or a
 * maximization diagram).
 */
template <typename T_Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T1, class T2>
#endif
          class T_Dcel = Envelope_pm_dcel>
class Envelope_diagram_2 :
  public Arrangement_2<T_Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
                       T_Dcel<T_Traits,
                              typename T_Traits::Xy_monotone_surface_3>
#else
                       typename T_Dcel::template Dcel<T_Traits,
                                                      typename T_Traits::Xy_monotone_surface_3>
#endif
                       >
{
public:
  typedef T_Traits                                      Traits_3;
  typedef typename Traits_3::Xy_monotone_surface_3      Xy_monotone_surface_3;

protected:
  typedef T_Dcel<Traits_3, Xy_monotone_surface_3>       Env_dcel;
  typedef Envelope_diagram_2<Traits_3, T_Dcel>          Self;
  friend class Arr_accessor<Self>;

public:
  typedef Arrangement_2<Traits_3, Env_dcel>             Base;
  typedef typename Env_dcel::Dcel_data_const_iterator   Surface_const_iterator;

  /*! Default constructor. */
  Envelope_diagram_2() :
    Base()
  {}

  /*! Constructor with a traits-class instance. */
  Envelope_diagram_2 (Traits_3* tr) :
    Base (tr)
  {}

};

/*!
 * Construct the lower envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <typename InputIterator, typename T_Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T1, class T2>
#endif
          class T_Dcel>
void lower_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<T_Traits, T_Dcel> & min_diagram)
{
  typedef T_Traits                                            Traits_3;
  typedef typename Envelope_diagram_2<Traits_3, T_Dcel>::Base Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.traits(), ENVELOPE_LOWER);
  env_alg.construct_lu_envelope (begin, end, min_diagram);
}

/*!
 * Construct the lower envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <typename InputIterator, typename T_Traits>
void lower_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<T_Traits> & min_diagram)
{
  typedef T_Traits                                            Traits_3;
  typedef typename Envelope_diagram_2<Traits_3>::Base         Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.traits(), ENVELOPE_LOWER);
  env_alg.construct_lu_envelope (begin, end, min_diagram);
}

/*!
 * Construct the upper envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T1, class T2>
#endif
          class T_Dcel>
void upper_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits, T_Dcel>& max_diagram)
{
  typedef Traits                                              Traits_3;
  typedef typename Envelope_diagram_2<Traits_3, T_Dcel>::Base Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.traits(), ENVELOPE_UPPER);
  env_alg.construct_lu_envelope (begin, end, max_diagram);
}

/*!
 * Construct the upper envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits>
void upper_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits>& max_diagram)
{
  typedef Traits                                        Traits_3;
  typedef typename Envelope_diagram_2<Traits_3>::Base   Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2>
                                                        Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.traits(), ENVELOPE_UPPER);
  env_alg.construct_lu_envelope (begin, end, max_diagram);
}

/*!
 * Construct the lower envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy-monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T1, class T2>
#endif
          class T_Dcel>
void
lower_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                              Envelope_diagram_2<Traits, T_Dcel>& min_diagram)
{
  typedef Traits                                              Traits_3;
  typedef typename Envelope_diagram_2<Traits_3, T_Dcel>::Base Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.traits(), ENVELOPE_LOWER);
  env_alg.construct_envelope_xy_monotone (begin, end, min_diagram);
}

/*!
 * Construct the lower envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy-monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits>
void lower_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                                   Envelope_diagram_2<Traits>& min_diagram)
{
  typedef Traits                                              Traits_3;
  typedef typename Envelope_diagram_2<Traits_3>::Base         Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.traits(), ENVELOPE_LOWER);
  env_alg.construct_envelope_xy_monotone (begin, end, min_diagram);
}

/*!
 * Construct the upper envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy_monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T1, class T2>
#endif
          class T_Dcel>
void
upper_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                              Envelope_diagram_2<Traits, T_Dcel>& max_diagram)
{
  typedef Traits                                              Traits_3;
  typedef typename Envelope_diagram_2<Traits_3, T_Dcel>::Base Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Base_arr_2> Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.traits(), ENVELOPE_UPPER);
  env_alg.construct_envelope_xy_monotone (begin, end, max_diagram);

  return;
}

/*!
 * Construct the upper envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy_monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class Traits>
void upper_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                                   Envelope_diagram_2<Traits>& max_diagram)
{
  typedef Traits                                       Traits_3;
  typedef typename Envelope_diagram_2<Traits_3>::Base  Base_arr_2;
  typedef Envelope_divide_and_conquer_3<Traits_3,
                                        Base_arr_2>    Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.traits(), ENVELOPE_UPPER);
  env_alg.construct_envelope_xy_monotone (begin, end, max_diagram);

  return;
}

CGAL_END_NAMESPACE

#endif

#endif // CGAL_ARR_SINGLE_CELL_2_H

