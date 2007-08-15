// Copyright (c) 2007 Max-Planck-Institute Saarbruecken (Germany).
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

#ifndef CGAL_ARR_TORUS_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_TORUS_TOPOLOGY_TRAITS_2_H

/*! \file
 * Definition of the Arr_torus_topology_traits_2<GeomTraits> class.
 */

#include <map>
#include <algorithm>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_default_dcel.h>
#if 1
// TASK replace point location strategy
#include <CGAL/Arr_naive_point_location.h>
#else
//#include <CGAL/Arr_walk_along_line_point_location.h>
#endif
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_observed_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_overlay_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_overlay_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_overlay_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_batched_pl_sl_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>

#include <CGAL/Arr_topology_traits/Arr_torus_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_torus_insertion_helper.h>
#include <CGAL/Arr_topology_traits/Arr_torus_overlay_helper.h>
#include <CGAL/Arr_topology_traits/Arr_torus_batched_pl_helper.h>
#include <CGAL/Arr_topology_traits/Arr_planar_inc_insertion_zone_visitor.h>

// TODO helpers

CGAL_BEGIN_NAMESPACE

// Forward declaration:
template <class GeomTraits_, class TopTraits_> 
class Arrangement_on_surface_2;

/*! \class Arr_torus_topology_traits_2
 * A topology-traits class that encapsulates the embedding of 2D arrangements
 * on a torus
 */
template <class GeomTraits_,
          class Dcel_ = Arr_default_dcel<GeomTraits_> >
class Arr_torus_topology_traits_2
{
    
public:
    ///! \name The geometry-traits types.
    //@{
    typedef GeomTraits_                                     Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2             Point_2;
    typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
    //@}
    
    ///! \name The DCEL types.
    //@{
    typedef Dcel_                                           Dcel;
    typedef typename Dcel::Size                             Size;
    typedef typename Dcel::Vertex                           Vertex;
    typedef typename Dcel::Halfedge                         Halfedge;
    typedef typename Dcel::Face                             Face;
    typedef typename Dcel::Outer_ccb                        Outer_ccb;
    typedef typename Face::Outer_ccb_const_iterator   Outer_ccb_const_iterator;
    typedef typename Dcel::Inner_ccb                        Inner_ccb;
    typedef typename Dcel::Isolated_vertex                  Isolated_vertex;
    //}@
    
    typedef Arr_torus_topology_traits_2<Geometry_traits_2, Dcel> Self;

    /*! \struct
     * An auxiliary structure for rebinding the topology traits with a new 
     * geometry-traits class and a new DCEL class.
     */
    template<typename T, typename D>
    struct rebind
    {
        typedef Arr_torus_topology_traits_2<T,D> other;
    };
    
protected:
    
    typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>    Traits_adaptor_2;

    enum Identification_crossing {
        AFTER_TO_BEFORE = 1,
        BEFORE_TO_AFTER = 2
    };

    struct Point_2_less_NS {
        /*! Construct default */
        Point_2_less_NS() : 
            _m_traits(NULL) {
        }
        
        /*! Construct */    
        Point_2_less_NS(Traits_adaptor_2 * traits) : 
            _m_traits(traits) {
        }

        Traits_adaptor_2 * _m_traits;
        
        bool operator()(const Point_2& p1, const Point_2& p2) const {
            return (_m_traits->compare_x_2_object()(
                            p1, p2
                    ) == CGAL::SMALLER);
        }
    };
    
    struct Point_2_less_WE {
        
        /*! Construct default */
        Point_2_less_WE() : 
            _m_traits(NULL) {
        }
        
        /*! Construct */    
        Point_2_less_WE(Traits_adaptor_2 * traits) : 
            _m_traits(traits) {
        }
        
        Traits_adaptor_2 * _m_traits;

        bool operator()(const Point_2& p1, const Point_2& p2) const {
            return (_m_traits->compare_y_at_x_2_object()(
                            p1, p2
                    ) == CGAL::SMALLER);
        }
    };
    
    friend class Point_2_less_WE;
    friend class Point_2_less_NS;
    
    //! type of line of discontinuity
    typedef std::map< Point_2, Vertex*, Point_2_less_NS > 
    Identification_NS;
    typedef std::map< Point_2, Vertex*, Point_2_less_WE >  
    Identification_WE;
    
    // TODO check Vertex_less
    struct Vertex_less {
        bool operator() (Vertex v1, Vertex v2) {
            return &v1 < &v2;
        }
    };

    typedef std::map< Vertex, typename Identification_NS::iterator, 
                      Vertex_less >
    Vertices_on_identification_NS;
    
    typedef std::map< Vertex, typename Identification_WE::iterator, 
                      Vertex_less >
    Vertices_on_identification_WE;

    
    // Data members:
    //! the DCEL
    Dcel _m_dcel;
    
    //! the geometry-traits adapter
    Traits_adaptor_2 *_m_traits;
    
    //! Inidicate whether we should evetually free the traits object.
    bool _m_own_traits;

    //! the top face
    mutable  Face *_m_f_top;
    
    //! used to locate curve-ends on the WE-identification curve
    mutable Identification_WE _m_identification_WE;
    
    //! used to locate curve-ends on the NS-identification curve
    mutable Identification_NS _m_identification_NS;
    
    //! used to locate vertices on the WE-identification curve
    mutable Vertices_on_identification_WE _m_vertices_on_identification_WE;
    
    //! used to locate vertices on the NS-identification curve
    mutable Vertices_on_identification_NS _m_vertices_on_identification_NS;

    // Copy constructor and assignment operator - not supported.
    Arr_torus_topology_traits_2 (const Self& );

    // assign operator
    Self& operator= (const Self& );
    
public:
    
    ///! \name Construction methods.
    //@{
    
    /*! Default constructor. */
    Arr_torus_topology_traits_2 ();
    
    
    /*! Constructor with a geometry-traits class. */
    Arr_torus_topology_traits_2 (Geometry_traits_2 *tr);
    
    /*! Assign the contents of another topology-traits class. */
    void assign (const Self& other);
    
    //@}

public:
    ///! \name Accessing the DCEL and constructing iterators.
    //@{

    /*! Get the DCEL (const version). */
    const Dcel& dcel () const
    {
        return (_m_dcel);
    }
    
    /*! Get the DCEL (non-const version). */
    Dcel& dcel ()
    {
        return (_m_dcel);
    }

    /*! Determine whether the DCEL reprsenets an empty structure. */
    bool is_empty_dcel () const
    {
        // An empty arrangement contains a single bounded face
        return (this->_m_dcel.size_of_faces() == 1 && 
                this->_m_dcel.size_of_vertices() == 0 &&
                this->_m_dcel.size_of_halfedges() == 0);
    }
    
    /*! Check if the given vertex is concrete (associated with a point). */
    bool is_concrete_vertex (const Vertex *v) const
    {
        //std::cout << "Arr_torus_topological_traits_2::is_concrete_vertex"
        //          << std::endl;
        return (! v->has_null_point());
    }
    
    /*! Get the number of concrete vertices. */
    Size number_of_concrete_vertices () const
    {
        //std::cout << "Arr_torus_topological_traits_2::" 
        //          << "number_of_concrete_vertices"
        //          << std::endl;
        // All vertices not lying at infinity are concrete.
        return (this->_m_dcel.size_of_vertices());
    }
    
    /*! Check if the given vertex is valid (not a fictitious one). */
    bool is_valid_vertex (const Vertex *v) const
    {
        //std::cout << "Arr_torus_topological_traits_2::is_valid_vertex"
        //          << std::endl;
        // all vertices are valid - even v_left/right lying at inf
        return (true);
    }
    
    /*! Get the number of valid vertices. */
    Size number_of_valid_vertices () const
    {
        //std::cout << "Arr_torus_topological_traits_2::number_of_valid_vertices"
        //          << std::endl;
        // all vertices are valid - even v_left/right lying at inf
        return (this->_m_dcel.size_of_vertices());
    }
    
    /*! Check if the given halfedge is valid (not a fictitious one). */
    bool is_valid_halfedge (const Halfedge *he) const
    {
        //std::cout << "Arr_torus_topological_traits_2::is_valid_halfedge"
        //          << std::endl;
        // all halfedges are valid
        return (true);
    }
    
    /*! Get the number of valid halfedges. */
    Size number_of_valid_halfedges () const
    {
        //std::cout << "Arr_torus_topological_traits_2::" 
        //          << "number_of_valid_halfedges"
        //          << std::endl;
        // all halfedges are valid
        return (this->_m_dcel.size_of_halfedges());
  }

    /*! Check if the given face is valid. */
    bool is_valid_face (const Face *f) const
    {
        //std::cout << "Arr_torus_topological_traits_2::is_valid_face"
        //          << std::endl;
        // all faces are valid
        return (true);
    }
    
    /*! Get the number of valid faces. */
    Size number_of_valid_faces () const
    {
        //std::cout << "Arr_torus_topological_traits_2::number_of_valid_faces"
        //          << std::endl;
        // all faces are valid
        return (this->_m_dcel.size_of_faces());
    }
    //@}

private:
    
    /// \name Auxiliary type definitions.
    //@{
    typedef Arrangement_on_surface_2<Geometry_traits_2, Self>    Arr;
    
    // Type definition for the constuction sweep-line visitor.
    typedef Arr_construction_subcurve<Geometry_traits_2>         CSubcurve; 
    typedef Arr_construction_event<Geometry_traits_2,
                                 CSubcurve,
                                 Arr>                          CEvent;
    typedef Arr_torus_construction_helper<Geometry_traits_2,
                                             Arr,
                                             CEvent,
                                             CSubcurve>        CHelper;
    
    // Type definition for the basic insertion sweep-line visitor.
    typedef Arr_basic_insertion_traits_2<Geometry_traits_2, Arr> BInsTraits;
    typedef Arr_construction_subcurve<BInsTraits>                BISubcurve; 
    typedef Arr_construction_event<BInsTraits,
                                 BISubcurve,
                                 Arr>                          BIEvent;
    typedef Arr_torus_insertion_helper<BInsTraits,
                                          Arr,
                                          BIEvent,
                                          BISubcurve>          BIHelper;
    
    // Type definition for the insertion sweep-line visitor.
    typedef Arr_insertion_traits_2<Geometry_traits_2, Arr>       InsTraits;
    typedef Arr_construction_subcurve<InsTraits>                 ISubcurve; 
    typedef Arr_construction_event<InsTraits,
                                 ISubcurve,
                                 Arr>                          IEvent;
    typedef Arr_torus_insertion_helper<InsTraits,
                                          Arr,
                                          IEvent,
                                          ISubcurve>           IHelper;
    
    // Type definition for the batched point-location sweep-line visitor.
    typedef Arr_batched_point_location_traits_2<Arr>             BplTraits;
    typedef Arr_torus_batched_pl_helper<BplTraits, Arr>     BplHelper;
    
    // Type definition for the overlay sweep-line visitor.
    template <class ExGeomTraits_, class ArrangementA_, class ArrangementB_>
    struct _Overlay_helper : public Arr_torus_overlay_helper
    <ExGeomTraits_, ArrangementA_, ArrangementB_, Arr,
       Arr_construction_event<ExGeomTraits_,
                              Arr_overlay_subcurve<ExGeomTraits_>,
                              Arr>,
       Arr_overlay_subcurve<ExGeomTraits_> >
    {
        typedef Arr_torus_overlay_helper
        <ExGeomTraits_, ArrangementA_, ArrangementB_, Arr,
               Arr_construction_event<ExGeomTraits_,
                                      Arr_overlay_subcurve<ExGeomTraits_>,
                                      Arr>,
               Arr_overlay_subcurve<ExGeomTraits_> >     Base;
        
        typedef typename Base::Traits_2                    Traits_2;
        typedef typename Base::Arrangement_red_2           Arrangement_red_2;
        typedef typename Base::Arrangement_blue_2          Arrangement_blue_2;
        typedef typename Base::Arrangement_2               Arrangement_2;
        typedef typename Base::Event                       Event;
        typedef typename Base::Subcurve                    Subcurve;
        typedef typename Base::Construction_helper         Construction_helper;

        _Overlay_helper (const ArrangementA_ *arrA,
                         const ArrangementB_ *arrB) :
            Base (arrA, arrB)
        {}
    };
    //@}
    
public:
    
    ///! \name Visitor types.
    //@{
    
    typedef Arr_construction_sl_visitor<CHelper>
    Sweep_line_construction_visitor;
    
    typedef Arr_insertion_sl_visitor<IHelper>
    Sweep_line_insertion_visitor;
    
    typedef Sweep_line_construction_visitor
    Sweep_line_non_intersecting_construction_visitor;
    
    typedef Arr_basic_insertion_sl_visitor<BIHelper>
    Sweep_line_non_intersecting_insertion_visitor;
    
    template <class OutputIterator_>
    struct Sweep_line_bacthed_point_location_visitor :
        public Arr_batched_pl_sl_visitor<BplHelper, OutputIterator_>
    {
        typedef OutputIterator_                               Output_iterator;
        typedef Arr_batched_pl_sl_visitor<BplHelper,
                                      Output_iterator>        Base;
        
        typedef typename Base::Traits_2                       Traits_2;
        typedef typename Base::Event                          Event;
        typedef typename Base::Subcurve                       Subcurve;
        
        Sweep_line_bacthed_point_location_visitor (const Arr *arr,
                                                   Output_iterator *oi) :
            Base (arr, oi)
        {}
    };
    
    template <class ArrangementA_, class ArrangementB_, class OverlayTraits_>
    struct Sweep_line_overlay_visitor :
        public Arr_overlay_sl_visitor
    <_Overlay_helper<Arr_overlay_traits_2<Geometry_traits_2,
                                              ArrangementA_,
                                              ArrangementB_>,
                         ArrangementA_, 
                         ArrangementB_>,
         OverlayTraits_>
    {
        typedef ArrangementA_                            ArrangementA_2;
        typedef ArrangementB_                            ArrangementB_2;
        typedef Arr                                      Arrangement_result_2;
        typedef OverlayTraits_                           Overlay_traits;
        
        typedef Arr_overlay_traits_2<Geometry_traits_2,
                                 ArrangementA_2,
                                 ArrangementB_2>     Geom_ovl_traits_2;
        
        typedef _Overlay_helper<Geom_ovl_traits_2,
                            ArrangementA_2,
                            ArrangementB_2>          Ovl_helper;
        
        typedef Arr_overlay_sl_visitor<Ovl_helper,
                                   Overlay_traits>   Base;
        
        typedef typename Base::Traits_2                  Traits_2;
        typedef typename Base::Event                     Event;
        typedef typename Base::Subcurve                  Subcurve;
        
        Sweep_line_overlay_visitor (const ArrangementA_2 *arrA,
                                    const ArrangementB_2 *arrB,
                                    Arrangement_result_2 *arr_res,
                                    Overlay_traits *overlay_tr) :
            Base (arrA, arrB, arr_res, overlay_tr)
        {}
    };
    
    typedef Arr_planar_inc_insertion_zone_visitor<Arr>
    Zone_insertion_visitor;
    

    //! the point location strategy
    typedef Arr_naive_point_location<Arr> Default_point_location_strategy;

    //@}
    
    ///! \name Topology-traits methods.
    //@{
    
    /*!
     * Initialize an empty DCEL structure.
     */
    void init_dcel ();
    
    /*!
     * Compare the relative y-position of the given point and the given edge
     * \param p The point.
     * \param he The edge (one of the pair of halfedges).
     * \pre p should lie in the x-range of the given edge.
     * \return The relative y-position of the point p and the edge.
     */
    Comparison_result compare_y_at_x (const Point_2& p,
                                      const Halfedge* he) const;
    

    /*!
     * Check if the given vertex is associated with the given curve end.
     * \param v The vertex.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param bound_x The boundary condition of the curve end in x.
     * \param bound_y The boundary condition of the curve end in y.
     * \pre The curve has a boundary condition in either x or y.
     * \return Whether v represents the given curve end.
     */
    bool are_equal (const Vertex *v,
                    const X_monotone_curve_2& cv, Curve_end ind,
                    Boundary_type bound_x, Boundary_type bound_y) const;

      /*!
     * Given a curve end with boundary conditions and a face that contains the
     * interior of the curve, find a place for a boundary vertex that will
     * represent the curve end along the face boundary.
     * \param f The face.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param bound_x The boundary condition of the curve end in x.
     * \param bound_y The boundary condition of the curve end in y.
     * \pre The curve has a boundary condition in either x or y.
     * \return An object that contains the curve end.
     *         In our case this object is either empty, or it may wrap a
     *         vertex with boundary conditions.
     */
    CGAL::Object place_boundary_vertex (Face *f,
                                        const X_monotone_curve_2& cv,
                                        Curve_end ind,
                                        Boundary_type bound_x,
                                        Boundary_type bound_y);

     /*!
      * Locate the predecessor halfedge for the given curve around a given
      * vertex with boundary conditions.
      * \param v The vertex.
      * \param cv The x-monotone curve.
      * \param ind The curve end.
      * \param bound_x The boundary condition of the curve end in x.
      * \param bound_y The boundary condition of the curve end in y.
      * \pre The curve has a boundary condition in either x or y, and should be
      *      incident to the vertex v.
      * \return An object that contains the curve end.
      */
    Halfedge* locate_around_boundary_vertex (Vertex *v,
                                             const X_monotone_curve_2& cv,
                                             Curve_end ind,
                                             Boundary_type bound_x,
                                             Boundary_type bound_y) const;
    
    /*!
     * Receive a notification on the creation of a new boundary vertex that
     * corresponds to the given curve end.
     * \param v The new boundary vertex.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param bound_x The boundary condition of the curve end in x.
     * \param bound_y The boundary condition of the curve end in y.
     */
    void notify_on_boundary_vertex_creation (Vertex *v,
                                             const X_monotone_curve_2& cv,
                                             Curve_end ind,
                                             Boundary_type bound_x,
                                             Boundary_type bound_y) const;
    
    /*!
     * Locate a DCEL feature that contains the given unbounded curve end.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param bound_x The boundary condition of the curve end in x.
     * \param bound_y The boundary condition of the curve end in y.
     * \pre The curve end is unbounded in either x or y.
     * \return An object that contains the curve end.
     *         In our case this object may either wrap an unbounded face,
     *         or an edge with an end-vertex at 
     *         infinity (in case of an overlap).
     */
    CGAL::Object locate_unbounded_curve_end (const X_monotone_curve_2& cv,
                                             Curve_end ind,
                                             Boundary_type bound_x,
                                             Boundary_type bound_y) {
        CGAL_assertion (false);
        return CGAL::Object();
    }
    
    /*!
     * Given two halfedges, determine if the path from the source vertex of the
     * first halfedge to the source vertex of the second haldedge (i.e. we go
     * from the first halfedge until reaching the second) is perimetric.
     * \param e1 The first halfedge.
     * \param e2 The second halfedge.
     * \return Whether the path from e1 to e2 (not inclusive) is perimetric.
     */
    bool is_perimetric_path (const Halfedge *e1,
                             const Halfedge *e2) const;
    
    /*!
     * Given two predecessor halfedges that will be used for inserting a
     * new halfedge pair (prev1 will be the predecessor of the halfedge he1,
     * and prev2 will be the predecessor of its twin he2), such that the
     * insertion will create a new perimetric face that forms a hole inside
     * an existing perimetric face, determine whether he1 will be incident to
     * this new face.
     * \param prev1 The first predecessor halfedge.
     * \param prev2 The second predecessor halfedge.
     * \param cv The x-monotone curve we use to connect prev1's target and
     *           prev2's target vertex.
     * \pre prev1 and prev2 belong to the same inner connected component.
     * \return true if he1 (and prev1) lies in the interior of the face we
     *         are about to create, false otherwise - in which case he2
     *         (and prev2) must be incident to this new face.
     */
    bool is_on_new_perimetric_face_boundary (const Halfedge *prev1,
                                             const Halfedge *prev2,
                                             const X_monotone_curve_2& cv) 
        const;

    /*!
     * Determine whether the two halfedges, belonging to different outer CCBs,
     * belong to the outer boundary of the same face.
     * \param e1 The first halfedge.
     * \param e2 The second halfedge.
     * \return Whether the two halfedge belong to the outer boundary of the 
     *         same face.
     */
    bool boundaries_of_same_face (const Halfedge *e1,
                                  const Halfedge *e2) const;
    
    
    /*!
     * Determine whether the given point lies in the interior of the given 
     * face.
     * \param f The face.
     * \param p The query point.
     * \param v The vertex associated with p (if exists).
     * \param f must not be fictitious, and v must not lie at infinity.
     * \return Whether p is contained in f's interior.
     */
    bool is_in_face (const Face *f, const Point_2& p, const Vertex *v) const;
    
    
    /*!
     * Split a fictitious edge using the given vertex.
     * \param e The edge to split (one of the pair of halfedges).
     * \param v The split vertex.
     * \pre e is a fictitious halfedge.
     * \return A halfedge whose direction is the 
     *         same as e's and whose target is
     *         the split vertex v.
     */
    Halfedge* split_fictitious_edge (Halfedge *e, Vertex *v) {
        // this topology never introduces fictious halfedges
        //std::cout << "Arr_torus_topology_traits_2::split_fictious_edge" 
        //          << std::endl;
        CGAL_assertion (false);
        return (0);
    }
    
    /*!
     * Determine whether the given face is unbounded.
     * \param f The face.
     * \return Whether f is unbounded.
     */
    bool is_unbounded (const Face *f) const {
        return false;
    }
    
    /*!
     * Determine whether the given boundary vertex is redundant.
     * \param v The vertex.
     * \return Whether v is redundant, and should be erased.
     */
    bool is_redundant (const Vertex *v) const;
        
    
    /*!
     * Erase the given redundant vertex by merging a fictitious edge.
     * The function does not free the vertex v itself.
     * \param v The vertex.
     * \pre v is a redundant vertex.
     * \return One of the pair of halfedges that form the merged edge.
     */
    Halfedge* erase_redundant_vertex (Vertex *v);
    
    //@}
    
    /// \name Additional accessors, specialized for this topology-traits class.
    //@{

    /*! Get top face (const version). */
    const Face* top_face () const
    {
        return (_m_f_top);
    }
    
    /*! Get top face */
    Face* top_face ()
    {
        return (_m_f_top);
    }

    /*! Get the vertex on line of identification associated with \c pt*/
    Vertex* vertex_NS(const Point_2& key) {
        typename Identification_NS::iterator it = 
            this->_m_identification_NS.find(key);
        if (it != this->_m_identification_NS.end()) {
            return it->second;
        }
        // else
        return NULL;
    }
    
    /*! Get the vertex on line of identification associated with \c pt*/
    Vertex* vertex_WE(const Point_2& key) {
        typename Identification_WE::iterator it = 
            this->_m_identification_WE.find(key);
        if (it != this->_m_identification_WE.end()) {
            return it->second;
        }
        // else
        return NULL;
    }
    
    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Identification_NS::iterator 
    curve_ends_and_vertices_on_identification_NS_begin() {
        return _m_identification_NS.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Identification_NS::iterator 
    curve_ends_and_vertices_on_identification_NS_end() {
        return _m_identification_NS.end();
    }

    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Identification_NS::const_iterator 
    curve_ends_and_vertices_on_identification_NS_begin() const {
        return _m_identification_NS.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Identification_NS::const_iterator 
    curve_ends_and_vertices_on_identification_NS_end() const {
        return _m_identification_NS.end();
    }

    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Identification_WE::iterator 
    curve_ends_and_vertices_on_identification_WE_begin() {
        return _m_identification_WE.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Identification_WE::iterator 
    curve_ends_and_vertices_on_identification_WE_end() {
        return _m_identification_WE.end();
    }

    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Identification_WE::const_iterator 
    curve_ends_and_vertices_on_identification_WE_begin() const {
        return _m_identification_WE.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Identification_WE::const_iterator 
    curve_ends_and_vertices_on_identification_WE_end() const {
        return _m_identification_WE.end();
    }


    /*! Get the geometry traits */
    Geometry_traits_2* geometry_traits ()
    {
        return (_m_traits);
    }
    
    //@}
    
protected:
    
    /// \name Auxiliary functions.
    //@{
    
    /*!
     * Computes the number of crossing of a path with the line of discontinuity
     * \param he1 Beginning of path
     * \param he2 End of path
     * \param leftmost OUTPUT: returns whether the leftmost intersection
     *                         is AFTER_TO_BEFORE or BEFORE_TO_AFTER
     * \param loop Indicates whether we should consider he1 as defining a 
     *        ccb-loop
     * \return a pair of values, crossings from AFTER_TO_BEFORE and 
               BEFORE_TO_AFTER
     */
    std::pair< unsigned int, unsigned int >
    _crossings_with_identification_NS(
            const Halfedge* he1, const Halfedge* he2,
            Identification_crossing& leftmost) const;

    /*!\brief
     * checks whether boundary condition in x and y is valid
     */
    inline 
    bool  _valid(CGAL::Boundary_type bound_x, CGAL::Boundary_type bound_y) 
        const {
        return (                        
                ((bound_x == CGAL::AFTER_DISCONTINUITY || 
                  bound_x == CGAL::BEFORE_DISCONTINUITY) &&
                 bound_y == CGAL::NO_BOUNDARY) 
                ||
                ((bound_y == CGAL::AFTER_DISCONTINUITY || 
                  bound_y == CGAL::BEFORE_DISCONTINUITY) &&
                 bound_x == CGAL::NO_BOUNDARY)
        );
    }       
    
    //@}
};

CGAL_END_NAMESPACE

#include <CGAL/Arr_topology_traits/Arr_torus_topology_traits_2_impl.h>

#endif
