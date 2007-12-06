// Copyright (c) 2006 Max-Planck-Institute Saarbruecken (Germany).
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

#ifndef CGAL_ARR_QDX_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_QDX_TOPOLOGY_TRAITS_2_H

/*! \file
 * Definition of the Arr_qdx_topology_traits_2<GeomTraits> class.
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
#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_overlay_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_overlay_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_overlay_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_batched_pl_sl_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>

#include <CGAL/Arr_topology_traits/Arr_qdx_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_qdx_insertion_helper.h>
#include <CGAL/Arr_topology_traits/Arr_qdx_overlay_helper.h>
#include <CGAL/Arr_topology_traits/Arr_qdx_batched_pl_helper.h>
#include <CGAL/Arr_topology_traits/Arr_inc_insertion_zone_visitor.h>

#if QdX_USE_AcX
#include <SoX/GAPS/Restricted_cad_3.h>
#endif

CGAL_BEGIN_NAMESPACE

// Forward declaration:
template <class GeomTraits_, class TopTraits_> 
class Arrangement_on_surface_2;

/*! \class Arr_qdx_topology_traits_2
 * A topology-traits class that encapsulates the embedding of 2D arrangements
 * of unbounded curves on the plane.
 */
template <class GeomTraits_,
          class Dcel_ = Arr_default_dcel<GeomTraits_> >
class Arr_qdx_topology_traits_2
{
    
public:
    ///! \name The geometry-traits types.
    //@{
    typedef GeomTraits_                                     Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2             Point_2;
    typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
    //@}
    
    ///! \name Embedded typed
    //@{
    typedef typename Geometry_traits_2::Surface_3           Quadric_3;

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
    
    typedef Arr_qdx_topology_traits_2<Geometry_traits_2,
                                           Dcel>          Self;

  /*! \struct
   * An auxiliary structure for rebinding the topology traits with a new 
   * geometry-traits class and a new DCEL class.
   */
  template<typename T, typename D>
  struct rebind
  {
    typedef Arr_qdx_topology_traits_2<T,D> other;
  };
    
protected:
    
    typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>    Traits_adaptor_2;

    enum Discontinuity_crossing {
        AFTER_TO_BEFORE = 1,
        BEFORE_TO_AFTER = 2
    };
    
    // less functor
    struct Point_2_less {
        Point_2_less() {
        }

        // constructor for traits
        Point_2_less(Traits_adaptor_2 * traits) : m_traits(traits) {
        }
        Traits_adaptor_2 * m_traits;

        bool operator() (const Point_2& p1, const Point_2& p2) {
            return 
                (m_traits->compare_x_on_identification_2_object()(p1, p2) == 
                 CGAL::SMALLER);
        }
    };

    struct Vertex_less {
        bool operator() (Vertex *v1, Vertex *v2) {
            return &(*v1) < &(*v2);
        }
    };

    friend class Point_2_less;

    //! type of line of discontinuity
    typedef std::map< Point_2, Vertex*, Point_2_less > Line_of_discontinuity;

    typedef std::map< Vertex*, typename Line_of_discontinuity::iterator,
                      Vertex_less >
    Vertices_on_line_of_discontinuity;
    
    // Data members:
    //! the DCEL
    Dcel m_dcel;
    
    //! the geometry-traits adapter
    Traits_adaptor_2 *m_traits;
    
    //! Inidicate whether we should evetually free the traits object.
    bool m_own_traits;

    mutable Quadric_3 m_quadric;
    
    //! indicates kind of left boundary
    mutable CGAL::Arr_parameter_space m_left;
    
    //! indicates kind of right boundary
    mutable CGAL::Arr_parameter_space m_right;
    
    // if non-concrete then inf, if concrete then singularity
    //! a vertex representing the left singularity/infinity of a quadric
    mutable Vertex *v_left; 

    //! a vertex representing the right singularity/infinity of a quadric
    mutable Vertex *v_right;

    //! the top face
    mutable  Face *f_top;

    //! used to locate curve-ends on the line of discontinuity
    mutable Line_of_discontinuity m_line_of_discontinuity;

    //! used to locate vertices on the line of discontinuity
    mutable Vertices_on_line_of_discontinuity 
    m_vertices_on_line_of_discontinuity;

    // Copy constructor and assignment operator - not supported.
    Arr_qdx_topology_traits_2 (const Self& );

    // assign operator
    Self& operator= (const Self& );
    
public:
    
    ///! \name Construction methods.
    //@{
    
    /*! Default constructor. */
    Arr_qdx_topology_traits_2 ();
    
    
    /*! Constructor with a geometry-traits class. */
    Arr_qdx_topology_traits_2 (Geometry_traits_2 *tr);
    
    /*! Assign the contents of another topology-traits class. */
    void assign (const Self& other);
    
    //@}

private:
    //!\name Initialising
    //!@{
    
    /*! initiales toptraits with a specific quadric \c base
     *
     * \pre \c base is elliposid, elliptic cylinder or elliptic paraboloid
     *      whose directrix is not parallel to the y-axis nor z-axis
     */
    void _initialize_with_quadric(const Quadric_3& base) {
        CGAL_precondition(base.is_ellipsoid() || base.is_elliptic_cylinder() ||
                          base.is_elliptic_paraboloid());
#if !QdX_USE_AcX        
        std::vector < typename Quadric_3::P_curve_2 > sils;
        base.silhouette(std::back_inserter(sils));
        int csil = static_cast< int >(sils.size());
        
        // check position of quadric
        if (base.is_ellipsoid()) {
            // extreme points of projected silhouette match singular points
            CGAL_precondition(csil == 1);
            m_left = CGAL::AFTER_SINGULARITY;
            m_right = CGAL::BEFORE_SINGULARITY;
        }
        if (base.is_elliptic_cylinder()) {
            // ensure to be a "non-vertical" cylinder:
            // test whether projected silhouette has two parallel lines
            // going to infinity
            
            if (csil == 1) {
                CGAL_precondition(sils[0].arcs_over_interval(0) == 2);
            } else {
                CGAL_precondition(csil == 2);
                CGAL_precondition(sils[0].arcs_over_interval(0) == 1 &&
                                  sils[1].arcs_over_interval(0) == 1);
            }
            m_left = CGAL::MINUS_INFINITY;
            m_right = CGAL::PLUS_INFINITY;
        }
        if (base.is_elliptic_paraboloid()) {
            CGAL_precondition(csil == 1);
            // check whether not "vertical"
            // TASK otherwise we have singularity in y
            CGAL_precondition(sils[0].f().degree() == 2);
            // ensure to be left- or right oriented
            // and extreme of projected silhouette matches to singularity
            CGAL_precondition(sils[0].arcs_over_interval(0) != 1); 
            // == 2 || == 0
            if (sils[0].arcs_over_interval(0) == 2) {
                m_left = CGAL::MINUS_INFINITY;
                m_right = CGAL::BEFORE_SINGULARITY;
            } else {
                m_left = CGAL::AFTER_SINGULARITY;
                m_right = CGAL::PLUS_INFINITY;
            }
        }
#else
        typedef QdX::Quadric_3_z_at_xy_isolator_traits< Quadric_3 > Traits;
        typedef SoX::Create_restricted_cad_3< Traits > Creator;
        typename Creator::Restricted_cad_3 cad = Creator()(base);
        
        if (base.is_ellipsoid()) {
            // extreme points of projected silhouette match singular points
            
            CGAL_precondition(cad.number_of_vertices() == 2);
            CGAL_precondition(cad.number_of_edges() == 2);
            CGAL_precondition(cad.number_of_faces() == 2);
            CGAL_precondition(cad.number_of_unbounded_faces() == 1);
            
            m_left = CGAL::AFTER_SINGULARITY;
            m_right = CGAL::BEFORE_SINGULARITY;
        }
        if (base.is_elliptic_cylinder()) {
            // ensure to be a "non-vertical" cylinder:
            // test whether projected silhouette has two parallel lines
            // going to infinity
            
            CGAL_precondition(cad.number_of_vertices() == 0);
            CGAL_precondition(cad.number_of_edges() == 2);
            CGAL_precondition(cad.number_of_faces() == 3);
            CGAL_precondition(cad.number_of_unbounded_faces() == 3);
            
            m_left = CGAL::MINUS_INFINITY;
            m_right = CGAL::PLUS_INFINITY;
        }
        if (base.is_elliptic_paraboloid()) {
            
            CGAL_precondition(cad.number_of_vertices() == 1);
            CGAL_precondition(cad.number_of_edges() == 2);
            CGAL_precondition(cad.number_of_faces() == 2);
            CGAL_precondition(cad.number_of_unbounded_faces() == 2);
            
            int number_of_vertices_at_minus_inf = 0;
            int number_of_vertices_at_plus_inf = 0;
            
            for (typename Creator::Restricted_cad_3::Edge_const_iterator
                     eit = cad.edges_begin(); eit != cad.edges_end();
                 eit++) {
#if 1 // TODO use traits instead
                CGAL_precondition(
                        eit->curve().get_parameter_space_in_y(CGAL::ARR_MIN_END)
                        == CGAL::ARR_INTERIOR
                ); 
                if (eit->curve().get_parameter_space_in_x(CGAL::ARR_MIN_END)
                    == CGAL::MINUS_INFINITY) {
                    number_of_vertices_at_minus_inf++;
                }
                CGAL_precondition(
                        eit->curve().get_parameter_space_in_y(CGAL::ARR_MAX_END)
                        == CGAL::ARR_INTERIOR
                ); 
                if (eit->curve().get_parameter_space_in_x(CGAL::ARR_MAX_END)
                    == CGAL::PLUS_INFINITY) {
                    number_of_vertices_at_plus_inf++;
                }                    
#else
                CGAL_precondition(
                        m_traits->parameter_space_in_y_2_object()(
                                eit->curve(), CGAL::ARR_MIN_END
                        ) == CGAL::ARR_INTERIOR
                ); 
                if (m_traits->parameter_space_in_x_2_object()(
                            eit->curve(), CGAL::ARR_MIN_END
                    ) == CGAL::MINUS_INFINITY) {
                    number_of_vertices_at_minus_inf++;
                }
                CGAL_precondition(
                        m_traits->parameter_space_in_y_2_object()(
                                eit->curve(), CGAL::ARR_MAX_END
                        ) == CGAL::ARR_INTERIOR
                ); 
                if (m_traits->parameter_space_in_x_2_object()(
                            eit->curve(), CGAL::ARR_MAX_END
                    ) == CGAL::PLUS_INFINITY) {
                    number_of_vertices_at_plus_inf++;
                }
#endif
            }

            CGAL_assertion(number_of_vertices_at_minus_inf + 
                           number_of_vertices_at_plus_inf == 4);
            CGAL_assertion(number_of_vertices_at_minus_inf != 1);
            CGAL_assertion(number_of_vertices_at_plus_inf != 1);
            
            // == 2 || == 0
            if (number_of_vertices_at_minus_inf > 0) {
                m_left = CGAL::MINUS_INFINITY;
                m_right = CGAL::BEFORE_SINGULARITY;
            } else {
                m_left = CGAL::AFTER_SINGULARITY;
                m_right = CGAL::PLUS_INFINITY;
            }
        }
#endif
        m_quadric = base;
        this->init_dcel();
    }
    
    //!@}
    
public:
    ///! \name Accessing the DCEL and constructing iterators.
    //@{

    /*! Get the DCEL (const version). */
    const Dcel& dcel () const
    {
        return (m_dcel);
    }
    
    /*! Get the DCEL (non-const version). */
    Dcel& dcel ()
    {
        return (m_dcel);
    }
    
    /*! Determine whether the DCEL reprsenets an empty structure. */
    bool is_empty_dcel () const
    {
        // An empty arrangement contains a single face, unbounded or bounded
        return (this->m_dcel.size_of_faces() == 1 && 
                this->m_dcel.size_of_vertices() == 0 &&
                this->m_dcel.size_of_halfedges() == 0);
    }
    
    /*! Check if the given vertex is concrete (associated with a point). */
    bool is_concrete_vertex (const Vertex *v) const
    {
        //std::cout << "Arr_qdx_topological_traits_2::is_concrete_vertex"
        //          << std::endl;
        return (! v->has_null_point());
    }
    
    /*! Get the number of concrete vertices. */
    Size number_of_concrete_vertices () const
    {
        //std::cout << "Arr_qdx_topological_traits_2::" 
        //          << "number_of_concrete_vertices"
        //          << std::endl;
        // All vertices not lying at infinity are concrete.
        return (this->m_dcel.size_of_vertices() - 
                (v_left != 0 && 
                 v_left->parameter_space_in_x() == CGAL::MINUS_INFINITY ?
                 1 : 0) -
                (v_right != 0 && 
                 v_right->parameter_space_in_x() == CGAL::PLUS_INFINITY ?
                 1 : 0)
        );
    }
    
    /*! Check if the given vertex is valid (not a fictitious one). */
    bool is_valid_vertex (const Vertex *v) const
    {
        //std::cout << "Arr_qdx_topological_traits_2::is_valid_vertex"
        //          << std::endl;
        // all vertices are valid - even v_left/right lying at inf
        return (true);
    }
    
    /*! Get the number of valid vertices. */
    Size number_of_valid_vertices () const
    {
        //std::cout << "Arr_qdx_topological_traits_2::number_of_valid_vertices"
        //          << std::endl;
        // all vertices are valid - even v_left/right lying at inf
        return (this->m_dcel.size_of_vertices());
    }
    
    /*! Check if the given halfedge is valid (not a fictitious one). */
    bool is_valid_halfedge (const Halfedge *he) const
    {
        //std::cout << "Arr_qdx_topological_traits_2::is_valid_halfedge"
        //          << std::endl;
        // all halfedges are valid
        return (true);
    }
    
    /*! Get the number of valid halfedges. */
    Size number_of_valid_halfedges () const
    {
        //std::cout << "Arr_qdx_topological_traits_2::" 
        //          << "number_of_valid_halfedges"
        //          << std::endl;
        // all halfedges are valid
        return (this->m_dcel.size_of_halfedges());
  }

    /*! Check if the given face is valid. */
    bool is_valid_face (const Face *f) const
    {
        //std::cout << "Arr_qdx_topological_traits_2::is_valid_face"
        //          << std::endl;
        // all faces are valid
        return (true);
    }
    
    /*! Get the number of valid faces. */
    Size number_of_valid_faces () const
    {
        //std::cout << "Arr_qdx_topological_traits_2::number_of_valid_faces"
        //          << std::endl;
        // all faces are valid
        return (this->m_dcel.size_of_faces());
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
    typedef Arr_qdx_construction_helper<Geometry_traits_2,
                                             Arr,
                                             CEvent,
                                             CSubcurve>        CHelper;
    
    // Type definition for the basic insertion sweep-line visitor.
    typedef Arr_basic_insertion_traits_2<Geometry_traits_2, Arr> BInsTraits;
    typedef Arr_construction_subcurve<BInsTraits>                BISubcurve; 
    typedef Arr_construction_event<BInsTraits,
                                 BISubcurve,
                                 Arr>                          BIEvent;
    typedef Arr_qdx_insertion_helper<BInsTraits,
                                          Arr,
                                          BIEvent,
                                          BISubcurve>          BIHelper;
    
    // Type definition for the insertion sweep-line visitor.
    typedef Arr_insertion_traits_2<Geometry_traits_2, Arr>       InsTraits;
    typedef Arr_construction_subcurve<InsTraits>                 ISubcurve; 
    typedef Arr_construction_event<InsTraits,
                                 ISubcurve,
                                 Arr>                          IEvent;
    typedef Arr_qdx_insertion_helper<InsTraits,
                                          Arr,
                                          IEvent,
                                          ISubcurve>           IHelper;
    
    // Type definition for the batched point-location sweep-line visitor.
    typedef Arr_batched_point_location_traits_2<Arr>             BplTraits;
    typedef Arr_qdx_batched_pl_helper<BplTraits, Arr>     BplHelper;
    
    // Type definition for the overlay sweep-line visitor.
    template <class ExGeomTraits_, class ArrangementA_, class ArrangementB_>
    struct _Overlay_helper : public Arr_qdx_overlay_helper
    <ExGeomTraits_, ArrangementA_, ArrangementB_, Arr,
       Arr_construction_event<ExGeomTraits_,
                              Arr_overlay_subcurve<ExGeomTraits_>,
                              Arr>,
       Arr_overlay_subcurve<ExGeomTraits_> >
    {
        typedef Arr_qdx_overlay_helper
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
    
    typedef Arr_inc_insertion_zone_visitor<Arr>
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
     * Make the necessary updates after the DCEL structure have been updated.
     */
    void dcel_updated ();

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
                    const X_monotone_curve_2& cv, Arr_curve_end ind,
                    Arr_parameter_space bound_x, Arr_parameter_space bound_y) const;

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
                                        Arr_curve_end ind,
                                        Arr_parameter_space bound_x,
                                        Arr_parameter_space bound_y);

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
                                             Arr_curve_end ind,
                                             Arr_parameter_space bound_x,
                                             Arr_parameter_space bound_y) const;
    
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
                                             Arr_curve_end ind,
                                             Arr_parameter_space bound_x,
                                             Arr_parameter_space bound_y) const;
    
    /*!
     * Locate a DCEL feature that contains the given curve end.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param bound_x The boundary condition of the curve end in x.
     * \param bound_y The boundary condition of the curve end in y.
     * \pre The curve end is incident to the boundary.
     * \return An object that contains the curve end.
     *         In our case this object may either wrap an unbounded face,
     *         or an edge with an end-vertex at 
     *         infinity (in case of an overlap).
     */
    CGAL::Object locate_curve_end (const X_monotone_curve_2& cv,
                                   Arr_curve_end ind,
                                   Arr_parameter_space bound_x,
                                   Arr_parameter_space bound_y);
    
    /*!
     * Given two predecessor halfedges that belong to the same inner CCB of
     * a face, determine what happens when we insert an edge connecting the
     * target vertices of the two edges.
     * \param prev1 The first predecessor halfedge.
     * \param prev2 The second predecessor halfedge.
     * \param cv The curve to be inserted
     * \pre The two halfedges belong to the same inner CCB.
     * \return A pair indicating whether the insertion will cause the face
     *         to split (the first flag), and if so - whether the split face
     *         will form a hole in the original face.
     */
    std::pair<bool, bool>
    face_split_after_edge_insertion (const Halfedge *prev1,
                                     const Halfedge *prev2,
                                     const X_monotone_curve_2& cv) const;

    /*!
     * Determine whether the removal of the given edge will cause the creation
     * of a hole.
     *  \param he The halfedge to be removed.
     * \pre Both he and its twin lie on an outer CCB of their incident faces.
     * \return Whether a new hole will be created.
     */
    bool hole_creation_after_edge_removal (const Halfedge *he) const;

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
        //std::cout << "Arr_qdx_topology_traits_2::split_fictious_edge" 
        //          << std::endl;
        CGAL_error();
        return (0);
    }

    /*!
     * Determine whether the given face is unbounded.
     * \param f The face.
     * \return Whether f is unbounded.
     */
    bool is_unbounded (const Face *f) const;

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

    /*! Get the quadric this class has been initialized with */
    const Quadric_3& quadric() const 
    {
        return (m_quadric);
    }

    /*! Get the type of the left boundary */
    CGAL::Arr_parameter_space left_boundary () const
    {
        return (m_left);
    }

    /*! Get the type of the right boundary */
    CGAL::Arr_parameter_space right_boundary () const
    {
        return (m_right);
    }

    /*! Get the leftmost vertex (const version). */
    const Vertex* leftmost_vertex () const
    {
        return (v_left);
    }

    /*! Get the leftmost vertex (non-const version). */
    Vertex* leftmost_vertex ()
    {
        return (v_left);
    }

    /*! Get the rightmost vertex (const version). */
    const Vertex* rightmost_vertex () const
    {
        return (v_right);
    }

    /*! Get the rightmost vertex (non-const version). */
    Vertex* rightmost_vertex ()
    {
        return (v_right);
    }

    /*! Get top face (const version). */
    const Face* top_face () const
    {
        return (f_top);
    }

    /*! Get top face */
    Face* top_face ()
    {
        return (f_top);
    }

    /*! Get bottom face (const version). */
    const Face* bottom_face () const
    {
        typename Line_of_discontinuity::const_iterator  it = 
            m_line_of_discontinuity.begin();
        if (it == m_line_of_discontinuity.end()) {
            return (f_top);
        } else {
            return (_face_before_vertex_on_discontinuity (it->second));
        }
    }

    /*! Get bottom face */
    Face* bottom_face ()
    {
        typename Line_of_discontinuity::const_iterator  it = 
            m_line_of_discontinuity.begin();
        if (it == m_line_of_discontinuity.end()) {
            return (f_top);
        } else {
            return (_face_before_vertex_on_discontinuity (it->second));
        }
    }
    
    /*! Get the vertex on line of discontinuity associated with \c pt*/
    Vertex* discontinuity_vertex(const Point_2& pt) {
        typename Line_of_discontinuity::iterator it = 
            this->m_line_of_discontinuity.find(pt);
        if (it != this->m_line_of_discontinuity.end()) {
            return it->second;
        }
        // else
        return NULL;
    }
    
    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Line_of_discontinuity::iterator 
    curve_ends_and_vertices_on_line_of_discontinuity_begin() {
        return m_line_of_discontinuity.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity
     */
    typename Line_of_discontinuity::iterator 
    curve_ends_and_vertices_on_line_of_discontinuity_end() {
        return m_line_of_discontinuity.end();
    }

    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Line_of_discontinuity::const_iterator 
    curve_ends_and_vertices_on_line_of_discontinuity_begin() const {
        return m_line_of_discontinuity.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the line of discontinuity (const version)
     */
    typename Line_of_discontinuity::const_iterator 
    curve_ends_and_vertices_on_line_of_discontinuity_end() const {
        return m_line_of_discontinuity.end();
    }

    /*! Get the geometry traits */
    Geometry_traits_2* geometry_traits ()
    {
        return (m_traits);
    }
    
    //@}
    
protected:
    
    /// \name Auxiliary functions.
    //@{
    
    /*!
     * Locate given curve end at a boundary vertex to the left end
     * or right end.
     * \param v The boundary vertex 
     * \param cv The curve
     * \param ind ARR_MIN_END if the vertex is induced by the minimal end;
     *            ARR_MAX_END if it is induced by the curve's maximal end.
     * \param equal will be set to true if equal curve found
     * \param allow_equal if set to true also checks for equal curve
     * \return The predecessing halfedge of the curve-end (cv, ind) 
     *         in the circular of incident halfedges around 
     *         the boundary vertex.
     */
    Halfedge* _locate_around_vertex_with_boundary_at_x(
            Vertex* v,
            const X_monotone_curve_2& cv, Arr_curve_end ind,
            bool& equal,
            bool allow_equal) const;
    
    /*!
     * Locate given curve end at a boundary vertex on the line of
     * discontinuity.
     * \param v The boundary vertex 
     * \param cv The curve
     * \param ind ARR_MIN_END if the vertex is induced by the minimal end;
     *            ARR_MAX_END if it is induced by the curve's maximal end.
     * \return The predecessing halfedge of the curve-end (cv, ind) 
     *         in the circular of incident halfedges around 
     *         the boundary vertex.
     */
    Halfedge* _locate_around_vertex_on_discontinuity(
            Vertex* v,
            const X_monotone_curve_2 & cv,
            Arr_curve_end ind) const;
    
#if 0
    /*!
     * Get the curve associated with a boundary vertex.
     * \param v The boundary vertex.
     * \param ind Output: ARR_MIN_END if the vertex is induced by the minimal end;
     *                    ARR_MAX_END if it is induced by the curve's maximal end.
     * \pre v is a valid boundary vertex.
     * \return The curve that induces v.
     */
    const X_monotone_curve_2& _curve (const Vertex *v,
                                          Arr_curve_end& ind) const;

    /*!
     * Compares two curve-ends around a point on the line of discontinuity
     * \param cv1 First curve
     * \param ind1 End of first curve
     * \param cv2 Second curve
     * \param ind2 End of second curve
     * \return the comparison result
     */
    CGAL::Comparison_result _cw_order_at_boundary_vertex (
            const X_monotone_curve_2& cv1,
            const CGAL::Arr_curve_end& ind1,
            const X_monotone_curve_2& cv2,
            const CGAL::Arr_curve_end& ind2) const;

     /*!
     * Get the predecessing halfedge of a curve-end at a boundary vertex.
     * \param v The boundary vertex 
     * \param cv The curve
     * \param ind ARR_MIN_END if the vertex is induced by the minimal end;
     *            ARR_MAX_END if it is induced by the curve's maximal end.
     * \param overlaps OUTPUT: is true in case (cv,ind) overlaps with an 
     *                         existing halfedge which is then returned
     * \pre v is a valid boundary vertex.
     * \return The predecessing halfedge of the curve-end (cv, ind) 
     *         in the circular of incident halfedges around 
     *         the boundary vertex.
     */
    Halfedge* _find_predecessor_at_boundary_vertex(
            Vertex *v,
            const X_monotone_curve_2& cv,
            CGAL::Arr_curve_end ind,
            bool& overlaps) const;

#endif

    /*! Return the face that lies before the given vertex, which lies
     * on the line of discontinuity.
     */
    Face * _face_before_vertex_on_discontinuity (Vertex * v) const;
    
    
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
    _crossings_with_line_of_discontinuity(
            const Halfedge* he1, const Halfedge* he2,
            Discontinuity_crossing& leftmost) const;

    /*!
     * Given two halfedges, determine if the path from the source vertex of the
     * first halfedge to the source vertex of the second haldedge (i.e. we go
     * from the first halfedge until reaching the second) is perimetric.
     * \param e1 The first halfedge.
     * \param e2 The second halfedge.
     * \return Whether the path from e1 to e2 (not inclusive) is perimetric.
     */
    bool _is_perimetric_path (const Halfedge *e1,
                              const Halfedge *e2) const;

    //@}
};

CGAL_END_NAMESPACE

#include <CGAL/Arr_topology_traits/Arr_qdx_topology_traits_2_impl.h>

#endif
