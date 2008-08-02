// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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

#ifndef CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE 
#define CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE 0
#endif

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

#include <CGAL/Arr_topology_traits/Sign_of_path.h>

#include <SoX/GAPS/Restricted_cad_3.h>
#include <QdX/Quadric_3_z_at_xy_isolator_traits.h>

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

    //! type of curve of identification
    typedef std::map< Point_2, Vertex*, Point_2_less > Identification;

    typedef std::map< Vertex*, typename Identification::iterator,
                      Vertex_less >
    Vertices_on_identification;
    
    // Data members:
    //! the DCEL
    Dcel _m_dcel;
    
    //! the geometry-traits adapter
    Traits_adaptor_2 *_m_traits;
    
    //! Inidicate whether we should evetually free the traits object.
    bool _m_own_traits;

    mutable Quadric_3 _m_quadric;
    
    //! indicates kind of left boundary
    mutable CGAL::Arr_boundary_type _m_left;
    
    //! indicates kind of right boundary
    mutable CGAL::Arr_boundary_type _m_right;
    
    // if non-concrete then inf, if concrete then singularity
    //! a vertex representing the left singularity/infinity of a quadric
    mutable Vertex *v_left; 

    //! a vertex representing the right singularity/infinity of a quadric
    mutable Vertex *v_right;

    //! the top face
    mutable  Face *f_top;

    //! used to locate curve-ends on the curve of identification
    mutable Identification _m_identification;

    //! used to locate vertices on the curve of identification
    mutable Vertices_on_identification 
    _m_vertices_on_identification;

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


 ///! \name Boundaries
    //@{
    
    /*! Obtain the boundary type for a given parameter space.
     * \param ps the parameter space.
     * \return the boundary type along ps.
     * \pre ps must not be ARR_INTERIOR.
     */
    CGAL::Arr_boundary_type boundary_type(
            const CGAL::Arr_parameter_space ps
    ) const {
        CGAL_precondition(ps != CGAL::ARR_INTERIOR);
        switch (ps) {
        case ARR_LEFT_BOUNDARY:
            return _m_left;
        case ARR_RIGHT_BOUNDARY: 
            return _m_right;
            
        case ARR_BOTTOM_BOUNDARY:
        case ARR_TOP_BOUNDARY: 
            return CGAL::ARR_IDENTIFICATION;
        default: CGAL_error();
        }
        // Cannot reach here!
        return CGAL::ARR_NUMBER_OF_BOUNDARY_TYPES;
    }
    
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

        typedef QdX::Quadric_3_z_at_xy_isolator_traits< 
            typename Geometry_traits_2::Curved_kernel_via_analysis_2, 
            Quadric_3 
        > Traits;
        typedef SoX::Create_restricted_cad_3< Traits > Creator;
        typename Creator::Restricted_cad_3 cad = Creator()(base);
        
        if (base.is_ellipsoid()) {
            // extreme points of projected silhouette match singular points
            
            CGAL_precondition(cad.number_of_vertices() == 2);
            CGAL_precondition(cad.number_of_edges() == 2);
            CGAL_precondition(cad.number_of_faces() == 2);
            CGAL_precondition(cad.number_of_unbounded_faces() == 1);
            
            _m_left = CGAL::ARR_CONTRACTION;
            _m_right = CGAL::ARR_CONTRACTION;
        }
        if (base.is_elliptic_cylinder()) {
            // ensure to be a "non-vertical" cylinder:
            // test whether projected silhouette has two parallel lines
            // going to infinity
            
            CGAL_precondition(cad.number_of_vertices() == 0);
            CGAL_precondition(cad.number_of_edges() == 2);
            CGAL_precondition(cad.number_of_faces() == 3);
            CGAL_precondition(cad.number_of_unbounded_faces() == 3);
            
            _m_left = CGAL::ARR_UNBOUNDED;
            _m_right = CGAL::ARR_UNBOUNDED;
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
                CGAL_precondition(
                        eit->curve().location(CGAL::ARR_MIN_END)
                        == CGAL::ARR_INTERIOR ||
                        eit->curve().location(CGAL::ARR_MIN_END)
                        == CGAL::ARR_LEFT_BOUNDARY ||
                        eit->curve().location(CGAL::ARR_MIN_END)
                        == CGAL::ARR_RIGHT_BOUNDARY
                ); 
                if (eit->curve().location(CGAL::ARR_MIN_END)
                    == CGAL::ARR_LEFT_BOUNDARY &&
                    _m_left == CGAL::ARR_UNBOUNDED) {
                    number_of_vertices_at_minus_inf++;
                }
                CGAL_precondition(
                        eit->curve().location(CGAL::ARR_MAX_END)
                        == CGAL::ARR_INTERIOR ||
                        eit->curve().location(CGAL::ARR_MAX_END)
                        == CGAL::ARR_LEFT_BOUNDARY ||
                        eit->curve().location(CGAL::ARR_MAX_END)
                        == CGAL::ARR_RIGHT_BOUNDARY
                ); 
                if (eit->curve().location(CGAL::ARR_MAX_END)
                    == CGAL::ARR_RIGHT_BOUNDARY &&
                    _m_right == CGAL::ARR_UNBOUNDED) {
                    number_of_vertices_at_plus_inf++;
                }                    
            }

            CGAL_assertion(number_of_vertices_at_minus_inf + 
                           number_of_vertices_at_plus_inf == 4);
            CGAL_assertion(number_of_vertices_at_minus_inf != 1);
            CGAL_assertion(number_of_vertices_at_plus_inf != 1);
            
            // == 2 || == 0
            if (number_of_vertices_at_minus_inf > 0) {
                _m_left = CGAL::ARR_UNBOUNDED;
                _m_right = CGAL::ARR_CONTRACTION;
            } else {
                _m_left = CGAL::ARR_CONTRACTION;
                _m_right = CGAL::ARR_UNBOUNDED;
            }
        }
        _m_quadric = base;
        this->init_dcel();
    }
    
    //!@}
    
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
        // An empty arrangement contains a single face, unbounded or bounded
        return (this->_m_dcel.size_of_faces() == 1 && 
                this->_m_dcel.size_of_vertices() == 0 &&
                this->_m_dcel.size_of_halfedges() == 0);
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
        return (this->_m_dcel.size_of_vertices() - 
                ((v_left != NULL && _m_left == CGAL::ARR_UNBOUNDED) ? 1 : 0) -
                ((v_right != NULL && _m_right == CGAL::ARR_UNBOUNDED) ? 1 : 0)
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
        return (this->_m_dcel.size_of_vertices());
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
        return (this->_m_dcel.size_of_halfedges());
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
     * \param ps_x The boundary condition of the curve end in x.
     * \param ps_y The boundary condition of the curve end in y.
     * \pre The curve has a boundary condition in either x or y.
     * \return Whether v represents the given curve end.
     */
    bool are_equal (const Vertex *v,
                    const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
                    CGAL::Arr_parameter_space ps_x, 
                    CGAL::Arr_parameter_space ps_y) const;

      /*!
     * Given a curve end with boundary conditions and a face that contains the
     * interior of the curve, find a place for a boundary vertex that will
     * represent the curve end along the face boundary.
     * \param f The face.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param ps_x The boundary condition of the curve end in x.
     * \param ps_y The boundary condition of the curve end in y.
     * \pre The curve has a boundary condition in either x or y.
     * \return An object that contains the curve end.
     *         In our case this object is either empty, or it may wrap a
     *         vertex with boundary conditions.
     */
    CGAL::Object place_boundary_vertex (Face *f,
                                        const X_monotone_curve_2& cv,
                                        CGAL::Arr_curve_end ind,
                                        CGAL::Arr_parameter_space ps_x,
                                        CGAL::Arr_parameter_space ps_y);

     /*!
      * Locate the predecessor halfedge for the given curve around a given
      * vertex with boundary conditions.
      * \param v The vertex.
      * \param cv The x-monotone curve.
      * \param ind The curve end.
      * \param ps_x The boundary condition of the curve end in x.
      * \param ps_y The boundary condition of the curve end in y.
      * \pre The curve has a boundary condition in either x or y, and should be
      *      incident to the vertex v.
      * \return An object that contains the curve end.
      */
    Halfedge* locate_around_boundary_vertex (Vertex *v,
                                             const X_monotone_curve_2& cv,
                                             CGAL::Arr_curve_end ind,
                                             CGAL::Arr_parameter_space ps_x,
                                             CGAL::Arr_parameter_space ps_y) 
        const;
    
    /*!
     * Receive a notification on the creation of a new boundary vertex that
     * corresponds to the given curve end.
     * \param v The new boundary vertex.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param ps_x The boundary condition of the curve end in x.
     * \param ps_y The boundary condition of the curve end in y.
     */
    void notify_on_boundary_vertex_creation (Vertex *v,
                                             const X_monotone_curve_2& cv,
                                             CGAL::Arr_curve_end ind,
                                             CGAL::Arr_parameter_space ps_x,
                                             CGAL::Arr_parameter_space ps_y) 
        const;
    
    /*!
     * Locate a DCEL feature that contains the given curve end.
     * \param cv The x-monotone curve.
     * \param ind The curve end.
     * \param ps_x The boundary condition of the curve end in x.
     * \param ps_y The boundary condition of the curve end in y.
     * \pre The curve end is incident to the boundary.
     * \return An object that contains the curve end.
     *         In our case this object may either wrap an unbounded face,
     *         or an edge with an end-vertex at 
     *         infinity (in case of an overlap).
     */
    CGAL::Object locate_curve_end (const X_monotone_curve_2& cv,
                                   CGAL::Arr_curve_end ind,
                                   CGAL::Arr_parameter_space ps_x,
                                   CGAL::Arr_parameter_space ps_y);
    
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

#if CGAL_NEW_FACE_SPLIT_STRATEGY
    /*!
     * Given two predecessor halfedges that belong to the same CCB of
     * a face, determine what happens when we insert an edge connecting the
     * target vertices of the two edges.
     * \param prev1 The first predecessor halfedge.
     * \param prev2 The second predecessor halfedge.
     * \param cv The curve to be inserted
     * \pre The two halfedges belong to the same inner CCB.
     * \return A pair indicating whether the insertion will cause the face
     *         to split (the first flag), and if so - whether the prev1 will be
     *         incident to the split face (second flag).
     */
    std::pair<bool, bool>
    face_update_upon_edge_insertion (const Halfedge *prev1,
                                     const Halfedge *prev2,
                                     const X_monotone_curve_2& cv) const;
#endif
    
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
        return (_m_quadric);
    }

    /*! Get the type of the left boundary */
    CGAL::Arr_parameter_space left_boundary () const
    {
        return (_m_left);
    }

    /*! Get the type of the right boundary */
    CGAL::Arr_parameter_space right_boundary () const
    {
        return (_m_right);
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
        typename Identification::const_iterator  it = 
            _m_identification.begin();
        if (it == _m_identification.end()) {
            return (f_top);
        } else {
            return (_face_before_vertex_on_identification(it->second));
        }
    }

    /*! Get bottom face */
    Face* bottom_face ()
    {
        typename Identification::const_iterator  it = 
            _m_identification.begin();
        if (it == _m_identification.end()) {
            return (f_top);
        } else {
            return (_face_before_vertex_on_identification(it->second));
        }
    }
    
    /*! Get the vertex on curve of identification associated with \c pt*/
    Vertex* vertex_on_identification(const Point_2& pt) {
        typename Identification::iterator it = 
            this->_m_identification.find(pt);
        if (it != this->_m_identification.end()) {
            return it->second;
        }
        // else
        return NULL;
    }
    
    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the curve of identification
     */
    typename Identification::iterator 
    curve_ends_and_vertices_on_identification_begin() {
        return _m_identification.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the curve of identification
     */
    typename Identification::iterator 
    curve_ends_and_vertices_on_identification_end() {
        return _m_identification.end();
    }

    /*! Get the beginning of all pairs of curve-end and its vertices
     *  along the curve of identification (const version)
     */
    typename Identification::const_iterator 
    curve_ends_and_vertices_on_identification_begin() const {
        return _m_identification.begin();
    }

    /*! Get the past-the-end value of all pairs of curve-end and its vertices
     *  along the curve of identification (const version)
     */
    typename Identification::const_iterator 
    curve_ends_and_vertices_on_identification_end() const {
        return _m_identification.end();
    }

    /*! Get the geometry traits */
    Geometry_traits_2* geometry_traits() { 
        return (_m_traits);
    }

    /*! Get the geometry traits (const version) */
    Geometry_traits_2* geometry_traits() const {
        return (_m_traits);
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
            const X_monotone_curve_2& cv, CGAL::Arr_curve_end ind,
            bool& equal,
            bool allow_equal) const;
    
    /*!
     * Locate given curve end at a boundary vertex on the curve of
     * identification
     * \param v The boundary vertex 
     * \param cv The curve
     * \param ind ARR_MIN_END if the vertex is induced by the minimal end;
     *            ARR_MAX_END if it is induced by the curve's maximal end.
     * \return The predecessing halfedge of the curve-end (cv, ind) 
     *         in the circular of incident halfedges around 
     *         the boundary vertex.
     */
    Halfedge* _locate_around_vertex_on_identification(
            Vertex* v,
            const X_monotone_curve_2 & cv,
            CGAL::Arr_curve_end ind) const;
    
    /*! Return the face that lies before the given vertex, which lies
     * on the curve of identification
     */
    Face * _face_before_vertex_on_identification(Vertex * v) const;
    
    /*!
     * Computes the sign of two halfedges approaching and leaving the
     * boundary
     * \param he1 The halfedge entering the boundary
     * \param he2 The halfedge leaving the boundary
     * \return the perimetricity of the subpath
     */
    CGAL::Sign _sign_of_subpath(const Halfedge* he1, const Halfedge* he2) 
        const;
    
    /*!
     * Computes the sign of a halfedge and a curve approaching and leaving the
     * boundary
     * \param he1 The halfedge entering the boundary
     * \param cv2 The curve leaving the boundary
     * \param end2 The end of the curve leaving the boundary
     * \return the perimetricity of the subpath
     */
    CGAL::Sign _sign_of_subpath(const Halfedge* he1, 
                                const X_monotone_curve_2& cv2,
                                const CGAL::Arr_curve_end& end2) const;
    
    //! type of Sign_of_path functor
    typedef CGALi::Sign_of_path< Geometry_traits_2, Self > Sign_of_path;
    
    friend class CGALi::Sign_of_path< Geometry_traits_2, Self >;
    
    //@}
};

CGAL_END_NAMESPACE

#include <CGAL/Arr_topology_traits/Arr_qdx_topology_traits_2_impl.h>

#endif
