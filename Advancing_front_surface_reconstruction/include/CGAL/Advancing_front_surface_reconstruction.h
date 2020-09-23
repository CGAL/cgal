// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_H
#define CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>

#include <CGAL/disable_warnings.h>

// In order to activate lazy evaluation:
// #define LAZY

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Cartesian_converter.h>

#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include <list>
#include <set>

#include <CGAL/Advancing_front_surface_reconstruction_vertex_base_3.h>
#include <CGAL/Advancing_front_surface_reconstruction_cell_base_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/internal/AFSR/Surface_vertex_base_2.h>
#include <CGAL/internal/AFSR/Surface_face_base_2.h>
#include <CGAL/internal/AFSR/construct_surface_2.h>
#include <CGAL/internal/AFSR/construct_polyhedron.h>
#include <CGAL/internal/AFSR/write_triple_indices.h>

namespace CGAL {

  // This iterator allows to visit all contours. It has the particularity
  // that it visits the entry point of the contour twice. This allows to
  // detect that the traversal of the border is finished. One more increment
  // brings us to the next vertex.



  template < class Surface>
  class Advancing_front_surface_reconstruction_boundary_iterator {
  private:
    Advancing_front_surface_reconstruction_boundary_iterator();

    typedef typename Surface::Finite_vertices_iterator Finite_vertices_iterator;
    typedef Advancing_front_surface_reconstruction_boundary_iterator<Surface>  Self;
    typedef typename Surface::Vertex_handle            Vertex_handle;
    typedef typename Surface::Vertex                   Vertex;

    const Surface& S;
    int mark;
    Finite_vertices_iterator first_vertex;
    Vertex_handle pos;
    bool first, last;

  public:
    Advancing_front_surface_reconstruction_boundary_iterator(const Surface& S_, int m)
      : S(S_), mark(m), first_vertex(S.triangulation_3().finite_vertices_begin()), pos(first_vertex)
    {
      if (pos->number_of_incident_border() <= 0){
        advance_to_next_boundary();
      }
      first = true;
      last = false;
    }

    Advancing_front_surface_reconstruction_boundary_iterator(const Surface& S_)
      : S(S_), pos(nullptr)
    {}

    Advancing_front_surface_reconstruction_boundary_iterator(const Self& s)
      : S(s.S), mark(s.mark), first_vertex(s.first_vertex), pos(s.pos), first(s.first), last(s.last)
    {}

    bool operator==(const Self &s) const
    {
      return pos == s.pos;
    }

    bool operator!=(const Self &s) const
    {
      return pos != s.pos;
    }


    Self operator++()
    {
      if(pos == nullptr) {
        return *this;
      }
      if(first){
        advance_on_boundary();
        first = false;
      } else if (last) {
        advance_to_next_boundary();
        first = true;
        last = false;
      } else {
        advance_on_boundary();
        if(&*pos == &*first_vertex){
          last = true;
        }
      }
      return *this;
    }

    Vertex_handle operator*()
    {
      return pos;
    }

    void advance_on_boundary()
    {
      if(pos == nullptr) {
        return;
      }
      pos = pos->first_incident()->first;
      pos->set_post_mark(mark);
    }

    void advance_to_next_boundary()
    {
      if(pos == nullptr) {
        return;
      }
      do {
        first_vertex++;
      } while((first_vertex != S.triangulation_3().finite_vertices_end()) &&
              (! ((first_vertex->is_on_border())
                  && ! first_vertex->is_post_marked(mark))));
      if(first_vertex != S.triangulation_3().finite_vertices_end()) {
        pos = first_vertex;
        pos->set_post_mark(mark);
        CGAL_assertion(pos->is_on_border());

      } else {
        pos = nullptr;
      }
    }
  };

  namespace AFSR{
    struct Default_priority {
      template <typename AdvancingFront, typename Cell_handle>
      double operator() (const AdvancingFront& adv, Cell_handle& c,
                         const int& index) const
      {
        return adv.smallest_radius_delaunay_sphere (c, index);
      }
    };
  } //end of namespace AFSR


  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstructionRef

  The class `Advancing_front_surface_reconstruction` enables advanced users to provide the unstructured
  point cloud in a 3D Delaunay triangulation. The reconstruction algorithm then marks vertices and faces
  in the triangulation as being on the 2D surface embedded in 3D space, and constructs a 2D triangulation
  data structure that describes the surface.  The vertices and facets of the 2D triangulation data structure
  store handles to the vertices and faces of the 3D triangulation, which enables the user to explore the
  2D as well as 3D neighborhood of vertices and facets of the surface.

  \tparam Dt must be a `Delaunay_triangulation_3` with
  `Advancing_front_surface_reconstruction_vertex_base_3` and `Advancing_front_surface_reconstruction_cell_base_3` blended into the vertex and cell type.
  The default uses the `Exact_predicates_inexact_constructions_kernel` as geometric traits class.

  \tparam P must be a functor with `double operator()(AdvancingFront,Cell_handle,int)` returning the
  priority of the facet `(Cell_handle,int)`. This functor enables the user to choose how candidate
  triangles are prioritized. If a facet should not appear in the output,
  `infinity()` must be returned. It defaults to a functor that returns the
  `smallest_radius_delaunay_sphere()`.

  */
  template <
    class Dt = Default,
    class P = Default>
  class Advancing_front_surface_reconstruction {

    typedef typename Default::Get<Dt,Delaunay_triangulation_3<Exact_predicates_inexact_constructions_kernel, Triangulation_data_structure_3<Advancing_front_surface_reconstruction_vertex_base_3<Exact_predicates_inexact_constructions_kernel>, Advancing_front_surface_reconstruction_cell_base_3<Exact_predicates_inexact_constructions_kernel> > > >::type Triangulation;
    typedef typename Default::Get<P,AFSR::Default_priority>::type Priority;
  public:

#ifdef DOXYGEN_RUNNING
  /// \name Types
  /// @{

  /*!
    The type of the 2D triangulation data structure describing the reconstructed surface, being a model of `TriangulationDataStructure_2`.
    - The type `Triangulation_data_structure_2::Vertex` is model of the concept `TriangulationDataStructure_2::Vertex` and has additionally the
    method `vertex_3()` that returns a `#Vertex_handle` to the associated 3D vertex.
    - The type `Triangulation_data_structure_2::Face` is model of the concept `TriangulationDataStructure_2::Face` and  has additionally the
    method `facet()` that returns the associated `#Facet`, and a method `bool is_on_surface()`
    for testing if a face is part of the reconstructed surface or a face incident to a boundary edge.

    In case the surface has boundaries, the 2D surface has one vertex which is associated to the infinite
    vertex of the 3D triangulation.
  */
    typedef unspecified_type Triangulation_data_structure_2;

  /*!
  The type of the 3D triangulation.
  */
    typedef unspecified_type Triangulation_3;

  /*!
  The type of the facet priority functor.
  */
    typedef unspecified_type Priority;

  /*!
  The point type.
  */
    typedef typename Triangulation_3::Point Point;

  /*!
  The vertex handle type of the 3D triangulation.
  */
    typedef typename Triangulation_3::Vertex_handle Vertex_handle;

  /*!
  The cell handle type of the 3D triangulation.
  */
    typedef typename Triangulation_3::Cell_handle Cell_handle;

  /*!
  The facet type of the 3D triangulation.
  */
    typedef typename Triangulation_3::Facet Facet;

  /*!
    A bidirectional iterator range which enables to enumerate all points that were removed
    from the 3D Delaunay triangulation during the surface reconstruction. The value type
    of the iterator is `#Point`.
  */
    typedef unspecified_type Outlier_range;

  /*!
    A bidirectional iterator range which enables to visit all boundaries.
    The value type of the iterator is `Vertex_on_boundary_range`.
  */
    typedef unspecified_type Boundary_range;

 /*!
    A bidirectional iterator range which enables to visit all vertices on a boundary.
    The value type of the iterator is  `#Vertex_handle`
  */
    typedef unspecified_type Vertex_on_boundary_range;
  /// @}
#endif

    typedef Triangulation Triangulation_3;
    typedef typename Triangulation_3::Geom_traits Kernel;
    typedef Advancing_front_surface_reconstruction<Dt,P> Extract;
    typedef typename Triangulation_3::Geom_traits Geom_traits;

    typedef typename Kernel::FT coord_type;

    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Segment_3  Segment;
    typedef typename Kernel::Triangle_3  Triangle;
    typedef typename Kernel::Sphere_3 Sphere;

    typedef typename Triangulation_3::Cell  Cell;
    typedef typename Triangulation_3::Vertex Vertex;
    typedef typename Triangulation_3::Edge Edge;
    typedef typename Triangulation_3::Facet Facet;
    typedef typename Triangulation_3::Cell_handle  Cell_handle;
    typedef typename Triangulation_3::Vertex_handle Vertex_handle;

    typedef typename Triangulation_3::Cell_circulator  Cell_circulator;
    typedef typename Triangulation_3::Facet_circulator Facet_circulator;

    typedef typename Triangulation_3::Locate_type Locate_type;

    typedef typename Triangulation_3::Finite_cells_iterator  Finite_cells_iterator;
    typedef typename Triangulation_3::Finite_facets_iterator Finite_facets_iterator;
    typedef typename Triangulation_3::Finite_vertices_iterator  Finite_vertices_iterator;
    typedef typename Triangulation_3::Finite_edges_iterator  Finite_edges_iterator;

    typedef typename Triangulation_3::All_cells_iterator  All_cells_iterator;
    typedef typename Triangulation_3::All_facets_iterator All_facets_iterator;
    typedef typename Triangulation_3::All_vertices_iterator  All_vertices_iterator;
    typedef typename Triangulation_3::All_edges_iterator  All_edges_iterator;

    typedef typename Triangulation_3::Vertex::Edge_incident_facet Edge_incident_facet;
    typedef typename Triangulation_3::Vertex::IO_edge_type IO_edge_type;
    typedef typename Triangulation_3::Vertex::criteria criteria;
    typedef typename Triangulation_3::Vertex::Radius_edge_type Radius_edge_type;
    typedef typename Triangulation_3::Vertex::Border_elt Border_elt;
    typedef typename Triangulation_3::Vertex::Next_border_elt Next_border_elt;
    typedef typename Triangulation_3::Vertex::Intern_successors_type Intern_successors_type;
    typedef typename Triangulation_3::Vertex::Radius_ptr_type Radius_ptr_type;

    typedef typename Triangulation_3::Vertex::Incidence_request_iterator Incidence_request_iterator;
    typedef typename Triangulation_3::Vertex::Incidence_request_elt Incidence_request_elt;

    typedef std::pair< Vertex_handle, Vertex_handle > Edge_like;
    typedef CGAL::Triple< Vertex_handle, Vertex_handle, Vertex_handle > Facet_like;

    typedef std::set<Radius_ptr_type> Ordered_border_type;

    typedef typename Ordered_border_type::iterator Ordered_border_iterator;

    enum Validation_case {NOT_VALID, NOT_VALID_CONNECTING_CASE, FINAL_CASE,
                          EAR_CASE, EXTERIOR_CASE, CONNECTING_CASE};

    //=====================================================================
    //=====================================================================


    typedef const std::list<std::list<Vertex_handle> > Boundary_range;
    typedef const std::list<Vertex_handle> Vertex_on_boundary_range;

  private:

    mutable std::list<std::list<Vertex_handle> >  m_boundaries;

    Timer postprocess_timer, extend_timer, extend2_timer, init_timer;

    Triangulation_3& T;

    Ordered_border_type _ordered_border;
    int _number_of_border;

    const coord_type COS_ALPHA_SLIVER;
    coord_type COS_BETA;
    const int NB_BORDER_MAX;
    coord_type DELTA; // = sampling quality of the border
    coord_type K, min_K;
    const coord_type eps;
    const coord_type inv_eps_2; // 1/(eps^2)
    const coord_type eps_3; // test de ^3 donc points tel 1e-7 soit petit
    const criteria STANDBY_CANDIDATE;
    const criteria STANDBY_CANDIDATE_BIS;
    const criteria NOT_VALID_CANDIDATE;

    //---------------------------------------------------------------------
    //Pour une visu correcte
    //pour retenir les facettes selectionnees
    int _vh_number;
    int _facet_number;

    //---------------------------------------------------------------------
    //Pour le post traitement
    mutable int _postprocessing_counter;
    int _size_before_postprocessing;

    std::list<Point> m_outliers;

    int _number_of_connected_components;

    Vertex_handle added_vertex;
    bool deal_with_2d;
    Priority priority;
    int max_connected_component;
    double K_init, K_step;
    std::list<Vertex_handle> interior_edges;
    std::list< Incidence_request_elt > incidence_requests;
    typename std::list< Incidence_request_elt >::iterator sentinel;
    typename std::list<Vertex_handle>::iterator ie_sentinel;
    std::list<Next_border_elt> nbe_pool;
    std::list<Intern_successors_type> ist_pool;


    Intern_successors_type* new_border()
    {
      nbe_pool.resize(nbe_pool.size()+1);

      Next_border_elt* p1 = & nbe_pool.back();
      nbe_pool.resize(nbe_pool.size()+1);
      Next_border_elt* p2 = & nbe_pool.back();

      Intern_successors_type ist(p1,p2);
      ist_pool.push_back(ist);

      Intern_successors_type* ret = &ist_pool.back();

      ret->first->first = nullptr;
      ret->second->first = nullptr;
      return ret;
    }


    inline bool is_on_border(Vertex_handle vh, const int& i) const
    {
      if (vh->m_incident_border == nullptr) return false; //vh is interior
      if (vh->m_incident_border->first->first != nullptr)
        {
          if (vh->m_incident_border->second->first != nullptr)
            return ((vh->m_incident_border->first->second.second == i)||
                    (vh->m_incident_border->second->second.second == i));
          return (vh->m_incident_border->first->second.second == i);
        }
      return false; //vh is still exterior
    }


    void remove_border_edge(Vertex_handle w, Vertex_handle v)
    {
      if (w->m_incident_border != nullptr)
        {
          if (w->m_incident_border->second->first == v)
            {
              w->m_incident_border->second->first = nullptr;
              set_interior_edge(w,v);
              return;
            }
          if (w->m_incident_border->first->first == v)
            {
              if (w->m_incident_border->second->first != nullptr)
                {
                  Next_border_elt* tmp = w->m_incident_border->first;
                  w->m_incident_border->first = w->m_incident_border->second;
                  w->m_incident_border->second = tmp;
                  w->m_incident_border->second->first = nullptr;
                  set_interior_edge(w,v);
                  return;
                }
              else
                {
                  w->m_incident_border->first->first = nullptr;
                  set_interior_edge(w,v);
                  return;
                }
            }
        }
    }


    inline bool is_interior_edge(Vertex_handle w, Vertex_handle v) const
    {

      bool r1;
      if(w->m_ie_first == ie_sentinel){
        r1 = false;
      }else {
        typename std::list<Vertex_handle>::iterator b(w->m_ie_first), e(w->m_ie_last);
        e++;
        typename std::list<Vertex_handle>::iterator r = std::find(b, e, v);
        r1 = ( r != e);
      }

      return r1;
    }

    //-------------------------------------------------------------------
    // pour gerer certaines aretes interieures: a savoir celle encore connectee au
    // bord (en fait seule, les aretes interieures reliant 2 bords nous
    // interressent...)

    inline void set_interior_edge(Vertex_handle w, Vertex_handle v)
    {
      if(w->m_ie_last == ie_sentinel){ // empty set
        CGAL_assertion(w->m_ie_first == w->m_ie_last);
        w->m_ie_last = interior_edges.insert(w->m_ie_last, v);
        w->m_ie_first = w->m_ie_last;
      } else {
        typename std::list<Vertex_handle>::iterator e(w->m_ie_last);
        e++;
#ifdef DEBUG
        typename std::list<Vertex_handle>::iterator r = std::find(w->m_ie_first, e, v);
        CGAL_assertion(r == e);
#endif
        w->m_ie_last = interior_edges.insert(e, v);
      }
    }


    inline void remove_interior_edge(Vertex_handle w, Vertex_handle v)
    {
      if(w->m_ie_first == ie_sentinel){
        CGAL_assertion(w->m_ie_last == w->m_ie_first);
      } else if(w->m_ie_first == w->m_ie_last){ // there is only one element
        if(*(w->m_ie_first) == v){
          interior_edges.erase(w->m_ie_first);
          w->m_ie_last = ie_sentinel;
          w->m_ie_first = w->m_ie_last;
        }
      } else {
        typename std::list<Vertex_handle>::iterator b(w->m_ie_first), e(w->m_ie_last);
        e++;
        typename std::list<Vertex_handle>::iterator r = std::find(b, e, v);
        if(r != e){
          if(r == w->m_ie_first){
            w->m_ie_first++;
          }
          if(r == w->m_ie_last){
            w->m_ie_last--;
          }
          interior_edges.erase(r);
        }
      }
    }


    //-------------------------------------------------------------------

    inline void set_incidence_request(Vertex_handle w, const Incidence_request_elt& ir)
    {
      if(w->m_ir_last == sentinel ){
        CGAL_assertion(w->m_ir_first == w->m_ir_last);
        w->m_ir_last = incidence_requests.insert(w->m_ir_last, ir);
        w->m_ir_first = w->m_ir_last;
      } else {
        typename std::list<Incidence_request_elt>::iterator e(w->m_ir_last);
        e++;
        w->m_ir_last = incidence_requests.insert(e, ir);
      }
    }

    inline bool is_incidence_requested(Vertex_handle w) const
    {
      if(w->m_ir_last == sentinel ){
        CGAL_assertion(w->m_ir_first == sentinel );
      }
      return (w->m_ir_last != sentinel );
    }

    inline Incidence_request_iterator incidence_request_begin(Vertex_handle w)
    {
      return w->m_ir_first;
    }

    inline Incidence_request_iterator incidence_request_end(Vertex_handle w)
    {
      if(w->m_ir_last != sentinel ){
        CGAL_assertion(w->m_ir_first != sentinel );
        Incidence_request_iterator it(w->m_ir_last);
        it++;
        return it;
      }
      return w->m_ir_last;
    }

    inline void erase_incidence_request(Vertex_handle w)
    {
      if(w->m_ir_last != sentinel ){
        CGAL_assertion(w->m_ir_first != sentinel );
        w->m_ir_last++;
        incidence_requests.erase(w->m_ir_first, w->m_ir_last);
        w->m_ir_first = sentinel ;
        w->m_ir_last = sentinel ;
      }
    }


    void re_init(Vertex_handle w)
    {
      if (w->m_incident_border != nullptr)
        {
          w->delete_border();
        }

      if(w->m_ir_first != sentinel ){
        CGAL_assertion(w->m_ir_last != sentinel );
        typename std::list< Incidence_request_elt >::iterator b(w->m_ir_first), e(w->m_ir_last);
        e++;
        incidence_requests.erase(b, e);
        w->m_ir_first = sentinel ;
        w->m_ir_last = sentinel ;
      }

      w->m_incident_border = new_border();
      w->m_mark = -1;
      w->m_post_mark = -1;
    }



    inline void dec_mark(Vertex_handle w)
    {
      w->m_mark--;
      if(w->m_mark == 0)
        {
          w->delete_border();
          erase_incidence_request(w);
        }
    }



    void initialize_vertices_and_cells()
    {

      Incidence_request_elt ire;
      incidence_requests.clear();
      incidence_requests.push_back(ire);
      sentinel = incidence_requests.begin();

      interior_edges.clear();
      interior_edges.push_back(Vertex_handle());
      ie_sentinel = interior_edges.begin();

      for(All_vertices_iterator fit = T.all_vertices_begin();
          fit != T.all_vertices_end();
          ++fit){
        fit->m_ie_first =  fit->m_ie_last = ie_sentinel;
        fit->m_ir_first = fit->m_ir_last = sentinel ;
        fit->m_incident_border = new_border();
      }
    }

    //-------------------- DESTRUCTOR -----------------------------------

    void clear_vertex(Vertex_handle w)
    {
      if (w->m_incident_border != nullptr)
        {
          w->delete_border();
        }
      if(w->m_ir_first != sentinel ){
        CGAL_assertion(w->m_ir_last != sentinel );
        typename std::list< Incidence_request_elt >::iterator b(w->m_ir_first), e(w->m_ir_last);
        e++;
        incidence_requests.erase(b, e);
      }

      if(w->m_ie_first != ie_sentinel){
        CGAL_assertion(w->m_ie_last != ie_sentinel);
        typename std::list<Vertex_handle>::iterator b(w->m_ie_first), e(w->m_ie_last);
        e++;
        interior_edges.erase(b, e);
      }
    }


    void clear_vertices()
    {
      for (All_vertices_iterator vit = T.all_vertices_begin();
           vit != T.all_vertices_end();
           ++vit){
        clear_vertex(vit);
      }
    }


  public:
    /// \name Creation
    /// @{

    /*!
    Constructor for the unstructured point cloud given as 3D Delaunay triangulation.
    */
    Advancing_front_surface_reconstruction(Triangulation_3& dt,
                                           Priority priority = Priority())
      : T(dt), _number_of_border(1), COS_ALPHA_SLIVER(-0.86),
        NB_BORDER_MAX(15), DELTA(.86), min_K(infinity()),
        eps(1e-7), inv_eps_2(coord_type(1)/(eps*eps)), eps_3(eps*eps*eps),
        STANDBY_CANDIDATE(3), STANDBY_CANDIDATE_BIS(STANDBY_CANDIDATE+1),
        NOT_VALID_CANDIDATE(STANDBY_CANDIDATE+2),
      _vh_number(static_cast<int>(T.number_of_vertices())), _facet_number(0),
      _postprocessing_counter(0), _size_before_postprocessing(0), _number_of_connected_components(0),
      deal_with_2d(false), priority(priority), max_connected_component(-1), K_init(1.1), K_step(.1)

    {
      if(T.dimension() == 2){
        deal_with_2d = true;
        Finite_vertices_iterator it = T.finite_vertices_begin();
        const Point& p = it->point();
        ++it;
        const Point& q = it->point();
        do{
          ++it;
        }while(collinear(p,q,it->point()));
        const Point& r = it->point();
        Vector u = q-r;
        Vector v = q-p;
        Vector w = r-p;
        Vector vw = cross_product(v,w);
        double len = (std::max)(u*u,(std::max)(v*v,w*w));
        Point s = p + 10* len * (vw/(vw*vw));
        added_vertex = T.insert(s);
      }
    }

    /// @}


    /*
      ~Advancing_front_surface_reconstruction()
      {

      std::cerr << "postprocessing" << postprocess_timer.time() << std::endl;
      std::cerr << "extend        " << extend_timer.time() << std::endl;
      std::cerr << "extend2       " << extend2_timer.time() << std::endl;
      std::cerr << "init          " << postprocess_timer.time() << std::endl;
      std::cerr << "#outliers     " << number_of_outliers() << std::endl;
      }
    */

    typedef Advancing_front_surface_reconstruction_boundary_iterator<Extract> Boundary_iterator;

    Boundary_iterator boundaries_begin() const
    {
      return Boundary_iterator(*this, next_mark());
    }


    Boundary_iterator boundaries_end() const
    {
      return Boundary_iterator(*this);
    }

    typedef std::list<Point> Outlier_range;

    typedef CGAL::Triangulation_data_structure_2<AFSR::Surface_vertex_base_2<Kernel,Vertex_handle>,
                                                 AFSR::Surface_face_base_2<Kernel, typename Triangulation_3::Facet> > TDS_2;

    typedef TDS_2 Triangulation_data_structure_2;

    mutable TDS_2 _tds_2;

    mutable typename TDS_2::Vertex_handle _tds_2_inf;

    /// \name Operations
    /// @{

    /*!
    runs the surface reconstruction function.

    \param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
           are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`.
           Described in Section \ref AFSR_Boundaries
    \param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates.
           Described in Section \ref AFSR_Selection

    */
    void run(double radius_ratio_bound=5, double beta= 0.52)
    {
      K = radius_ratio_bound;
      COS_BETA = cos(beta);
      if(T.dimension() < 3){
        return;
      }
      initialize_vertices_and_cells();
      bool re_init = false;
      do
        {
          _number_of_connected_components++;

          if ( (re_init = init(re_init)) )
            {
              //std::cerr << "Growing connected component " << _number_of_connected_components << std::endl;
              extend_timer.start();
              extend();
              extend_timer.stop();

              if ((number_of_facets() > static_cast<int>(T.number_of_vertices()))&&
                  (NB_BORDER_MAX > 0))
                // en principe 2*nb_sommets = nb_facettes: y a encore de la marge!!!
                {
                  while(postprocessing()){
                    extend2_timer.start();
                    extend();

                    extend2_timer.stop();
                  }
                }
            }
        }while(re_init &&
               ((_number_of_connected_components < max_connected_component)||
                (max_connected_component < 0)));

      _tds_2_inf = AFSR::construct_surface(_tds_2, *this);

      boundaries();
      clear_vertices();
    }

    /*!
    returns the reconstructed surface.
    */
    const Triangulation_data_structure_2& triangulation_data_structure_2() const
    {
      return _tds_2;
    }

    /*!
    returns the underlying 3D Delaunay triangulation.
    */
    Triangulation_3&
    triangulation_3() const
    {
      return T;
    }

    /*!
    returns an iterator range over the outliers.
    */
    const Outlier_range& outliers() const
    {
      return m_outliers;
    }

    /*!
    returns an iterator range over the boundaries.
    */
    const Boundary_range& boundaries() const
    {
      if(has_boundaries() && m_boundaries.empty()){
        Boundary_iterator b =  boundaries_begin();
        Boundary_iterator e =  boundaries_end();
        for(; b!= e; ++b){
          Vertex_handle v = *b;
          std::list<Vertex_handle> border;
          m_boundaries.push_back(border);
          do {
            m_boundaries.back().push_back(*b);
            ++b;
          }while(*b != v);
        }
      }

      return  m_boundaries;
    }
    /// @}

/// \name Predicates
/// @{

    /*!
    returns `true` if the reconstructed surface has boundaries.
    */
    bool
    has_boundaries() const
    {
      return _tds_2_inf != typename TDS_2::Vertex_handle();
    }

    /*!
    returns `true` if the facet is on the surface.
    */
    bool
    has_on_surface(Facet f) const
    {
      return f.first->has_facet_on_surface(f.second);
    }

    /*!
    returns `true` if the facet `f2` is on the surface.
    */
    bool
    has_on_surface(typename Triangulation_data_structure_2::Face_handle f2) const
    {
      return f2->is_on_surface();
    }

    /*!
    returns `true` if the vertex `v2` is on the surface.
    */
    bool
    has_on_surface(typename Triangulation_data_structure_2::Vertex_handle v2) const
    {
      return v2 != _tds_2_inf;
    }
    /// @}

    int number_of_connected_components() const
    {
      return _number_of_connected_components;
    }

    int number_of_facets() const
    {
      return _facet_number;
    }


    int number_of_vertices() const
    {
      return _vh_number;
    }


    int number_of_outliers() const
    {
      return static_cast<int>(m_outliers.size());
    }

    typedef typename std::list<Point>::const_iterator Outlier_iterator;


    Outlier_iterator outliers_begin() const
    {
      return m_outliers.begin();
    }


    Outlier_iterator m_outliers_end() const
    {
      return m_outliers.end();
    }

    int next_mark() const
    {
      _postprocessing_counter++;
      return _postprocessing_counter;
    }


    Next_border_elt* border_elt(const Vertex_handle& v1, const Vertex_handle& v2) const
    {
      return v1->border_elt(v2);
    }


    //public

    IO_edge_type* border_IO_elt(const Vertex_handle& v1, const Vertex_handle& v2)
    {
      return &border_elt(v1,v2)->second.first.second;
    }


    IO_edge_type* set_border_elt(const Vertex_handle& v1, const Vertex_handle& v2,
                                 const Border_elt& e)
    {
      v1->set_next_border_elt(Next_border_elt (v2, e));
      return border_IO_elt(v1, v2);
    }


    IO_edge_type* set_again_border_elt(const Vertex_handle& v1, const Vertex_handle& v2,
                                       const Border_elt& e)
    {
      border_elt(v1,v2)->second = e;
      return border_IO_elt(v1, v2);
    }

    //---------------------------------------------------------------------
    bool is_border_elt(Edge_like& key, Border_elt& result) const
    {
      Next_border_elt* it12 = border_elt(key.first, key.second);
      if (it12 != nullptr)
        {
          result = it12->second;
          return true;
        }

      Next_border_elt* it21 =  border_elt(key.second, key.first);
      if (it21 != nullptr)
        {
          result = it21->second;
          std::swap(key.first, key.second);
          return true;
        }
      return false;
    }

    //---------------------------------------------------------------------
    bool is_border_elt(Edge_like& key) const {
      Next_border_elt* it12 =  border_elt(key.first, key.second);
      if (it12 != nullptr)
        {
          return true;
        }

      Next_border_elt* it21 =  border_elt(key.second, key.first);
      if (it21 != nullptr)
        {
          std::swap(key.first, key.second);
          return true;
        }
      return false;
    }

    //---------------------------------------------------------------------
    bool is_ordered_border_elt(const Edge_like& key, Border_elt& result) const
    {
      Next_border_elt* it12 =  border_elt(key.first, key.second);
      if (it12 != nullptr)
        {
          result = it12->second;
          return true;
        }
      return false;
    }

    //---------------------------------------------------------------------
    void
    remove_border_elt(const Edge_like& ordered_key)
    {
      remove_border_edge(ordered_key.first, ordered_key.second);
    }

    //---------------------------------------------------------------------
    bool is_ordered_border_elt(const Edge_like& e,
                               IO_edge_type* &ptr) const
    {
      Vertex_handle v1 = e.first;

      Next_border_elt* it12 =  border_elt(v1, e.second);
      if (it12 != nullptr)
        {
          ptr = &it12->second.first.second;
          return true;
        }
      return false;
    }

    //---------------------------------------------------------------------
    void set_incidence_request(const Vertex_handle& v,
                               const criteria& value,
                               const Edge_like& e)
    {
      Incidence_request_elt incident_elt(value, e);
      set_incidence_request(v, incident_elt);
    }

    //---------------------------------------------------------------------
    bool is_interior_edge(const Edge_like& key) const
    // pour gerer certaines aretes interieures: a savoir celle encore connectee au
    // bord (en fait seule, les aretes interieures reliant 2 bords nous
    // interressent...)
    {
      return (is_interior_edge(key.first, key.second)||
              is_interior_edge(key.second, key.first));
    }

    //---------------------------------------------------------------------
#ifdef AFSR_LAZY

    coord_type lazy_squared_radius(const Cell_handle& c)
    {
      if (c->lazy_squared_radius() != nullptr)
        return *(c->lazy_squared_radius());

      c->set_lazy_squared_radius
        (CGAL::squared_radius(c->vertex(0)->point(),
                              c->vertex(1)->point(),
                              c->vertex(2)->point(),
                              c->vertex(3)->point()));
      return *(c->lazy_squared_radius());
    }

    Point lazy_circumcenter(const Cell_handle& c)
    {
      if (c->lazy_circumcenter() != nullptr)
        return *(c->lazy_circumcenter());

      c->set_lazy_circumcenter
        (CGAL::circumcenter(c->vertex(0)->point(),
                      c->vertex(1)->point(),
                      c->vertex(2)->point(),
                      c->vertex(3)->point()));
      return *(c->lazy_circumcenter());
    }

#endif //NOLAZY

    //---------------------------------------------------------------------
    Edge_incident_facet next(const Edge_incident_facet& e) const
    {
      Cell_handle c = e.first.first;
      int i = e.second;
      int i1 = e.first.second, i2 = e.first.third;
      int i3 = (6 - e.second - i1 - i2);

      Cell_handle n = c->neighbor(i);
      int j1 = n->index(c->vertex(i1)), j2 = n->index(c->vertex(i2));
      int j =  n->index(c->vertex(i3));
      return Edge_incident_facet(Edge(n, j1, j2), j);
    }

    //---------------------------------------------------------------------
    Edge_incident_facet previous(const Edge_incident_facet& e) const
    {
      Cell_handle c = e.first.first;
      int i = e.second;
      int i1 = e.first.second, i2 = e.first.third;
      int i3 = (6 - e.second - i1 - i2);

      Cell_handle n = c->neighbor(i3);
      int j1 = n->index(c->vertex(i1)), j2 = n->index(c->vertex(i2));
      int j =  n->index(c->vertex(i));
      return Edge_incident_facet(Edge(n, j1, j2), j);
    }

    //---------------------------------------------------------------------
    bool
    my_coplanar(const Point& p, const Point& q,
                const Point& r, const Point& s) const
    {
      coord_type qpx = q.x()-p.x();
      coord_type qpy = q.y()-p.y();
      coord_type qpz = q.z()-p.z();
      coord_type rpx = r.x()-p.x();
      coord_type rpy = r.y()-p.y();
      coord_type rpz = r.z()-p.z();
      coord_type spx = s.x()-p.x();
      coord_type spy = s.y()-p.y();
      coord_type spz = s.z()-p.z();

      coord_type den = CGAL::determinant(qpx,qpy,qpz,
                                         rpx,rpy,rpz,
                                         spx,spy,spz);
      return (CGAL_NTS abs(den) < eps_3);
    }

    //---------------------------------------------------------------------
    bool
    my_collinear(const Point& p, const Point& q, const Point& s) const
    {
      coord_type psx = p.x()-s.x();
      coord_type psy = p.y()-s.y();
      coord_type psz = p.z()-s.z();
      coord_type qsx = q.x()-s.x();
      coord_type qsy = q.y()-s.y();
      coord_type qsz = q.z()-s.z();
      coord_type rsx = psy*qsz-psz*qsy;
      coord_type rsy = psz*qsx-psx*qsz;
      coord_type rsz = psx*qsy-psy*qsx;

      coord_type den = CGAL::determinant(psx,psy,psz,
                                         qsx,qsy,qsz,
                                         rsx,rsy,rsz);

      return (CGAL_NTS abs(den) < eps_3);
    }

    //---------------------------------------------------------------------
    void
    select_facet(const Cell_handle& c, const int& i)
    {
      c->select_facet(i);
      _facet_number++;
      c->set_facet_number(i, _facet_number);
    }


    void
    unselect_facet(const Cell_handle& c, const int& i)
    {
      c->unselect_facet(i);
    }


    int
    number_of_border_edges()
    {
      int _border_count(0);
      for(Finite_edges_iterator e_it=T.finite_edges_begin();
          e_it!=T.finite_edges_end();
          e_it++)
        {
          Cell_handle c = (*e_it).first;
          int i1 = (*e_it).second, i2 = (*e_it).third;
          Edge_like key(c->vertex(i1), c->vertex(i2));

          if (is_border_elt(key))
            _border_count++;
        }
      return _border_count;
    }


    //=====================================================================


    /// \name Priority values
    /// @{

    /*!

      computes the priority of the facet `(c,index)` such that the
      facet with the smallest radius of Delaunay sphere has the
      highest priority.

    \param c handle to the cell containing the facet
    \param index index of the facet in `c`

    */
    coord_type
    smallest_radius_delaunay_sphere(const Cell_handle& c,
                                    const int& index) const
    {
      int i1, i2, i3;

      if(deal_with_2d && ( (c->vertex((index+1) & 3) == added_vertex)
                           || (c->vertex((index+2) & 3) == added_vertex)
                           || (c->vertex((index+3) & 3) == added_vertex) ))
        {
          return infinity();
        }
      Cell_handle n = c->neighbor(index);
      // lazy evaluation ...
      coord_type value = c->smallest_radius(index);
      if ((value >= 0)&&(n->smallest_radius(n->index(c)) == value))
        return value;

      const Point& cp0 = c->vertex(index)->point();
      const Point& cp1 = c->vertex((index+1) & 3)->point();
      const Point& cp2 = c->vertex((index+2) & 3)->point();
      const Point& cp3 = c->vertex((index+3) & 3)->point();

      const Point& np0 = n->vertex(0)->point();
      const Point& np1 = n->vertex(1)->point();
      const Point& np2 = n->vertex(2)->point();
      const Point& np3 = n->vertex(3)->point();

      bool c_is_plane(my_coplanar(cp0, cp1, cp2, cp3));
      bool n_is_plane(my_coplanar(np0, np1, np2, np3));

      bool c_is_infinite(T.is_infinite(c));
      bool n_is_infinite(T.is_infinite(n));
      if ((c_is_plane && n_is_plane)||
          (c_is_plane && n_is_infinite)||
          (n_is_plane && c_is_infinite)||
          my_collinear(cp1, cp2, cp3))
        value = infinity();
      else
        {
          if (c_is_infinite||n_is_infinite||c_is_plane||n_is_plane)
            {
              int ind;
              Cell_handle cc;
              if(c_is_infinite||c_is_plane)
                {
                  cc = n;
                  ind = n->index(c);
                }
              else
                {
                  cc = c;
                  ind = index;
                }
              i1 = (ind+1) & 3;
              i2 = (ind+2) & 3;
              i3 = (ind+3) & 3;

              const Point& pp0 = cc->vertex(ind)->point();
              const Point& pp1 = cc->vertex(i1)->point();
              const Point& pp2 = cc->vertex(i2)->point();
              const Point& pp3 = cc->vertex(i3)->point();

              Sphere facet_sphere(pp1, pp2, pp3);
              if (squared_distance(facet_sphere.center(), pp0) <
                  facet_sphere.squared_radius())
                {
#ifdef AFSR_LAZY
                  value = lazy_squared_radius(cc);
#else
                  // qualified with CGAL, to avoid a compilation error with clang
                  if(volume(pp0, pp1, pp2, pp3) != 0){
                    value = CGAL::squared_radius(pp0, pp1, pp2, pp3);
                  } else {
                    typedef Exact_predicates_exact_constructions_kernel EK;
                    Cartesian_converter<Kernel, EK> to_exact;
                    Cartesian_converter<EK, Kernel> back_from_exact;
                    value = back_from_exact(CGAL::squared_radius(to_exact(pp0),
                                                                 to_exact(pp1),
                                                                 to_exact(pp2),
                                                                 to_exact(pp3)));
                  }
#endif
                }
              else
                value = facet_sphere.squared_radius();
            }
          else
            {
              Point cc, cn;
#ifdef AFSR_LAZY
              cc = lazy_circumcenter(c);
              cn = lazy_circumcenter(n);
#else
              cc = CGAL::circumcenter(cp0, cp1, cp2, cp3);
              cn = CGAL::circumcenter(np0, np1, np2, np3);
#endif
              // computation of the distance of  cp1  to the  dual segment cc, cn...
              Vector V(cc - cn), Vc(cc - cp1), Vn(cp1 - cn);
              coord_type ac(V * Vc), an(V * Vn), norm_V(V * V);
              if ((ac > 0) && (an > 0))
                {
                  value = (Vc*Vc) - ac*ac/norm_V;
                  if ((value < 0)||(norm_V > inv_eps_2)){
                    // qualified with CGAL, to avoid a compilation error with clang
                    value = CGAL::squared_radius(cp1, cp2, cp3);
                  }
                }
              else
                {
                  if (ac <= 0)
                    value = squared_distance(cc, cp1);
                  else // (an <= 0)
                    value = squared_distance(cn, cp1);
                }
            }
        }
      // cache computed values
      c->set_smallest_radius(index, value);
      n->set_smallest_radius(n->index(c), value);

      return value;
    }

    /*!

      returns the infinite floating value that prevents a facet to be used.
    */
    coord_type infinity() const { return std::numeric_limits<coord_type>::infinity(); }
    /// @}

    //---------------------------------------------------------------------
    // For a border edge e we determine the incident facet which has the highest
    // chance to be a natural extension of the surface

    Radius_edge_type
    compute_value(const Edge_incident_facet& e)
    {
      Cell_handle c = e.first.first;
      int i = e.second;
      int i1 = e.first.second, i2 = e.first.third;
      int i3 = 6 - e.second - i1 - i2;

      Edge_incident_facet e_it = e;

      coord_type min_valueP = NOT_VALID_CANDIDATE,
        min_valueA = infinity();
      Facet min_facet, min_facetA;
      bool border_facet(false);

      coord_type pscal;
      const Point& p1 = c->vertex(i1)->point();
      const Point& p2 = c->vertex(i2)->point();
      const Point& pc = c->vertex(i3)->point();

      Vector P2P1 = p1-p2, P2Pn, PnP1;

      Vector v2, v1 = cross_product(pc-p2, P2P1);

      coord_type norm, norm1 = v1*v1;
      coord_type norm12 = P2P1*P2P1;

      e_it = next(e_it);

      do
        {
          Cell_handle neigh = e_it.first.first;
          Facet facet_it(neigh, e_it.second);

          if (!T.is_infinite(facet_it))
            {
              int n_ind = facet_it.second;
              int n_i1 = e_it.first.second;
              int n_i2 = e_it.first.third;
              int n_i3 = 6 - n_ind - n_i1 - n_i2;

              coord_type tmp=0;

              tmp = priority (*this, neigh, n_ind);

              Edge_like el1(neigh->vertex(n_i1),neigh->vertex(n_i3)),
                el2(neigh->vertex(n_i2),neigh->vertex(n_i3));

              if ((tmp != infinity())&&
                  neigh->vertex(n_i3)->not_interior()&&
                  (!is_interior_edge(el1))&&(!is_interior_edge(el2)))
                {
                  const Point& pn = neigh->vertex(n_i3)->point();

                  P2Pn = pn-p2;
                  v2 = cross_product(P2P1,P2Pn);

                  //pas necessaire de normer pour un bon echantillon:
                  //            on peut alors tester v1*v2 >= 0
                  norm =  sqrt(norm1 * (v2*v2));
                  pscal = v1*v2;
                  // check if the triangle will produce a sliver on the surface
                  bool sliver_facet = (pscal <= COS_ALPHA_SLIVER*norm);

                  if (!sliver_facet)
                    {
                      if (tmp < min_valueA)
                        {
                          PnP1 = p1-pn;
                          // DELTA represente la qualite d'echantillonnage du bord
                          // We skip triangles having an internal angle along e
                          // whose cosinus is smaller than -DELTA
                          // that is the angle is larger than arcos(-DELTA)
                          border_facet = !((P2P1*P2Pn >=
                                            -DELTA*sqrt(norm12*(P2Pn*P2Pn)))&&
                                           (P2P1*PnP1 >=
                                            -DELTA*sqrt(norm12*(PnP1*PnP1))));
                          // \todo investigate why we simply do not skip this triangle
                          // but continue looking for a better candidate
                          // if (!border_facet){
                          min_facetA = facet_it;
                          min_valueA = tmp;
                          min_valueP = pscal/norm;
                        }
                    }
                }
            }
          e_it = next(e_it);
        }
      while(e_it.first.first != c);

      criteria value;

      if ((min_valueA == infinity()) || border_facet) // bad facets case
        {
          min_facet = Facet(c, i); // !!! sans aucune signification....
          value = NOT_VALID_CANDIDATE; // Attention a ne pas inserer dans PQ
        }
      else
        {
          min_facet = min_facetA;

          //si on considere seulement la pliure value appartient a [0, 2]
          //value = coord_type(1) - min_valueP;

          // si la pliure est bonne on note suivant le alpha sinon on prend en compte la
          // pliure seule... pour discriminer entre les bons slivers...
          // si on veut discriminer les facettes de bonnes pliures plus finement
          // alors -(1+1/min_valueA) app a [-inf, -1]
          // -min_valueP app a [-1, 1]

          if (min_valueP > COS_BETA)
            value = -(coord_type(1) + coord_type(1)/min_valueA);
          else
            {
              //on refuse une trop grande non-uniformite
              coord_type tmp = priority (*this, c, i);
              if (min_valueA <= K * tmp)
                value = - min_valueP;
              else
                {
                  value = STANDBY_CANDIDATE; // tres mauvais candidat mauvaise pliure
                  // + grand alpha... a traiter plus tard....
                  min_K =
                    (std::min)(min_K,
                               min_valueA/tmp);
                }
            }
        }

      Cell_handle n = min_facet.first;
      int ni1 = n->index(c->vertex(i1)), ni2 = n->index(c->vertex(i2));

      return
        Radius_edge_type(value, IO_edge_type(e, Edge_incident_facet
                                             (Edge(n, ni1, ni2),
                                              min_facet.second)));
    }

    //=====================================================================
    // The parameter re_init is false the first time only
    // Returns true, iff it found a face where the next surface can grow
    bool
    init(const bool& re_init)
    {
      init_timer.start();
      Facet min_facet;
      coord_type min_value = infinity();
      int i1, i2, i3;

      if (!re_init){
        Finite_facets_iterator end = T.finite_facets_end();
        for(Finite_facets_iterator facet_it = T.finite_facets_begin();
            facet_it != end;
            ++facet_it)
          {
            coord_type value = priority (*this, (*facet_it).first,
                                         (*facet_it).second);
            if (value < min_value)
              {
                min_facet = *facet_it;
                min_value = value;
              }
          }
      }else{ //if (re_init)
        Finite_facets_iterator end = T.finite_facets_end();
        for(Finite_facets_iterator facet_it = T.finite_facets_begin();
            facet_it != end;
            ++facet_it)
          {
            Cell_handle c = (*facet_it).first;
            int index = (*facet_it).second;
            if (c->vertex((index+1) & 3)->is_exterior())
              if (c->vertex((index+2) & 3)->is_exterior())
                if (c->vertex((index+3) & 3)->is_exterior())
                  {
                    coord_type value = priority (*this, c, index);

                    // we might not want the triangle, for example because it is too large
                    if(value == infinity()){
                      value = min_value;
                    }

                    if (value < min_value)
                      {
                        min_facet = *facet_it;
                        min_value = value;
                      }
                  }
          }
      }

      if (min_value != infinity())
        {
          Cell_handle c_min = min_facet.first;

          int ind = min_facet.second;
          i1 = (ind+1) & 3;
          i2 = (ind+2) & 3;
          i3 = (ind+3) & 3;

          Radius_edge_type e12, e23, e31;

          e12 = compute_value(Edge_incident_facet(Edge(c_min, i1, i2), ind));
          e23 = compute_value(Edge_incident_facet(Edge(c_min, i2, i3), ind));
          e31 = compute_value(Edge_incident_facet(Edge(c_min, i3, i1), ind));

          IO_edge_type* p12 = set_border_elt(c_min->vertex(i1), c_min->vertex(i2),
                                             Border_elt(e12, _number_of_border));
          IO_edge_type* p23 = set_border_elt(c_min->vertex(i2), c_min->vertex(i3),
                                             Border_elt(e23, _number_of_border));
          IO_edge_type* p31 = set_border_elt(c_min->vertex(i3), c_min->vertex(i1),
                                             Border_elt(e31, _number_of_border));

          c_min->vertex(i1)->inc_mark();
          c_min->vertex(i2)->inc_mark();
          c_min->vertex(i3)->inc_mark();
          _ordered_border.insert(Radius_ptr_type (e12.first, p12));
          _ordered_border.insert(Radius_ptr_type (e23.first, p23));
          _ordered_border.insert(Radius_ptr_type (e31.first, p31));

          select_facet(c_min, ind);
          init_timer.stop();
          return true;
        }
      init_timer.stop();
      return false;
    }

    //---------------------------------------------------------------------
    // test de reciprocite avant de recoller une oreille anti-singularite
    int
    test_merge(const Edge_like& ordered_key, const Border_elt& result,
               const Vertex_handle& v, const coord_type& ear_alpha)
    {
      Edge_incident_facet Ifacet = result.first.second.first;

      const Point& p1 = (ordered_key.first)->point();
      const Point& p2 = (ordered_key.second)->point();
      const Point& pc = v->point();

      Cell_handle neigh = Ifacet.first.first;
      int n_ind = Ifacet.second;
      int n_i1 = Ifacet.first.second;
      int n_i2 = Ifacet.first.third;
      int n_i3 = (6 - n_ind - n_i1 - n_i2);

      const Point& pn = neigh->vertex(n_i3)->point();
      Vector v1 = cross_product(pc-p2,p1-p2),
        v2 = cross_product(p1-p2,pn-p2);
      coord_type norm = sqrt((v1*v1)*(v2*v2));

      if (v1*v2 > COS_BETA*norm)
        return 1; // label bonne pliure sinon:

      if (ear_alpha <= K * priority(*this, neigh, n_ind))
        return 2; // label alpha coherent...

      return 0; //sinon oreille a rejeter...
    }


    //---------------------------------------------------------------------
    void
    ordered_map_erase(const criteria& value, const IO_edge_type* pkey)
    {
      _ordered_border.erase(Radius_ptr_type(value,(IO_edge_type*)pkey));
    }

    //---------------------------------------------------------------------
    void
    force_merge(const Edge_like& ordered_key, const Border_elt& result)
    {
      criteria value = result.first.first;
      IO_edge_type* pkey = border_IO_elt(ordered_key.first, ordered_key.second);

      ordered_map_erase(value, pkey);

      remove_border_elt(ordered_key);
    }

    //---------------------------------------------------------------------
    void
    dequeue_incidence_request(const Vertex_handle& v)
    {
      if (is_incidence_requested(v))
        {
          for(Incidence_request_iterator v_it = incidence_request_begin(v);
              v_it != incidence_request_end(v);
              v_it++)
            {
              IO_edge_type* ptr;

              if (is_ordered_border_elt(v_it->second, ptr))
                _ordered_border.insert(Radius_ptr_type(v_it->first, ptr));
            }
          erase_incidence_request(v);
        }
    }


    //---------------------------------------------------------------------
    void
    merge_ear(const Edge_like& ordered_el1, const Border_elt& result1,
              const Edge_like& ordered_key,
              const Vertex_handle& v1, const Vertex_handle& v2,
              const Edge_incident_facet& edge_Ifacet_2)
    {
      remove_border_elt(ordered_key);
      force_merge(ordered_el1, result1);
      Radius_edge_type e2 = compute_value(edge_Ifacet_2);
      IO_edge_type* p2;
      if (ordered_el1.first == v1)
        p2 = set_border_elt(v2, ordered_el1.second,
                            Border_elt(e2,result1.second));
      else
        p2 = set_border_elt(ordered_el1.first, v2,
                            Border_elt(e2,result1.second));
      dec_mark(v1);

      _ordered_border.insert(Radius_ptr_type(e2.first, p2));

      //depiler les eventuelles requettes de connections avortees... zones etoilees,
      //en effet le bord a change donc on peut peut etre maintenant.
      dequeue_incidence_request(v2);
      if (ordered_el1.first == v1)
        dequeue_incidence_request(ordered_el1.second);
      else
        dequeue_incidence_request(ordered_el1.first);
    }

    //---------------------------------------------------------------------
    void
    border_extend(const Edge_like& ordered_key, const Border_elt& result12,
                  const Vertex_handle& v1, const Vertex_handle& v2,
                  const Vertex_handle& v3,
                  const Radius_edge_type& e1, const Radius_edge_type& e2,
                  IO_edge_type* &p1, IO_edge_type* &p2)
    {
      remove_border_elt(ordered_key);

      //depiler v3 avant de le mettre a jour... pour reperer s'il est sur un bord
      if (v3->is_on_border())
        dequeue_incidence_request(v3);

      if (ordered_key.first == v1)
        {
          p1 = set_border_elt(v1, v3, Border_elt(e1,result12.second));
          p2 = set_border_elt(v3, v2, Border_elt(e2,result12.second));
        }
      else
        {
          p2 = set_border_elt(v2, v3, Border_elt(e2,result12.second));
          p1 = set_border_elt(v3, v1, Border_elt(e1,result12.second));
        }

      v3->inc_mark();

      //depiler les eventuelles requettes de connections avortees... zones etoilees,
      //en effet le bord a change donc on peut peut etre maintenant.
      dequeue_incidence_request(v1);
      dequeue_incidence_request(v2);
    }

    //=====================================================================
    Validation_case
    validate(const Edge_incident_facet& edge_Efacet,
             const criteria& value)
    {
      int i = (6 - edge_Efacet.second
               - edge_Efacet.first.second
               - edge_Efacet.first.third);
      Cell_handle c =  edge_Efacet.first.first;

      Vertex_handle v1 = c->vertex(edge_Efacet.first.second),
        v2 = c->vertex(edge_Efacet.first.third);

      Edge_like ordered_el1(c->vertex(i), v1);
      Edge_like ordered_el2(c->vertex(i), v2);
      Border_elt result1, result2, result12;

      Edge_like ordered_key(v1,v2);

      if (!is_border_elt(ordered_key, result12))
        std::cerr << "+++probleme coherence bord <validate>" << std::endl;

      bool is_border_el1 = is_border_elt(ordered_el1, result1),
        is_border_el2 = is_border_elt(ordered_el2, result2);

      Radius_edge_type e1, e2;
      if (c->vertex(i)->not_interior())
        {
          if ((!is_interior_edge(ordered_el1))&&
              (!is_interior_edge(ordered_el2)))
            {
              //toujours utile meme avec l'essai de try_to_close_border avant
              //validate pour la resolution de singularite par oreille qui elle
              //doit etre dans Delaunay.
              if (is_border_el1&&is_border_el2)
                {
                  remove_border_elt(ordered_key);
                  force_merge(ordered_el1, result1);
                  force_merge(ordered_el2, result2);

                  dec_mark(v1);
                  dec_mark(v2);
                  dec_mark(c->vertex(i));

                  select_facet(c, edge_Efacet.second);

                  return FINAL_CASE;
                }
              //---------------------------------------------------------------------
              //on peut alors marquer v1 et on pourrait essayer de merger
              //sans faire de calcul inutile???
              if (is_border_el1)
                {
                  Edge_incident_facet edge_Ifacet_2(Edge(c, i, edge_Efacet.first.third),
                                                    edge_Efacet.second);
                  merge_ear(ordered_el1, result1,
                            ordered_key, v1, v2, edge_Ifacet_2);

                  select_facet(c, edge_Efacet.second);

                  return EAR_CASE;
                }
              //---------------------------------------------------------------------
              //idem pour v2
              if (is_border_el2)
                {
                  Edge_incident_facet edge_Ifacet_1(Edge(c, i, edge_Efacet.first.second),
                                                    edge_Efacet.second);
                  merge_ear(ordered_el2, result2,
                            ordered_key, v2, v1, edge_Ifacet_1);
                  select_facet(c, edge_Efacet.second);

                  return EAR_CASE;
                }

              //---------------------------------------------------------------------
              if ((!is_border_el1)&&(!is_border_el2))
                {
                  // si on veut s'interdir de spliter un bord (pelure d'orange....)
                  // seulement c->vertex(i)->is_exterior()
                  // pour s'autoriser des split de bord surface a bord->sphere ou Moebius...
                  // alors || is_on_same_border:
                  //       if (c->vertex(i)->is_exterior() || is_on_same_border)
                  // pour passer au tore (changementde type de topologie)
                  // recoller deux bord different...
                  //       if (c->vertex(i)->not_interior() deja teste en haut

                  if(c->vertex(i)->is_exterior())
                    {
                      Edge_incident_facet edge_Ifacet_1(Edge(c, i, edge_Efacet.first.second),
                                                        edge_Efacet.second);

                      Edge_incident_facet edge_Ifacet_2(Edge(c, i, edge_Efacet.first.third),
                                                        edge_Efacet.second);
                      e1 = compute_value(edge_Ifacet_1);
                      e2 = compute_value(edge_Ifacet_2);

                      IO_edge_type* p1;
                      IO_edge_type* p2;

                      border_extend(ordered_key, result12,
                                    v1, v2, c->vertex(i),
                                    e1, e2, p1, p2);

                      // if e1 contain infinity there is no candidates to
                      // continue: compute_value is not valid...

                      _ordered_border.insert(Radius_ptr_type(e1.first, p1));

                      _ordered_border.insert(Radius_ptr_type(e2.first, p2));

                      select_facet(c, edge_Efacet.second);

                      return EXTERIOR_CASE;
                    }
                  else // c->vertex(i) is a border point (and now there's only 1
                  // border incident to a point... _mark<1 even if th orientation
                  // may be such as one vh has 2 successorson the same border...
                  {
                    // a ce niveau on peut tester si le recollement se fait en
                    // maintenant la compatibilite d'orientation des bords (pour
                    // surface orientable...) ou si elle est brisee...
                    Edge_incident_facet edge_Ifacet_1(Edge(c, i, edge_Efacet.first.second),
                                                      edge_Efacet.second);
                    Edge_incident_facet edge_Ifacet_2(Edge(c, i, edge_Efacet.first.third),
                                                      edge_Efacet.second);

                    e1 = compute_value(edge_Ifacet_1);
                    e2 = compute_value(edge_Ifacet_2);

                    if ((e1.first >= STANDBY_CANDIDATE)&&(e2.first >= STANDBY_CANDIDATE))
                      return NOT_VALID_CONNECTING_CASE;

                    // vu compute value: les candidats oreilles fournis sont sans
                    // aretes interieures et le sommet oppose n'est pas non plus interieur
                    Edge_incident_facet ear1 = e1.second.second;
                    Edge_incident_facet ear2 = e2.second.second;

                    int ear1_i = (6 - ear1.second
                                  - ear1.first.second
                                  - ear1.first.third);
                    Cell_handle ear1_c =  ear1.first.first;
                    Border_elt result_ear1;

                    int ear2_i = (6 - ear2.second
                                  - ear2.first.second
                                  - ear2.first.third);
                    Cell_handle ear2_c =  ear2.first.first;
                    Border_elt result_ear2;

                    Edge_like ear1_e, ear2_e;
                    // pour maintenir la reconstruction d'une surface orientable :
                    // on verifie que les bords se recollent dans des sens opposes
                    if (ordered_key.first==v1)
                      {
                        ear1_e = Edge_like(c->vertex(i), ear1_c ->vertex(ear1_i));
                        ear2_e = Edge_like(ear2_c ->vertex(ear2_i), c->vertex(i));
                      }
                    else
                      {
                        ear1_e = Edge_like(ear1_c ->vertex(ear1_i), c->vertex(i));
                        ear2_e = Edge_like(c->vertex(i), ear2_c ->vertex(ear2_i));
                      }

                    //maintient la surface orientable
                    bool is_border_ear1 = is_ordered_border_elt(ear1_e, result_ear1);
                    bool is_border_ear2 = is_ordered_border_elt(ear2_e, result_ear2);
                    bool ear1_valid(false), ear2_valid(false);
                    if (is_border_ear1&&(e1.first < STANDBY_CANDIDATE)&&
                        (e1.first <=  value)&&
                        (result12.second==result_ear1.second))
                      {
                        ear1_valid = test_merge(ear1_e, result_ear1, v1,
                                                priority(*this, ear1_c, ear1.second)) != 0;
                      }
                    if (is_border_ear2&&(e2.first < STANDBY_CANDIDATE)&&
                        (e2.first <= value)&&
                        (result12.second==result_ear2.second))
                      {
                        ear2_valid = test_merge(ear2_e, result_ear2, v2,
                                                priority(*this, ear2_c, ear2.second)) != 0;
                      }
                    if ((!ear1_valid)&&(!ear2_valid))
                      return NOT_VALID_CONNECTING_CASE;

                    IO_edge_type* p1;
                    IO_edge_type* p2;

                    border_extend(ordered_key, result12,
                                  v1, v2, c->vertex(i),
                                  e1, e2, p1, p2);

                    if (ear1_valid&&ear2_valid&&(ear1_e==ear2_e))
                      {
                        if (e1.first < e2.first)
                          {
                            Validation_case res = validate(ear1, e1.first);
                            if (!((res == EAR_CASE)||(res == FINAL_CASE)))
                              std::cerr << "+++probleme de recollement : cas "
                                        << res << std::endl;
                            e2 = compute_value(edge_Ifacet_2);

                            if (ordered_key.first == v1)
                              p2 = set_again_border_elt(c->vertex(i), v2,
                                                        Border_elt(e2, result2.second));
                            else
                              p2 = set_again_border_elt(v2, c->vertex(i),
                                                        Border_elt(e2, result2.second));

                            _ordered_border.insert(Radius_ptr_type(e2.first, p2));
                          }
                        else
                          {
                            Validation_case res = validate(ear2, e2.first);
                            if (!((res == EAR_CASE)||(res == FINAL_CASE)))
                              std::cerr << "+++probleme de recollement : cas "
                                        << res << std::endl;
                            e1 = compute_value(edge_Ifacet_1);

                            if (ordered_key.first == v1)
                              p1 = set_again_border_elt(v1, c->vertex(i),
                                                        Border_elt(e1, result1.second));
                            else
                              p1 = set_again_border_elt(c->vertex(i), v1,
                                                        Border_elt(e1, result1.second));

                            _ordered_border.insert(Radius_ptr_type(e1.first, p1));
                          }
                      }
                    else// les deux oreilles ne se recollent pas sur la meme arete...
                      {
                        // on resoud la singularite.
                        if (ear1_valid)
                          {
                            Validation_case res = validate(ear1, e1.first);
                            if (!((res == EAR_CASE)||(res == FINAL_CASE)))
                              std::cerr << "+++probleme de recollement : cas "
                                        << res << std::endl;
                          }
                        if (ear2_valid)
                          {
                            Validation_case res = validate(ear2, e2.first);
                            if (!((res == EAR_CASE)||(res == FINAL_CASE)))
                              std::cerr << "+++probleme de recollement : cas "
                                        << res << std::endl;
                          }
                        // on met a jour la PQ s'il y a lieu... mais surtout pas
                        // avant la resolution de la singularite
                        if (!ear1_valid)
                          {
                            _ordered_border.insert(Radius_ptr_type(e1.first, p1));
                          }
                        if (!ear2_valid)
                          {
                            _ordered_border.insert(Radius_ptr_type(e2.first, p2));
                          }
                      }
                    select_facet(c, edge_Efacet.second);
                    return CONNECTING_CASE;
                  }
                }
            }
        }
      return NOT_VALID;
    }

    //=====================================================================
    void
    re_compute_values()
    {
      if(!_ordered_border.empty())
        {
          Ordered_border_type _ordered_border_tmp;
          do
            {
              Ordered_border_iterator e_it = _ordered_border.begin();
              Edge_incident_facet mem_Ifacet =  e_it->second->first;
              Cell_handle c_tmp = mem_Ifacet.first.first;
              _ordered_border.erase(e_it);
              Vertex_handle v1 = c_tmp->vertex(mem_Ifacet.first.second);
              Vertex_handle v2 = c_tmp->vertex(mem_Ifacet.first.third);

              Radius_edge_type new_candidate;
              new_candidate = compute_value(mem_Ifacet);

              if (new_candidate.first == STANDBY_CANDIDATE)
                {
                  // a garder pour un K un peu plus grand...
                  new_candidate.first = STANDBY_CANDIDATE_BIS;
                }

              Border_elt result;
              Edge_like key_tmp(v1,v2);
              is_border_elt(key_tmp, result);
              IO_edge_type* pnew =
                set_again_border_elt(key_tmp.first, key_tmp.second,
                                     Border_elt (new_candidate, result.second));
              _ordered_border_tmp.insert(Radius_ptr_type(new_candidate.first, pnew));
            }
          while(!_ordered_border.empty());

          _ordered_border.swap(_ordered_border_tmp);
        }
    }

    //---------------------------------------------------------------------
    void
    extend()
    {
      // initilisation de la variable globale K: qualite d'echantillonnage requise
      K = K_init; // valeur d'initialisation de K pour commencer prudemment...
      coord_type K_prev = K;

      Vertex_handle v1, v2;
      if (_ordered_border.empty()){
        return;
      }
      do
        {
          min_K = infinity(); // pour retenir le prochain K necessaire pour progresser...
          do
            {

              Ordered_border_iterator e_it = _ordered_border.begin();

              criteria value = e_it->first;
              if (value >= STANDBY_CANDIDATE)
                re_compute_values();
              else
                {
                  Edge_incident_facet candidate = e_it->second->second;
                  Cell_handle c_ext =  candidate.first.first;
                  int i1, i2 , i3;
                  i1 = candidate.first.second;
                  i2 = candidate.first.third;
                  i3 = (6 - i1- i2 - candidate.second);

                  Edge_incident_facet mem_Ifacet =  e_it->second->first;
                  Cell_handle c_tmp =  mem_Ifacet.first.first;

                  v1 = c_tmp->vertex(mem_Ifacet.first.second);
                  v2 = c_tmp->vertex(mem_Ifacet.first.third);

                  Radius_edge_type mem_e_it(e_it->first, *e_it->second);

                  _ordered_border.erase(e_it);
                  Validation_case validate_result = validate(candidate, value);
                  if ((validate_result == NOT_VALID)||
                      (validate_result == NOT_VALID_CONNECTING_CASE))
                    {
                      Radius_edge_type new_candidate;
                      Border_elt result;
                      Edge_like key_tmp(v1,v2);
                      is_border_elt(key_tmp, result);

                      if (validate_result == NOT_VALID_CONNECTING_CASE)
                        set_incidence_request(c_ext->vertex(i3), value, key_tmp);

                      if (validate_result == NOT_VALID)
                        {
                          new_candidate = compute_value(mem_Ifacet);
                          if ((new_candidate != mem_e_it))
                            //                               &&(new_candidate.first < NOT_VALID_CANDIDATE))
                            {
                              IO_edge_type* pnew =
                                set_again_border_elt(key_tmp.first, key_tmp.second,
                                                     Border_elt (new_candidate, result.second));
                              _ordered_border.insert(Radius_ptr_type(new_candidate.first,
                                                                     pnew));
                            }
                        }
                    }
                }
            }
          while((!_ordered_border.empty())&&
                (_ordered_border.begin()->first < STANDBY_CANDIDATE_BIS));
          K_prev = K;
          K += (std::max)(K_step, min_K - K + eps);
          // on augmente progressivement le K mais on a deja rempli sans
          // faire des betises auparavant...
        }
      while((!_ordered_border.empty())&&(K <= K)&&(min_K != infinity())&&(K!=K_prev));

#ifdef VERBOSE
      if ((min_K < infinity())&&(!_ordered_border.empty())) {
        std::cout << "   [ next K required = " << min_K << " ]" << std::endl;
      }
#endif // VERBOSE
    }


    //---------------------------------------------------------------------
    // En principe, si l'allocateur de cellules etait bien fait on aurait pas besoin
    // de mettre a jour les valeurs rajoutees pour les cellules a  la main...

    void
    re_init_for_free_cells_cache(const Vertex_handle& vh)
    {
      std::list<Cell_handle> ch_set;
      T.incident_cells(vh, std::back_inserter(ch_set));
      for (typename std::list<Cell_handle>::iterator c_it = ch_set.begin();
           c_it != ch_set.end();
           c_it++)
        (*c_it)->clear();
    }

    //---------------------------------------------------------------------
    void
    swap_selected_facets_on_conflict_boundary(const Vertex_handle& vh)
    {
      std::list<Cell_handle> ch_set;
      T.incident_cells(vh, std::back_inserter(ch_set));
      for (typename std::list<Cell_handle>::iterator c_it = ch_set.begin();
           c_it != ch_set.end(); c_it++)
        {
          Cell_handle c = *c_it;
          int index = c->index(vh);
          Cell_handle neigh = c->neighbor(index);
          int n_ind = neigh->index(c);
          neigh->set_smallest_radius(n_ind, -1); // pour obliger le recalcul
          // si c est selectionnee c'est qu'elle est aussi le mem_IFacet renvoye par
          // compute_value... donc a swapper aussi
          if (c->is_selected_facet(index))
            {
              int fn = c->facet_number(index);
              unselect_facet(c, index);
              neigh->select_facet(n_ind);
              neigh->set_facet_number(n_ind, fn);
              int i1 = (n_ind+1) & 3;
              int i2 = (n_ind+2) & 3;
              int i3 = (n_ind+3) & 3;
              Edge_like key(neigh->vertex(i1), neigh->vertex(i2));

              if (is_border_elt(key))
                {
                  Edge_incident_facet ei_facet(Edge(neigh, i1, i2),
                                               n_ind);
                  *border_IO_elt(key.first, key.second) =
                    IO_edge_type(ei_facet, ei_facet);
                }
              key = Edge_like(neigh->vertex(i1), neigh->vertex(i3));
              if (is_border_elt(key))
                {
                  Edge_incident_facet ei_facet(Edge(neigh, i1, i3),
                                               n_ind);
                  *border_IO_elt(key.first, key.second) =
                    IO_edge_type(ei_facet, ei_facet);
                }
              key = Edge_like(neigh->vertex(i3), neigh->vertex(i2));
              if (is_border_elt(key))
                {
                  Edge_incident_facet ei_facet(Edge(neigh, i3, i2),
                                               n_ind);
                  *border_IO_elt(key.first, key.second) =
                    IO_edge_type(ei_facet, ei_facet);
                }
            }
        }
    }

    //---------------------------------------------------------------------
    Facet
    next_surface_facet(const Edge_incident_facet& start)
    {
      Edge_incident_facet circ = next(start);
      Cell_handle c =  start.first.first;
      do
        {
          Cell_handle ch =  circ.first.first;
          int ind = circ.second;
          Cell_handle neigh = ch->neighbor(ind);
          int n_ind = neigh->index(ch);
          if (ch->is_selected_facet(ind)){
            return Facet(ch, ind);
          }
          if (neigh->is_selected_facet(n_ind)){
            return Facet(neigh, n_ind);
          }
          circ = next(circ);
        }
      while(circ.first.first != c);
      // si on passe par la, alors y a eu un probleme....
      std::cerr << "+++probleme dans la MAJ avant remove..." << std::endl;
      return Facet(c, start.second);
    }

    //---------------------------------------------------------------------
    void
    retract_border_for_incident_facets(const Vertex_handle& vh)
    {
      Next_border_elt border_elt =  *(vh->first_incident());
      int border_index = border_elt.second.second;
      Vertex_handle vh_succ = border_elt.first;
      IO_edge_type io_edge = border_elt.second.first.second;
      Edge_incident_facet i_facet = io_edge.first;
      Cell_handle c = i_facet.first.first;
      int i1 = c->index(vh);
      int i2 = c->index(vh_succ);
      int index = i_facet.second;
      int i3 = 6 - index - i1 - i2;
      Vertex_handle vh_int = c->vertex(i3);
      ordered_map_erase(border_elt.second.first.first,
                        border_IO_elt(vh, vh_succ));
      remove_border_edge(vh, vh_succ);
      // 1- a virer au cas ou car vh va etre detruit
      remove_interior_edge(vh_succ, vh);
      bool while_cond(true);
      do
        {
          _facet_number--;

          CGAL_assertion(c->is_selected_facet(index));
          unselect_facet(c, index);

          Facet f32 =
            next_surface_facet(Edge_incident_facet(Edge(c, i3, i2),
                                                   index));

          if (!vh_int->is_on_border())
            {
              re_init(vh_int);
              vh_int->inc_mark();
            }

          Edge_incident_facet e32(Edge(f32.first,
                                       f32.first->index(vh_int),
                                       f32.first->index(vh_succ)), f32.second);
          Radius_edge_type rad_elt_32(STANDBY_CANDIDATE, IO_edge_type(e32, e32));
          Border_elt result;
          if (is_ordered_border_elt(Edge_like(vh_int, vh), result))
            {
              ordered_map_erase(result.first.first, border_IO_elt(vh_int, vh));
              remove_border_edge(vh_int, vh);
              // 1- a virer au cas ou car vh va etre detruit
              remove_interior_edge(vh_int, vh);
              while_cond = false;
            }
          // a titre  preventif... on essaye de s'assurer de marquer les aretes
          // interieures au sens large...

          // 2- a virer a tout pris pour que maintenir le sens de interior edge
          remove_interior_edge(vh_int, vh_succ);
          remove_interior_edge(vh_succ, vh_int);

          IO_edge_type* p32 = set_border_elt(vh_int, vh_succ,
                                             Border_elt(rad_elt_32, border_index));
          _ordered_border.insert(Radius_ptr_type (STANDBY_CANDIDATE, p32));

          // incrementation...
          if (while_cond)
            {
              Facet f31 =
                next_surface_facet(Edge_incident_facet(Edge(c, i3, i1),
                                                       index));

              c = f31.first;
              index = f31.second;
              i1 = c->index(vh);
              vh_succ = vh_int;
              i2 = c->index(vh_int);
              i3 = 6 - index - i1 - i2;
              vh_int = c->vertex(i3);
            }
        }
      while(while_cond);
    }

    //---------------------------------------------------------------------
    bool
    create_singularity(const Vertex_handle& vh)
    {
      // Pour reperer le cas de triangle isole
      if (vh->is_on_border())
        {
          // vh sommet 0
          Next_border_elt border_elt =  *(vh->first_incident());
          Vertex_handle vh_1 = border_elt.first;// sommet 1
          border_elt =  *(vh_1->first_incident());
          Vertex_handle vh_2 = border_elt.first;// sommet 2
          border_elt =  *(vh_2->first_incident());
          Vertex_handle vh_3 = border_elt.first;// sommet 0 ???
          Cell_handle c;
          int i, j, k;
          if ((vh_3 == vh)&&(T.is_facet(vh, vh_1, vh_2, c, i ,j ,k)))
            {
              int l = 6-i-j-k;
              Cell_handle neigh = c->neighbor(l);

              if
                (c->is_selected_facet(l)||neigh->is_selected_facet(neigh->index(c)))
                return true;
            }
        }


      // Reperer le cas d'aretes interieures...
      std::list<Vertex_handle> vh_list;
      T.incident_vertices(vh, std::back_inserter(vh_list));

      for (typename std::list<Vertex_handle>::iterator v_it = vh_list.begin();
           v_it != vh_list.end(); v_it++)
        if ((*v_it)->is_on_border() && is_interior_edge(Edge_like(vh, *v_it)))
          return true;
      return false;
    }


    //---------------------------------------------------------------------
    void
    store_outlier(const Point& p){
      m_outliers.push_back(p);
    }


    void dec_vh_number()
    {
      _vh_number--;
    }


    struct Remove : public CGAL::cpp98::unary_function<Vertex_handle, bool>
    {

      Extract& E;
      Triangulation_3& T;

      Remove(Extract& E_, Triangulation_3& T_) : E(E_), T(T_) {}

      bool operator()(Vertex_handle vh) {
        if (vh->is_exterior())
          {
            E.swap_selected_facets_on_conflict_boundary(vh);
            E.re_init_for_free_cells_cache(vh);
            Point p = vh->point();
            T.remove(vh);
            E.dec_vh_number();
            E.store_outlier(p);

            return true;
          }
        else if (vh->is_on_border()&&(!E.create_singularity(vh)))
          {
            E.swap_selected_facets_on_conflict_boundary(vh);
            E.retract_border_for_incident_facets(vh);
            E.re_init_for_free_cells_cache(vh);
            Point p = vh->point();
            T.remove(vh);
            E.dec_vh_number();
            E.store_outlier(p);

            return true;
          }
        else
          { }
        return false;
      }
    };


    //---------------------------------------------------------------------
    bool
    postprocessing()
    {
      postprocess_timer.start();

      _postprocessing_counter++;

      std::list<Vertex_handle> L_v;

      //  Pour controler les sommets choisis sur le bord...

      // nombre d'aretes a partir duquel on considere que c'est irrecuperable NB_BORDER_MAX

      int vh_on_border_inserted(0);
      for(Finite_vertices_iterator v_it = T.finite_vertices_begin();
          v_it != T.finite_vertices_end();
          v_it++)
        {
          erase_incidence_request(v_it);
          if ((v_it->is_on_border())&&
              (!v_it->is_post_marked(_postprocessing_counter)))
            {
              std::list<Vertex_handle> L_v_tmp;
              Vertex_handle vprev_it(v_it), done(vprev_it), vh_it;
              //           Vertex_handle vsucc_it;
              int v_count(0);
              // collect all vertices on the border
              do
                {
                  vh_it =  vprev_it->first_incident()->first;
                  L_v_tmp.push_back(vh_it);
                  vh_it->set_post_mark(_postprocessing_counter);
                  vprev_it = vh_it;
                  v_count++;
                }
              while((vprev_it != done)&&(v_count < NB_BORDER_MAX));
              // we stopped either because we did a complete tour, or because
              // the border was so long that we consider it as too big to close
              // e.g., if it is a terrain with only one real border at the exterior
              if (v_count < NB_BORDER_MAX)
                {
                  L_v.insert(L_v.begin(), L_v_tmp.begin(), L_v_tmp.end());
                  vh_on_border_inserted += v_count;
                }

            }
          if (v_it->is_exterior())
            L_v.push_back(v_it);
        }

      std::size_t itmp, L_v_size_mem;
      L_v_size_mem = L_v.size();
      if ((vh_on_border_inserted != 0)&& // pour ne post-traiter que les bords
          (L_v.size() < .1 * _size_before_postprocessing))
        {
          {
            do
              {
                itmp = L_v.size();
                typename std::list<Vertex_handle>::iterator new_end =
                  std::remove_if(L_v.begin(), L_v.end(), Remove(*this,T));
                L_v.erase(new_end, L_v.end());
              }
            while (!L_v.empty() && (L_v.size() < itmp));
          }
#ifdef VERBOSE
          if(L_v.size() > 0){
            std::cout << "   " << L_v.size() << " non regular points." << std::endl;
          }
#endif // VERBOSE
          re_compute_values();
        }
      else{
        postprocess_timer.stop();
        return false;
      }
      // we stop if we removed more than 10% of points or after 20 rounds
      if ((L_v_size_mem == L_v.size())||
          ((_size_before_postprocessing - T.number_of_vertices()) >
           .1 * _size_before_postprocessing)||
          (_postprocessing_counter > 20)){
        postprocess_timer.stop();
        return false;
      }

      min_K = infinity();
      // fin--
      //   if (_postprocessing_counter < 5)
      //     return true;
      postprocess_timer.stop();
      return true;
    }

  }; // class Advancing_front_surface_reconstruction

  namespace AFSR {

    template <typename T>
    struct Auto_count : public CGAL::cpp98::unary_function<const T&,std::pair<T,std::size_t> >{
      mutable std::size_t i;

      Auto_count()
        : i(0)
      {}

      std::pair<T,std::size_t> operator()(const T& p) const {
        return std::make_pair(p,i++);
      }
    };

    template <typename T, typename CC>
    struct Auto_count_cc : public CGAL::cpp98::unary_function<const T&,std::pair<T,std::size_t> >{
      mutable std::size_t i;
      CC cc;

      Auto_count_cc(CC cc)
        : i(0), cc(cc)
      {}

      template <typename T2>
      std::pair<T,std::size_t> operator()(const T2& p) const {
        return std::make_pair(cc(p),i++);
      }
    };
  }

  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstructionRef

  For a sequence of points computes a sequence of index triples
  describing the faces of the reconstructed surface.

  \tparam PointInputIterator must be an input iterator with 3D points as value type.  This point type must
  be convertible to `Exact_predicates_inexact_constructions_kernel::Point_3` with the `Cartesian_converter`.
  \tparam IndicesOutputIterator must be an output iterator to which
  `std::array<std::size_t, 3>` can be assigned.

  \param b iterator on the first point of the sequence
  \param e past the end iterator of the point sequence
  \param out output iterator
  \param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
         are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`.
         Described in Section \ref AFSR_Boundaries
  \param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates.
         Described in Section \ref AFSR_Selection

  */
  template <typename PointInputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  advancing_front_surface_reconstruction(PointInputIterator b,
                                         PointInputIterator e,
                                         IndicesOutputIterator out,
                                         double radius_ratio_bound = 5,
                                         double beta = 0.52 )
  {
    typedef Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
    typedef Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

    typedef Triangulation_data_structure_3<LVb,LCb> Tds;
    typedef Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

    typedef Advancing_front_surface_reconstruction<Triangulation_3> Reconstruction;
    typedef typename std::iterator_traits<PointInputIterator>::value_type InputPoint;
    typedef typename Kernel_traits<InputPoint>::Kernel InputKernel;
    typedef Cartesian_converter<InputKernel,Kernel> CC;
    typedef Kernel::Point_3 Point_3;

    CC cc=CC();
    Triangulation_3 dt( boost::make_transform_iterator(b, AFSR::Auto_count_cc<Point_3,CC>(cc)),
                        boost::make_transform_iterator(e, AFSR::Auto_count_cc<Point_3,CC>(cc) )  );

    Reconstruction R(dt);
    R.run(radius_ratio_bound, beta);
    write_triple_indices(out, R);
    return out;
  }

  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstructionRef

  For a sequence of points computes a sequence of index triples
  describing the faces of the reconstructed surface.

  \tparam PointInputIterator must be an input iterator with 3D points as value type.  This point type must
  be convertible to `Exact_predicates_inexact_constructions_kernel::Point_3` with the `Cartesian_converter`.
  \tparam IndicesOutputIterator must be an output iterator to which
  `std::array<std::size_t, 3>` can be assigned.
  \tparam Priority must be a functor with `double operator()(AdvancingFront,Cell_handle,int)` returning the
  priority of the facet `(Cell_handle,int)`.

  \param b iterator on the first point of the sequence
  \param e past the end iterator of the point sequence
  \param out output iterator
  \param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
         are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`.
         Described in Section \ref AFSR_Boundaries
  \param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates.
         Described in Section \ref AFSR_Selection
  \param priority allows the user to choose how candidate triangles are prioritized.

  */
  template <typename PointInputIterator, typename IndicesOutputIterator, typename Priority>
  IndicesOutputIterator
  advancing_front_surface_reconstruction(PointInputIterator b,
                                         PointInputIterator e,
                                         IndicesOutputIterator out,
                                         Priority priority,
                                         double radius_ratio_bound = 5,
                                         double beta = 0.52 )
  {
    typedef Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
    typedef Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

    typedef Triangulation_data_structure_3<LVb,LCb> Tds;
    typedef Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

    typedef Advancing_front_surface_reconstruction<Triangulation_3,Priority> Reconstruction;
    typedef typename std::iterator_traits<PointInputIterator>::value_type InputPoint;
    typedef typename Kernel_traits<InputPoint>::Kernel InputKernel;
    typedef Cartesian_converter<InputKernel,Kernel> CC;
    typedef Kernel::Point_3 Point_3;

    CC cc=CC();
    Triangulation_3 dt( boost::make_transform_iterator(b, AFSR::Auto_count_cc<Point_3,CC>(cc)),
                        boost::make_transform_iterator(e, AFSR::Auto_count_cc<Point_3,CC>(cc) )  );

    Reconstruction R(dt, priority);
    R.run(radius_ratio_bound, beta);
    write_triple_indices(out, R);
    return out;
  }


  template <typename PointInputIterator, typename Kernel, typename Items,  template < class T, class I, class A> class HDS, typename Alloc,typename Priority>
  void
  advancing_front_surface_reconstruction(PointInputIterator b,
                                         PointInputIterator e,
                                         Polyhedron_3<Kernel,Items,HDS,Alloc>& polyhedron,
                                         Priority priority,
                                         double radius_ratio_bound = 5,
                                         double beta = 0.52)
  {
    typedef Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
    typedef Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

    typedef Triangulation_data_structure_3<LVb,LCb> Tds;
    typedef Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

    typedef Advancing_front_surface_reconstruction<Triangulation_3,Priority> Reconstruction;
    typedef typename std::iterator_traits<PointInputIterator>::value_type InputPoint;
    typedef typename Kernel_traits<InputPoint>::Kernel InputKernel;
    typedef Cartesian_converter<InputKernel,Kernel> CC;
    typedef typename Kernel::Point_3 Point_3;

    CC cc=CC();
    Triangulation_3 dt( boost::make_transform_iterator(b, AFSR::Auto_count_cc<Point_3,CC>(cc)),
                        boost::make_transform_iterator(e, AFSR::Auto_count_cc<Point_3,CC>(cc) )  );

    Reconstruction R(dt, priority);
    R.run(radius_ratio_bound, beta);
    AFSR::construct_polyhedron(polyhedron, R);
  }


  template <typename PointInputIterator, typename Kernel, typename Items, template < class T, class I, class A> class HDS, typename Alloc>
  void
  advancing_front_surface_reconstruction(PointInputIterator b,
                                         PointInputIterator e,
                                         Polyhedron_3<Kernel,Items,HDS,Alloc>& polyhedron,
                                         double radius_ratio_bound = 5,
                                         double beta = 0.52)
  {
    typedef Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
    typedef Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

    typedef Triangulation_data_structure_3<LVb,LCb> Tds;
    typedef Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

    typedef Advancing_front_surface_reconstruction<Triangulation_3> Reconstruction;
    typedef typename std::iterator_traits<PointInputIterator>::value_type InputPoint;
    typedef typename Kernel_traits<InputPoint>::Kernel InputKernel;
    typedef Cartesian_converter<InputKernel,Kernel> CC;
    typedef typename Kernel::Point_3 Point_3;
    CC cc=CC();
    Triangulation_3 dt( boost::make_transform_iterator(b, AFSR::Auto_count_cc<Point_3,CC>(cc)),
                        boost::make_transform_iterator(e, AFSR::Auto_count_cc<Point_3,CC>(cc) )  );

    Reconstruction R(dt);
    R.run(radius_ratio_bound, beta);
    AFSR::construct_polyhedron(polyhedron, R);
  }



} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_H
