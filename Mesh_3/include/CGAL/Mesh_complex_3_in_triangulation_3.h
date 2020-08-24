// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb, Clement Jamin
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H
#define CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/Mesh_3/Mesh_complex_3_in_triangulation_3_fwd.h>
#include <CGAL/disable_warnings.h>
#include <CGAL/iterator.h>
#include <CGAL/Mesh_3/utilities.h>
#include <CGAL/Mesh_3/Mesh_complex_3_in_triangulation_3_base.h>
#include <CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h>
#include <CGAL/Time_stamper.h>

#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/mpl/if.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {


template <typename Tr,
          typename CornerIndex,
          typename CurveIndex>
class Mesh_complex_3_in_triangulation_3 :
  public Mesh_3::Mesh_complex_3_in_triangulation_3_base<
    Tr, typename Tr::Concurrency_tag>
  , public CGAL::Mesh_3::internal::Debug_messages_tools
{
public:
  typedef typename Tr::Concurrency_tag                   Concurrency_tag;

private:
  typedef Mesh_complex_3_in_triangulation_3<
    Tr,CornerIndex,CurveIndex>                                    Self;
  typedef Mesh_3::Mesh_complex_3_in_triangulation_3_base<
                                          Tr,Concurrency_tag>     Base;

public:
  typedef typename Base::size_type                        size_type;

  typedef typename Tr::Point                              Point;
  typedef typename Base::Edge                             Edge;
  typedef typename Base::Facet                            Facet;
  typedef typename Base::Vertex_handle                    Vertex_handle;
  typedef typename Base::Cell_handle                      Cell_handle;
  typedef CornerIndex                                     Corner_index;
  typedef CurveIndex                                      Curve_index;

  typedef CGAL::Hash_handles_with_or_without_timestamps   Hash_fct;

#ifndef CGAL_NO_DEPRECATED_CODE
  typedef CurveIndex Curve_segment_index;
#endif

  typedef typename Base::Triangulation                    Triangulation;
  typedef typename Base::Subdomain_index                  Subdomain_index;

  using Base::surface_patch_index;

private:
  // Type to store the edges:
  //  - a set of std::pair<Vertex_handle,Vertex_handle> (ordered at insertion)
  //  - which allows fast lookup from one Vertex_handle
  //  - each element of the set has an associated info (Curve_index) value
  typedef boost::bimaps::bimap<
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::set_of_relation<>,
    boost::bimaps::with_info<Curve_index> >           Edge_map;

  typedef typename Edge_map::value_type               Internal_edge;

  // Type to store the corners
  typedef boost::unordered_map<Vertex_handle,
                               Corner_index,
                               Hash_fct>              Corner_map;

  // Type to store far vertices
  typedef std::vector<Vertex_handle>                  Far_vertices_vec;

public:
  /**
   * Constructor
   */
  Mesh_complex_3_in_triangulation_3() = default;

  /**
   * Copy constructor
   */
  Mesh_complex_3_in_triangulation_3(const Self& rhs);

  /**
   * Move constructor
   */
  Mesh_complex_3_in_triangulation_3(Self&& rhs)
    : Base(std::move(rhs))
    , edges_(std::move(rhs.edges_))
    , corners_(std::move(rhs.corners_))
    , far_vertices_(std::move(rhs.far_vertices_))
  {}

  /**
   * Assignement operator, also serves as move-assignement
   */
  Self& operator=(Self rhs)
  {
    swap(rhs);
    return *this;
  }

  /**
   * Swaps this & rhs
   */
  void swap(Self& rhs)
  {
    Base::swap(rhs);
    edges_.swap(rhs.edges_);
    corners_.swap(rhs.corners_);
    far_vertices_.swap(rhs.far_vertices_);
  }

  /**
   * Clears data of c3t3
   */
  void clear()
  {
    Base::clear();
    edges_.clear();
    corners_.clear();
    far_vertices_.clear();
  }

  /// Import Base functions
  using Base::is_in_complex;
  using Base::add_to_complex;
  using Base::remove_from_complex;
  using Base::triangulation;
  using Base::set_surface_patch_index;



  /**
   * Add edge e to complex, with Curve_index index
   */
  void add_to_complex(const Edge& e,
                      const Curve_index& index)
  {
    add_to_complex(e.first->vertex(e.second),
                   e.first->vertex(e.third),
                   index);
  }

  /**
   * Add edge (v1,v2) to complex, with Curve_index index
   */
  void add_to_complex(const Vertex_handle& v1,
                      const Vertex_handle& v2,
                      const Curve_index& index)
  {
    add_to_complex(make_internal_edge(v1,v2), index);
  }

  /**
   * Mark vertex \c v as a corner of the complex
   */
  void add_to_complex(const Vertex_handle& v, const Corner_index& index)
  {
    v->set_dimension(0);
    corners_.insert(std::make_pair(v,index));
  }

  /**
   * Remove edge \c e from complex
   */
  void remove_from_complex(const Edge& e)
  {
    remove_from_complex(e.first->vertex(e.second), e.first->vertex(e.third));
  }

  /**
   * Remove edge (v1,v2) from complex
   */
  void remove_from_complex(const Vertex_handle& v1, const Vertex_handle& v2)
  {
    remove_from_complex(make_internal_edge(v1,v2));
  }

  /**
   * Remove vertex \c v from complex
   */
  void remove_from_complex(const Vertex_handle& v)
  {
    corners_.erase(v);
    v->set_dimension(-1);
  }

  std::size_t number_of_far_points() const
  {
    return far_vertices_.size();
  }

  void add_far_point(const Point &p)
  {
    far_vertices_.push_back(triangulation().insert(p));
  }

  void add_far_point(Vertex_handle vh)
  {
    far_vertices_.push_back(vh);
  }

  void remove_far_points()
  {
    Triangulation &tr = triangulation();
    //triangulation().remove(far_vertices_.begin(), far_vertices_.end());
    typename Far_vertices_vec::const_iterator it = far_vertices_.begin();
    typename Far_vertices_vec::const_iterator it_end = far_vertices_.end();
    for ( ; it != it_end ; ++it)
    {
      std::vector<Cell_handle> new_cells;
      new_cells.reserve(32);
      tr.remove_and_give_new_cells(*it, std::back_inserter(new_cells));

      typename std::vector<Cell_handle>::iterator nc_it = new_cells.begin();
      typename std::vector<Cell_handle>::iterator nc_it_end = new_cells.end();
      for ( ; nc_it != nc_it_end ; ++nc_it)
      {
        Cell_handle c = *nc_it;
        for (int i = 0 ; i < 4 ; ++i)
        {
          Facet mirror_facet = tr.mirror_facet(std::make_pair(c, i));
          if (is_in_complex(mirror_facet))
          {
            set_surface_patch_index(c, i,
                                    surface_patch_index(mirror_facet));
            c->set_facet_surface_center(i,
              mirror_facet.first->get_facet_surface_center(mirror_facet.second));
          }
        }
        /*int i_inf;
        if (c->has_vertex(tr.infinite_vertex(), i_inf))
        {
          Facet mirror_facet = tr.mirror_facet(std::make_pair(c, i_inf));
          if (is_in_complex(mirror_facet))
          {
            set_surface_patch_index(c, i_inf,
                                    surface_patch_index(mirror_facet));
          }
        }*/
      }
    }
    far_vertices_.clear();
  }

  /**
   * Returns the number of edges of c3t3
   */
  size_type number_of_edges_in_complex() const
  {
    return edges_.size();
  }
  size_type number_of_edges() const
  {
    return edges_.size();
  }

  /**
   * Returns the number of corners of c3t3
   */
  size_type number_of_vertices_in_complex() const
  {
    return corners_.size();
  }
  size_type number_of_corners() const
  {
    return corners_.size();
  }

  void rescan_after_load_of_triangulation();

  /**
   * Returns true if edge \c e is in complex
   */
  bool is_in_complex(const Edge& e) const
  {
    return is_in_complex(e.first->vertex(e.second), e.first->vertex(e.third));
  }

  /**
   * Returns true if edge (v1,v2) is in C3T3
   */
  bool is_in_complex(const Vertex_handle& v1, const Vertex_handle& v2) const
  {
    return is_in_complex(make_internal_edge(v1,v2));
  }

  /**
   * Returns true if \c v is a 0-dimensionnal feature in the c3t3
   */
  bool is_in_complex(const Vertex_handle& v) const
  {
    return (corners_.find(v) != corners_.end());
  }

  /**
   * Returns Curve_index of edge \c e
   */
  Curve_index curve_index(const Edge& e) const
  {
    return curve_index(e.first->vertex(e.second),
                       e.first->vertex(e.third));
  }

  Curve_index curve_index(const Vertex_handle& v1,
                          const Vertex_handle& v2) const
  {
    return curve_index(make_internal_edge(v1,v2));
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  CGAL_DEPRECATED
  Curve_index curve_segment_index(const Edge& e) const
  {
    return curve_index(e);
  }

  CGAL_DEPRECATED
  Curve_index curve_segment_index(const Vertex_handle& v1,
                                  const Vertex_handle& v2) const
  {
    return curve_index(v1, v2);
  }
#endif // CGAL_NO_DEPRECATED_CODE

  /**
   * Returns Corner_index of vertex \c v
   */
  Corner_index corner_index(const Vertex_handle& v) const
  {
    typename Corner_map::const_iterator it = corners_.find(v);
    if ( corners_.end() != it ) { return it->second; }
    return Corner_index();
  }

  /**
   * Outputs the outer boundary of the entire domain with facets oriented outward.
   */
  std::ostream& output_boundary_to_off(std::ostream& out) const
  {
    internal::output_boundary_of_c3t3_to_off(*this, 0, out, false);
    return out;
  }

  /**
   * Outputs the outer boundary of the selected subdomain with facets oriented outward.
   */
  std::ostream& output_boundary_to_off(std::ostream& out, Subdomain_index subdomain) const
  {
    output_boundary_of_c3t3_to_off(*this, subdomain, out);
    return out;
  }

  /**
   * Outputs the surface facets with a consistent orientation at the interface of two subdomains.
   */
  std::ostream& output_facets_in_complex_to_off(std::ostream& out) const
  {
    internal::output_facets_in_complex_to_off(*this, out);
    return out;
  }

  /**
   * Fills \c out with incident edges (1-dimensional features of \c v.
   * OutputIterator value type is std::pair<Vertex_handle,Curve_index>
   * \pre v->in_dimension() < 2
   */
  template <typename OutputIterator>
  OutputIterator
  adjacent_vertices_in_complex(const Vertex_handle& v, OutputIterator out) const;

  // -----------------------------------
  // Undocumented
  // -----------------------------------

  /**
   * Returns true if c3t3 is valid
   */
  bool is_valid(bool verbose = false) const;

  // -----------------------------------
  // Complex traversal
  // -----------------------------------
private:
  class Edge_iterator_not_in_complex
  {
    const Self& c3t3_;
    const Curve_index index_;
  public:
    Edge_iterator_not_in_complex(const Self& c3t3,
                                 const Curve_index& index = Curve_index())
    : c3t3_(c3t3)
    , index_(index) { }

    template <typename Iterator>
    bool operator()(Iterator it) const
    {
      if ( index_ == Curve_index() ) { return ! c3t3_.is_in_complex(*it); }
      else { return c3t3_.curve_index(*it) != index_;  }
    }
  };

  class Vertex_iterator_not_in_complex
  {
    const Self& c3t3_;
    const Corner_index index_;
  public:
    Vertex_iterator_not_in_complex(const Self& c3t3,
                                   const Corner_index& index = Corner_index())
    : c3t3_(c3t3)
    , index_(index) { }

    template <typename ItMap>
    bool operator()(const ItMap it) const
    {
      if ( index_ == Corner_index() ) { return false; }
      else { return it->second != index_;  }
    }
  };

  // Filtered iterator
  typedef Filter_iterator<
    typename Corner_map::const_iterator,
    Vertex_iterator_not_in_complex >            Vertex_map_filter_iterator;

  // Iterator type to get the first element of pair
  typedef boost::transform_iterator <
    Mesh_3::internal::First_of<typename Vertex_map_filter_iterator::value_type>,
    Vertex_map_filter_iterator >                Vertex_map_iterator_first;

  // Iterator type to remove a level of referencing
  class Vertex_map_iterator_first_dereference
    : public boost::iterator_adaptor <
        Vertex_map_iterator_first_dereference,
        Vertex_map_iterator_first,
        typename Vertex_map_iterator_first::value_type::value_type,
        boost::use_default,
        typename Vertex_map_iterator_first::value_type::reference >
  {
    typedef Vertex_map_iterator_first_dereference Self;
    typedef boost::iterator_adaptor <
        Vertex_map_iterator_first_dereference,
        Vertex_map_iterator_first,
        typename Vertex_map_iterator_first::value_type::value_type,
        boost::use_default,
        typename Vertex_map_iterator_first::value_type::reference > iterator_adaptor_;
  public:
    typedef typename  Vertex_map_iterator_first::reference  pointer;
    typedef typename iterator_adaptor_::reference           reference;

    Vertex_map_iterator_first_dereference() : Self::iterator_adaptor_() { }

    template < typename Iterator >
    Vertex_map_iterator_first_dereference(Iterator i)
      : Self::iterator_adaptor_(typename Self::iterator_adaptor_::base_type(i))
    { }

    pointer operator->() const { return *(this->base()); }
    reference operator*() const { return **(this->base()); }

    operator Vertex_handle() { return Vertex_handle(*(this->base())); }
  };

public:
  /// Iterator type to visit the edges of the 1D complex.
  typedef Filter_iterator<
    typename Triangulation::Finite_edges_iterator,
    Edge_iterator_not_in_complex >          Edges_in_complex_iterator;

  /// Returns a Facets_in_complex_iterator to the first facet of the 1D complex
  Edges_in_complex_iterator edges_in_complex_begin() const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this),
                                 this->triangulation().finite_edges_begin());
  }

  /// Returns a Facets_in_complex_iterator to the first facet of the 1D complex
  Edges_in_complex_iterator
  edges_in_complex_begin(const Curve_index& index) const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this,index),
                                 this->triangulation().finite_edges_begin());
  }

  /// Returns past-the-end iterator on facet of the 1D complex
  Edges_in_complex_iterator edges_in_complex_end(const Curve_index& = Curve_index()) const
  {
    return CGAL::filter_iterator(this->triangulation().finite_edges_end(),
                                 Edge_iterator_not_in_complex(*this));
  }

  /// Iterator type to visit the edges of the 0D complex.
  typedef Vertex_map_iterator_first_dereference Vertices_in_complex_iterator;

  /// Returns a Vertices_in_complex_iterator to the first vertex of the 0D complex
  Vertices_in_complex_iterator vertices_in_complex_begin() const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this),
                                 corners_.begin());
  }

  /// Returns a Vertices_in_complex_iterator to the first vertex of the 0D complex
  Vertices_in_complex_iterator
  vertices_in_complex_begin(const Corner_index& index) const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this,index),
                                 corners_.begin());
  }

  /// Returns past-the-end iterator on facet of the 0D complex
  Vertices_in_complex_iterator vertices_in_complex_end() const
  {
    return CGAL::filter_iterator(corners_.end(),
                                 Vertex_iterator_not_in_complex(*this));
  }


private:
  /**
   * Creates an Internal_edge object (i.e a pair of ordered Vertex_handle)
   */
  Internal_edge make_internal_edge(const Vertex_handle& v1,
                                   const Vertex_handle& v2) const
  {
    if ( v1 < v2 ) { return Internal_edge(v1,v2); }
    else { return Internal_edge(v2,v1); }
  }

  /**
   * Returns true if \c edge is in C3T3
   */
  bool is_in_complex(const Internal_edge& edge) const
  {
    return (curve_index(edge) != Curve_index() );
  }

  /**
   * Add edge \c edge to complex, with Curve_index index
   */
  void add_to_complex(const Internal_edge& edge, const Curve_index& index)
  {
    CGAL_precondition(!is_in_complex(edge));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Add edge ( " << disp_vert(edge.left)
              << " , " << disp_vert(edge.right) << " ), curve_index=" << index
              << " to c3t3.\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG
    std::pair<typename Edge_map::iterator, bool> it = edges_.insert(edge);
    it.first->info = index;
  }

  /**
   * Remove edge \c edge from complex
   */
  void remove_from_complex(const Internal_edge& edge)
  {
    edges_.erase(edge);
  }

  /**
   * Returns Curve_index of edge \c edge
   */
  Curve_index curve_index(const Internal_edge& edge) const
  {
    typename Edge_map::const_iterator it = edges_.find(edge);
    if ( edges_.end() != it ) { return it->info; }
    return Curve_index();
  }

private:
  Edge_map edges_;
  Corner_map corners_;
  Far_vertices_vec far_vertices_;
};


template <typename Tr, typename CI_, typename CSI_>
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
Mesh_complex_3_in_triangulation_3(const Self& rhs)
  : Base(rhs)
  , edges_()
  , corners_()
{
  // Copy edges
  for ( typename Edge_map::const_iterator it = rhs.edges_.begin(),
       end = rhs.edges_.end() ; it != end ; ++it )
  {
    const Vertex_handle& va = it->right;
    const Vertex_handle& vb = it->left;

    Vertex_handle new_va;
    this->triangulation().is_vertex(rhs.triangulation().point(va), new_va);

    Vertex_handle new_vb;
    this->triangulation().is_vertex(rhs.triangulation().point(vb), new_vb);

    this->add_to_complex(make_internal_edge(new_va,new_vb), it->info);
  }

  // Copy corners
  for ( typename Corner_map::const_iterator it = rhs.corners_.begin(),
       end = rhs.corners_.end() ; it != end ; ++it )
  {
    Vertex_handle new_v;
    this->triangulation().is_vertex(rhs.triangulation().point(it->first), new_v);
    this->add_to_complex(new_v, it->second);
  }

  // Parse vertices to identify far vertices
  if (rhs.far_vertices_.size() > 0)
  {
    Triangulation &tr = triangulation();
    typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
    for(typename Tr::Finite_vertices_iterator end = tr.finite_vertices_end();
        vit != end ; ++vit)
    {
      if (vit->in_dimension() == -1)
        far_vertices_.push_back(vit);
    }
    CGAL_assertion(far_vertices_.size() == rhs.far_vertices_.size());
  }
}


template <typename Tr, typename CI_, typename CSI_>
template <typename OutputIterator>
OutputIterator
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
adjacent_vertices_in_complex(const Vertex_handle& v, OutputIterator out) const
{
  CGAL_precondition(v->in_dimension() < 2);

  typedef typename Edge_map::right_const_iterator Rcit;
  typedef typename Edge_map::left_const_iterator Lcit;

  // Add edges containing v is on the left
  std::pair<Rcit,Rcit> range_right = edges_.right.equal_range(v);
  for ( Rcit rit = range_right.first ; rit != range_right.second ; ++rit )
  {
    *out++ = std::make_pair(rit->second, rit->info);
  }

  // Add edges containing v on the right
  std::pair<Lcit,Lcit> range_left = edges_.left.equal_range(v);
  for ( Lcit lit = range_left.first ; lit != range_left.second ; ++lit )
  {
    *out++ = std::make_pair(lit->second, lit->info);
  }

  return out;
}


template <typename Tr, typename CI_, typename CSI_>
bool
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
is_valid(bool verbose) const
{
  typedef typename Tr::Weighted_point                         Weighted_point;
  typedef boost::unordered_map<Vertex_handle, int, Hash_fct>  Vertex_map;

  Vertex_map vertex_map;

  // Fill map counting neighbor number for each vertex of an edge
  for ( typename Edge_map::const_iterator it = edges_.begin(),
       end = edges_.end() ; it != end ; ++it )
  {
    const Vertex_handle& v1 = it->right;
    if ( vertex_map.find(v1) == vertex_map.end() ) { vertex_map[v1] = 1; }
    else { vertex_map[v1] += 1; }

    const Vertex_handle& v2 = it->left;
    if ( vertex_map.find(v2) == vertex_map.end() ) { vertex_map[v2] = 1; }
    else { vertex_map[v2] += 1; }
  }

  // Verify that each vertex has 2 neighbors if it's not a corner
  for ( typename Vertex_map::iterator vit = vertex_map.begin(),
       vend = vertex_map.end() ; vit != vend ; ++vit )
  {
    if ( vit->first->in_dimension() != 0 && vit->second != 2 )
    {
      if(verbose)
        std::cerr << "Validity error: vertex " << (void*)(&*vit->first)
                  << " (" << this->triangulation().point(vit->first) << ") "
                  << "is not a corner (dimension " << vit->first->in_dimension()
                  << ") but has " << vit->second << " neighbor(s)!\n";
      return false;
    }
  }

  // Verify that balls of each edge intersect
  for ( typename Edge_map::const_iterator it = edges_.begin(),
       end = edges_.end() ; it != end ; ++it )
  {
    typename Tr::Geom_traits::Compute_weight_3 cw =
      this->triangulation().geom_traits().compute_weight_3_object();
    typename Tr::Geom_traits::Construct_point_3 cp =
      this->triangulation().geom_traits().construct_point_3_object();
    typename Tr::Geom_traits::Construct_sphere_3 sphere =
      this->triangulation().geom_traits().construct_sphere_3_object();
    typename Tr::Geom_traits::Do_intersect_3 do_intersect =
      this->triangulation().geom_traits().do_intersect_3_object();

    const Weighted_point& itrwp = this->triangulation().point(it->right);
    const Weighted_point& itlwp = this->triangulation().point(it->left);

    if ( ! do_intersect(sphere(cp(itrwp), cw(itrwp)), sphere(cp(itlwp), cw(itlwp))) )
    {
      std::cerr << "Points p[" << disp_vert(it->right) << "], dim=" << it->right->in_dimension()
                << " and q[" << disp_vert(it->left) << "], dim=" << it->left->in_dimension()
                << " form an edge but do not intersect !\n";
      return false;
    }
  }

  return true;
}

template <typename Tr, typename CI_, typename CSI_>
void
Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::
rescan_after_load_of_triangulation() {
  corners_.clear();
  for(typename Tr::Finite_vertices_iterator
        vit = this->triangulation().finite_vertices_begin(),
        end = this->triangulation().finite_vertices_end();
      vit != end; ++vit)
  {
    if ( vit->in_dimension() == 0 ) {
      add_to_complex(vit, Corner_index(1));
    }
  }
  Base::rescan_after_load_of_triangulation();
}

template <typename Tr, typename CI_, typename CSI_>
std::ostream &
operator<< (std::ostream& os,
            const Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_> &c3t3)
{
  // TODO: implement edge saving
  typedef typename Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::Concurrency_tag Concurrency_tag;
  return os << static_cast<
    const Mesh_3::Mesh_complex_3_in_triangulation_3_base<Tr, Concurrency_tag>&>(c3t3);
}


template <typename Tr, typename CI_, typename CSI_>
std::istream &
operator>> (std::istream& is,
            Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_> &c3t3)
{
  // TODO: implement edge loading
  typedef typename Mesh_complex_3_in_triangulation_3<Tr,CI_,CSI_>::Concurrency_tag Concurrency_tag;
  is >> static_cast<
    Mesh_3::Mesh_complex_3_in_triangulation_3_base<Tr, Concurrency_tag>&>(c3t3);
  c3t3.rescan_after_load_of_triangulation();
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_COMPLEX_3_IN_TRIANGULATION_3_H
