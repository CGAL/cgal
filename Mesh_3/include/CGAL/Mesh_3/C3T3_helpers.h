// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_C3T3_HELPERS_H
#define CGAL_MESH_3_C3T3_HELPERS_H

#include <CGAL/Mesh_3/config.h>
#include <CGAL/use.h>

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>

#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/type_traits/is_same.hpp>

#include <functional>
#include <vector>
#include <set>

namespace CGAL {
namespace Mesh_3 {

  
#ifdef CGAL_INTRUSIVE_LIST
template <typename Type>
class Intrusive_list {
public:
  
  typedef Type Type_handle;
  typedef Type_handle& reference;
  typedef const Type_handle& const_reference;
  typedef Type_handle value_type;

  Intrusive_list()
    : f(), b(), n(0)
  {}

  ~Intrusive_list()
  {
    clear();
  }


  Intrusive_list(const Intrusive_list& rhs)
  {
    CGAL_assertion(false);
  }

#ifdef CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR
  template <typename IT>
  Intrusive_list(IT first, IT last)
	: f(), b(), n(0)
  {
    if(first == last){
      return;
    }
   
    f = *first;
    Type_handle ch = f;
    ++n;
    ++first;
	while(first != last){
      if((ch != Type(*first)) && ((*first)->next_intrusive()==Type_handle())){
        // not yet inserted
        ch->set_next_intrusive(*first);
        (*first)->set_previous_intrusive(ch);
        ch = *first;
        ++n;
      }
      ++first;
    } 
    b = ch;
    b->set_next_intrusive(f);
    f->set_previous_intrusive(b);
  }
#endif

  bool
  is_valid() const
  {
    if(n < 0){
      std::cerr << "n < 0" << std::endl;
      return false;
    }
    if(n == 0){
      if (f != Type_handle()){
        std::cerr << "n==0, but f!= Type_handle()" << std::endl;
        return false;
      }
      if (b != Type_handle()){
        std::cerr << "n==0, but b!= Type_handle()" << std::endl;
        return false;
      }
    }else{
      if(f->previous_intrusive() != b){
        std::cerr << "f->previous_intrusive() != b" << std::endl;
        return false;
      }
      if(b->next_intrusive() != f){
        std::cerr << "b->next_intrusive() != f" << std::endl;
      return false;
      }


      Type_handle ch = f;
      for(std::size_t i = 1; i < n; i++){
        if(ch->next_intrusive()->previous_intrusive() != ch){
          std::cerr << "ch->next_intrusive()->previous_intrusive() != ch" << std::endl;
          return false;
        }
        ch = ch->next_intrusive();
      }
      if(ch != b){
        std::cerr << "ch!= b)" << std::endl;
        return false;
      }
    }
    return true;
  }

 
  void clear() 
  {
    if(!empty()){
      while( f!= b ){
        Type_handle h = f;
        f=f->next_intrusive();
        h->set_previous_intrusive(Type_handle());
        h->set_next_intrusive(Type_handle());
      }
      b->set_previous_intrusive(Type_handle());
      b->set_next_intrusive(Type_handle());
      f = b = Type_handle();
    }
    n = 0;
  }

  std::size_t size() const
  {
    return n;
  }


  struct iterator {
    Type_handle pos, b;

    typedef Type_handle                      value_type;
    typedef const Type_handle*                     pointer;
    typedef const Type_handle&                     reference;
    typedef std::size_t                      size_type;
    typedef std::ptrdiff_t                   difference_type;
    typedef std::forward_iterator_tag  iterator_category;

    iterator(Type_handle f, Type_handle b)
      : pos(f), b(b)
    {}

    iterator()
      : pos()
    {}

    iterator& operator++()
    {
      if(pos != Type_handle()){
        if(pos == b){
          pos = Type_handle(); // past the end
		    }else {
          pos = pos->next_intrusive();
		    }
      }
      return *this;
    }

    iterator operator++(int)
    {
      iterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const iterator& i) const
    { 
      return pos == i.pos;
    }
    
    bool operator!=(const iterator& i) const
    {
      return !(*this == i);
    }

    reference operator*() const
    {
      return pos;
    }

    pointer operator->() const
    {
      return pos;
    }
  }; // struct iterator


  iterator begin()
  {
    return iterator(f,b);
  }

  iterator end()
  {
    return iterator();
  }


  Type_handle front() const
  {
    return f;
  }

  Type_handle& front()
  {
    return f;
  }


  Type_handle back() const
  {
    return b;
  }

  Type_handle& back()
  {
    return b;
  }

  iterator insert(iterator /* position */,
                  const Type_handle& ch)
  {
    CGAL_assertion( (ch->next_intrusive() == Type_handle() && ch->previous_intrusive() == Type_handle()) || 
            (ch->next_intrusive() != Type_handle() && ch->previous_intrusive() != Type_handle()) );
    CGAL_expensive_assertion(is_valid());
    
    if(ch->next_intrusive() != Type_handle()){
      return iterator(ch->next_intrusive()/*first*/, ch/*last*/);
    }
    else{
      insert(ch);
      return iterator(ch->next_intrusive()/*first*/, ch/*last*/);
    }
  }

  void insert(Type_handle ch)
  {
    CGAL_assertion( (ch->next_intrusive() == Type_handle() && ch->previous_intrusive() == Type_handle()) || 
            (ch->next_intrusive() != Type_handle() && ch->previous_intrusive() != Type_handle()) );
    CGAL_expensive_assertion(is_valid());
    
    if(ch->next_intrusive() != Type_handle()){
      return;
    }
    if(empty()){
      f = b = ch;
      ch->set_next_intrusive(ch);
      ch->set_previous_intrusive(ch);
    } else {
      ch->set_next_intrusive(f);
      ch->set_previous_intrusive(b);
      f->set_previous_intrusive(ch);
      b->set_next_intrusive(ch);
      b = ch;
    }
    n++;
  }

  void erase(Type_handle ch)
  {
    CGAL_assertion( (ch->next_intrusive() == Type_handle() && ch->previous_intrusive() == Type_handle()) || 
            (ch->next_intrusive() != Type_handle() && ch->previous_intrusive() != Type_handle()) );
    CGAL_expensive_assertion(is_valid());
    if(ch->next_intrusive() == Type_handle()){
      return;
    }
    if(f == b){ // only 1 element in the list
      CGAL_assertion(f == ch);
      CGAL_assertion(n == 1);
      
      f = b = Type_handle();
    } else {
      if(f == ch){
        f = f->next_intrusive();
      }
      if(b == ch){
        b = b->previous_intrusive();
      }
      Type_handle p = ch->previous_intrusive(), n = ch->next_intrusive();
      p->set_next_intrusive(n);
      n->set_previous_intrusive(p);
    }
    ch->set_next_intrusive(Type_handle());
    ch->set_previous_intrusive(Type_handle());
    CGAL_assertion(ch->next_intrusive() == Type_handle());
    CGAL_assertion(ch->previous_intrusive() == Type_handle());
    n--;
  }

  bool empty() const
  {
    if(f == Type_handle()){
      CGAL_assertion(b == Type_handle());
      CGAL_assertion(n == 0);
    }
    return f == Type_handle();
  }

  bool contains(Type_handle th) const
  {
    if(th->next_intrusive() == Type_handle())
    {
      CGAL_assertion(th->previous_intrusive() == Type_handle());
      return true;
    }
    else return false;
  }
  
  void push_back(Type_handle ch)
  {
    insert(ch);
  }

private:
  Type_handle f,b; 
  std::size_t n;
};
#endif // #ifdef CGAL_INTRUSIVE_LIST

template <typename C3T3, typename MeshDomain>
class C3T3_helpers
{
  // -----------------------------------
  // Private types
  // -----------------------------------
  typedef typename C3T3::Triangulation  Tr;
  typedef Tr                            Triangulation;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Gt::Vector_3         Vector_3;
  typedef typename Gt::Point_3          Point_3;
  typedef typename Gt::Plane_3          Plane_3;
  typedef typename Gt::FT               FT;
  
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Cell             Cell;
  typedef typename Tr::Facet            Facet;
  
  typedef typename C3T3::Surface_patch_index  Surface_patch_index;
  typedef typename C3T3::Subdomain_index      Subdomain_index;
  typedef typename C3T3::Index                Index;
  
  typedef std::vector<Cell_handle>      Cell_vector;
  typedef std::set<Cell_handle>         Cell_set;

  typedef std::vector<Facet>            Facet_vector;
  typedef std::vector<Vertex_handle>    Vertex_vector;
  typedef std::set<Vertex_handle>       Vertex_set;
  
#ifdef CGAL_INTRUSIVE_LIST
  typedef Intrusive_list<Cell_handle>   Outdated_cell_set;
#else 
  typedef Cell_set  Outdated_cell_set;
#endif //CGAL_INTRUSIVE_LIST

#ifdef CGAL_INTRUSIVE_LIST
  typedef Intrusive_list<Vertex_handle>  Moving_vertices_set;
#else
  typedef Vertex_set Moving_vertices_set;
#endif //CGAL_INTRUSIVE_LIST

  
private:
  // Facet_boundary stores the boundary of surface facets
  typedef std::pair<Vertex_handle,Vertex_handle> Ordered_edge;
  typedef std::pair<Surface_patch_index, std::pair<int, Index> > Facet_topology_description;
  typedef std::map<Ordered_edge,Facet_topology_description>  Facet_boundary;
  
  typedef Triangulation_helpers<Tr> Th;
  
public:
  // -----------------------------------
  // Public interface
  // -----------------------------------
  typedef boost::optional<Vertex_handle> Update_mesh;
  
  /**
   * Constructor
   */
  C3T3_helpers(C3T3& c3t3, const MeshDomain& domain)
    : c3t3_(c3t3)
    , tr_(c3t3.triangulation())
    , domain_(domain) { }
 
  /**
   * @brief tries to move \c old_vertex to \c new_position in the mesh
   * @param new_position the new position of \c old_vertex
   * @param old_vertex the old vertex
   * @param criterion the criterion which will be used to verify the new 
   *    position is ok. c3t3 minimal value of new criterion shall not decrease.
   * @param modified_vertices contains the vertices incident to cells which 
   *    may have been impacted by relocation
   * @return a pair which contains:
   *    - a bool which is \c true if the move has been done.
   *    - a Vertex_handle which is always filled and may be the new vertex (if
   *      the move is a success), or the vertex which lies at \c v's position in
   *      the updated c3t3.
   */
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh(const Point_3& new_position,
              const Vertex_handle& old_vertex,
              const SliverCriterion& criterion,
              OutputIterator modified_vertices);

  /** @brief tries to move \c old_vertex to \c new_position in the mesh
   *
   * Same as update_mesh, but with the precondition that 
   * Th().no_topological_change(tr_, old_vertex, new_position,
   * incident_cells) return false.
   */
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh_topo_change(const Point_3& new_position,
                          const Vertex_handle& old_vertex,
                          const SliverCriterion& criterion,
                          OutputIterator modified_vertices);
  
  /**
   * Updates mesh moving vertex \c old_vertex to \c new_position. Returns the
   * new vertex of the triangulation.
   *
   * Insert into modified vertices the vertices which are impacted by to move.
   */
  template <typename OutputIterator>
  Vertex_handle update_mesh(const Point_3& new_position,
                            const Vertex_handle& old_vertex,
                            OutputIterator modified_vertices,
                            bool fill_modified_vertices = true);
  
  /**
   * Updates mesh moving vertex \c old_vertex to \c new_position. Returns the
   * new vertex of the triangulation.
   */
  Vertex_handle update_mesh(const Point_3& new_position,
                            const Vertex_handle& old_vertex)
  {
    return update_mesh(new_position, old_vertex, Emptyset_iterator(), false);
  }
  
  /**
   * Rebuilds restricted Delaunay
   */
  template <typename ForwardIterator>
  void rebuild_restricted_delaunay(ForwardIterator first_cell,
                                   ForwardIterator last_cell,
                                   Moving_vertices_set& moving_vertices);
  
#ifdef CGAL_INTRUSIVE_LIST
  template <typename OutdatedCells>
  void rebuild_restricted_delaunay(OutdatedCells& outdated_cells,
                                   Moving_vertices_set& moving_vertices);
#endif

  /**
   * @brief Project \c p on surface, using incident facets of \c v
   * @param p The point to project
   * @param v The vertex from which p was moved
   * @param index The index of the surface patch where v lies, if known.
   * @return the projected point
   *
   * \c p is projected as follows using normal of least square fitting plane
   * on \c v incident surface points. If \c index is specified, only
   * surface points that are on the same surface patch are used to compute
   * the fitting plane.
   */
  Point_3
  project_on_surface(const Point_3& p, const Vertex_handle& v, 
                     Surface_patch_index index = Surface_patch_index()) const;

  /**
   * Returns the minimum value for criterion for incident cells of \c vh 
   */
  template <typename SliverCriterion> 
  FT min_incident_value(const Vertex_handle& vh,
                        const SliverCriterion& criterion) const;
  
  /**
   * Moves \c old_vertex to \c new_position
   * Stores the cells which have to be updated in \c outdated_cells
   * Updates the Vertex_handle old_vertex to its new value in \c moving_vertices
   */
  Vertex_handle move_point(const Vertex_handle& old_vertex,
                           const Point_3& new_position,
                           Outdated_cell_set& outdated_cells_set,
                           Moving_vertices_set& moving_vertices);

  /**
   * Outputs to out the sliver (wrt \c criterion and \c sliver_bound) incident
   * to \c v
   */
  template <typename SliverCriterion, typename OutputIterator>
  OutputIterator
  incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound,
                   OutputIterator out) const;


  template <typename SliverCriterion, typename OutputIterator>
  OutputIterator
  new_incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound,
                   OutputIterator out) const;

  /**
   * Returns the sliver (wrt \c criterion and \c sliver_bound) incident to \c v
   */ 
  template <typename SliverCriterion> 
  Cell_vector
  incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound) const
  {
    Cell_vector slivers;
    incident_slivers(v, criterion, sliver_bound, std::back_inserter(slivers));
    return slivers;
  }

  template <typename SliverCriterion> 
  Cell_vector
  new_incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound) const
  {
    Cell_vector slivers;
    new_incident_slivers(v, criterion, sliver_bound, std::back_inserter(slivers));
    return slivers;
  }


  /**
   * Returns the number of slivers incident to \c v
   */
  template <typename SliverCriterion> 
  std::size_t
  number_of_incident_slivers(const Vertex_handle& v,
                             const SliverCriterion& criterion,
                             const FT& sliver_bound) const;
 
  template <typename SliverCriterion> 
  bool
  is_sliver(const Cell_handle& ch,
            const SliverCriterion& criterion,
            const FT& sliver_bound) const;

  /**
   * Returns the minimum criterion value of cells contained in \c cells
   * Precondition: cells of \c cells must not be infinite.
   * Warning: Here we don't check if cells are in c3t3
   */
  template <typename SliverCriterion>
  FT min_sliver_value(const Cell_vector& cells,
                      const SliverCriterion& criterion,
                      const bool use_cache = true) const;  
  
  /**
   * Reset cache validity of all cells of c3t3_
   */
  void reset_cache() const
  {
    namespace bl = boost::lambda;
    std::for_each(c3t3_.cells_in_complex_begin(),c3t3_.cells_in_complex_end(),
                  bl::bind(&Cell::reset_cache_validity, bl::_1) );
  }
  
private:  
  // -----------------------------------
  // Usefull Functors
  // -----------------------------------
  /**
   * @class Get_all_facets
   *
   * A functor which adds to an output iterator canonical facets of a cell
   */
  template <typename OutputIterator>
  class Get_all_facets
  {
  public:
    Get_all_facets(const Triangulation& tr, OutputIterator out)
      : tr_(tr)
      , out_(out) {}
    
    void operator()(const Cell_handle& cell)
    {
#ifndef CGAL_MESH_3_NEW_GET_FACETS
      for ( int i=0 ; i<4 ; ++i )
        if ( !tr_.is_infinite(cell,i) )
          *out_++ = canonical_facet(cell,i);
#else
      // Instead of iterating over the facets we iterate over the vertices
      // If a vertex is infinite we report only the facet opposite to it and return
      // If all vertices are finite we report all facets
      // This approach makes less tests if vertices are infinite
      int i=0;
      for ( ; i<4 ; ++i ){
        if ( tr_.is_infinite(cell->vertex(i)) ){         
          *out_++ = canonical_facet(cell,i);
          return;
        }
      }
      for ( i=0 ; i<4 ; ++i ){
        *out_++ = canonical_facet(cell,i);
      }
#endif
    }
    
  private:
    Facet canonical_facet(const Cell_handle& c, const int i) const
    {
#ifndef CGAL_MESH_3_NEW_GET_FACETS
      Facet facet(c,i);
      Facet mirror = tr_.mirror_facet(facet);
      return ( (mirror<facet)?mirror:facet );
#else
  	  Cell_handle n = c->neighbor(i);
      if(c < n){
        return Facet(c,i);
      }else{
        return Facet(n,n->index(c));
      }
#endif
    }
    
  private:
    const Triangulation& tr_;
    OutputIterator out_;
  };
  
  
  /**
   * @class Is_in_c3t3
   *
   * A functor which returns true if a given handle is in c3t3
   */
  template <typename Handle>
  class Is_in_c3t3 : public std::unary_function<Handle, bool>
  {
  public:
    Is_in_c3t3(const C3T3& c3t3) : c3t3_(c3t3) { }
    bool operator()(const Handle& h) const { return c3t3_.is_in_complex(h); }
    
  private:
    const C3T3& c3t3_;
  };
  

  /**
   * @class Is_sliver
   *
   * A functor which answers true if a Cell_handle is a sliver
   */
  template <typename SliverCriterion>
  struct Is_sliver : public std::unary_function<Cell_handle,bool>
  {
    Is_sliver(const C3T3& c3t3,
              const SliverCriterion& criterion,
              const FT& bound)
      : c3t3_(c3t3)
      , criterion_(criterion)
      , bound_(bound) { }
    
    bool operator()(const Cell_handle& c) const
    {
      if ( c3t3_.is_in_complex(c) )
      {
        CGAL_assertion(!c3t3_.triangulation().is_infinite(c));
        
        if ( ! c->is_cache_valid() )
        {
          FT sliver_value = criterion_(c3t3_.triangulation().tetrahedron(c));
          c->set_sliver_value(sliver_value);
        }
        return ( c->sliver_value() <= bound_ );
      }
      else
        return false;
    }
    
  private:
    const C3T3& c3t3_;
    const SliverCriterion& criterion_;
    const FT bound_;
  };
  
  
  /**
   * @class Update_c3t3
   *
   * A functor which updates c3t3 w.r.t the domain.
   */
  class Update_c3t3
  {
  public:
    Update_c3t3(const MeshDomain& domain, C3T3& c3t3)
      : domain_(domain)
      , c3t3_(c3t3) {}
    
    /**
     * @brief Updates facet \c facet in c3t3
     * @param facet the facet to update
     * @param update if set to \c false, checking only is done
     * @return true if \c facet is in c3t3
     */
    bool operator()(const Facet& facet, const bool update = true) const
    {
      typedef typename C3T3::Triangulation::Geom_traits Gt;
      typedef typename Gt::Segment_3 Segment_3;
      typedef typename Gt::Ray_3 Ray_3;
      typedef typename Gt::Line_3 Line_3;
      
      // Nothing to do for infinite facets
      if ( c3t3_.triangulation().is_infinite(facet) )
        return false;
      
      // Functors
      typename Gt::Is_degenerate_3 is_degenerate = 
        Gt().is_degenerate_3_object();
      
      // Get dual of facet
      Object dual = c3t3_.triangulation().dual(facet);

      // The dual is a segment, a ray or a line
      if ( const Segment_3* p_segment = object_cast<Segment_3>(&dual) )
      {
        if (is_degenerate(*p_segment)) 
          return false;
        
        return dual_intersect(*p_segment,facet,update);
      }
      else if ( const Ray_3* p_ray = object_cast<Ray_3>(&dual) )
      {
        if (is_degenerate(*p_ray))
          return false;
        
        return dual_intersect(*p_ray,facet,update);
      }
      else if ( const Line_3* p_line = object_cast<Line_3>(&dual) )
      {
        return dual_intersect(*p_line,facet,update);
      }
      
      // Should not happen
      CGAL_assertion(false);
      return false;
    }
    
    /**
     * @brief Updates cell \c ch in c3t3
     * @param ch the cell to update
     * @param update if set to \c false, checking only is done
     * @return true if \c ch is in c3t3
     */
    bool operator()(const Cell_handle& ch, const bool update = true) const
    {
      typedef boost::optional<typename MeshDomain::Subdomain_index> Subdomain;
      
      if ( c3t3_.triangulation().is_infinite(ch) )
        return false;
      
      // treat cell
      const Subdomain subdomain =
        domain_.is_in_domain_object()(c3t3_.triangulation().dual(ch));
      
      if ( subdomain && update )
      {
        c3t3_.add_to_complex(ch,*subdomain);
      }
      
      return subdomain;
    }
    
  private:
    
    // Returns true if query intersects the surface.
    template <typename Query>
    bool dual_intersect(const Query& dual,
                        const Facet& facet,
                        const bool update) const
    {
      typedef boost::optional<typename MeshDomain::Surface_patch_index> Surface_patch;
      typedef typename MeshDomain::Intersection Intersection;
      
      typename MeshDomain::Construct_intersection construct_intersection =
        domain_.construct_intersection_object();

#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3     
      
      typename MeshDomain::Do_intersect_surface do_intersect_surface =
        domain_.do_intersect_surface_object();
      Surface_patch surface = do_intersect_surface(dual);

#else // CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3

      Intersection intersection = construct_intersection(dual);
      Surface_patch surface =  
        (CGAL::cpp0x::get<2>(intersection) == 0) ? Surface_patch() : 
        domain_.surface_patch_index(CGAL::cpp0x::get<1>(intersection));

#endif // CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
     
      // Update if needed
      if ( surface && update )
      {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
        Intersection intersection = construct_intersection(dual);
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3

        // Update facet surface center
        Point_3 surface_center = CGAL::cpp0x::get<0>(intersection);
        facet.first->set_facet_surface_center(facet.second,surface_center);
        
        // Update status in c3t3 
        c3t3_.add_to_complex(facet,*surface);          
      }
      
      return surface;
    }
    
    
  private:
    const MeshDomain& domain_;
    C3T3& c3t3_;
  }; //end class Update_c3t3

  class Facet_updater {

    std::set<Vertex_handle>& vertex_to_proj;
    C3T3& c3t3_;
    Update_c3t3& c3t3_updater_;

  public:
    typedef Facet& reference;
    typedef const Facet& const_reference;

    Facet_updater(C3T3& c3t3,
      std::set<Vertex_handle>& vertex_to_proj,
      Update_c3t3& c3t3_updater_)
      : vertex_to_proj(vertex_to_proj), c3t3_(c3t3), c3t3_updater_(c3t3_updater_)
    {}

    void
      operator()(const Facet& f)
    {
      // Update facet
      c3t3_.remove_from_complex(f);
      c3t3_updater_(f);

      // Update vertex_to_proj
      if ( c3t3_.is_in_complex(f) )
      {
        // Iterate on vertices
        int k = f.second;
        for ( int i=1 ; i<4 ; ++i )
        {
          const Vertex_handle& v = f.first->vertex((k+i)&3);
          if ( c3t3_.in_dimension(v) > 2 )
          { 
            vertex_to_proj.insert(v);
          }
        }
      }
    }

  }; // end class Facet_updater
  
  
  /**
   * @class Sliver_criterion_value
   *
   * A functor which returns sliver criterion value for a Cell_handle 
   */
  template <typename SliverCriterion>
  class Sliver_criterion_value
    : public std::unary_function<Cell_handle, FT>
  {
  public:
    Sliver_criterion_value(const Tr& tr,
                           const SliverCriterion& criterion)
      : p_tr_(&tr)
      , criterion_(criterion) {}
    
    FT operator()(const Cell_handle& ch) const
    {
      CGAL_precondition(!p_tr_->is_infinite(ch));
      
      if ( ! ch->is_cache_valid() )
      {
        FT sliver_value = criterion_(p_tr_->tetrahedron(ch));
        ch->set_sliver_value(sliver_value);
      }
      return ch->sliver_value();
    }
    
  private:
    // '=' is used, so p_tr_ must be a pointer ...
    const Tr* p_tr_;
    SliverCriterion criterion_;
  };
  
private:
  // -----------------------------------
  // Private methods
  // -----------------------------------
  /**
   * Returns the minimum criterion value of c3t3 cells contained in \c cells.
   */
  template <typename SliverCriterion>
  FT min_sliver_in_c3t3_value(const Cell_vector& cells,
                              const SliverCriterion& criterion,
                              const bool use_cache = true) const
  {
    // Get complex cells only
    Cell_vector c3t3_cells;
    std::remove_copy_if(cells.begin(),
                        cells.end(),
                        std::back_inserter(c3t3_cells),
                        std::not1(Is_in_c3t3<Cell_handle>(c3t3_)) );
    
    return min_sliver_value(c3t3_cells,criterion,use_cache);
  }
  
  /**
   * Removes objects of [begin,end[ range from \c c3t3_
   */
  template<typename ForwardIterator>
  void remove_from_c3t3(ForwardIterator begin, ForwardIterator end)
  {
    while ( begin != end )
      c3t3_.remove_from_complex(*begin++);
  }
  
  /**
   * Remove cells and facets of \c cells from c3t3
   */
  template < typename ForwardIterator >
  void remove_cells_and_facets_from_c3t3(ForwardIterator cells_begin,
                                         ForwardIterator cells_end)
  {
    Facet_vector facets = get_facets_not_inplace(cells_begin,cells_end);
    remove_from_c3t3(facets.begin(), facets.end());
    remove_from_c3t3(cells_begin, cells_end);    
  }
  
  /**
   * Insert into \c out the vertices of range [cells_begin,cells_end[
   */
  template <typename InputIterator, typename OutputIterator>
  void fill_modified_vertices(InputIterator cells_begin,
                              InputIterator cells_end,
                              const Vertex_handle& vertex,
                              OutputIterator out) const;
  
  
  /**
   * Update mesh iff sliver criterion value does not decrease.
   */
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh_no_topo_change(const Point_3& new_position,
                             const Vertex_handle& old_vertex,
                             const SliverCriterion& criterion,
                             OutputIterator modified_vertices,
                             const Cell_vector& conflict_cells);
   
  /**
   * Move point and returns the set of cells that are not valid anymore, and
   * the set of cells which have been deleted by the move process.
   */
  template < typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle move_point(const Vertex_handle& old_vertex,
                           const Point_3& new_position,
                           OutdatedCellsOutputIterator outdated_cells,
                           DeletedCellsOutputIterator deleted_cells);

  Vertex_handle 
  move_point_topo_change(const Vertex_handle& old_vertex,
                         const Point_3& new_position,
                         Outdated_cell_set& outdated_cells_set);

  template < typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle
  move_point_topo_change(const Vertex_handle& old_vertex,
                         const Point_3& new_position,
                         OutdatedCellsOutputIterator outdated_cells,
                         DeletedCellsOutputIterator deleted_cells);
  
  template < typename ConflictCellsInputIterator,
             typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle 
  move_point_topo_change_conflict_zone_known(
     const Vertex_handle& old_vertex,
     const Point_3& new_position,
     ConflictCellsInputIterator conflict_cells_begin,
     ConflictCellsInputIterator conflict_cells_end,
     OutdatedCellsOutputIterator outdated_cells,
     DeletedCellsOutputIterator deleted_cells);

  Vertex_handle move_point_topo_change(const Vertex_handle& old_vertex,
                                       const Point_3& new_position);
  
  template < typename OutdatedCellsOutputIterator >
  Vertex_handle
  move_point_no_topo_change(const Vertex_handle& old_vertex,
                            const Point_3& new_position,
                            OutdatedCellsOutputIterator outdated_cells);

  Vertex_handle
  move_point_no_topo_change(const Vertex_handle& old_vertex,
                            const Point_3& new_position);
  
  /**
   * Returns the least square plane from v, using adjacent surface points
   */
  Plane_3 get_least_square_surface_plane(const Vertex_handle& v,
                                         Point_3& ref_point,
                                         Surface_patch_index index = Surface_patch_index()) const;
  
  /**
   * @brief Returns the projection of \c p, using direction of 
   * \c projection_vector
   */
  Point_3
  project_on_surface_aux(const Point_3& p,
                         const Point_3& ref_point,
                         const Vector_3& projection_vector) const;
  
  /**
   * Reverts the move from \c old_point to \c new_vertex. Returns the inserted
   * vertex located at \c old_point.
   */
  Vertex_handle revert_move(const Vertex_handle& new_vertex,
                            const Point_3& old_point)
  {
    Cell_set outdated_cells;
       
    // Move vertex
    Vertex_handle revert_vertex = 
      move_point_topo_change(new_vertex, 
                             old_point,
                             std::inserter(outdated_cells, outdated_cells.end()), 
                             CGAL::Emptyset_iterator()); //deleted cells
    CGAL_assertion(Vertex_handle() != revert_vertex);
    
    // Restore cell & facet attributes
    restore_mesh(outdated_cells.begin(), outdated_cells.end());
    
    return revert_vertex;
  }
  
  /**
   * Returns the boundary of facets of \c facets
   */
  Facet_boundary get_surface_boundary(const Facet_vector& facets) const;
  
  /**
   * Returns the boundary of facets of \c cells
   */
  Facet_boundary get_surface_boundary(const Cell_vector& cells) const
  {
    return get_surface_boundary(get_facets(cells));
  }
  
  /**
   * Returns false if there is a vertex belonging to one facet of \c facets 
   * which has not his dimension < 3
   */
  bool check_no_inside_vertices(const Facet_vector& facets) const;
  
  /**
   * Returns the impacted cells when moving \c vertex to \c conflict_point
   */
  template <typename OutputIterator>
  OutputIterator
  get_conflict_zone_no_topo_change(const Vertex_handle& vertex,
                                   OutputIterator conflict_cells) const;
  
  template <typename OutputIterator>
  OutputIterator
  get_conflict_zone_topo_change(const Vertex_handle& vertex,
                                const Point_3& conflict_point,
                                OutputIterator conflict_cells) const;
  
  template <typename CellsOutputIterator,   
            typename FacetsOutputIterator>
  void
  get_conflict_zone_topo_change(const Vertex_handle& v,
                                const Point_3& conflict_point,
                                CellsOutputIterator insertion_conflict_cells,
                                FacetsOutputIterator insertion_conflict_boundary,
                                CellsOutputIterator removal_conflict_cells) const;

  
  template < typename ConflictCellsInputIterator,
             typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle 
  move_point_topo_change_conflict_zone_known(const Vertex_handle& old_vertex,
                                             const Point_3& new_position,
                                             const Facet& insertion_boundary_facet,
                                             ConflictCellsInputIterator insertion_conflict_cells_begin,
                                             ConflictCellsInputIterator insertion_conflict_cells_end,
                                             ConflictCellsInputIterator removal_conflict_cells_begin,
                                             ConflictCellsInputIterator removal_conflict_cells_end,    
                                             OutdatedCellsOutputIterator outdated_cells,
                                             DeletedCellsOutputIterator deleted_cells);

  /**
   * Updates \c boundary wrt \c edge: if edge is already in boundary we remove
   * it, else we add it.
   */
  void update_boundary(Facet_boundary& boundary,
                       const Ordered_edge& edge,
                       const Vertex_handle third_vertex,
                       const Surface_patch_index& surface_index) const
  {
    const typename Facet_boundary::value_type x = 
      std::make_pair(edge,
                     std::make_pair(surface_index,
                                    std::make_pair(c3t3_.in_dimension(third_vertex),
                                                   c3t3_.index(third_vertex)
                                                   )
                                    )
                     );
    typename Facet_boundary::iterator boundary_it =
      boundary.find(edge);
    
    if ( boundary_it != boundary.end() )
      boundary.erase(boundary_it);
    else
      boundary.insert(x);
  }
  
  /**
   * Returns the facets of \c cells (returns each facet only once i.e. use
   * canonical facet)
   */
  Facet_vector get_facets(const Cell_vector& cells) const
  {
    return get_facets(cells.begin(),cells.end());
  }
  
  //  TODO: write get_facets so that it uses update_facets with a FacetUpdater that calls push_back

#if defined(CGAL_MESH_3_GET_FACETS_USING_INTRUSIVE_LIST) && defined(CGAL_INTRUSIVE_LIST)
  template <typename ForwardIterator>
  Facet_vector get_facets(ForwardIterator first_cell,
                          ForwardIterator last_cell) const
  {
    Facet_vector result; // AF: todo: resize?
#ifdef CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR
    Intrusive_list<Cell_handle> outdated_cells(first_cell, last_cell);
#else     
    Intrusive_list<Cell_handle> outdated_cells;
    for( ;first_cell!= last_cell; ++first_cell){
      outdated_cells.insert(*first_cell); 
    }
#endif
    
    for(typename Intrusive_list<Cell_handle>::iterator it = outdated_cells.begin();
        it != outdated_cells.end();
        ++it){
      Cell_handle cell = *it;
      int i=0;
      bool inf = false;
      for ( ; i<4 && (!inf) ; ++i ){
        if ( tr_.is_infinite(cell->vertex(i)) ){
          inf = true;
          Cell_handle n = cell->neighbor(i);
          if(n->next_intrusive() != Cell_handle()){// the neighbor is also outdated
            if(cell < n){ // otherwise n will report it later
              result.push_back(Facet(cell,i));
            }
          } else { // report it now or never
            if(cell < n){ 
              result.push_back(Facet(cell,i));
            }else {
              result.push_back(Facet(n,n->index(cell)));
            }
          }
        }
      }
      if(! inf){
        for ( i=0 ; i<4 ; ++i ){
          Cell_handle n = cell->neighbor(i);
          if(n->next_intrusive() != Cell_handle()){// the neighbor is also outdated
            if(cell < n){ // otherwise n will report it later
              result.push_back(Facet(cell,i));
            }
          } else { // report it now or never
            if(cell < n){ 
              result.push_back(Facet(cell,i));
            }else {
              result.push_back(Facet(n,n->index(cell)));
            }
          }
        }
      }
    }
    return result;
  }
#else
  /**
   * Returns the facets of \c cells (returns each facet only once i.e. use
   * canonical facet)
   */
  template <typename ForwardIterator>
  Facet_vector get_facets(ForwardIterator first_cell,
                          ForwardIterator last_cell) const
  {
    // Get all facets
    typedef Get_all_facets<std::back_insert_iterator<Facet_vector> > Get_facets;
    
    Facet_vector all_facets;
    all_facets.reserve(64);
    std::for_each(first_cell,
                  last_cell,
                  Get_facets(tr_,std::back_inserter(all_facets)));
    
    std::sort(all_facets.begin(), all_facets.end());
    
    // Keep one copy of each facet (maybe copy could be avoided)
    //    typename Facet_vector::iterator all_facets_end =
    //      std::unique(all_facets.begin(), all_facets.end());
    Facet_vector facets;
    facets.reserve(64);
    std::unique_copy(all_facets.begin(),
                     all_facets.end(),
                     std::back_inserter(facets));
    
    return facets;
  }
#endif

#ifdef CGAL_INTRUSIVE_LIST
  template <typename FacetUpdater>
  void update_facets(Intrusive_list<Cell_handle>& outdated_cells, FacetUpdater updater)
  {
    typename Intrusive_list<Cell_handle>::iterator it;
    for(it = outdated_cells.begin();
        it != outdated_cells.end();
        ++it)
    {
      Cell_handle cell = *it;

      int i=0;
      bool inf = false;
      for ( ; i<4 && (!inf) ; ++i ){
        if ( tr_.is_infinite(cell->vertex(i)) ){
          inf = true;
          Cell_handle n = cell->neighbor(i);
          if(n->next_intrusive() != Cell_handle()){// the neighbor is also outdated
            if(cell < n){ // otherwise n will report it later
              updater(Facet(cell,i));
            }
          } else { // report it now or never
            if(cell < n){ 
              updater(Facet(cell,i));
            }else {
              updater(Facet(n,n->index(cell)));
            }
          }
        }
      }
      if(! inf){
        for ( i=0 ; i<4 ; ++i ){
          Cell_handle n = cell->neighbor(i);
          if(n->next_intrusive() != Cell_handle()){// the neighbor is also outdated
            if(cell < n){ // otherwise n will report it later
              updater(Facet(cell,i));
            }
          } else { // report it now or never
            if(cell < n){ 
              updater(Facet(cell,i));
            }else {
              updater(Facet(n,n->index(cell)));
            }
          }
        }
      }
    }
  }
#endif //CGAL_INTRUSIVE_LIST


  /**
   * Returns the facets of \c cells (returns each facet only once i.e. use
   * canonical facet)
   */
  template <typename ForwardIterator>
  Facet_vector get_facets_not_inplace(ForwardIterator first_cell,
                                      ForwardIterator last_cell) const
  {
    typedef Get_all_facets<std::back_insert_iterator<Facet_vector> > Get_facets;
    
    Facet_vector all_facets;
    all_facets.reserve(64);
    std::for_each(first_cell,
                  last_cell,
                  Get_facets(tr_,std::back_inserter(all_facets)));
    
    std::sort(all_facets.begin(), all_facets.end());
    
    // Keep one copy of each facet (maybe copy could be avoided)
    //    typename Facet_vector::iterator all_facets_end =
    //      std::unique(all_facets.begin(), all_facets.end());
    Facet_vector facets;
    facets.reserve(64);
    std::unique_copy(all_facets.begin(),
                     all_facets.end(),
                     std::back_inserter(facets));
    CGAL_HISTOGRAM_PROFILER("|facets|", facets.size());
    return facets;
  }
 

  /**
   * Returns true if all surface facets of cells are really in restricted
   * Delaunay.
   */
  bool verify_surface(const Cell_vector& cells) const
  {
    // Naive implementation.
    // Todo: improve this (maybe we don't have to check if no facet is on surface)
    Facet_vector facets = get_facets(cells);
    Facet_vector surface_facets;
    
    // Check that nothing changed
    Update_c3t3 checker(domain_,c3t3_);
    for ( typename Facet_vector::iterator fit = facets.begin() ;
          fit != facets.end() ;
          ++fit )
    {
      if ( c3t3_.is_in_complex(*fit) )
      {
        surface_facets.push_back(*fit);
      }
      
      if ( c3t3_.is_in_complex(*fit) != checker(*fit,false) )
        return false;
    }
    
    // Facet surface center must be updated if verify_surface is ok
    std::for_each(surface_facets.begin(),surface_facets.end(),checker);
    
    return true;
  }
  
  
  /**
   * Restore mesh for cells and facets of \c cells, using domain_
   */ 
  void restore_mesh(const Cell_vector& cells)
  {
    restore_mesh(cells.begin(), cells.end());
  }
  
  /**
   * Restore mesh for cells and facets of \c cells, using domain_
   */ 
  template <typename ForwardIterator>
  void restore_mesh(ForwardIterator first_cell, ForwardIterator last_cell)
  {
    Facet_vector facets = get_facets(first_cell, last_cell);
    restore_mesh(first_cell, last_cell, facets.begin(), facets.end());
  }
  
  /**
   * Restore mesh for cells of \c cells and facets of \c facets, using domain_
   */
  template <typename CellForwardIterator, typename FacetForwardIterator>
  void restore_mesh(CellForwardIterator first_cell,
                    CellForwardIterator last_cell,
                    FacetForwardIterator first_facet,
                    FacetForwardIterator last_facet)
  {
    // Update mesh
    Update_c3t3 updater(domain_,c3t3_);
    std::for_each(first_facet, last_facet, updater);
    std::for_each(first_cell, last_cell, updater);
  }
  
  /**
   * Returns true if facets of \c facets have the same boundary as 
   * \c old_boundary
   */
  bool check_surface_mesh(const Facet_vector& facets,
                          const Facet_boundary& old_boundary) const
  {
    Facet_boundary new_boundary = get_surface_boundary(facets);
    return ( old_boundary.size() == new_boundary.size()
            && std::equal(new_boundary.begin(),
                          new_boundary.end(),
                          old_boundary.begin()) );
  }
  
  /**
   * Restore mesh for cells and facets of \c cells, then check that the new
   * boundary of facets of \c cells is the same as \c old_boundary.
   */
  bool restore_and_check_mesh(const Cell_vector& cells,
                              const Facet_boundary& old_boundary)
  {
    Facet_vector facets = get_facets(cells);
    restore_mesh(cells.begin(), cells.end(), facets.begin(), facets.end());
    return check_mesh(facets, old_boundary);
  }
  
  void set_facet_visited(const Facet& facet)
  {
    facet.first->set_facet_visited(facet.second);
    const Facet mirror_facet = tr_.mirror_facet(facet);
    mirror_facet.first->set_facet_visited(mirror_facet.second);    
  }
  
  /**
   * Orders handles \c h1, \c h2 & \c h3
   */ 
  template <typename Handle>
  void order_handles(Handle& h1, Handle& h2, Handle& h3) const
  {
    if ( h2 < h1 )
      std::swap(h1,h2);
    
    if ( h3 < h2 )
    {
      std::swap(h2,h3);
      
      if ( h2 < h1 ) // don't need to compare h2 & h1 if h2 didn't change
        std::swap(h1,h2);
    }
  }

  template < typename ForwardIterator >
  void reset_cache_validity(ForwardIterator cells_begin,
                            ForwardIterator cells_end) const
  {
    namespace bl = boost::lambda;
    std::for_each(cells_begin, cells_end,
                  bl::bind(&Cell::reset_cache_validity, *bl::_1) );
  }
  

private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const MeshDomain& domain_;
};
  
  
template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh(const Point_3& new_position,
            const Vertex_handle& old_vertex,
            const SliverCriterion& criterion,
            OutputIterator modified_vertices)
{
  // std::cerr << "\nupdate_mesh[v1](" << new_position << ",\n"
  //           << "                " << (void*)(&*old_vertex) << "=" << old_vertex->point()
  //           << ")\n";

  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(old_vertex, std::back_inserter(incident_cells));
  if ( Th().no_topological_change(tr_, old_vertex, new_position, incident_cells) )
  {
    BOOST_FOREACH(Cell_handle& ch, std::make_pair(incident_cells.begin(), 
                                                  incident_cells.end()))
    {
      ch->invalidate_circumcenter();
    }
    return update_mesh_no_topo_change(new_position,
                                      old_vertex,
                                      criterion,
                                      modified_vertices,
                                      incident_cells);
  }
  else
  {
    return update_mesh_topo_change(new_position,
                                   old_vertex,
                                   criterion,
                                   modified_vertices);
  }
}

template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh_no_topo_change(const Point_3& new_position,
                           const Vertex_handle& vertex,
                           const SliverCriterion& criterion,
                           OutputIterator modified_vertices,
                           const Cell_vector& conflict_cells )
{
  // std::cerr << "update_mesh_no_topo_change(\n"
  //           << new_position << ",\n"
  //           << "                " << (void*)(&*vertex) << "=" << vertex->point()
  //           << ")\n";

  // Get old values
  FT old_sliver_value = min_sliver_in_c3t3_value(conflict_cells, criterion);
  Point_3 old_position = vertex->point();
  
  // Move point
  move_point_no_topo_change(vertex,new_position);
    
  // Get new criterion value (conflict_zone did not change)
  const FT new_sliver_value = 
    min_sliver_in_c3t3_value(conflict_cells, criterion, false);
  
  // Check that mesh is still valid
  if ( new_sliver_value > old_sliver_value && verify_surface(conflict_cells) )
  {
    fill_modified_vertices(conflict_cells.begin(), conflict_cells.end(),
                           vertex, modified_vertices);
    return std::make_pair(true,vertex);
  }
  else
  {
    // std::cerr << "update_mesh_no_topo_change: revert move to "
    //           << old_position << "\n";
    // revert move
    move_point_no_topo_change(vertex,old_position);
    reset_cache_validity(conflict_cells.begin(), conflict_cells.end());
    return std::make_pair(false,vertex);
  }
}


template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh_topo_change(const Point_3& new_position,
                        const Vertex_handle& old_vertex,
                        const SliverCriterion& criterion,
                        OutputIterator modified_vertices)
{
  // check_c3t3(c3t3_);
  // std::cerr << "\n"
  //           << "update_mesh_topo_change("<< new_position << ",\n"
  //           << "                        " << (void*)(&*old_vertex) << "=" << old_vertex->point()
  //           << ")\n";
  Cell_set insertion_conflict_cells;
  Cell_set removal_conflict_cells;
  Facet_vector insertion_conflict_boundary;
  insertion_conflict_boundary.reserve(64);
  get_conflict_zone_topo_change(old_vertex, new_position,
                                std::inserter(insertion_conflict_cells,insertion_conflict_cells.end()),
                                std::back_inserter(insertion_conflict_boundary),
                                std::inserter(removal_conflict_cells, removal_conflict_cells.end()));

  if(insertion_conflict_boundary.empty())
    return std::make_pair(false,old_vertex); //new_location is a vertex already

  Cell_vector conflict_cells;
  conflict_cells.reserve(insertion_conflict_cells.size()+removal_conflict_cells.size());
  std::set_union(insertion_conflict_cells.begin(), insertion_conflict_cells.end(),
                 removal_conflict_cells.begin(), removal_conflict_cells.end(),
                 std::back_inserter(conflict_cells)); 
  
  FT old_sliver_value = min_sliver_in_c3t3_value(conflict_cells, criterion);
  Point_3 old_position = old_vertex->point();
  
  // Keep old boundary
  Facet_boundary old_surface_boundary = get_surface_boundary(conflict_cells);
  
  Cell_vector outdated_cells;
  outdated_cells.reserve(64);
    Vertex_handle new_vertex = 
        move_point_topo_change_conflict_zone_known(old_vertex, new_position,
                                                    insertion_conflict_boundary[0],
                                                    insertion_conflict_cells.begin(),
                                                    insertion_conflict_cells.end(),
                                                    removal_conflict_cells.begin(),
                                                    removal_conflict_cells.end(),
                                                    std::back_inserter(outdated_cells),
                                                    CGAL::Emptyset_iterator());
  // If nothing changed, return
  if ( old_position == new_vertex->point() ) 
  {
    // std::cerr << "update_mesh_topo_change: no move!\n";
    // check_c3t3(c3t3_);
    return std::make_pair(false,old_vertex);
  }
  
  restore_mesh(outdated_cells.begin(),outdated_cells.end());
  FT new_sliver_value = min_sliver_in_c3t3_value(outdated_cells, criterion);
  
  // Check that surface boundary does not change.
  // This check ensures that vertices which are inside c3t3 stay inside. 
  if ( new_sliver_value > old_sliver_value
      && check_surface_mesh(get_facets(outdated_cells), old_surface_boundary) )
  {
    fill_modified_vertices(outdated_cells.begin(), outdated_cells.end(),
                           new_vertex, modified_vertices);
    // check_c3t3(c3t3_);
    return std::make_pair(true,new_vertex);
  }
  else
  {
    // Remove from c3t3 cells which will be destroyed by revert_move
    remove_cells_and_facets_from_c3t3(outdated_cells.begin(),
                                      outdated_cells.end());
    
    // std::cerr << "update_mesh_topo_change: revert move to "
    //           << old_position << "\n";
    // Revert move
    Vertex_handle revert_vertex = revert_move(new_vertex, old_position);
    
    // check_c3t3(c3t3_);
    return std::make_pair(false,revert_vertex);
  }
}
  
template <typename C3T3, typename MD>
template <typename OutputIterator>  
typename C3T3_helpers<C3T3,MD>::Vertex_handle
C3T3_helpers<C3T3,MD>::  
update_mesh(const Point_3& new_position,
            const Vertex_handle& old_vertex,
            OutputIterator modified_vertices,
            bool fill_vertices)
{
  // std::cerr << "\nupdate_mesh[v2](" << new_position << ",\n"
  //           << "                " << (void*)(&*old_vertex) << "=" << old_vertex->point()
  //           << ")\n";
  Cell_vector outdated_cells;
  Vertex_handle new_vertex = move_point(old_vertex,
                                        new_position,
                                        std::back_inserter(outdated_cells),
                                        CGAL::Emptyset_iterator());
  
  restore_mesh(outdated_cells.begin(),outdated_cells.end());
  
  // Fill modified vertices
  if ( fill_vertices 
        && !(boost::is_same<OutputIterator,CGAL::Emptyset_iterator>::value))
  {
    fill_modified_vertices(outdated_cells.begin(), outdated_cells.end(),
                           new_vertex, modified_vertices);        
  }
  
  return new_vertex;  
}
  
#ifdef CGAL_INTRUSIVE_LIST
template <typename C3T3, typename MD>
template <typename OutdatedCells>
void
C3T3_helpers<C3T3,MD>:: 
rebuild_restricted_delaunay(OutdatedCells& outdated_cells,
                            Moving_vertices_set& moving_vertices)
{
  typename OutdatedCells::iterator first_cell = outdated_cells.begin();
  typename OutdatedCells::iterator last_cell = outdated_cells.end();
  Update_c3t3 updater(domain_,c3t3_);
  
  // Updates cells
  while ( first_cell != last_cell )
  {
    const Cell_handle cell = *first_cell++;
    c3t3_.remove_from_complex(cell);
    updater(cell);
  }

  // Get facets (returns each canonical facet only once)
  //  Facet_vector facets;
  std::set<Vertex_handle> vertex_to_proj;
  Facet_updater facet_updater(c3t3_,vertex_to_proj, updater);
  update_facets(outdated_cells, facet_updater);
  
  // now we can clear
  outdated_cells.clear();

    CGAL_HISTOGRAM_PROFILER("|vertex_to_proj|=", vertex_to_proj.size());
  // Project interior vertices
  // TODO : iterate to be sure no interior vertice become on the surface
  // because of move ?
  for ( typename std::set<Vertex_handle>::iterator it = vertex_to_proj.begin() ;
       it != vertex_to_proj.end() ;
       ++it )
  {
    Point_3 new_pos = project_on_surface((*it)->point(),*it);

    if ( new_pos != Point_3() )
    {   
      //freezing needs 'erase' to be done before the vertex is actually destroyed
      // Update moving vertices (it becomes new_vertex)
      moving_vertices.erase(*it);
    
      Vertex_handle new_vertex = update_mesh(new_pos,*it);
      c3t3_.set_dimension(new_vertex,2);

      moving_vertices.insert(new_vertex);
    }
  }
}
#endif //CGAL_INTRUSIVE_LIST

template <typename C3T3, typename MD>
template <typename ForwardIterator>
void
C3T3_helpers<C3T3,MD>:: 
rebuild_restricted_delaunay(ForwardIterator first_cell,
                            ForwardIterator last_cell,
                            Moving_vertices_set& moving_vertices)
{
  Update_c3t3 updater(domain_,c3t3_);
  
  // Get facets (returns each canonical facet only once)
  Facet_vector facets = get_facets(first_cell, last_cell);
  
  // Updates cells
  while ( first_cell != last_cell )
  {
    const Cell_handle& cell = *first_cell++;
    c3t3_.remove_from_complex(cell);
    updater(cell);
  }
  
  // Updates facets
  std::set<std::pair<Vertex_handle, Surface_patch_index> > vertex_to_proj;
  for ( typename Facet_vector::iterator fit = facets.begin() ;
       fit != facets.end() ;
       ++fit )
  {
    // Update facet
    c3t3_.remove_from_complex(*fit);
    updater(*fit);

    // Update vertex_to_proj
    if ( c3t3_.is_in_complex(*fit) )
    {
      // Iterate on vertices
      int k = fit->second;
      for ( int i=1 ; i<4 ; ++i )
      {
        const Vertex_handle& v = fit->first->vertex((k+i)&3);
        if ( c3t3_.in_dimension(v) > 2 )
        { 
          vertex_to_proj.insert
            (std::make_pair(v, c3t3_.surface_patch_index(*fit)));
        }
      }
    }
  }
  
  // Project interior vertices
  // TODO : iterate to be sure no interior vertice become on the surface
  // because of move ?
  for ( typename std::set<std::pair<Vertex_handle, Surface_patch_index> >
          ::iterator it = vertex_to_proj.begin() ;
        it != vertex_to_proj.end() ;
        ++it )
  {
    Point_3 new_pos = project_on_surface((it->first)->point(),it->first,it->second);

    if ( new_pos != Point_3() )
    {
      //freezing needs 'erase' to be done before the vertex is actually destroyed
      // Update moving vertices (it becomes new_vertex)
      moving_vertices.erase(it->first);

      Vertex_handle new_vertex = update_mesh(new_pos,it->first);
      c3t3_.set_dimension(new_vertex,2);
      
      moving_vertices.insert(new_vertex);
    }
  }
}


template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator,
          typename DeletedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point(const Vertex_handle& old_vertex,
           const Point_3& new_position,
           OutdatedCellsOutputIterator outdated_cells,
           DeletedCellsOutputIterator deleted_cells)
{
  // std::cerr << "C3T3_helpers::move_point[v2](" 
  //           << (void*)(&*old_vertex) << " = " << old_vertex->point()
  //           << " , " << new_position << ")\n";
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(old_vertex, std::back_inserter(incident_cells));
  if ( Th().no_topological_change(tr_, old_vertex, new_position, incident_cells) )
  {
    BOOST_FOREACH(Cell_handle& ch, std::make_pair(incident_cells.begin(), 
                                                  incident_cells.end()))
    {
      ch->invalidate_circumcenter();
    }
    std::copy(incident_cells.begin(),incident_cells.end(), outdated_cells);
    return move_point_no_topo_change(old_vertex,
                                     new_position);
  }
  else
  {
    return move_point_topo_change(old_vertex,
                                  new_position,
                                  outdated_cells,
                                  deleted_cells);
  }
}

template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point(const Vertex_handle& old_vertex,
           const Point_3& new_position,
           Outdated_cell_set& outdated_cells_set, 
           Moving_vertices_set& moving_vertices)
{
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(old_vertex, std::back_inserter(incident_cells));
  if ( Th().no_topological_change(tr_, old_vertex, new_position, incident_cells) )
  {
    BOOST_FOREACH(Cell_handle& ch, std::make_pair(incident_cells.begin(), 
                                                  incident_cells.end()))
    {
      ch->invalidate_circumcenter();
    }
    std::copy(incident_cells.begin(),incident_cells.end(), 
      std::inserter(outdated_cells_set, outdated_cells_set.end()));
    return move_point_no_topo_change(old_vertex, new_position);
  }
  else
  {
    moving_vertices.erase(old_vertex);

    Vertex_handle new_vertex = move_point_topo_change(old_vertex, new_position, outdated_cells_set);

    moving_vertices.insert(new_vertex);
    return new_vertex;
  }  
}  

template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change(const Vertex_handle& old_vertex,
                       const Point_3& new_position,
                       Outdated_cell_set& outdated_cells_set)
{
  Cell_set insertion_conflict_cells;
  Cell_set removal_conflict_cells;
  Facet_vector insertion_conflict_boundary;
  insertion_conflict_boundary.reserve(64);
  
  get_conflict_zone_topo_change(old_vertex, new_position,
                                std::inserter(insertion_conflict_cells,insertion_conflict_cells.end()),
                                std::back_inserter(insertion_conflict_boundary),
                                std::inserter(removal_conflict_cells, removal_conflict_cells.end()));

  for(typename Cell_set::iterator it = insertion_conflict_cells.begin();
      it != insertion_conflict_cells.end(); ++it)
      outdated_cells_set.erase(*it);
  for(typename Cell_set::iterator it = removal_conflict_cells.begin();
      it != removal_conflict_cells.end(); ++it)
      outdated_cells_set.erase(*it);

  Cell_vector outdated_cells;
  Vertex_handle nv = move_point_topo_change_conflict_zone_known(old_vertex, new_position,
                                insertion_conflict_boundary[0],
                                insertion_conflict_cells.begin(),
                                insertion_conflict_cells.end(),
                                removal_conflict_cells.begin(),
                                removal_conflict_cells.end(),
                                std::back_inserter(outdated_cells),
                                CGAL::Emptyset_iterator()); // deleted_cells

  for(typename Cell_vector::iterator it = outdated_cells.begin();
      it != outdated_cells.end(); ++it)
      outdated_cells_set.insert(*it);
  
  return nv;
}

template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator,
          typename DeletedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change(const Vertex_handle& old_vertex,
                       const Point_3& new_position,
                       OutdatedCellsOutputIterator outdated_cells,
                       DeletedCellsOutputIterator deleted_cells)
{
  Cell_set insertion_conflict_cells;
  Cell_set removal_conflict_cells;
  Facet_vector insertion_conflict_boundary;
  insertion_conflict_boundary.reserve(64);
  
  get_conflict_zone_topo_change(old_vertex, new_position,
                                std::inserter(insertion_conflict_cells,insertion_conflict_cells.end()),
                                std::back_inserter(insertion_conflict_boundary),
                                std::inserter(removal_conflict_cells, removal_conflict_cells.end()));
  
  Vertex_handle nv = move_point_topo_change_conflict_zone_known(old_vertex, new_position,
                                insertion_conflict_boundary[0],
                                insertion_conflict_cells.begin(),
                                insertion_conflict_cells.end(),
                                removal_conflict_cells.begin(),
                                removal_conflict_cells.end(),
                                outdated_cells,
                                deleted_cells);
  return nv;
}


template <typename C3T3, typename MD>
template < typename ConflictCellsInputIterator,
           typename OutdatedCellsOutputIterator,
           typename DeletedCellsOutputIterator >
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change_conflict_zone_known(
    const Vertex_handle& old_vertex,
    const Point_3& new_position,
    const Facet& insertion_boundary_facet,
    ConflictCellsInputIterator insertion_conflict_cells_begin,//ordered
    ConflictCellsInputIterator insertion_conflict_cells_end,
    ConflictCellsInputIterator removal_conflict_cells_begin,//ordered
    ConflictCellsInputIterator removal_conflict_cells_end,    
    OutdatedCellsOutputIterator outdated_cells,
    DeletedCellsOutputIterator deleted_cells)//warning : this should not be an iterator to Intrusive_list
                                             //o.w. deleted_cells will point to null pointer or so and crash
{
  Point_3 old_position = old_vertex->point();
  // make one set with conflict zone  
  Cell_set conflict_zone;
  std::set_union(insertion_conflict_cells_begin, insertion_conflict_cells_end,
                 removal_conflict_cells_begin, removal_conflict_cells_end,
                 std::inserter(conflict_zone, conflict_zone.end()));

  // Remove conflict zone cells from c3t3 (they will be deleted by insert/remove)
  remove_cells_and_facets_from_c3t3(conflict_zone.begin(), conflict_zone.end());
 
// Start Move point // Insert new_vertex, remove old_vertex
  int dimension = c3t3_.in_dimension(old_vertex);
  Index vertex_index = c3t3_.index(old_vertex);
  FT meshing_info = old_vertex->meshing_info();

  // insert new point
  Vertex_handle new_vertex = tr_.insert_in_hole(new_position, 
                                                insertion_conflict_cells_begin,
                                                insertion_conflict_cells_end,
                                                insertion_boundary_facet.first, 
                                                insertion_boundary_facet.second);

  // If new_position is hidden, update what should be and return default constructed handle
  if ( Vertex_handle() == new_vertex ) 
  { 
    std::copy(conflict_zone.begin(), conflict_zone.end(), outdated_cells);
    return old_vertex; 
  }
  // remove old point
  tr_.remove(old_vertex);
  
  c3t3_.set_dimension(new_vertex,dimension);
  c3t3_.set_index(new_vertex,vertex_index);
  new_vertex->set_meshing_info(meshing_info);
  // End Move point

  //// Fill outdated_cells
  // Get conflict zone in new triangulation and set cells outdated
  Cell_vector new_conflict_cells;
  new_conflict_cells.reserve(64);
  get_conflict_zone_topo_change(new_vertex, old_position,
                                std::back_inserter(new_conflict_cells)); 
  std::copy(new_conflict_cells.begin(),new_conflict_cells.end(),outdated_cells);

  // Fill deleted_cells
  if(! boost::is_same<DeletedCellsOutputIterator,CGAL::Emptyset_iterator>::value)
    std::copy(conflict_zone.begin(), conflict_zone.end(), deleted_cells);

  return new_vertex;
}

template <typename C3T3, typename MD>
template < typename ConflictCellsInputIterator,
           typename OutdatedCellsOutputIterator,
           typename DeletedCellsOutputIterator >
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change_conflict_zone_known(
    const Vertex_handle& old_vertex,
    const Point_3& new_position,
    ConflictCellsInputIterator conflict_cells_begin,
    ConflictCellsInputIterator conflict_cells_end,
    OutdatedCellsOutputIterator outdated_cells,
    DeletedCellsOutputIterator deleted_cells)
{
  Point_3 old_position = old_vertex->point();
  
  // Remove conflict zone cells from c3t3 (cells will be destroyed)  
  remove_cells_and_facets_from_c3t3(conflict_cells_begin, conflict_cells_end);
  
#ifdef CGAL_INTRUSIVE_LIST
  // AF: moved here from below, because the cells still exist
  //     and as we want to remove on the fly from the inplace list
  std::copy(conflict_cells_begin, conflict_cells_end, deleted_cells);
#endif

  // Move point
  Vertex_handle new_vertex = move_point_topo_change(old_vertex,new_position);
  
  // If nothing changed, return
  if ( Vertex_handle() == new_vertex )
  {
    std::copy(conflict_cells_begin,conflict_cells_end,outdated_cells);
    return old_vertex;
  }
  
  // Get conflict zone in new triangulation and set cells outdated
  Cell_vector new_conflict_cells;
  new_conflict_cells.reserve(64);
  get_conflict_zone_topo_change(new_vertex, old_position,
                                std::back_inserter(new_conflict_cells));
  
  std::copy(new_conflict_cells.begin(),new_conflict_cells.end(),outdated_cells);

   // Fill deleted_cells
#ifndef CGAL_INTRUSIVE_LIST
  //AF: move this higher up so that we can remove in the inplace list 
  std::copy(conflict_cells_begin, conflict_cells_end, deleted_cells);
#endif

  return new_vertex;
}


template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>::
move_point_topo_change(const Vertex_handle& old_vertex,
                       const Point_3& new_position)
{
  // Insert new_vertex, remove old_vertex
  int dimension = c3t3_.in_dimension(old_vertex);
  Index vertex_index = c3t3_.index(old_vertex);
  FT meshing_info = old_vertex->meshing_info();

  // insert new point
  Vertex_handle new_vertex = tr_.insert(new_position,old_vertex->cell());
  // If new_position is hidden, return default constructed handle
  if ( Vertex_handle() == new_vertex ) { return Vertex_handle(); }
  // remove old point
  tr_.remove(old_vertex);
  
  c3t3_.set_dimension(new_vertex,dimension);
  c3t3_.set_index(new_vertex,vertex_index);
  new_vertex->set_meshing_info(meshing_info);

  return new_vertex;
}
  

template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_no_topo_change(const Vertex_handle& old_vertex,
                          const Point_3& new_position,
                          OutdatedCellsOutputIterator outdated_cells)
{
  get_conflict_zone_no_topo_change(old_vertex, outdated_cells);
  return move_point_no_topo_change(old_vertex, new_position);
}


template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_no_topo_change(const Vertex_handle& old_vertex,
                          const Point_3& new_position)
{  
  // Change vertex position
  old_vertex->set_point(new_position);
  return old_vertex;  
}


/**
 * @brief Returns the projection of \c p, using direction of 
 * \c projection_vector
 */
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Point_3
C3T3_helpers<C3T3,MD>:: 
project_on_surface_aux(const Point_3& p,
                       const Point_3& ref_point,
                       const Vector_3& projection_vector) const
{
  typedef typename Gt::Segment_3 Segment_3;
  typedef typename MD::Intersection Intersection;

  // Build a segment directed as projection_direction,
  typename Gt::Compute_squared_distance_3 sq_distance =
    Gt().compute_squared_distance_3_object();
  
  typename Gt::Compute_squared_length_3 sq_length =
    Gt().compute_squared_length_3_object();
  
  typename Gt::Construct_scaled_vector_3 scale =
    Gt().construct_scaled_vector_3_object();
    
  typename Gt::Is_degenerate_3 is_degenerate =
    Gt().is_degenerate_3_object();

  typename MD::Construct_intersection construct_intersection =
    domain_.construct_intersection_object();
  
  const FT sq_dist = sq_distance(p,ref_point);
  const FT sq_proj_length = sq_length(projection_vector);
  
  if ( CGAL_NTS is_zero(sq_proj_length) )
    return ref_point;
  
  const Vector_3 projection_scaled_vector =
    scale(projection_vector, CGAL::sqrt(sq_dist/sq_proj_length));
  
  const Point_3 source = p + projection_scaled_vector;
  const Point_3 target = p - projection_scaled_vector;
  
  const Segment_3 proj_segment(source,target);
  
  if ( is_degenerate(proj_segment) )
    return ref_point;

#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3

  typename MD::Do_intersect_surface do_intersect =
    domain_.do_intersect_surface_object();

  if ( do_intersect(proj_segment) )
    return CGAL::cpp0x::get<0>(construct_intersection(proj_segment));
  else
    return ref_point;  

#else // CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3

  Intersection intersection = construct_intersection(proj_segment);
  if(CGAL::cpp0x::get<2>(intersection) == 2)
    return CGAL::cpp0x::get<0>(intersection);
  else 
    return ref_point;

#endif // CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
}

  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Plane_3
C3T3_helpers<C3T3,MD>:: 
get_least_square_surface_plane(const Vertex_handle& v,
                               Point_3& reference_point,
                               Surface_patch_index patch_index) const
{
  // Get incident facets
  Facet_vector facets;
  tr_.finite_incident_facets(v,std::back_inserter(facets));

  // Get adjacent surface points
  std::vector<Point_3> surface_point_vector;
  for ( typename Facet_vector::iterator fit = facets.begin() ;
       fit != facets.end() ;
       ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) && 
         (patch_index == Surface_patch_index() || 
          c3t3_.surface_patch_index(*fit) == patch_index) )
    {
      const Cell_handle& cell = fit->first;
      const int& i = fit->second;
      
      surface_point_vector.push_back(cell->get_facet_surface_center(i));
    }
  }

  // In some cases point is not a real surface point
  if ( surface_point_vector.empty() )
    return Plane_3();
    
  // Compute least square fitting plane
  Plane_3 plane;
  CGAL::linear_least_squares_fitting_3(surface_point_vector.begin(),
                                       surface_point_vector.end(),
                                       plane,
                                       Dimension_tag<0>());
  
  reference_point = surface_point_vector.front();

  return plane;
}
  
  
  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Point_3
C3T3_helpers<C3T3,MD>:: 
project_on_surface(const Point_3& p,
                   const Vertex_handle& v,
                   Surface_patch_index index) const
{
  // return domain_.project_on_surface(p);
  // Get plane
  Point_3 reference_point(CGAL::ORIGIN);
  Plane_3 plane = get_least_square_surface_plane(v,reference_point, index);
  
  if ( reference_point == CGAL::ORIGIN )
    return p;
  
  // Project
  if ( p != v->point() )
    return project_on_surface_aux(p,
                                  v->point(),
                                  plane.orthogonal_vector());
  else
    return project_on_surface_aux(p,
                                  reference_point,
                                  plane.orthogonal_vector());
}

  
  
template <typename C3T3, typename MD>
template <typename SliverCriterion> 
typename C3T3_helpers<C3T3,MD>::FT
C3T3_helpers<C3T3,MD>::
min_incident_value(const Vertex_handle& vh,
                   const SliverCriterion& criterion) const
{
  Cell_vector incident_cells;
  tr_.finite_incident_cells(vh,std::back_inserter(incident_cells));
  
  return min_sliver_in_c3t3_value(incident_cells, criterion);
}

template <typename OutputIterator, typename CH, typename Fct>
struct Filter {

  mutable OutputIterator out;
  const Fct& fct;

  Filter(OutputIterator out, const Fct& fct) 
    : out(out), fct(fct)
  {}

  void operator()(CH cell_handle) const
  {
    if(fct(cell_handle)){
      *out++ = cell_handle;
    }
  }

};
  
template <typename CH, typename Fct>
struct Counter {

  const Fct& fct;
  std::size_t& count;

  Counter(const Fct& fct, std::size_t& count) 
    : fct(fct), count(count)
  {}

  void operator()(CH cell_handle)
  {
    if(fct(cell_handle)){
      ++count;
    }
  }

};


template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
incident_slivers(const Vertex_handle& v,
                 const SliverCriterion& criterion,
                 const FT& sliver_bound,
                 OutputIterator out) const
{
  typedef SliverCriterion Sc;
  
  std::vector<Cell_handle> incident_cells;
  tr_.incident_cells(v, std::back_inserter(incident_cells));
  
  std::remove_copy_if(incident_cells.begin(),
                      incident_cells.end(),
                      out,
                      std::not1(Is_sliver<Sc>(c3t3_,criterion,sliver_bound)));
  
  return out;
}

template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
new_incident_slivers(const Vertex_handle& v,
                     const SliverCriterion& criterion,
                     const FT& sliver_bound,
                     OutputIterator out) const
{
  typedef SliverCriterion Sc;
  typedef Filter<OutputIterator,Cell_handle,Is_sliver<Sc> > F;
  
  Is_sliver<Sc> i_s(c3t3_, criterion, sliver_bound);
  F f(out, i_s);
  tr_.incident_cells(v,boost::make_function_output_iterator(f));

  return f.out;
}

template <typename C3T3, typename MD>
template <typename SliverCriterion>
bool
C3T3_helpers<C3T3,MD>::
is_sliver(const Cell_handle& ch,
          const SliverCriterion& criterion,
          const FT& sliver_bound) const
{
  Is_sliver<SliverCriterion> iss(c3t3_,criterion,sliver_bound);
  return iss(ch);
}

template <typename C3T3, typename MD>
template <typename SliverCriterion>
std::size_t
C3T3_helpers<C3T3,MD>::
number_of_incident_slivers(const Vertex_handle& v,
                           const SliverCriterion& criterion,
                           const FT& sliver_bound) const
{
  typedef SliverCriterion Sc;
  typedef Counter<Cell_handle,Is_sliver<Sc> > C;

  std::size_t count = 0;
  Is_sliver<Sc> is_sliver(c3t3_,criterion,sliver_bound);
  C c(is_sliver, count);
  tr_.incident_cells(v, boost::make_function_output_iterator(c));

  return count;
}

  
template <typename C3T3, typename MD>
template <typename SliverCriterion>
typename C3T3_helpers<C3T3,MD>::FT
C3T3_helpers<C3T3,MD>::
min_sliver_value(const Cell_vector& cells,
                 const SliverCriterion& criterion,
                 const bool use_cache) const
{
  using boost::make_transform_iterator;
  
  if ( cells.empty() )
    return SliverCriterion::max_value;
  
  if ( ! use_cache )
  { 
    reset_cache_validity(cells.begin(),cells.end());
  }
  
  // Return min dihedral angle
  Sliver_criterion_value<SliverCriterion> sc_value(tr_,criterion);
  
  return *(std::min_element(make_transform_iterator(cells.begin(),sc_value),
                            make_transform_iterator(cells.end(),sc_value)));
}
  
  
template <typename C3T3, typename MD>
template <typename InputIterator, typename OutputIterator>
void
C3T3_helpers<C3T3,MD>::
fill_modified_vertices(InputIterator cells_begin,
                       InputIterator cells_end,
                       const Vertex_handle& vertex,
                       OutputIterator out) const
{
  std::set<Vertex_handle> already_inserted_vertices;
  // Dont insert vertex in out
  already_inserted_vertices.insert(vertex);
  
  for ( InputIterator it = cells_begin ; it != cells_end ; ++it )
  {
    for ( int i=0 ; i<4 ; ++i )
    {
	  // Insert vertices if not already inserted
      const Vertex_handle& current_vertex = (*it)->vertex(i);
      if ( !tr_.is_infinite(current_vertex)
          && already_inserted_vertices.insert(current_vertex).second )
      {
        *out++ = current_vertex;
      }
    }
  }
}
  
  
template <typename C3T3, typename MD>
template <typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
get_conflict_zone_no_topo_change(const Vertex_handle& vertex,
                                 OutputIterator conflict_cells) const
{  
  return tr_.incident_cells(vertex,conflict_cells);
}

template <typename C3T3, typename MD>
template <typename CellsOutputIterator,
          typename FacetsOutputIterator>
void
C3T3_helpers<C3T3,MD>::
get_conflict_zone_topo_change(const Vertex_handle& v,
                              const Point_3& conflict_point,
                              CellsOutputIterator insertion_conflict_cells,
                              FacetsOutputIterator insertion_conflict_boundary,
                              CellsOutputIterator removal_conflict_cells) const
{
  // Get triangulation_vertex incident cells : removal conflict zone
  tr_.incident_cells(v, removal_conflict_cells);

  // Get conflict_point conflict zone 
  int li=0;
  int lj=0;
  typename Tr::Locate_type lt;
  Cell_handle cell = tr_.locate(conflict_point, lt, li, lj, v->cell());
  
  if ( lt == Tr::VERTEX ) // Vertex removal is forbidden 
    return;
  
  // Find conflict zone
  tr_.find_conflicts(conflict_point,
                     cell,
                     insertion_conflict_boundary,
                     insertion_conflict_cells);
}
  
template <typename C3T3, typename MD>
template <typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
get_conflict_zone_topo_change(const Vertex_handle& vertex,
                              const Point_3& conflict_point,
                              OutputIterator conflict_cells) const
{
  // Get triangulation_vertex incident cells
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(vertex, std::back_inserter(incident_cells));
  
  // Get conflict_point conflict zone
  Cell_vector deleted_cells;
  deleted_cells.reserve(64);
  
  // Vertex removal is forbidden 
  int li=0;
  int lj=0;
  typename Tr::Locate_type locate_type;
  Cell_handle cell = tr_.locate(conflict_point,
                                locate_type,
                                li,
                                lj,
                                vertex->cell());
  
  if ( Tr::VERTEX == locate_type )
    return conflict_cells;
  
  // Find conflict zone
  tr_.find_conflicts(conflict_point,
                     cell,
                     CGAL::Emptyset_iterator(),
                     std::back_inserter(deleted_cells),
                     CGAL::Emptyset_iterator());
  
  // Compute union of conflict_point conflict zone and triangulation_vertex
  // incident cells
  std::sort(deleted_cells.begin(),deleted_cells.end());
  std::sort(incident_cells.begin(),incident_cells.end());
  
  std::set_union(deleted_cells.begin(), deleted_cells.end(),
                 incident_cells.begin(), incident_cells.end(),
                 conflict_cells);
  
  return conflict_cells;
}
  
  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Facet_boundary
C3T3_helpers<C3T3,MD>::
get_surface_boundary(const Facet_vector& facets) const
{
  Facet_boundary boundary;
  typename Facet_vector::const_iterator fit = facets.begin();
  for ( ; fit != facets.end() ; ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      const Surface_patch_index surface_index = c3t3_.surface_patch_index(*fit);
      const int k = fit->second;
      Vertex_handle v1 = fit->first->vertex((k+1)&3);
      Vertex_handle v2 = fit->first->vertex((k+2)&3);
      Vertex_handle v3 = fit->first->vertex((k+3)&3);
      
      // Check that each vertex is a surface one
      // This is a trick to ensure that in_domain vertices stay inside
      if ( c3t3_.in_dimension(v1) > 2
          || c3t3_.in_dimension(v2) > 2
          || c3t3_.in_dimension(v3) > 2 )
      {
        // Ordered_edge(tr_.infinite_vertex(),v1) can't be in boundary
        // So if there is a boundary facets which is not built on 3 boundary
        // vertices, check of boundary equality before and after the move will
        // fail (we know that this is not the case before)
        update_boundary(boundary,
                        Ordered_edge(Vertex_handle(),Vertex_handle()),
                        v3,
                        Surface_patch_index());
        return boundary;
      }
      
      order_handles(v1,v2,v3);
      
      CGAL_assertion(v1<v2);
      CGAL_assertion(v2<v3);
      
      update_boundary(boundary, Ordered_edge(v1,v2), v3, surface_index);
      update_boundary(boundary, Ordered_edge(v1,v3), v2, surface_index);
      update_boundary(boundary, Ordered_edge(v2,v3), v1, surface_index);
    }
  }

  // std::cerr.precision(17);
  // std::cerr << "boundary { ";
  // BOOST_FOREACH(const typename Facet_boundary::value_type& v,
  //               boundary)
  // {
  //   std::cerr << "(" << v.first.first->point() << ", " << v.first.second->point() << ", " << v.second.first << ") ";
  // }
  // std::cerr << "}\n";
  return boundary;
}
  
template <typename C3T3, typename MD>
bool
C3T3_helpers<C3T3,MD>::
check_no_inside_vertices(const Facet_vector& facets) const
{
  typename Facet_vector::const_iterator fit = facets.begin();
  for ( ; fit != facets.end() ; ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      const int k = fit->second;
      const Vertex_handle& v1 = fit->first->vertex((k+1)&3);
      const Vertex_handle& v2 = fit->first->vertex((k+2)&3);
      const Vertex_handle& v3 = fit->first->vertex((k+3)&3);
      
      // Check that each vertex is a surface one
      if ( c3t3_.in_dimension(v1) > 2
          || c3t3_.in_dimension(v2) > 2
          || c3t3_.in_dimension(v3) > 2 )
      {
        return false;
      }
    }
  }
  
  return true;
}
  
  
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_C3T3_HELPERS_H
