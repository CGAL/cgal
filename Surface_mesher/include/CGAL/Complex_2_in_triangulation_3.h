// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_3_H

#include <CGAL/license/Surface_mesher.h>


// TODO: add the iterators
// TODO: document the output/input function of C2T3?

#include <CGAL/circulator.h>
#include <CGAL/iterator.h>
#include <CGAL/Union_find.h>
#include <set>
#include <map>
#include <vector>
#include <boost/format.hpp>

namespace CGAL {

  namespace details {

    template <typename Tr, typename Edge_info>
    class C2t3_helper_class
    {
    protected:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef std::pair<Vertex_handle, Vertex_handle> Pair_of_vertices;
      // computes and return an ordered pair of Vertex
      Pair_of_vertices
      make_ordered_pair(const Vertex_handle vh1, const Vertex_handle vh2) const {
	if (vh1 < vh2) {
	  return std::make_pair(vh1, vh2);
	}
	else {
	  return std::make_pair(vh2, vh1);
	}
      }
    };

    template <typename Tr, typename Edge_info>
    class C2t3_mark_edges_helper_class : public C2t3_helper_class<Tr,Edge_info>
    {
    protected:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Edge Edge;
      typedef typename C2t3_helper_class<Tr,Edge_info>::Pair_of_vertices Pair_of_vertices;
      typedef std::map<Pair_of_vertices, Edge_info> Marked_edges;
      Marked_edges marked_edges;
    public:
      void mark(const Cell_handle& c, const int& i, const int& j,
		const Edge_info& info = Edge_info())
      {
	marked_edges.insert(std::make_pair(this->make_ordered_pair(c->vertex(i),
								   c->vertex(j)),
					   info));
      }
      void mark(const Edge& e,
		const Edge_info& info = Edge_info())
      {
	mark(e.first, e.second, e.third, info);
      }
      bool unmark(const Cell_handle& c, const int& i, const int& j)
      {
	return marked_edges.erase(this->make_ordered_pair(c->vertex(i),
							  c->vertex(j))) > 0;
      }
      bool unmark(const Edge& e)
      {
	return unmark(e.first, e.second, e.third);
      }
      Edge_info& 
      get_info(const Edge& e)
      {
	return get_info(e.first, e.second, e.third);
      }
      Edge_info& 
      get_info(const Cell_handle& c, const int& i, const int& j)
      {
	return this->marked_edges[this->make_ordered_pair(c->vertex(i),
							  c->vertex(j))];
      }

    }; // end of class C2t3_mark_edges_helper_class<Tr, Edge_info>

    template <typename Tr>
    class C2t3_mark_edges_helper_class<Tr,void> 
      : public C2t3_helper_class<Tr,void>
    {
    protected:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Edge Edge;
      typedef typename C2t3_helper_class<Tr,void>::Pair_of_vertices Pair_of_vertices;
      typedef std::set<Pair_of_vertices> Marked_edges;
      Marked_edges marked_edges;
    public:
      void mark(const Cell_handle& c, const int& i, const int& j)
      {
	marked_edges.insert(this->make_ordered_pair(c->vertex(i),
						    c->vertex(j)));
      }
      void mark(const Edge& e)
      {
	mark(e.first, e.second, e.third);
      }
      bool unmark(const Cell_handle& c, const int& i, const int& j)
      {
	return marked_edges.erase(this->make_ordered_pair(c->vertex(i), 
							  c->vertex(j))) > 0;
      }
      bool unmark(const Edge& e)
      {
	return unmark(e.first, e.second, e.third);
      }
    }; // end of specialiazation C2t3_mark_edges_helper_class<Tr,void>

    
  } // end nested-namespace details (in CGAL::)

namespace Surface_mesher {

template < class Tr>
typename Tr::size_type number_of_facets_on_surface(const Tr& T) {
  typename Tr::size_type result=0;
  for (typename Tr::Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end(); ++fit)
    if (fit->first->is_facet_on_surface (fit->second))
      ++result;
  return result;
}

} // end nested-namespace Surface_mesher (in CGAL)

template < class Tr, typename Edge_info_ = void >
class Complex_2_in_triangulation_3 :
    public details::C2t3_mark_edges_helper_class<Tr, Edge_info_>
{
public:
  typedef Complex_2_in_triangulation_3 <Tr, Edge_info_> Self;
  typedef details::C2t3_mark_edges_helper_class<Tr, Edge_info_> Base;
  typedef typename Base::Marked_edges Marked_edges;

  typedef Tr Triangulation;
  typedef Edge_info_ Edge_info;

  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Triangulation::Facet Facet;
  typedef typename Triangulation::Edge Edge;

  typedef std::set<Facet> Facets;

  typedef std::size_t size_type;

  typedef Const_circulator_from_container<Facets> Facet_circulator;

  typedef typename Base::Pair_of_vertices Pair_of_vertices;
  typedef std::map <Pair_of_vertices,
		    std::pair<int, Facets > >
                                                  Edge_facet_counter;
  enum Face_status{ NOT_IN_COMPLEX = 0,
                    ISOLATED = 1, // - An ISOLATED edge is a marked edge,
                                  //   without any incident facets.
                    BOUNDARY,     // - An edge is on BOUNDARY if it has only
                                  //   one incident facet.
                                  // - A vertex is on BOUNDARY if all its
                                  //   incident edges are REGULAR or on
                                  //   BOUNDARY, at least one is on
                                  //   BOUNDARY, and the incident facets
                                  //   form only one connected component.
                    REGULAR,      // - A facet that is in the complex is
                                  //   REGULAR.
                                  // - An edge is REGULAR if it has
                                  //   exactly two incident facets.
                                  // - A vertex is REGULAR if all it
                                  //   incident edges are REGULAR, and the
                                  //   incident facets form only one
                                  //   connected component.
                    SINGULAR};    // - SINGULAR is for all other cases.

  class Iterator_not_in_complex {
    const Self* self;
  public:
    Iterator_not_in_complex(const Self* self = 0) : self(self) 
    {
    }
    
    template <typename Iterator> // Facet or Edges iterators
    bool operator()(Iterator it) const {
      if(self)
	return ! self->is_in_complex(*it);
      else
	return true;
    }
  }; // end struct Iterator_not_in_complex

  class Vertex_not_in_complex {
    Self* self;
  public:
    Vertex_not_in_complex(){} //added for SWIG wrapping
    Vertex_not_in_complex(Self* self) : self(self) 
    {
    }
    
    bool operator()(Vertex_handle v) const { // Takes as argument an iterator to a
                                             // Vertex, convertible to Vertex_handle.
      return ! self->is_in_complex(v);
    }
  }; // end struct Vertex_not_in_complex

  class Facet_not_in_complex {
    Self* self;
  public:
    Facet_not_in_complex(Self* self) : self(self) 
    {
    }
    
    bool operator()(Facet f) const {
      return ! self->is_in_complex(f);
    }
  }; // end struct Facet_not_in_complex

  class Iterator_not_on_boundary {
    Self* self;
  public:
    Iterator_not_on_boundary(){}  //added for SWIG wrapping
    Iterator_not_on_boundary(Self* self) : self(self) 
    {
    }

    template <class Edge_iterator>
    bool operator()(Edge_iterator eit) const {
      return self->face_status(*eit)!= BOUNDARY;
    }

  };

  typedef Filter_iterator<typename Triangulation::Finite_facets_iterator,
                          Iterator_not_in_complex> Facet_iterator;
  typedef Filter_iterator<typename Triangulation::Finite_edges_iterator,
                          Iterator_not_in_complex> Edge_iterator;

  // class to ensure that Vertex_iterator is convertible to Vertex_handle
  class Vertex_iterator : 
    public Filter_iterator<typename Triangulation::Finite_vertices_iterator,
                           Vertex_not_in_complex>
  {
    typedef typename Triangulation::Finite_vertices_iterator Tr_iterator;
    typedef Filter_iterator<typename Triangulation::Finite_vertices_iterator,
                            Vertex_not_in_complex> Base;
    typedef typename Base::Predicate Predicate;
    typedef Vertex_iterator Self;
  public:
    Vertex_iterator(){} //added for SWIG wrapping
    Vertex_iterator(Base i) : Base(i)
    {
    }

    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
  
    operator Vertex_handle() const // const added for SWIG wrapping
    {
      return Vertex_handle(this->base());
    }
  };

  typedef Filter_iterator<typename Triangulation::Finite_edges_iterator,
                          Iterator_not_on_boundary> Boundary_edges_iterator;

protected:
  Triangulation& tr;
  Edge_facet_counter  edge_facet_counter;
  size_type m_number_of_facets;

public:
  Facet canonical_facet(Cell_handle c, int i) const {
    Cell_handle c2 = c->neighbor(i);
    return (c2 < c) ? std::make_pair(c2,c2->index(c)) : std::make_pair(c,i);
  }

  Facet opposite_facet(Facet f) const {
    Cell_handle c2 = f.first->neighbor(f.second);
    return std::make_pair(c2,c2->index(f.first));
  }
 public:

  // Constructors

  Complex_2_in_triangulation_3 (Triangulation& t) 
    : tr(t), m_number_of_facets(0)
  {
  }

  void clear()
  {
    m_number_of_facets = 0;
    edge_facet_counter.clear();
    this->marked_edges.clear();
  }

  // Access functions

  Triangulation& triangulation()
  {
    return tr;
  }

  const Triangulation& triangulation() const
  {
    return tr;
  }

  Face_status face_status (const Facet& f) const {
    return face_status (f.first, f.second);
  }

  Face_status face_status (const Cell_handle c, const int i) const {
    return (c->is_facet_on_surface(i)) ? REGULAR : NOT_IN_COMPLEX;
  }

  Face_status face_status (const Edge& e) const {
    return face_status(e.first->vertex(e.second), e.first->vertex(e.third));
  }

  Face_status face_status (const Cell_handle c, const int i, const int j) const {
    return face_status(c->vertex(i), c->vertex(j));
  }

  Face_status face_status (const Vertex_handle& va,
                           const Vertex_handle& vb) const
  {
    typename Edge_facet_counter::const_iterator it =
      edge_facet_counter.find(this->make_ordered_pair(va, vb));

    if (it == edge_facet_counter.end())
    {
      if(is_marked(va, vb))
	return ISOLATED;
      else
	return NOT_IN_COMPLEX;
    }

    switch (it->second.first)
    {
    case 0 : return ISOLATED;
    case 1 : return BOUNDARY;
    case 2 : return REGULAR;
    default : return SINGULAR;
    }
  } // end face_status(const Vertex_handle&, const Vertex_handle&)

  Face_status face_status (const Vertex_handle& v)
  {
    if(v->is_c2t3_cache_valid() && v->cached_number_of_incident_facets() == 0)
      return NOT_IN_COMPLEX;

    //test incident edges for REGULARITY and count BOUNDARY edges
    typename std::vector<Vertex_handle> vertices;
    vertices.reserve(64);
    tr.incident_vertices(v, std::back_inserter(vertices));
    int number_of_boundary_incident_edges = 0; //COULD BE a bool
    for (typename std::vector<Vertex_handle>::iterator vit=vertices.begin();
	 vit != vertices.end();
	 vit++ ) 
    {
      switch( face_status(v, *vit) )
      {
      case NOT_IN_COMPLEX: case REGULAR: break;
      case BOUNDARY: ++number_of_boundary_incident_edges; break;
      default : return SINGULAR;
      }
    }

    // from now on incident edges (in complex) are REGULAR or BOUNDARY

    int nb_incident_facets, nb_components;
    union_find_of_incident_facets(v, nb_incident_facets, nb_components);

    if ( nb_incident_facets == 0 ) 
      return NOT_IN_COMPLEX;
    else if ( nb_components > 1 ) 
      return SINGULAR;
    else // REGULAR OR BOUNDARY
    {
      if (number_of_boundary_incident_edges != 0)
        return BOUNDARY;
      else
        return REGULAR;
    }
  } //end of face_status(Vertex_handle)

  bool is_marked(const Vertex_handle& va, const Vertex_handle& vb) const
  {
    typename Marked_edges::const_iterator it = 
      this->marked_edges.find(this->make_ordered_pair(va, vb));
    return it != this->marked_edges.end();
  }

  bool is_marked(const Cell_handle& c, const int& i, const int& j) const
  {
    return is_marked(c->vertex(i),c->vertex(j));
  }

  bool is_marked(const Edge& e) const
  {
    return is_marked(e.first, e.second, e.third);
  }

  // This function should be called only when incident edges
  // are known to be REGULAR OR BOUNDARY
  bool is_regular_or_boundary_for_vertices(Vertex_handle v)  {
    int i,j;
    union_find_of_incident_facets(v,i,j);
    return (j == 1);
  }

   bool is_in_complex (Vertex_handle v) {
    int i,j;
    union_find_of_incident_facets(v,i,j);
    return ( i != 0);
  }

  // extract the subset F of facets of the complex incident to v
  // set i to the number of facets in F
  // set j to the number of connected component of the adjacency graph
  //     of F
  void union_find_of_incident_facets(const Vertex_handle v, int& i, int& j) {
    if( v->is_c2t3_cache_valid() )
    {
      i = v->cached_number_of_incident_facets();
      j = v->cached_number_of_components();
      return;
    }

    Union_find<Facet> facets;
    incident_facets(v, std::back_inserter(facets));

    typedef std::map<Vertex_handle, 
      typename Union_find<Facet>::handle>  Vertex_Set_map;
    typedef typename Vertex_Set_map::iterator Vertex_Set_map_iterator;

    Vertex_Set_map vsmap;

    for(typename Union_find<Facet>::iterator it = facets.begin();
	it != facets.end();
	++it){
      const Cell_handle& ch = (*it).first;
      const int& i = (*it).second;
      for(int j=0; j < 3; ++j){
	const Vertex_handle w = ch->vertex(tr.vertex_triple_index(i,j));
	if(w != v){
	  Vertex_Set_map_iterator vsm_it = vsmap.find(w);
	  if(vsm_it != vsmap.end()){
	    facets.unify_sets(vsm_it->second, it);
	  } else {
	    vsmap.insert(std::make_pair(w, it));
	  }
	}
      }
    }
    
    i = static_cast<int>(facets.size());  // we cast as it cannot be too many
    j = static_cast<int>(facets.number_of_sets());
    v->set_c2t3_cache(i, j);
    return;
  }
  
  bool is_in_complex (const Facet& f) const {
    return is_in_complex (f.first, f.second);
  } 

  bool is_in_complex (const Cell_handle c, const int i) const {
    return  face_status(c,i) != NOT_IN_COMPLEX;
  }

  bool is_in_complex (const Cell_handle& c,const int i, const int j) const {
    return  face_status(c,i,j) != NOT_IN_COMPLEX;
  }  
  
  bool is_in_complex (const Edge& e) const {
    return  face_status(e) != NOT_IN_COMPLEX;
  }

  size_type number_of_facets() const
  {
    return m_number_of_facets;
  }

  size_type number_of_edges() const
  {
    return edge_facet_counter.size();
  }

  size_type number_of_marked_edges() const
  {
    return this->marked_edges.size();
  }

  Facet_circulator incident_facets (const Edge& e) {
    typename Edge_facet_counter::iterator it = 
      edge_facet_counter.find(this->make_ordered_pair(e.first->vertex(e.second),
						      e.first->vertex(e.third)));
    
    if( it == edge_facet_counter.end() )
      return Facet_circulator();
    else
    {
      // position the circulator on the first element of the facets set
      Facets& lof = it->second.second;
      return Facet_circulator(&lof);
    }
  }

  /** @TODO: document this class in the
      SurfaceMeshComplex_2InTriangulation_3 concept.
   */
  template <typename OutputIterator>
  OutputIterator incident_facets(const Vertex_handle v, OutputIterator it)
  {
    // TODO: review this function (Laurent Rineau)

    // We assume that for the generated facets the Cell_handle is smaller than the opposite one
    tr.incident_facets(v,
		       CGAL::filter_output_iterator(it, 
						    Facet_not_in_complex(this)));
    return it;
  }

  /** This function assumes that the edge is regular. */
  Facet neighbor(Cell_handle ch, int index, int j) const 
  {
    const int i1  = tr.vertex_triple_index(index, tr. cw(j));
    const int i2  = tr.vertex_triple_index(index, tr.ccw(j));

    Edge edge = Edge(ch, i1, i2);
    CGAL_assertion(face_status(edge) == REGULAR);

    typename Tr::Facet_circulator facet_circ = 
      tr.incident_facets(edge, ch,index);
    do { 
      ++facet_circ;
    } while(! is_in_complex(*facet_circ) );
    return opposite_facet(*facet_circ);
  }
  
    /** This function assumes that the edge is regular. */
  Facet neighbor(Facet f, int j) const 
  {
    return neighbor(f.first,f.second,j);
  }

  // Setting functions

  void add_to_complex (const Facet& f) {
    add_to_complex (f.first, f.second);
  }

  void add_to_complex (const Cell_handle c, const int i) {
    change_in_complex_status<true, false>(c, i);
  }

  // backward compatibility with implementation of CGAL-3.2
  void set_in_complex (const Facet& f) {
    add_to_complex(f);
  }
  // backward compatibility with implementation of CGAL-3.2
  void set_in_complex (const Cell_handle c, const int i) {
    add_to_complex(c, i);
  }

  template <bool in_complex, bool force_modification>
  void change_in_complex_status(const Cell_handle c, const int i)
  {
    // if not already in the complex
    if ( force_modification || 
         (in_complex ? 
          face_status (c, i) == NOT_IN_COMPLEX
          : face_status (c, i) != NOT_IN_COMPLEX) )
    {
      if(in_complex) 
        ++m_number_of_facets;
      else
        --m_number_of_facets;

      const Facet f = canonical_facet(c, i);

      c->set_facet_on_surface(i, in_complex);

      switch( tr.dimension() )
      {
      case 3:
        {
          const Cell_handle& c2 = c->neighbor(i);
          const int& i2 = c2->index(c);
          c2->set_facet_on_surface(i2, in_complex);
        }
        break;
      case 2:
        break;
      default:
        CGAL_error();
      }

      const int dimension_plus_1 = tr.dimension() + 1;

      // update c2t3 for edges of f
      // We consider only pairs made by vertices without i
      for (int j = 0; j < dimension_plus_1; j++) {
        for (int k = j + 1; k < dimension_plus_1; k++) {
          if ( (i != j) && (i != k) ){

            const Pair_of_vertices e = 
              this->make_ordered_pair(c->vertex(j),
				      c->vertex(k));

            if(in_complex)
            {
              ++(edge_facet_counter[e].first);
              edge_facet_counter[e].second.insert(f); // @TODO: beurk.
                                                           // Recode this!
            }
            else 
            {
              typename Edge_facet_counter::iterator it =
                edge_facet_counter.find(e);

              CGAL_assertion( it != edge_facet_counter.end() );
              
	      it->second.second.erase(f);
	      --(it->second.first);
	      CGAL_assertion(it->second.first >= 0);
	      if(it->second.first == 0)
	      {
		// if the edge is marked, leave it ISOLATED.
		if(!is_marked(e.first, e.second))
		  edge_facet_counter.erase(it);
	      }
            }
          }
        }
      }

      // update c2t3 for vertices of f
      for (int j = 0; j < dimension_plus_1; j++) {
        if (j != i)
          c->vertex(j)->invalidate_c2t3_cache();
      }
    }
  }

 
  void remove_from_complex (const Facet& f) {
    remove_from_complex (f.first, f.second);
  }

  void remove_from_complex (const Cell_handle c, const int i) {
    change_in_complex_status<false, false>(c, i);
  }

  Facet_iterator facets_begin() const {
    return CGAL::filter_iterator(tr.finite_facets_end(),
                                 Iterator_not_in_complex(this),
                                 tr.finite_facets_begin());
  }

  Facet_iterator facets_end() const {
    return CGAL::filter_iterator(tr.finite_facets_end(),
                                 Iterator_not_in_complex(this));
  }

  
  Edge_iterator edges_begin(){
    return CGAL::filter_iterator(tr.finite_edges_end(),
                                 Iterator_not_in_complex(this),
                                 tr.finite_edges_begin());
  }

  Edge_iterator edges_end(){
    return CGAL::filter_iterator(tr.finite_edges_end(),
                                 Iterator_not_in_complex(this));
  }

  Vertex_iterator vertices_begin(){
    return CGAL::filter_iterator(tr.finite_vertices_end(),
                                 Vertex_not_in_complex(this),
                                 tr.finite_vertices_begin());
  }

  Vertex_iterator vertices_end(){
    return CGAL::filter_iterator(tr.finite_vertices_end(),
                                 Vertex_not_in_complex(this));
  }

  Boundary_edges_iterator boundary_edges_begin() {
    return CGAL::filter_iterator(tr.finite_edges_end(),
                                 Iterator_not_on_boundary(this),
                                 tr.finite_edges_begin()); 
  }

  Boundary_edges_iterator boundary_edges_end() {
    return CGAL::filter_iterator(tr.finite_edges_end(),
                                 Iterator_not_on_boundary(this)); 
  }

  bool is_valid(bool verbose = false)
  {
    const typename Tr::size_type nb = number_of_facets_on_surface(tr);
    if(number_of_facets() != nb)
    {
      if(verbose) {
        std::cerr << boost::format("C2t3: Invalid number of facet: %1% (should be %2%)!\n")
          % number_of_facets() % nb;
      }
      return false;
    }
    for(Facet_iterator it = facets_begin(),
          end = facets_end();
        it != end; ++it)
    {
      CGAL_assertion(it->first->is_facet_on_surface(it->second));
      const Facet& f = tr.mirror_facet(*it);
      if(!f.first->is_facet_on_surface(f.second))
      {
        if(verbose) {
          std::cerr << 
            boost::format("C2t3: facet (%1%, %2%) is marked on surface"
                          "will its mirror facet (%3, %4) is not!\n")
            % &*it->first % it->second
            % &*f.first % f.second;
        }
        return false;
      }
    }
    return true;
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Tr>()();
  }
#endif
}; // end Complex_2_in_triangulation_3

template < class Tr, typename Edge_info>
std::istream & 
operator>> (std::istream& is, Complex_2_in_triangulation_3<Tr, Edge_info>& c2t3)
{
  c2t3.clear();
  is >> c2t3.triangulation();

  // restore datas of c2t3
  for(typename Tr::Finite_facets_iterator fit = 
        c2t3.triangulation().finite_facets_begin();
      fit != c2t3.triangulation().finite_facets_end();
      ++fit)
    if(fit->first->is_facet_on_surface(fit->second))
      c2t3.template change_in_complex_status<true, true>(fit->first, fit->second);

  return is;
}

template < class Tr, typename Edge_info>
std::ostream & 
operator<< (std::ostream& os, const Complex_2_in_triangulation_3<Tr, Edge_info> &c2t3)
{
  return os << c2t3.triangulation();
}

} // end namespace CGAL

#endif // CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
