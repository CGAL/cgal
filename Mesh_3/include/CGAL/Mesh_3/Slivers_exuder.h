// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_SLIVERS_EXUDER_H
#define CGAL_MESH_3_SLIVERS_EXUDER_H

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Double_map.h>
#include <CGAL/iterator.h>
#include <map>
#include <vector>
#include <set>
#include <iomanip> // std::setprecision
#include <iostream> // std::cerr/cout
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>

#include <CGAL/radius_ratio.h>

#include <CGAL/Implicit_surfaces_mesher_3.h> // for Mesh_3::check_c2t3

namespace CGAL {

namespace Mesh_3 {

template <typename K>
typename K::FT
radius_ratio(const typename K::Tetrahedron_3& t, K k = K())
{
  return radius_ratio(t[0], t[1], t[2], t[3], k);
}

  namespace details { // various function objects
    template <typename Map1, typename Map2>
    class Map_value_of : 
      std::unary_function<typename Map1::Direct_entry,
                          typename Map2::mapped_type>
    {
      Map2* m;
    public:
      Map_value_of(Map2& map) : m(&map) 
      {
      }

      typedef typename std::unary_function<typename Map1::Direct_entry,
					   typename Map2::mapped_type> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;

      result_type&
      operator()(const argument_type& p) const
      { // returns (*m)[p.first]
        typename Map2::iterator it = m->find(p.first);
        CGAL_assertion( it != m->end() );
	return it->second;
      }
    }; // end class Map_value_of

    template <typename Map>
    struct First_of : 
      public std::unary_function<typename Map::value_type,
                                 const typename Map::key_type&>
    {
      typedef std::unary_function<typename Map::value_type,
				  const typename Map::key_type&> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;
      
      const typename Map::key_type&
      operator()(const typename Map::value_type& p) const
      {
	return p.first;
      }
    }; // end class First_of

    template <typename Map>
    struct Second_of : 
      public std::unary_function<typename Map::value_type,
                                 const typename Map::mapped_type&>
    {
      typedef std::unary_function<typename Map::value_type,
				  const typename Map::mapped_type&> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;
      
      const typename Map::mapped_type&
      operator()(const typename Map::value_type& p) const
      {
	return p.second;
      }
    }; // end class Second_of

    template <typename Map>
    class Value_of :
      std::unary_function<typename Map::data_type,
			  typename Map::mapped_type>
    {
      Map* m;
    public:
      typedef std::unary_function<typename Map::key_type,
				  typename Map::mapped_type&> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;
      
      Value_of(Map* map) : m(map) 
      {
      }
      
      typename Map::mapped_type&
      operator()(const typename Map::key_type& k)
      {
	return (*m)[k];
      }
    }; // end class Value_of

    template <typename Gt, typename Vertex_handle>
    class Distance_from_v :
      public std::unary_function<Vertex_handle, typename Gt::RT>
    {
      Vertex_handle * v;
      Gt gt;
    public:
      Distance_from_v(Vertex_handle& vh,
		   Gt geom_traits = Gt()) 
	: v(&vh), gt(geom_traits)
      {
      }
      
      typename Gt::RT
      operator()(const Vertex_handle& vh) const
      {
        typedef typename Gt::Compute_squared_distance_3
            Compute_squared_distance_3;
        Compute_squared_distance_3 distance =
            gt.compute_squared_distance_3_object();

	return distance((*v)->point(), vh->point());
      }
    }; // end class Distance_from_v

  } // end namespace details

template <
  typename C2T3, 
  typename FT = typename C2T3::Triangulation::Geom_traits::FT
  >
class Slivers_exuder
{
public: //types
  typedef typename C2T3::Triangulation Tr;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Vertex_handle Vertex_handle;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
  typedef typename Weighted_point::Weight Weight;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Geom_traits::Tetrahedron_3 Tetrahedron_3;
  
private: //types
  /** Umbrella will represent the link of the umbrella of a point in the
   *  restricted Delaunay. */
  typedef std::set<std::pair<Vertex_handle, Vertex_handle> > Umbrella;

  
  /** Pre_star will represent the pre-star of a point. It is a
   *  (double)-map of Facet (viewed from cells inside the star),
   *  ordered by the critial_radius of the point with the cell lies on
   *  the facet, at the exterior of the pre-star. */
  typedef CGAL::Double_map<Facet, double> Pre_star;

public: // methods
  Slivers_exuder(C2T3 c2t3, double d = 0.45)
    : c2t3(c2t3), tr(c2t3.triangulation()),
      sq_delta(d*d), num_of_pumped_vertices(0),
      num_of_ignored_vertices(0), total_pumping(0.0),
      initialized(false)
  {
  }

  template <typename Handle>
  static
  void order_two_handles(Handle& h1, Handle& h2)
  {
    if( h2 < h1 )
      std::swap(h1, h2);
  }

  static 
  void find_the_two_other_indices(int i1, int i2,
				  int& i3, int& i4)
    // given i1 and i2 (i1<>i2), returns i3 and i4 so that
    // i3 < i4 and i1+i2+i3+i4=6
  {
    CGAL_precondition( 0 <= i1 && i1 <= 3 );
    CGAL_precondition( 0 <= i2 && i2 <= 3 );
    CGAL_precondition( i1 != i2 );
    
    for(i3 = 0; i3 < 4; ++i3)
      if( i3 != i1 && i3 != i2 )
	break;
    i4 = 6 - i1 - i2 - i3;

    CGAL_postcondition( 0 <= i3 && i3 < i4 && i4 <= 3 );
    CGAL_postcondition( i3 != i1 && i3 != i2 );
    CGAL_postcondition( i4 != i1 && i4 != i2 );
  }

  /** This fonction verifies that every vertex of the umbrella is found
      twice (which means almost that the umbrella is a circle). */
  static
  bool check_umbrella(Umbrella& umbrella)
  {
    typedef std::map<Vertex_handle, int> My_map;
    My_map my_map;

    for(typename Umbrella::const_iterator p_it = umbrella.begin();
	p_it != umbrella.end();
	++p_it)
    {
      ++my_map[p_it->first];
      ++my_map[p_it->second];
    }
    for(typename My_map::const_iterator it = my_map.begin();
	it != my_map.end(); ++it)
      if( it->second != 2 )
	return false;

    return true;
  } // end check_umbrella

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
  /** This function verifies that the pre_star contains exactly the set of
      facets given by the sequence [begin, end[. 

      If v!=0, it also fills another Pre_star object, from the sequence [begin,
      end[, and checks that is in the same order as pre_star.
  */
  template <class Input_facet_it>
  bool check_pre_star(const Pre_star& pre_star, 
                      Input_facet_it begin,
                      Input_facet_it end,
                      const Vertex_handle v = Vertex_handle())
  {
    Pre_star pre_star_copy = pre_star;
    Pre_star pre_star2;

    // fill pre_star2
    for(Input_facet_it fit = begin;
        fit != end;
        ++fit)
    {
      if(v != Vertex_handle())
      {
        const Facet opposite_facet = tr.mirror_facet(*fit);
        pre_star2.insert(*fit, compute_critical_radius(v,
                                                       opposite_facet.first));
      }
    }

    if(v != Vertex_handle())
    {
      while(!pre_star_copy.empty())
      {
        if(pre_star2.empty()) {
          std::cerr << "pre_star is too big!\n";
          return false;
        }
        if(pre_star_copy.front()->first != pre_star2.front()->first) {
          std::cerr << "bad order\n";
          return false;
        }
        pre_star2.pop_front();
        pre_star_copy.pop_front();
      }
    }

    if( ! pre_star2.empty() ) {
      std::cerr << "pre_star is too small!\n";
      return false;
    }

    pre_star_copy = pre_star;

    for(Input_facet_it fit = begin;
        fit != end;
        ++fit)
    {
      if(!pre_star_copy.erase(*fit))
        return false;
      if(v != Vertex_handle())
      {
        const Facet opposite_facet = tr.mirror_facet(*fit);
        pre_star2.insert(*fit, compute_critical_radius(v,
                                                       opposite_facet.first));
      }
    }
    if( !pre_star_copy.empty() )
      return false;

    return true;
  }

  /** This function verifies that the pre_star contains exactly the set of
      facets on the boundary of the conflict zone of the weighted point wp.
      The vertex handle vh is an hint for the location of wp.
      
      It also fills another Pre_star object, and checks that is in the same
      order as pre_star.
  */
  template <typename Cell_handle_output_iterator>
  bool check_pre_star(const Pre_star& pre_star, 
                      const Weighted_point& wp,
                      const Vertex_handle& vh,
                      Cell_handle_output_iterator cells_out = Emptyset_iterator())
  {
    std::vector<Facet> internal_facets;
    std::vector<Facet> boundary_facets;
    internal_facets.reserve(64);
    boundary_facets.reserve(64);
    
    tr.find_conflicts(wp,
                      vh->cell(),
                      std::back_inserter(boundary_facets),
                      cells_out,
                      std::back_inserter(internal_facets));

    const bool result =  check_pre_star(pre_star,
                                        boundary_facets.begin(),
                                        boundary_facets.end(),
                                        vh);
    if( ! result )
      std::cerr << "boundary_facets.size()=" << boundary_facets.size()
                << "\npre_star.size()=" << pre_star.size()
                << "\ntested wp=" << wp 
                << "\n";
    return result;
  }

  // ad hoc debug function
  std::string display_point_with_pre_star(Vertex_handle v, 
                                          const Pre_star& pre_star)
  {
    std::stringstream stream;
    stream <<  &*v << " (" << v->point() << ")";
    for(typename Pre_star::const_iterator it = pre_star.begin(); 
        it != pre_star.end();
        ++it)
    {
      const Facet& f = it->first;
      const Cell_handle& c = f.first;
      const int& index = f.second;

      for(int i = 0; i < 4 ; ++i)
        if(i != index && c->vertex(i) == v)
        {
          stream << " on pre_star";
          return stream.str();
        }
    }
    return  stream.str();
  }
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

  void init(double radius_ratio_limit = 1.)
  {
    CGAL_assertion_code(
                        std::cerr << "checking c2t3 before exudation...\n";
                        CGAL::Mesh_3::check_c2t3(c2t3);
                        std::cerr << "done.\n";
                        )
 
    stop_limit_on_radius_ratio = radius_ratio_limit;

    cells_queue.clear();
    compute_aspect_ratio_priority_map();
    std::cerr << "slivers: " << cells_queue.size() << "\n";
    initialized = true;
  }

  void pump_vertices(double radius_ratio_limit = 1.)
  {
    if( ! initialized )
      init();

    // store radius_ratio_limit in the member stop_limit_on_radius_ratio
    stop_limit_on_radius_ratio = radius_ratio_limit;

    while( !cells_queue.empty() )
    {
      std::pair<double, Cell_handle> front = *(cells_queue.front());
      Cell_handle c = front.second;
	  
      int i;
      for( i = 0; i < 4; ++i )
      {
        // do not pump surface vertices
        if( true )//c->vertex(i)->point().surface_index() == 0 )
          if( pump_vertex(c->vertex(i)) )
          {
            std::cout << "P"; // vertex has been pumped
            break;
          }
          else
            std::cout << "."; // cannot increase the weight of this vertex
        else
          std::cout << "s"; // vertex is on a surface
      }

      // if the tet could not be deleted
      if (i==4)
        cells_queue.pop_front();

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      std::cerr << "\nchecking c2t3...\n";
      CGAL::Mesh_3::check_c2t3(c2t3);
      std::cerr << "done.\n";
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

    }
  } // end function pump_vertices

  bool pump_vertex(Vertex_handle v)
  {
    CGAL_expensive_assertion(tr.is_valid());
    using boost::make_transform_iterator;

    std::vector<Cell_handle> incident_cells;
    std::vector<Vertex_handle> incident_vertices;
    
    incident_cells.reserve(64);
    incident_vertices.reserve(128);
    
    tr.incident_cells(v, std::back_inserter(incident_cells));
    tr.incident_vertices(v, std::back_inserter(incident_vertices));
    
    typedef std::map<Facet, double> Ratios;
    Ratios ratios;
    
    Pre_star pre_star;
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
    std::set<Cell_handle> pre_star_interior;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

    double best_weight = 0;
    double worst_radius_ratio = 1;

    details::Distance_from_v<Geom_traits, Vertex_handle> 
      distance_from_v(v, tr.geom_traits());

    const double sq_d_v = 
      *(std::min_element(make_transform_iterator(incident_vertices.begin(),
						 distance_from_v),
			 make_transform_iterator(incident_vertices.end(),
						 distance_from_v)));

    // This boolean value will be set to true if one of the incident cells
    // is in the domain.
    CGAL_assertion_code(bool is_incident_to_a_cell_in_domain = false);

    // filling 'pre_star' and 'ratios' with initial values
    for(typename std::vector<Cell_handle>::const_iterator cit = 
          incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
    {
      CGAL_assertion( v->point().surface_index() > 0 || (*cit)->is_in_domain() );

      const int index = (*cit)->index(v);
      const Facet f = Facet(*cit, index);
      if( (*cit)->is_in_domain() )
      {
        CGAL_assertion_code(is_incident_to_a_cell_in_domain = true);
        const double r = CGAL::to_double(radius_ratio(tr.tetrahedron(*cit), 
                                                      Geom_traits()));
        ratios[f] = r;
        if( r < worst_radius_ratio ) worst_radius_ratio = r;
      }
      const Facet opposite_facet = tr.mirror_facet(f);
      if( !tr.is_infinite(opposite_facet.first) )
      {
        pre_star.insert(f, compute_critical_radius(v,
                                                   opposite_facet.first));
        // if the facet is on the convex hull, do not insert it in pre_star
        // however, the corresponding cell *cit is inserted in
        // pre_star_interior
      }
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      pre_star_interior.insert(*cit);
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
    }

    // every point in the mesh should be incident to a cell in domain:
    // if it is on a surface, at least one side of the surface should be in
    // domain.
    CGAL_assertion(is_incident_to_a_cell_in_domain);

    double first_radius_ratio = worst_radius_ratio;

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
    CGAL_assertion_code(typename Tr::size_type size_of_pre_star = pre_star.size());

    CGAL_assertion_code(Pre_star best_pre_star = pre_star);
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      
//     pre_star.pop_front();

//     std::cerr << "v=" << &*v << std::endl;
    
    while( !pre_star.empty() && 
           pre_star.front()->first < (sq_delta * sq_d_v) &&
           // La condition suivante teste que la facet qui va etre flippee
           // n'est pas contrainte.
           ! c2t3.is_in_complex(pre_star.front()->second) )
    {
        std::pair<double, Facet> facet = *(pre_star.front());
        Facet link = facet.second;

        // first critial radius
        double critical_r = facet.first;

// 	std::cerr << "critical_r(" << &*link.first
// 		  << "," << link.second
// 		  << ", s=" << pre_star.size() << ")="
// 		  << critical_r << " ";

        // link.first n'est pas contrainte donc on est dans le domaine.
//         CGAL_assertion( !tr.is_infinite(link.first) );
//         CGAL_assertion( link.first->is_in_domain() );
//         CGAL_assertion( ! link.first->
// 			is_facet_on_surface(link.second) );
 
        Facet opposite_facet = tr.mirror_facet(link);
        const Cell_handle& opposite_cell = opposite_facet.first;
        
//         CGAL_assertion( ! opposite_cell->
//                             is_facet_on_surface(opposite_facet.second) );
// 	CGAL_assertion( !tr.is_infinite(opposite_cell) );
//         CGAL_assertion(	opposite_cell->is_in_domain() );
 
	
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
	int number_of_erased_facets = 0;
	int number_of_new_facets = 0;
        Pre_star old_pre_star = pre_star;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
	int index_of_new_facet = -1;


        for(int i = 0; i<4 ; ++i)
          {
//             std::cerr << critical_r << " " 
//                       << power_product(v->point(),
//                                        opposite_cell->vertex(i)->point())
//                       << " " << sq_d_v << std::endl;
            
            const Facet new_facet = Facet(opposite_cell, i);
	    const Facet temp_facet = tr.mirror_facet(new_facet);
            if(!pre_star.erase(temp_facet))
	      {
		// in this case temp_facet is a new link facet
                // viewed from outside the prestar
                CGAL_assertion( opposite_cell->neighbor(i) != link.first);
                pre_star.insert(new_facet,
				compute_critical_radius(v,
							temp_facet.first));
		const Cell_handle& c = new_facet.first;
		const int& k = new_facet.second;
  
                CGAL_assertion( c != link.first );

// 		std::cerr << &*(c->vertex((k+1)&3)) << "\n"
// 			  << &*(c->vertex((k+2)&3)) << "\n"
// 			  << &*(c->vertex((k+3)&3)) << std::endl;
		ratios[new_facet] = 
		  CGAL::to_double(radius_ratio(v->point(),
                                               c->vertex((k+1)&3)->point(),
                                               c->vertex((k+2)&3)->point(),
                                               c->vertex((k+3)&3)->point(),
                                               Geom_traits()));
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
		++ number_of_new_facets;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
		index_of_new_facet = i;
	      }
	    else
            {
              if(temp_facet.first->is_in_domain())
              {
                
                CGAL_assertion_code(typename Ratios::size_type number_of_erased = )
  
                  ratios.erase(temp_facet);

                CGAL_assertion( number_of_erased == 1);
              }

		// in this case temp_facet is an old  link facet
                // viewed from inside the prestar
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
                ++number_of_erased_facets;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
	      }
          }

// 	std::cerr << "n(" << number_of_erased_facets << ")";
// 	std::cerr << &*(opposite_cell) << "\n";
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
	    std::cerr << "\nFlip " << number_of_erased_facets 
		      << "-" << number_of_new_facets << "\n"
                      << "wp=" << display_point_with_pre_star(v, old_pre_star) << "\n"
		      << "critical_r=" << critical_r << "\n"
		      << "power_product=" 
		      << power_product(
                           opposite_cell->vertex(index_of_new_facet)->point(),
			   v->point().point()) << "\n";
	CGAL_assertion((number_of_erased_facets == 1 && 
                        number_of_new_facets == 3 ) || // flip 2-3
                       (number_of_erased_facets == 2 && 
                        number_of_new_facets == 2 ));  // flip 3-2

        CGAL_assertion_code(size_of_pre_star = size_of_pre_star + number_of_new_facets - number_of_erased_facets);
        CGAL_assertion( pre_star.size() == size_of_pre_star );
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

//         // DOCUMENT THIS FUNCTOR
//         details::Map_value_of<Pre_star,
//           std::map<Facet, double> > value_of_ratios(ratios);

        

//         details::Value_of<std::map<Facet, double> > apply_ratios(&ratios);
//         details::First_of<Radius_radius_priority_queue> get_first;

        using details::Second_of;

	double min_of_pre_star = 
	  *(std::min_element(make_transform_iterator(ratios.begin(),
						     Second_of<Ratios>()),
			     make_transform_iterator(ratios.end(),
						     Second_of<Ratios>())));
        CGAL_assertion( min_of_pre_star > 0 );

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        std::set<Cell_handle> old_pre_star_interior = pre_star_interior;

        pre_star_interior.clear();

        const bool check_star = check_pre_star(pre_star,
                                               Weighted_point(v->point().point(),
                                                              (pre_star.empty() ? critical_r : (critical_r + pre_star.front()->first) / 2),
                                                              v->point().surface_index()),
                                               v,
                                               inserter(pre_star_interior));

        // first check that old_pre_star_interior in
        // included in pre_star_interior

        CGAL_assertion(std::includes(pre_star_interior.begin(), pre_star_interior.end(),
                                     old_pre_star_interior.begin(), old_pre_star_interior.end()));

        std::vector<Cell_handle> difference;
        std::set_difference(pre_star_interior.begin(), pre_star_interior.end(),
                            old_pre_star_interior.begin(), old_pre_star_interior.end(),
                            std::back_inserter(difference));
        
        CGAL_assertion(difference.size() == pre_star_interior.size() - old_pre_star_interior.size() );

        for(typename std::vector<Cell_handle>::const_iterator cit = difference.begin();
            cit != difference.end();
            ++cit)
        {
          std::cerr << "added cell " << &**cit << " " << (*cit)->is_in_domain();
          if( *cit != opposite_cell )
            std::cerr << " (!)";
          std::cerr << " : \n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(0), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(1), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(2), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(3), old_pre_star) << "\n";
        }
        if( ! check_star )
          std::cerr << "critical_r=" << critical_r 
                    << "\npre_star.front()->first=" << pre_star.front()->first
                    << "\n";
        std::cerr << "difference=" << pre_star_interior.size() - old_pre_star_interior.size() << "\n";
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        if( min_of_pre_star > worst_radius_ratio )
	  {
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
            CGAL_assertion_code(best_pre_star = pre_star;)
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER


	    worst_radius_ratio = min_of_pre_star;
	    if( !pre_star.empty() )
            {
//               CGAL_assertion( critical_r < pre_star.front()->first );
	      best_weight = (critical_r + pre_star.front()->first) / 2;
            }
	    else
            {
              CGAL_assertion(false);
              best_weight = critical_r;
            }
	  }


        facet = *(pre_star.front());
        link = facet.second;
        critical_r = facet.first;
// 	pre_star.pop_front();
      }

//    if( link.first->is_facet_on_surface(link.second) )
//       std::cerr << "X";   
    Bare_point p = v->point().point();
    int index = v->point().surface_index();
    //    tr.remove(v);
    if(best_weight > v->point().weight())
    {
      ++num_of_pumped_vertices;
//       std::cerr << "Best weight for point(" << p << ") is "
//           << best_weight << "\nand worst_radius_ratio is "
//           << worst_radius_ratio << " and was " << first_radius_ratio << std::endl;;
      
      total_pumping += (worst_radius_ratio - first_radius_ratio);
      
      Weighted_point wp = Weighted_point(p, best_weight, index);

      std::vector<Cell_handle> deleted_cells;
      std::vector<Facet> internal_facets;
      std::vector<Facet> boundary_facets;
      deleted_cells.reserve(64);
      internal_facets.reserve(64);
      boundary_facets.reserve(64);
      
      tr.find_conflicts(wp,
			v->cell(),
                        std::back_inserter(boundary_facets),
			std::back_inserter(deleted_cells),
			std::back_inserter(internal_facets));

     
      const bool seed_domain_marker = boundary_facets[0].first->is_in_domain();
      bool other_seed_domain_marker = seed_domain_marker;

      CGAL_assertion(!tr.is_infinite(boundary_facets[0].first) ||
		     !seed_domain_marker);
      

      // delete cells from the list of bad cells
      for(typename std::vector<Cell_handle>::iterator
	    cit = deleted_cells.begin(); 
	  cit != deleted_cells.end();
	  ++cit)
      {
        if((*cit)->is_in_domain() != seed_domain_marker)
          other_seed_domain_marker = ! seed_domain_marker;

        CGAL_assertion( v->point().surface_index() > 0 || (*cit)->is_in_domain() );

	cells_queue.erase(*cit);
      }

      const Facet one_boundary_facet_from_outside =
	tr.mirror_facet(boundary_facets[0]);
      // Here, one_boundary_facet_from_outside is a facet on the
      // boundary of the conflicts zone, seen from outside the
      // conflicts zone.

      Umbrella umbrella;

      {
	typename std::vector<Facet>::iterator ifit = internal_facets.begin();
	for ( ; ifit != internal_facets.end() ; ++ifit) {
          const int index = ifit->second;
          const Cell_handle& cell = ifit->first;
	  if ( cell->is_facet_on_surface(index) ) {
            // for each edges of all marked internal facets
            for(int i1 = 0; i1 < 4; ++i1)
              for(int i2 = i1 + 1; i2 < 4; ++i2)
                if( i1 != index && i2 != index )
                {
                  Vertex_handle v1 = cell->vertex(i1);
                  Vertex_handle v2 = cell->vertex(i2);
                  order_two_handles ( v1, v2 );
                  CGAL_assertion( v1 < v2 );
                  const std::pair<Vertex_handle, Vertex_handle> edge = 
                    std::make_pair(v1, v2);
                  typename Umbrella::iterator it = umbrella.find(edge);
                  if( it != umbrella.end() )
                    umbrella.erase(it);
                  else
                    // @TODO: remark: umbrella.find()/umbrella.insert() can
                    // be improved using lower_bound()!
                    // See the trick in include/CGAL/Double_map.h
                    umbrella.insert(edge);
                }
            // remove already present edges, so that only external
            // edges of the umbrella are in umbrella at the end
	  }
	}
      }

      // The next two assertions mean
      //    index == 0 <==> umbrella.empty()
      CGAL_assertion( index > 0 || umbrella.empty());

      CGAL_assertion(check_umbrella(umbrella));

     //  std::set<Facet> deja_vu_cells;
//       for(typename std::vector<Cell_handle>::iterator
// 	    cit = deleted_cells.begin();
// 	  cit != deleted_cells.end();
// 	  ++cit)
// 	{
// 	  cells_queue.erase(*cit);
// 	  for(int i = 0; i < 4; ++i)
// 	    if( boundary_facets_of_conflicts_zone.find(Facet(*cit, i)) ==
//                   boundary_facets_of_conflicts_zone.end()
//                   // (*cit,i) is not on the boundary of the zone
//                 &&
// 		deja_vu_cells.find(tr.mirror_facet(Facet(*cit, i))) ==
// 		  deja_vu_cells.end() )
// 	      {
// 		deja_vu_cells.insert(Facet(*cit, i));
// 		if( (*cit)->is_facet_on_surface(i) )
// 		  {
// 		    int index_of_v_in_cit = (*cit)->index(v);
// 		    int i1;
// 		    int i2;
// 		    find_the_two_other_indices(i, index_of_v_in_cit,
// 					       i1, i2);
// 		    Vertex_handle v1 = (*cit)->vertex(i1);
// 		    Vertex_handle v2 = (*cit)->vertex(i2);
// 		    if( v2 < v1 ) std::swap(v1, v2);
// 		    umbrella.insert(std::make_pair(v1, v2));
// 		  }
// 	      }
// 	}
//      }


      Vertex_handle new_vertex = tr.insert(wp);

      std::vector<Cell_handle> new_cells;
      new_cells.reserve(64);

      tr.incident_cells(new_vertex, std::back_inserter(new_cells));

      /* Now we will restore markers in cells:     *
       *    - the ones that mark in-domain cells   *
       *    - the ones that mark on-surface facets *
       *                                           */
     //  // first, we start with the on-surface markers of the facets of
//       // the boundary of the star
//       for(typename std::vector<Cell_handle>::iterator cit = new_cells.begin();
// 	  cit != new_cells.end();
// 	  ++cit)
// 	if( !tr.is_infinite(*cit) &&
// 	    (*cit)->is_in_domain() )
// 	{ // ** re-marks facets on the boundary of the star **
// 	  const int index = (*cit)->index(new_vertex);
// 	  const Facet f = Facet(*cit, index);
// 	  const Facet opposite_facet = tr.mirror_facet(f);
// 	  bool marker = 
// 	    opposite_facet.first->is_facet_on_surface(opposite_facet.second);
// 	  (*cit)->set_surface_facet(index, marker);
// 	  // ** re-inserts cells of the star in cells_queue **
// 	  cells_queue.insert(*cit, radius_ratio(tr.tetrahedron(*cit)));
// 	}

      // then we walk in the star, to check for internal facets that are on 
      // the surface, using "umbrella".

      // before the walk, restores the in_conflict_flag of cells
      for(typename std::vector<Cell_handle>::iterator cit = new_cells.begin();
	  cit != new_cells.end(); ++cit)
      {
	(*cit)->set_in_conflict_flag(0);
      }

      walk_on_the_boundary_of_the_star(
        tr.mirror_facet(one_boundary_facet_from_outside),
	seed_domain_marker,
        other_seed_domain_marker,
	umbrella,
	new_vertex);

      // after the walk, restores the in_conflict_flag of cells
      for(typename std::vector<Cell_handle>::iterator cit = new_cells.begin();
	  cit != new_cells.end(); ++cit)
      {
        CGAL_assertion((*cit)->get_in_conflict_flag()==1);
        (*cit)->set_in_conflict_flag(0);
        CGAL_assertion( v->point().surface_index() > 0 || (*cit)->is_in_domain() );
      }

      //CGAL_assertion (umbrella.empty());

      return true; // The vertex v has been pumped into new_vertex.
    }
    else
      {
	++num_of_ignored_vertices;
	return false;
      }
  } // end pump_vertex

  void walk_on_the_boundary_of_the_star(Facet f,
					bool in_domain_marker,
					bool in_domain_marker_2,
					Umbrella& umbrella,
					Vertex_handle v)
  {
    Cell_handle& c = f.first;
    const int& index = f.second;

    CGAL_assertion(!tr.is_infinite(c) || !in_domain_marker);

    CGAL_assertion(v == c->vertex(index));
    
    if (c->get_in_conflict_flag() != 0) 
      return;

    c->set_in_conflict_flag(1);
    c->set_in_domain(in_domain_marker);

    // ** re-inserts cells of the star in cells_queue **
    if (in_domain_marker)
    {
      const double ratio = CGAL::to_double(radius_ratio(tr.tetrahedron(c),
                                                        Geom_traits()));
      if( ratio < stop_limit_on_radius_ratio )
        cells_queue.insert(c, ratio);
    }

    for (int i=0; i<4; ++i) {
      if (i == index) {
	// ** re-marks facets on the boundary of the star **
	const Facet f = Facet(c, i);
	const Facet opposite_facet = tr.mirror_facet(f);
	bool marker = 
	  opposite_facet.first->is_facet_on_surface(opposite_facet.second);
	c->set_facet_on_surface(i, marker);
      }
      else {

	Cell_handle next_c = c->neighbor(i);

	int next_index = next_c->index(v);
      
	int i1, i2;
	find_the_two_other_indices(i, index, i1, i2);
	Vertex_handle v1 = c->vertex(i1);
	Vertex_handle v2 = c->vertex(i2);
	order_two_handles(v1, v2);
	CGAL_assertion(v1 < v2);
	if( umbrella.find(std::make_pair(v1, v2)) !=
	    umbrella.end() ) // the facet (c, i) is on the surface
	  {
	    c->set_facet_on_surface(i, true);
	    next_c->set_facet_on_surface(next_c->index(c), true);
	    walk_on_the_boundary_of_the_star(Facet(next_c, next_index),
					     in_domain_marker_2,
                                             in_domain_marker,
					     umbrella, v);
	  }
	else
	  {
	    CGAL_assertion(  umbrella.find(std::make_pair(v2, v1)) ==
			     umbrella.end() );
	    walk_on_the_boundary_of_the_star(Facet(next_c, next_index),
					     in_domain_marker,
                                             in_domain_marker_2,
					     umbrella, v);
	  }
      }
    }
  } // end walk_on_the_boundary_of_the_star
  
private: // data
  C2T3 c2t3;
  Tr& tr;
  double sq_delta;
  double stop_limit_on_radius_ratio;

  int num_of_pumped_vertices;
  int num_of_ignored_vertices;
  double total_pumping;

  bool initialized;
    
  typedef CGAL::Double_map<Cell_handle, double> Radius_radius_priority_queue;
  
  Radius_radius_priority_queue cells_queue;

private: // methods
  void compute_aspect_ratio_priority_map()
  {
    for(Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end();
        ++cit)
      if(cit->is_in_domain())
      {
        const double ratio = CGAL::to_double(radius_ratio(tr.tetrahedron(cit),
                                                          Geom_traits()));
        if( ratio < stop_limit_on_radius_ratio )
          cells_queue.insert(cit, ratio);
      }
  }

  typename Geom_traits::FT
  power_product(const Weighted_point& p1, const Weighted_point& p2) const
  {
    typedef typename Geom_traits::Compute_power_product_3
        Compute_power_product_3;
    Compute_power_product_3 product =
        tr.geom_traits().compute_power_product_3_object();
    return product(p1, p2);
  }

  double compute_critical_radius(const Vertex_handle& v,
				 const Cell_handle& c) const
  {
    typedef typename Geom_traits::Compute_critical_squared_radius_3
      Critical_radius;
    Critical_radius critical_radius = 
      tr.geom_traits().compute_critical_squared_radius_3_object();
    
    const FT result = critical_radius(c->vertex(0)->point(),
                                      c->vertex(1)->point(),
                                      c->vertex(2)->point(),
                                      c->vertex(3)->point(),
                                      v->point());
//     CGAL_assertion_code(
//       typedef typename Geom_traits::Construct_weighted_circumcenter_3
//         Construct_weighted_circumcenter_3;
//       typedef typename Geom_traits::Compute_power_product_3
//         Compute_power_product_3;
//       Compute_power_product_3 power_product =
//         tr.geom_traits().compute_power_product_3_object();
//       Construct_weighted_circumcenter_3 weighted_circumcenter =
//         tr.geom_traits().construct_weighted_circumcenter_3_object();

//       Weighted_point wc(weighted_circumcenter(c->vertex(0)->point(),
//                                               c->vertex(1)->point(),
//                                               c->vertex(2)->point(),
//                                               c->vertex(3)->point()))
//       );
    
//     CGAL_assertion( result == power_product(v->point(), wc) );
    
    return CGAL::to_double(result);
  }

//   void restore_markers_around(Vertex_handle v)
//   {
//     std::vector<Cell_handle> incident_cells;
//     incident_cells.reserve(64);
//     tr.incident_cells(v, std::back_inserter(incident_cells));

// //     std::cerr << "S(" << incident_cells.size() << ")";
//     bool test = false;

//     for(typename std::vector<Cell_handle>::const_iterator cit = 
//           incident_cells.begin();
//         cit != incident_cells.end();
//         ++cit)
//       {	
// 	const int index = (*cit)->index(v);
// 	const Facet f = Facet(*cit, index);
// 	const Facet opposite_facet = tr.mirror_facet(f);
// 	bool marker = 
// 	  opposite_facet.first->is_facet_on_surface(opposite_facet.second);
// 	(*cit)->set_facet_on_surface(index, marker);
// 	test = test || marker;
// 	(*cit)->set_in_domain(true);
//       }
// //     if( test ) 
// //       std::cerr << "Y";
//   }
}; // end class Slivers_exuder

} // end namespace Mesh_3

template < class Tr>
void
output_slivers_to_off (std::ostream& os, const Tr & tr, 
		       const double sliver_bound = 0.25) 
{
  typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Point Point;
  
  typedef std::set<Cell_handle> Slivers;
  Slivers slivers;
  
  for( Finite_cells_iterator cit = tr.finite_cells_begin(); 
       cit != tr.finite_cells_end(); ++cit)
  {
    if( cit->is_in_domain() &&
	CGAL::to_double(Mesh_3::radius_ratio(tr.tetrahedron(cit),
                                             typename Tr::Geom_traits()))
        < sliver_bound )
      slivers.insert(cit);
  }
  
  // Header.
  os << "OFF \n"
      << tr.number_of_vertices() << " " <<
      slivers.size() * 4 << 
      " " << 0 << "\n";

  os << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
  {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }
  
  int surface_slivers = 0;
  int flat_slivers = 0;
  
  // Finite cells indices.
  for( typename Slivers::iterator cit = slivers.begin();
       cit != slivers.end(); ++cit)
  {
    const Cell_handle& c = *cit;
    
    if( c->vertex(0)->point().surface_index() > 0 &&
        c->vertex(1)->point().surface_index() > 0 &&
        c->vertex(2)->point().surface_index() > 0 &&
        c->vertex(3)->point().surface_index() > 0 )
      ++flat_slivers;
    
    if( c->vertex(0)->point().surface_index() > 0 ||
        c->vertex(1)->point().surface_index() > 0 ||
        c->vertex(2)->point().surface_index() > 0 ||
        c->vertex(3)->point().surface_index() > 0 )
      ++surface_slivers;
    
    for(int i = 0; i < 4; ++i)
      {
        os << "3 ";
        for (int k=0; k<3; k++)
          os << V[c->vertex((i+k)&3)] << " ";
        
	os << "\n"; // without color.
      }
   }
   std::cerr << "Number of slivers: " << slivers.size()
             << "\nNumber of surface slivers: " << surface_slivers
             << "\nNumber of flat slivers: " << flat_slivers << std::endl;
}
} // end namespace CGAL


#endif // end CGAL_MESH_3_SLIVERS_EXUDER_H
