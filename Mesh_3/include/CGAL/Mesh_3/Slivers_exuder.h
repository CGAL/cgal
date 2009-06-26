// Copyright (c) 2004-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
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
// Author(s)     : Laurent Rineau

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
#include <boost/format.hpp>

#include <CGAL/radius_ratio.h>

#include <CGAL/Mesh_3/Slivers_exuder_cell_attributes_traits.h>

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

    // That functor Second_of takes a pair as input (the value type of a
    // map), and returns the ".second" member of that pair. It is used in
    // Sliver_exuder, to constructor a transform iterator.

    // It should be doable using STL bind operators, but i am not sure how
    // to use them. -- Laurent Rineau, 2007/07/27
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

    // That function is constructed with a vertex handle v1.
    // Then, its operator() takes an other vertex handle v2 as input, and
    // returns the distance d(v1, v2).
    // It is used in Sliver_exuder, to constructor a transform iterator.
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

  
  /** Pre_star will represent the pre-star of a point. It is a (double)-map
   *  of Facet (viewed from cells inside the star), ordered by the
   *  critial_radius of the point with the cell that lies on the facet, at
   *  the exterior of the pre-star. */
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

//   static 
//   void find_the_two_other_indices(int i1, int i2,
// 				  int& i3, int& i4)
//     // given i1 and i2 (i1<>i2), returns i3 and i4 so that
//     // i3 < i4 and i1+i2+i3+i4=6
//   {
//     CGAL_precondition( 0 <= i1 && i1 <= 3 );
//     CGAL_precondition( 0 <= i2 && i2 <= 3 );
//     CGAL_precondition( i1 != i2 );
    
//     for(i3 = 0; i3 < 4; ++i3)
//       if( i3 != i1 && i3 != i2 )
// 	break;
//     i4 = 6 - i1 - i2 - i3;

//     CGAL_postcondition( 0 <= i3 && i3 < i4 && i4 <= 3 );
//     CGAL_postcondition( i3 != i1 && i3 != i2 );
//     CGAL_postcondition( i4 != i1 && i4 != i2 );
//   }


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
    if(v != Vertex_handle())
    {
      Pre_star pre_star2;

      // fill pre_star2
      for(Input_facet_it fit = begin;
          fit != end;
          ++fit)
      {
        const Facet opposite_facet = tr.mirror_facet(*fit);
        if(! tr.is_infinite(opposite_facet.first) ) 
        {
          pre_star2.insert(*fit, compute_critical_radius(v,
                                                         opposite_facet.first));
        }
      }

      while(!pre_star_copy.empty() && !pre_star2.empty())
      {
        if(pre_star_copy.front()->first != pre_star2.front()->first) {
          std::cerr << "bad order\n";
          std::cerr << boost::format("pre_star.front()->first=%1%, should be %2%\n")
            % pre_star_copy.front()->first % pre_star2.front()->first;
          return false;
        }
        pre_star2.pop_front();
        pre_star_copy.pop_front();
      }

      if(pre_star2.empty() && ! pre_star_copy.empty()) {
        std::cerr << "pre_star is too big!\n";
        while(!pre_star_copy.empty())
        {
          const Facet f = pre_star_copy.front()->second;
          const double r = pre_star_copy.front()->first;
          pre_star_copy.pop_front();
          std::cerr << boost::format("extra facet (%1%,%2%) (infinite: %3%, opposite infinite: %4%), critical radius: %5%\n")
            % &*f.first % f.second % tr.is_infinite(f.first) % tr.is_infinite(f.first->neighbor(f.second))
            % r;
        }
        return false;
      }
      
      if( pre_star_copy.empty() && ! pre_star2.empty() ) {
        std::cerr << "pre_star is too small!\n";
        while(!pre_star2.empty())
        {
          const Facet f = pre_star2.front()->second;
          pre_star2.pop_front();
          std::cerr << boost::format("missing facet (%1%,%2%) (infinite: %3%, opposite infinite: %4%)\n")
            % &*f.first % f.second % tr.is_infinite(f.first) % tr.is_infinite(f.first->neighbor(f.second));
        }
        return false;
      }
    }

    pre_star_copy = pre_star;

    for(Input_facet_it fit = begin;
        fit != end;
        ++fit)
    {
      const Facet opposite_facet = tr.mirror_facet(*fit);
      if(!tr.is_infinite(opposite_facet.first) && !pre_star_copy.erase(*fit))
        return false;
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
    std::vector<Facet> boundary_facets;
    boundary_facets.reserve(64);
    
    tr.find_conflicts(wp,
                      vh->cell(),
                      std::back_inserter(boundary_facets),
                      cells_out,
                      CGAL::Emptyset_iterator());

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
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
    CGAL_assertion(c2t3.is_valid());
 
    std::cerr << std::boolalpha;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
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
    pump_vertices<true>(radius_ratio_limit);
  }

  template <
    bool pump_vertices_on_surfaces
    >
  void pump_vertices(double radius_ratio_limit = 1.)
  {
    if( ! initialized )
      init();

    // store radius_ratio_limit in the member stop_limit_on_radius_ratio
    stop_limit_on_radius_ratio = radius_ratio_limit;

    while( !cells_queue.empty() )
    {
      typename Radius_radius_priority_queue::Reverse_entry 
	front = *(cells_queue.front());
      Cell_handle c = front.second;

      int i;
      for( i = 0; i < 4; ++i )
      {
        // pump_vertices_on_surfaces is a boolean template parameter.  The
        // following condition is pruned at compiled time, if
        // pump_vertices_on_surfaces==false.
        if( pump_vertices_on_surfaces || 
	    c->vertex(i)->point().surface_index() == 0 )
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
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
    std::cerr << 
      boost::format("Pumping vertex %1% (point: %2%, on surface: %3%)\n")
      % &*v % v->point() % (v->point().surface_index() != 0);
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
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

    // That variable is a counter of the facet that have not been added to
    // the pre_start, because they lie in the convex hull.
    int number_of_facets_on_convex_hull = 0;

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
      else
        ++ number_of_facets_on_convex_hull;
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
    CGAL_assertion_code(typename Tr::size_type size_of_pre_star = 
                          pre_star.size() + number_of_facets_on_convex_hull);

    CGAL_assertion_code(Pre_star best_pre_star = pre_star);
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      
//     pre_star.pop_front();

//     std::cerr << "v=" << &*v << std::endl;
    
    bool can_flip = true; // If that boolean is set to false, it means that
                          // a facet in the complex is about to be
                          // flipped. In that case, the pumping is stopped.
    while( can_flip &&
           !pre_star.empty() && 
           pre_star.front()->first < (sq_delta * sq_d_v) &&
           // La condition suivante teste que la facet qui va etre flippee
           // n'est pas contrainte.
           ! c2t3.is_in_complex(pre_star.front()->second) )
    {
      typename Pre_star::reverse_iterator facet_handle = pre_star.front();

      Facet link = facet_handle->second;

      // first critial radius
      double critical_r = facet_handle->first;
 
      Facet opposite_facet = tr.mirror_facet(link);
      const Cell_handle& opposite_cell = opposite_facet.first;
        	
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      int number_of_erased_facets = 0;
      int number_of_new_facets = 0;
      Pre_star old_pre_star = pre_star;
      int index_of_new_facet = -1;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

      can_flip = true;
      for(int i = 0; i<4 ; ++i)
      {
        const Facet new_facet = Facet(opposite_cell, i);
        const Facet temp_facet = tr.mirror_facet(new_facet);
        if(!pre_star.erase(temp_facet))
        {
          if(! tr.is_infinite(temp_facet.first) )
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

            if(opposite_cell->is_in_domain())
            {
              ratios[new_facet] = 
                CGAL::to_double(radius_ratio(v->point(),
                                             c->vertex((k+1)&3)->point(),
                                             c->vertex((k+2)&3)->point(),
                                             c->vertex((k+3)&3)->point(),
                                             Geom_traits()));
            }
          }
          else
            ++number_of_facets_on_convex_hull;
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
          ++ number_of_new_facets;
          index_of_new_facet = i;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        }
        else
        {
          if(c2t3.is_in_complex(temp_facet)) {
            can_flip = false;
          }
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

      if(can_flip) 
      {
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        std::cerr << "\nFlip " << number_of_erased_facets 
                  << "-" << number_of_new_facets << "\n"
                  << "wp=" << display_point_with_pre_star(v, old_pre_star) << "\n"
                  << "critical_r=" << critical_r << "\n"
                  << "pre_star.front()->first=" << pre_star.front()->first
                  << "\npower_product=" 
                  << power_product(
                                   opposite_cell->vertex(index_of_new_facet)->point(),
                                   v->point().point()) << "\n";
        CGAL_assertion((number_of_erased_facets == 1 && 
                        number_of_new_facets == 3 ) || // flip 2-3
                       (number_of_erased_facets == 2 && 
                        number_of_new_facets == 2 ));  // flip 3-2

        CGAL_assertion_code(size_of_pre_star = size_of_pre_star + number_of_new_facets - number_of_erased_facets);
        CGAL_assertion( pre_star.size() + number_of_facets_on_convex_hull == size_of_pre_star );
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

        using details::Second_of;

        double min_of_pre_star = 
          *(std::min_element(make_transform_iterator(ratios.begin(),
                                                     Second_of<Ratios>()),
                             make_transform_iterator(ratios.end(),
                                                     Second_of<Ratios>())));
        //         CGAL_assertion( min_of_pre_star > 0 );

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        std::set<Cell_handle> old_pre_star_interior = pre_star_interior;

        pre_star_interior.clear();

        const bool check_star = check_pre_star(pre_star,
                                               Weighted_point(v->point().point(),
                                                              (pre_star.empty() ? critical_r : (critical_r + pre_star.front()->first) / 2),
                                                              v->point().surface_index()),
                                               v,
                                               inserter(pre_star_interior));
        if( ! check_star )
          std::cerr << "critical_r=" << critical_r 
                    << "\npre_star.front()->first=" << pre_star.front()->first
                    << "\n";
        CGAL_assertion(check_star);
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
          std::cerr << "added cell " << &**cit << " "
                    << ( (*cit)->is_in_domain() ? "(in domain)" : "" );
          if( *cit != opposite_cell )
            std::cerr << " (!)";
          std::cerr << " : \n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(0), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(1), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(2), old_pre_star) << "\n"
                    << "  " << display_point_with_pre_star((*cit)->vertex(3), old_pre_star) << "\n";
        }
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
            best_weight = (critical_r + pre_star.front()->first) / 2;
          }
          else
          {
            CGAL_assertion(false);
            best_weight = critical_r;
          }
        }


        facet_handle = pre_star.front();
        link = facet_handle->second;
        critical_r = facet_handle->first;
      } // end if(can_flip)
    } // end while(... can pump...)

    Bare_point p = v->point().point();
    int index = v->point().surface_index();
    //    tr.remove(v);
    if(best_weight > v->point().weight())
    {
      ++num_of_pumped_vertices;
      
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


      typedef Mesh_3::Slivers_exuder_cell_attributes_traits<typename Tr::Cell>
        Cell_attributes_traits;
      typedef typename Cell_attributes_traits::Cell_attributes
        Cell_attributes;
      Cell_attributes_traits attributes_traits;
      // That map boundary_facets_from_outside stores the facets (c, i) on
      // the boundary of the conflict zone, so that the cell handles c are
      // at the exterior of the conflict zone.
      // It also stores if the cell c->neighbor(i) is in the domain.
      typedef std::map<Facet, Cell_attributes> Boundary_facets_from_outside;
      Boundary_facets_from_outside boundary_facets_from_outside;

      for(typename std::vector<Facet>::const_iterator fit =
	    boundary_facets.begin(), end = boundary_facets.end();
	  fit != end;
	  ++fit) 
      {
	boundary_facets_from_outside.
	  insert(std::make_pair(tr.mirror_facet(*fit),
            attributes_traits.get_attributes(&*fit->first)));
      }

      // delete cells from the list of bad cells
      for(typename std::vector<Cell_handle>::iterator
	    cit = deleted_cells.begin(); 
	  cit != deleted_cells.end();
	  ++cit)
      {
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
            CGAL_assertion(tr.has_vertex(*ifit, v));
            // for each edges of all marked internal facets
            for(int i1 = 0; i1 < 4; ++i1)
              for(int i2 = i1 + 1; i2 < 4; ++i2)
                if( i1 != index && i2 != index)
                {
                  Vertex_handle v1 = cell->vertex(i1);
                  Vertex_handle v2 = cell->vertex(i2);
                  if(v1 == v || v2 == v) continue;
                  order_two_handles ( v1, v2 );
                  CGAL_assertion( v1 < v2 );
                  const std::pair<Vertex_handle, Vertex_handle> edge = 
                    std::make_pair(v1, v2);
//                   typename Umbrella::iterator it = umbrella.find(edge);
//                   if( it != umbrella.end() )
//                     umbrella.erase(it);
//                   else
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

      // The next assertion means
      //    index == 0 ==> umbrella.empty()
//       CGAL_assertion( index > 0 || umbrella.empty());

      // The following assertion cannot be checked, if the surface coded by
      // c2t3 is not manifold.

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      CGAL_assertion(c2t3.is_valid());

      int i,j;
      c2t3.union_find_of_incident_facets(v,i,j);

      std::cerr << boost::format("Vertex (%1%, %2%) is about to be pumped.\n"
                                 "  status: %3%\n"
                                 "  number_of_indicent_facets: %4%\n"
                                 "  number_of_components: %5%\n"
                                 "  size of umbrella: %6% (should be equal to %4%)\n")
        % &*v % v->point()
        % c2t3.face_status(v)
        % i
        % j
        % umbrella.size();
#endif

      Vertex_handle new_vertex = tr.insert(wp);

      std::vector<Cell_handle> new_cells;
      new_cells.reserve(64);

      tr.incident_cells(new_vertex, std::back_inserter(new_cells));

      // We walk in the star, to check for facets that are on 
      // the surface, using "umbrella".

      // Before the walk, reset the in_conflict_flag of cells (to 0),
      // restore the boolen is_in_domain(), and update the "bad cells"
      // (badly shaped) priority queue.
      for(typename std::vector<Cell_handle>::iterator cit = new_cells.begin();
	  cit != new_cells.end(); ++cit)
      {
	(*cit)->set_in_conflict_flag(0);
	const int index = (*cit)->index(new_vertex);
	const Facet f = std::make_pair(*cit, index); // the facet opposite to
						     // new_vertex, in *cit

	// Search the mirror facet of f in boundary_facets_from_outside.
	// That search cannot fail.
	const typename Boundary_facets_from_outside::const_iterator it =
	  boundary_facets_from_outside.find(tr.mirror_facet(f));
	CGAL_assertion(it != boundary_facets_from_outside.end());

	// And use the map entry to restore (*cit)'s attributes.
	attributes_traits.restore_attributes(&**cit, it->second);

	if((*cit)->is_in_domain()) {
	  // if the new cell is in the domain, and it ratio is less that
	  // the maximum, push it in the cells queue.
	  const double ratio = CGAL::to_double(radius_ratio(tr.tetrahedron(*cit),
							    Geom_traits()));
	  if( ratio < stop_limit_on_radius_ratio )
	    cells_queue.insert(*cit, ratio);
	}
      }

      walk_on_the_boundary_of_the_star(
        tr.mirror_facet(one_boundary_facet_from_outside),
	umbrella,
	new_vertex);

      // after the walk, reset the in_conflict_flag of cells (to 0).
      for(typename std::vector<Cell_handle>::iterator cit = new_cells.begin();
	  cit != new_cells.end(); ++cit)
      {
        CGAL_assertion((*cit)->get_in_conflict_flag()==1);
        (*cit)->set_in_conflict_flag(0);
        CGAL_assertion( new_vertex->point().surface_index() > 0 || (*cit)->is_in_domain() );
      }

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      c2t3.union_find_of_incident_facets(new_vertex,i,j);

      std::cerr << boost::format("Vertex (%1%, %2%) is has been pumped.\n"
                                 "  status: %3%\n"
                                 "  number_of_indicent_facets: %4%\n"
                                 "  number_of_components: %5%\n")
        % &*new_vertex % new_vertex->point()
        % c2t3.face_status(new_vertex)
        % i
        % j;

      CGAL_assertion(c2t3.is_valid(true));
#endif

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
					Umbrella& umbrella,
					Vertex_handle v)
  {
    Cell_handle& c = f.first;
    const int& index = f.second;

    CGAL_assertion(v == c->vertex(index));
    
    if (c->get_in_conflict_flag() != 0) 
      return;

    c->set_in_conflict_flag(1);

    for (int i=0; i<4; ++i) {
      if (i == index) {
	// ** re-marks facets on the boundary of the star **
	const Facet f = Facet(c, i);
	const Facet opposite_facet = tr.mirror_facet(f);
	const bool marker = 
	  opposite_facet.first->is_facet_on_surface(opposite_facet.second);
	c->set_facet_on_surface(i, marker);
      }
      else {
	const Cell_handle& next_c = c->neighbor(i);

	const int next_index = next_c->index(v);
      
	const int i1 = tr.next_around_edge(i, index);
	const int i2 = tr.next_around_edge(index, i);
	Vertex_handle v1 = c->vertex(i1);
	Vertex_handle v2 = c->vertex(i2);
	order_two_handles(v1, v2);
	CGAL_assertion(v1 < v2);

        const typename Umbrella::iterator um_it = 
          umbrella.find(std::make_pair(v1, v2));
	if( um_it != umbrella.end() ) // the facet (c, i) is on the surface
	{
          c->set_facet_on_surface(i, true);
//           next_c->set_facet_on_surface(next_c->index(c), true);
	}
        else 
          c->set_facet_on_surface(i, false);
	walk_on_the_boundary_of_the_star(Facet(next_c, next_index),
					 umbrella, v);
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
    return CGAL::to_double(result);
  }

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
