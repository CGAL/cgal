// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_SLIVER_EXUDER_H
#define CGAL_MESH_3_SLIVER_EXUDER_H

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Double_map.h>
#include <map>
#include <vector>
#include <boost/iterator/transform_iterator.hpp>
// #include <boost/bind.hpp>

#ifndef PI
#define PI 3.1415926535897932
#endif

template <typename K>
double
area(const CGAL::Tetrahedron_3<K>& t)
{// This function is (c) 2005 Pierre Alliez
  typedef typename K::Triangle_3 Triangle;
    
  Triangle t1 = Triangle(t[0],t[1],t[2]);
  Triangle t2 = Triangle(t[0],t[1],t[3]);
  Triangle t3 = Triangle(t[2],t[1],t[3]);
  Triangle t4 = Triangle(t[2],t[0],t[3]);
  double a1 = std::sqrt(CGAL_NTS to_double(t1.squared_area()));
  double a2 = std::sqrt(CGAL_NTS to_double(t2.squared_area()));
  double a3 = std::sqrt(CGAL_NTS to_double(t3.squared_area()));
  double a4 = std::sqrt(CGAL_NTS to_double(t4.squared_area()));
  return a1 + a2 + a3 + a4;
}

template <typename K>
double
circumradius(const CGAL::Tetrahedron_3<K>& t)
{ // This function is (c) 2005 Pierre Alliez
  typename K::Point_3 center = circumcenter( t.vertex(0),
                                             t.vertex(1),
                                             t.vertex(2),
                                             t.vertex(3));
  return CGAL::sqrt(CGAL::to_double( squared_distance( center, t.vertex(0))));
}

template <typename K>
double
len(const CGAL::Vector_3<K> &v)
{ // This function is (c) 2005 Pierre Alliez
  return std::sqrt(CGAL::to_double(v*v));
}

double fix_sine(double sine)
{ // This function is (c) 2005 Pierre Alliezx
  if(sine >= 1)
    return 1;
  else
    if(sine <= -1)
      return -1;
    else
      return sine;
}

template <typename K>
double 
angle_rad(const CGAL::Vector_3<K> &u,
          const CGAL::Vector_3<K> &v)
{ // This function is (c) 2005 Pierre Alliez
  // check
  double product = len(u)*len(v);
  if(product == 0.)
    return 0.0;

  // cosine
  double dot = CGAL::to_double(u*v);
  double cosine = dot / product;

  // sine
  typename K::Vector_3 w = CGAL::cross_product(u,v);
  double AbsSine = len(w) / product;

  if(cosine >= 0)
    return std::asin(fix_sine(AbsSine));
  else
    return PI-std::asin(fix_sine(AbsSine));
}

template <typename K>
double angle_deg(const CGAL::Vector_3<K> &u,
                 const CGAL::Vector_3<K> &v)
{
  static const double conv = 1.0/PI*180.0;

  return conv*angle_rad(u,v);
}

template <typename K>
typename CGAL::Vector_3<K>
normal(const CGAL::Point_3<K>& a,
       const CGAL::Point_3<K>& b,
       const CGAL::Point_3<K>& c)
{
  return CGAL::cross_product(b-a,c-a);
}

// dihedral angle at an edge [vi,vj]
template <typename K>
double dihedral_angle(const CGAL::Tetrahedron_3<K>& t,
                      const int i,
                      const int j)
{
  const CGAL::Vector_3<K> vi = normal(t[(i+1)&3],
                                      t[(i+2)&3],
                                      t[(i+3)&3]);
  const CGAL::Vector_3<K> vj = normal(t[(j+1)&3],
                                      t[(j+2)&3],
                                      t[(j+3)&3]);
  return 180-angle_deg(vi,vj);
}

template <typename K>
double
radius_ratio(const typename CGAL::Tetrahedron_3<K>& t)
{ // This function is (c) 2005 Pierre Alliez
  typename K::FT inradius = 3 * std::abs(t.volume()) / area(t);
  double circ = circumradius(t);
  if(circ == 0)
    return 0;
  else
    return (3 * CGAL::to_double(inradius) / circ);
}

namespace CGAL {

namespace Mesh_3 {

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
      {
	return (*m)[p.first];
      }
    }; // end class Map_value_of

    template <typename Map>
    struct First_of : 
      std::unary_function<typename Map::Direct_entry,
                          const typename Map::Key&>
    {
      typedef std::unary_function<typename Map::Direct_entry,
				  const typename Map::Key&> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;
      
      const typename Map::Key&
      operator()(const typename Map::Direct_entry& p) const
      {
	return p.first;
      }
    }; // end class First_of

    template <typename Map>
    class Value_of :
      std::unary_function<typename Map::key_type,
			  typename Map::mapped_type>
    {
      Map& m;
    public:
      typedef std::unary_function<typename Map::key_type,
				  typename Map::mapped_type&> Base;
      typedef typename Base::result_type result_type;
      typedef typename Base::argument_type argument_type;
      
      Value_of(Map& map) : m(map) 
      {
      }
      
      typename Map::mapped_type&
      operator()(const typename Map::key_type& k)
      {
	return m[k];
      }
    }; // end class Value_of

    template <typename Gt, typename Vertex_handle>
    class Power_from_v :
      public std::unary_function<Vertex_handle, typename Gt::RT>
    {
      Vertex_handle * v;
      Gt gt;
    public:
      Power_from_v(Vertex_handle& vh,
		   Gt geom_traits = Gt()) 
	: v(&vh), gt(geom_traits)
      {
      }
      
      typename Gt::RT
      operator()(const Vertex_handle& vh) const
      {
        typedef typename Gt::Compute_power_product_3
            Compute_power_product_3;
        Compute_power_product_3 product =
            gt.compute_power_product_3_object();

	return product((*v)->point(), vh->point());
      }
    }; // end class Power_from_v

  } // end namespace details

template <typename Tr, 
          typename FT = typename Tr::Geom_traits::FT >
class Slivers_exuder
{
public: //types
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Vertex_handle Vertex_handle;

  typedef typename Tr::Finite_cells_iterator Finite_cell_iterator;
  typedef typename Weighted_point::Weight Weight;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Geom_traits::Tetrahedron_3 Tetrahedron_3;
  
public: // methods
  Slivers_exuder(Tr& t, double d = 0.45)
    : tr(t), sq_delta(d*d), num_of_pumped_vertices(0),
      num_of_ignored_vertices(0), total_pumping(0.0)
  {
  }

  void init() 
  {
  }
  
  void pump_vertices()
  {
    num_of_pumped_vertices = 0;
    num_of_ignored_vertices = 0;
    total_pumping = 0.0;
    
    for(typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
	vit != tr.finite_vertices_end();
	++vit)
      if(vit->info() != true)
	pump_vertex(vit);
    
    std::cerr << "Pumped vertices:  " << num_of_pumped_vertices << std::endl;
    std::cerr << "Ignored vertices: " << num_of_ignored_vertices << std::endl;
    std::cerr << "Average pumping:  " << (total_pumping / num_of_pumped_vertices) << std::endl;
  }
//   void init()
//   {
//     cell_queue.clear();
//     compute_aspect_ratio_priority_map();
//   }

//   void pump_vertices()
//   {
//     typename Cell_radius_priority_queue::iterator it = cell_queue.begin();
//     while ( it != cell_queue.end() )
//       {
//         My_cell cell = it->second;
//         int i;
//         Cell_handle ch;
//         if( tr.is_cell(cell.V[0],
//                        cell.V[1],
//                        cell.V[2],
//                        cell.V[3],
//                        ch, i, i, i, i) )
//           {
//             i = 0;
//             while( (i < 4) && cell.V[i]->info() == true )
//               ++i;
//             if( i < 4 )
//               {
//                 std::cerr << ".";
//                 pump_vertex(cell.V[i]);
//               }
//             else
//               std::cerr << "!";
//           }
//         cell_queue.erase(it);
//         it = cell_queue.begin();
//       }
//     std::cerr << "\n";
//   }

  void pump_vertex(Vertex_handle v)
  {
    using boost::make_transform_iterator;

    std::vector<Cell_handle> incident_cells;
    std::vector<Vertex_handle> incident_vertices;
    
    incident_cells.reserve(64);
    incident_vertices.reserve(128);
    
    tr.incident_cells(v, std::back_inserter(incident_cells));
    tr.incident_vertices(v, std::back_inserter(incident_vertices));
    
    std::map<Facet, double> ratios;
    
    typedef CGAL::Double_map<Facet, double> Pre_star;
    Pre_star pre_star;

    double best_weight = 0;
    double worst_radius_ratio = 1;

    details::Power_from_v<Geom_traits, Vertex_handle> 
      distance_from_v(v, tr.geom_traits());

    const double sq_d_v = 
      *(std::min_element(make_transform_iterator(incident_vertices.begin(),
						 distance_from_v),
			 make_transform_iterator(incident_vertices.end(),
						 distance_from_v)));

    for(typename std::vector<Cell_handle>::const_iterator cit = 
          incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
      {
	const int index = (*cit)->index(v);
	const Facet f = Facet(*cit, index);
	const double r = radius_ratio(tr.tetrahedron(*cit));
	ratios[f] = r;
	if( r < worst_radius_ratio ) worst_radius_ratio = r;
	const Facet opposite_facet = compute_opposite_facet(f);
        if(opposite_facet.first->is_in_domain())
        {
          CGAL_assertion( ! f.first->is_facet_on_surface(f.second) );
          CGAL_assertion( ! opposite_facet.first->is_facet_on_surface(opposite_facet.second) );
          pre_star.insert(f, compute_critical_radius(v,
						     opposite_facet.first));
        }
      }

    double first_radius_ratio = worst_radius_ratio;
    
    details::Map_value_of<Pre_star,
      std::map<Facet, double> > value_of_ratios(ratios);

//     std::cerr << "worst_radius_ratio=" << worst_radius_ratio << std::endl;
//     std::cerr << "=" << 
// 	  *(std::min_element(make_transform_iterator(pre_star.begin(),
// 						     value_of_ratios),
// 			     make_transform_iterator(pre_star.end(),
// 						     value_of_ratios)))
// 	      << std::endl;

    CGAL_assertion(!pre_star.empty());

    std::pair<double, Facet> facet = *(pre_star.front());
    Facet link = facet.second;
    double critical_r = facet.first;
//     pre_star.pop_front();

//     std::cerr << "v=" << &*v << std::endl;
    
    while( (critical_r < (sq_delta * sq_d_v)) &&
           // La condition suivante teste que la facet qui va etre flippee
           // n'est pas contrainte.
           ! link.first->is_facet_on_surface(link.second) )
      {
// 	std::cerr << "critical_r(" << &*link.first
// 		  << "," << link.second
// 		  << ", s=" << pre_star.size() << ")="
// 		  << critical_r << " ";

        // link.first n'est pas contrainte donc on est dans le domaine.
        CGAL_assertion( !tr.is_infinite(link.first) );
	CGAL_assertion( link.first->is_in_domain() );
	CGAL_assertion( ! link.first->
			is_facet_on_surface(link.second) );
 
        Facet opposite_facet = compute_opposite_facet(link);
        const Cell_handle& opposite_cell = opposite_facet.first;
        const int& opposite_index = opposite_facet.second;
        
        CGAL_assertion( ! opposite_cell->
                            is_facet_on_surface(opposite_index) );
	CGAL_assertion( !tr.is_infinite(opposite_cell) );
	CGAL_assertion(	opposite_cell->is_in_domain() );
 
	
	int number_of_erased_facets = 0;

        for(int i = 0; i<4 ; ++i)
          {
            CGAL_assertion( critical_r < power_product(v->point(), opposite_cell->vertex(i)->point()) );
            const Facet new_facet = Facet(opposite_cell, i);
	    const Facet temp_facet = compute_opposite_facet(new_facet);
            if(!pre_star.erase(temp_facet))
	      {
                CGAL_assertion( opposite_cell->neighbor(i) != link.first);
                pre_star.insert(new_facet,
				compute_critical_radius(v,
							new_facet.first));
		const Cell_handle& c = new_facet.first;
		const int& k = new_facet.second;
  
                CGAL_assertion( c != link.first );

// 		std::cerr << &*(c->vertex((k+1)&3)) << "\n"
// 			  << &*(c->vertex((k+2)&3)) << "\n"
// 			  << &*(c->vertex((k+3)&3)) << std::endl;
		CGAL_assertion( tr.geom_traits().orientation_3_object()
				(v->point(),
				 c->vertex((k+1)&3)->point(),
				 c->vertex((k+2)&3)->point(),
				 c->vertex((k+3)&3)->point())
				!= COPLANAR );
		ratios[new_facet] = 
		  radius_ratio(Tetrahedron_3(v->point(),
					     c->vertex((k+1)&3)->point(),
					     c->vertex((k+2)&3)->point(),
					     c->vertex((k+3)&3)->point()));
	      }
	    else
	      {
                ++number_of_erased_facets;
// 		if(number_of_erased_facets > 1)
//  		  std::cerr << "[3-2]";
	      }
          }
	
	// 	std::cerr << "n(" << number_of_erased_facets << ")";
	CGAL_assertion( number_of_erased_facets > 0 );
	CGAL_assertion( number_of_erased_facets < 3 );

	double min_of_pre_star = 
	  *(std::min_element(make_transform_iterator(pre_star.begin(),
						     value_of_ratios),
			     make_transform_iterator(pre_star.end(),
						     value_of_ratios)));
	
// 	std::cerr << "min_of_pre_star=" << min_of_pre_star << std::endl;
        if( min_of_pre_star > worst_radius_ratio )
	  {
	    worst_radius_ratio = min_of_pre_star;
	    if( !pre_star.empty() )
	      best_weight = (critical_r + pre_star.front()->first) / 2;
	    else
	      best_weight = critical_r;
	  }

        facet = *(pre_star.front());
        link = facet.second;
        critical_r = facet.first;
// 	pre_star.pop_front();
      }

/*    if( link.first->is_facet_on_surface(link.second) )
      std::cerr << "X";
 */   
    Bare_point p = v->point().point();
    //    tr.remove(v);
    if(best_weight > v->point().weight())
    {
      ++num_of_pumped_vertices;
//       std::cerr << "Best weight for point(" << p << ") is "
//           << best_weight << "\nand worst_radius_ratio is "
//           << worst_radius_ratio << " and was " << first_radius_ratio << std::endl;;
      
      total_pumping += (worst_radius_ratio - first_radius_ratio);
      
      Weighted_point wp = Weighted_point(p, best_weight);
      Vertex_handle new_vertex = tr.insert(wp);
      restore_markers_around(new_vertex);
    }
    else
      ++num_of_ignored_vertices;
    
  } // end pump_vertex
  
private: // data
  Tr& tr;
  double sq_delta;

  int num_of_pumped_vertices;
  int num_of_ignored_vertices;
  double total_pumping;

  struct My_cell
  {
    My_cell(const Vertex_handle& v0, const Vertex_handle& v1,
            const Vertex_handle& v2, const Vertex_handle& v3)
    {
      V[0] = v0;
      V[1] = v1;
      V[2] = v2;
      V[3] = v3;
    }
    
    Vertex_handle V[4];
  }; // end struct My_cell
  
  typedef std::map<double, My_cell> Cell_radius_priority_queue;
  
  Cell_radius_priority_queue cell_queue;

private: // methods
  void compute_aspect_ratio_priority_map()
  {
    for(Finite_cell_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end();
        ++cit)
      if(cit->is_in_domain())
        cell_queue.insert(std::make_pair(radius_ratio(tr.tetrahedron(cit)),
                                         My_cell(cit->vertex(0),
                                                 cit->vertex(1),
                                                 cit->vertex(2),
                                                 cit->vertex(3))));

    // DEBUG
//     typename Cell_radius_priority_queue::iterator q_it = cell_queue.begin();
//     for(int i = 0; i<20; ++i)
//       std::cerr << "Radius ratio(" << i << ")="
//              << (q_it++)->first << std::endl;
    // end DEBUG
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
/*    typedef typename Geom_traits::Construct_weighted_circumcenter_3
      Construct_weighted_circumcenter_3;
    typedef typename Geom_traits::
      Compute_squared_radius_smallest_orthogonal_sphere_3
      Compute_squared_radius_smallest_orthogonal_sphere_3;

    Construct_weighted_circumcenter_3 weighted_circumcenter =
      tr.geom_traits().construct_weighted_circumcenter_3_object();
    Compute_squared_radius_smallest_orthogonal_sphere_3
      squared_radius_smallest_orthogonal_sphere =
      tr.geom_traits().
      compute_squared_radius_smallest_orthogonal_sphere_3_object();

    const Weighted_point& wp0 = c->vertex(0)->point();
    const Weighted_point& wp1 = c->vertex(1)->point();
    const Weighted_point& wp2 = c->vertex(2)->point();
    const Weighted_point& wp3 = c->vertex(3)->point();
    
    Weighted_point wc(weighted_circumcenter(wp0, wp1, wp2, wp3),
		      squared_radius_smallest_orthogonal_sphere(wp0,
								wp1,
								wp2,
								wp3));

    return CGAL::to_double(power_product(wc,v->point()));*/
    
//     std::cerr << "c(";
    typedef typename Geom_traits::Compute_critical_squared_radius_3
      Critical_radius;
    Critical_radius critical_radius = 
      tr.geom_traits().compute_critical_squared_radius_3_object();
    
    double result = critical_radius(c->vertex(0)->point(),
			            c->vertex(1)->point(),
				    c->vertex(2)->point(),
				    c->vertex(3)->point(),
				    v->point());
//     std::cerr << result << ")";
    return result;
  }

  static Facet compute_opposite_facet(Facet f)
  {
    const Cell_handle& neighbor_cell = f.first->neighbor(f.second);
    const int opposite_index = neighbor_cell->index(f.first);
    return Facet(neighbor_cell, opposite_index);
  }

  void restore_markers_around(Vertex_handle v)
  {
    std::vector<Cell_handle> incident_cells;
    incident_cells.reserve(64);
    tr.incident_cells(v, std::back_inserter(incident_cells));

//     std::cerr << "S(" << incident_cells.size() << ")";
    bool test = false;

    for(typename std::vector<Cell_handle>::const_iterator cit = 
          incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
      {	
	const int index = (*cit)->index(v);
	const Facet f = Facet(*cit, index);
	const Facet opposite_facet = compute_opposite_facet(f);
	bool marker = 
	  opposite_facet.first->is_facet_on_surface(opposite_facet.second);
	(*cit)->set_surface_facet(index, marker);
	test = test || marker;
	(*cit)->set_in_domain(true);
      }
//     if( test ) 
//       std::cerr << "Y";
  }
}; // end class Slivers_exuder

} // end namespace Mesh_3

template < class Tr>
void
output_slivers_to_off (std::ostream& os, const Tr & T, const double sliver_bound = 0.25) 
{
  typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Point Point;
  
  typedef std::set<Cell_handle> Slivers;
  Slivers slivers;
  
  for( Finite_cells_iterator cit = T.finite_cells_begin(); 
       cit != T.finite_cells_end(); ++cit)
  {
    if( cit->is_in_domain() &&
        radius_ratio(T.tetrahedron(cit)) < sliver_bound )
      slivers.insert(cit);
  }
  
  // Header.
  os << "OFF \n"
      << T.number_of_vertices() << " " <<
      slivers.size() * 4 << 
      " " << 0 << "\n";

  os << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  
  int inum = 0;
  for( Finite_vertices_iterator
       vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
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
    
    if( c->vertex(0)->info() &&
        c->vertex(1)->info() &&
        c->vertex(2)->info() &&
        c->vertex(3)->info() )
      ++flat_slivers;
    
    if( c->vertex(0)->info() ||
        c->vertex(1)->info() ||
        c->vertex(2)->info() ||
        c->vertex(3)->info() )
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


#endif // end CGAL_MESH_3_SLIVER_EXUDER_H
