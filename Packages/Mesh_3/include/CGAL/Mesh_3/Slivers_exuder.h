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

template <typename Tr, 
	  typename FT = typename Tr::Geom_traits::FT >
class Slivers_exuder
{
public: //types
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Vertex_handle Vertex_handle;

  typedef typename Tr::Finite_cells_iterator Finite_cell_iterator;
  typedef typename Weighted_point::Weight Weight;
  
public: // methods
  Slivers_exuder(Tr& t, double d = 0.45) : tr(t), delta(d)
  {
  }

  void init()
  {
    cell_queue.clear();
    compute_aspect_ratio_priority_map();
  }

  void pump_vertices()
  {
    typename Cell_radius_priority_queue::iterator it = cell_queue.begin();
    while ( it != cell_queue.end() )
      {
	My_cell cell = it->second;
	int i;
	Cell_handle ch;
	if( tr.is_cell(cell.V[0],
		       cell.V[1],
		       cell.V[2],
		       cell.V[3],
		       ch, i, i, i, i) )
	  {
	    i = 0;
	    while( (i < 4) && cell.V[i]->info() == true )
	      ++i;
	    if( i < 4 )
	      {
		std::cerr << ".";
		pump_vertex(cell.V[i]);
	      }
	    else
	      std::cerr << "!";
	  }
	cell_queue.erase(it);
	it = cell_queue.begin();
      }
    std::cerr << "\n";
  }

  void pump_vertex(Vertex_handle v)
  {
    std::vector<Cell_handle> incident_cells;
    incident_cells.reserve(64);
    tr.incident_cells(v, std::back_inserter(incident_cells));
    std::map<Facet, double> ratios;
    CGAL::Double_map<Facet, double> pre_start;

    double best_weight = 0;
    double best_weight_next = 0;
    double worst_radius_ratio = 1;

    for(typename std::vector<Cell_handle>::const_iterator cit = 
	  incident_cells.begin();
	cit != incident_cells.end();
	++cit)
      {
	int index = (*cit)->index(v);
	Facet f = Facet(*cit, index);
	const double r = radius_ratio(tr.tetrahedron(*cit));
	ratios[f] = r;
	if( r < worst_radius_ratio ) worst_radius_ratio = r;
	pre_start.insert(f, compute_critical_radius(v, f));
      }

    CGAL_assertion(!pre_start.empty());

    std::pair<double, Facet> facet = *(pre_start.front());
    Facet link = facet.second;
    double critical_r = facet.first;
    pre_start.pop_front();

    while( (critical_r < delta) &&
	   // La condition suivante teste que la facet qui va etre flippee
	   // n'est pas contrainte.
	   ! link.first->is_facet_on_surface(link.second) )
      {
	
	facet = *(pre_start.front());
	link = facet.second;
	critical_r = facet.first;
	pre_start.pop_front();
      }
    

  } // end pump_vertex
  
private: // data
  Tr& tr;
  double delta;
    
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
// 		<< (q_it++)->first << std::endl;
    // end DEBUG
  }

  double compute_critical_radius(Vertex_handle, Facet)
  {
    return 0.;
  }
  
}; // end class Slivers_exuder

} // end namespace CGAL


#endif // end CGAL_MESH_3_SLIVER_EXUDER_H
