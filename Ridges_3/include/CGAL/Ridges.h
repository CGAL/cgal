// Copyright (c) 2007  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Marc Pouget and Frédéric Cazals
#ifndef CGAL_RIDGE_3_H_
#define CGAL_RIDGE_3_H_

#include <utility>
#include <list>
#include <map>

#include <CGAL/basic.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>

#include <CGAL/property_map.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
 
enum Ridge_interrogation_type {MAX_RIDGE, MIN_RIDGE, CREST_RIDGE};

enum Ridge_type {NO_RIDGE=0, 
		 MAX_ELLIPTIC_RIDGE, MAX_HYPERBOLIC_RIDGE, MAX_CREST_RIDGE, 
		 MIN_ELLIPTIC_RIDGE, MIN_HYPERBOLIC_RIDGE, MIN_CREST_RIDGE};
 
//are ridges tagged as elliptic or hyperbolic using 3rd or 4th order
//differential quantitities?
//with Ridge_order_3 P1 and P2 are not used and the sharpness is not defined.
enum Ridge_order {Ridge_order_3 = 3, Ridge_order_4 = 4};
  
//---------------------------------------------------------------------------
//Ridge_line : a connected sequence of edges of a
//TriangularPolyhedralSurface crossed by a
//ridge (with a barycentric coordinate to compute the crossing point),
//with a Ridge_type and weights : strength and sharpness. Note
//sharpness is only available (more precisely only meaningful)
//if the Ridge_approximation has
//been computed with the Ridge_order Ridge_order_4.
//(else, if it is computed with Ridge_order_3 it keeps its initial
//value 0)
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh > class Ridge_line
{
public:
  typedef typename TriangulatedSurfaceMesh::Traits::FT         FT;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3   Vector_3;
  typedef typename TriangulatedSurfaceMesh::Traits::Point_3    Point_3;
  typedef typename TriangulatedSurfaceMesh::Halfedge_const_handle Halfedge_const_handle;
  typedef std::pair< Halfedge_const_handle, FT> ridge_halfhedge; 

  Ridge_type line_type() const {return m_line_type;}
  Ridge_type& line_type() {return m_line_type;}

  const FT strength() const {return m_strength;}
  FT& strength() {return m_strength;}

  const FT sharpness() const {return m_sharpness;}
  FT& sharpness() {return m_sharpness;}

  const std::list<ridge_halfhedge>* line() const { return &m_line;}
  std::list<ridge_halfhedge>* line() { return &m_line;}

  //constructor
  Ridge_line();
  
  /* The output is : line_type, strength, sharpness, list of points of
     the polyline. An insert operator << is also available.
   */
  void dump_4ogl(std::ostream& out_stream) const ;
  void dump_verbose(std::ostream& out_stream) const ;

protected:
  //one of MAX_ELLIPTIC_RIDGE, MAX_HYPERBOLIC_RIDGE, MAX_CREST_RIDGE,
  //MIN_ELLIPTIC_RIDGE, MIN_HYPERBOLIC_RIDGE or MIN_CREST_RIDGE
  Ridge_type m_line_type;  
  std::list<ridge_halfhedge> m_line;
  FT m_strength;// = integral of ppal curvature along the line
  FT m_sharpness;// = (integral of second derivative of curvature
		 // along the line) multiplied by the squared of 
                 // the size of the model
		 // (which is the radius of the smallest enclosing
		 // ball)
};

//--------------------------------------------------------------------------
// IMPLEMENTATION OF Ridge_line members
//--------------------------------------------------------------------------

 //constructor
template < class TriangulatedSurfaceMesh >
Ridge_line<TriangulatedSurfaceMesh>::
Ridge_line() : m_strength(0.), m_sharpness(0.)  {}
   

template < class TriangulatedSurfaceMesh >
void Ridge_line<TriangulatedSurfaceMesh>::
dump_4ogl(std::ostream& out_stream) const
{
  out_stream << line_type() << " "
	     << strength() << " "
	     << sharpness() << " ";

  typename std::list<ridge_halfhedge >::const_iterator
    iter = line()->begin(), 
    ite =  line()->end();
  for (;iter!=ite;iter++){
    //he: p->q, r is the crossing point
    Point_3 p = iter->first->opposite()->vertex()->point(),
            q = iter->first->vertex()->point();
    Point_3 r = CGAL::barycenter(p, iter->second, q);
    out_stream << " " << r ;	
  }
  out_stream  << std::endl;  
}

//verbose output
template < class TriangulatedSurfaceMesh >
void Ridge_line<TriangulatedSurfaceMesh>::
dump_verbose(std::ostream& out_stream) const
{
  out_stream << "Line type is : " << line_type() << std::endl
	     << "Strength is :  " << strength() << std::endl
	     << "Sharpness is : " << sharpness() << std::endl
	     << "Polyline point coordinates are : " << std::endl;

  typename std::list<ridge_halfhedge>::const_iterator
    iter = line()->begin(), 
    ite =  line()->end();
  for (;iter!=ite;iter++){
    //he: p->q, r is the crossing point
    Point_3 p = iter->first->opposite()->vertex()->point(),
            q = iter->first->vertex()->point();
    Point_3 r = CGAL::barycenter(p, iter->second, q);
    out_stream << r << std::endl;	
  }
}

template <class TriangulatedSurfaceMesh>
std::ostream& 
operator<<(std::ostream& out_stream, const Ridge_line<TriangulatedSurfaceMesh>& ridge_line)
{
  ridge_line.dump_verbose(out_stream);
  return out_stream;
}

//---------------------------------------------------------------------------
//Vertex2Data_Property_Map_with_std_map
// defines models for Vertex2FTPropertyMap and Vertex2VectorPropertyMap
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh >
class Vertex2Data_Property_Map_with_std_map 
{
 public:
  typedef typename TriangulatedSurfaceMesh::Traits::FT        FT;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3  Vector_3;
  typedef typename TriangulatedSurfaceMesh::Vertex_const_handle Vertex_const_handle;

  struct Vertex_cmp{
    bool operator()(Vertex_const_handle a,  Vertex_const_handle b) const{
      return &*a < &*b;
    }
  };

  typedef std::map<Vertex_const_handle, FT, Vertex_cmp> Vertex2FT_map;
  typedef boost::associative_property_map< Vertex2FT_map > Vertex2FT_property_map;

  typedef std::map<Vertex_const_handle, Vector_3, Vertex_cmp> Vertex2Vector_map;
  typedef boost::associative_property_map< Vertex2Vector_map > Vertex2Vector_property_map;
};


//---------------------------------------------------------------------------
//Ridge_approximation
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
class Ridge_approximation
{
 public:  
  typedef typename TriangulatedSurfaceMesh::Traits::FT        FT;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3  Vector_3;
  typedef typename TriangulatedSurfaceMesh::Vertex_const_handle     Vertex_const_handle;
  typedef typename TriangulatedSurfaceMesh::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename TriangulatedSurfaceMesh::Facet_const_handle      Facet_const_handle;
  typedef typename TriangulatedSurfaceMesh::Facet_const_iterator    Facet_const_iterator;

  //requirements for the templates TriangulatedSurfaceMesh and Vertex2FTPropertyMap or Vertex2VectorPropertyMap
  CGAL_static_assertion((boost::is_same<Vertex_const_handle, typename Vertex2FTPropertyMap::key_type>::value));
  CGAL_static_assertion((boost::is_same<Vertex_const_handle, typename Vertex2VectorPropertyMap::key_type>::value));
  CGAL_static_assertion((boost::is_same<FT, typename Vertex2FTPropertyMap::value_type>::value));
  CGAL_static_assertion((boost::is_same<Vector_3, typename Vertex2VectorPropertyMap::value_type>::value));

  typedef std::pair< Halfedge_const_handle, FT>    Ridge_halfhedge;
  typedef CGAL::Ridge_line<TriangulatedSurfaceMesh>  Ridge_line;

  Ridge_approximation(const TriangulatedSurfaceMesh &P,
		      const Vertex2FTPropertyMap& vertex2k1_pm, 
		      const Vertex2FTPropertyMap& vertex2k2_pm,
		      const Vertex2FTPropertyMap& vertex2b0_pm, 
		      const Vertex2FTPropertyMap& vertex2b3_pm,
		      const Vertex2VectorPropertyMap& vertex2d1_pm, 
		      const Vertex2VectorPropertyMap& vertex2d2_pm,
		      const Vertex2FTPropertyMap& vertex2P1_pm, 
		      const Vertex2FTPropertyMap& vertex2P2_pm);
 
  template <class OutputIterator>
  OutputIterator compute_max_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);
  template <class OutputIterator>
  OutputIterator compute_min_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);
  template <class OutputIterator>
  OutputIterator compute_crest_ridges(OutputIterator it, Ridge_order ord = Ridge_order_3);
  
  //Find MAX_RIDGE, MIN_RIDGE or CREST_RIDGE ridges iterate on P facets,
  //find a non-visited, regular (i.e. if there is a coherent
  //orientation of ppal dir at the facet vertices), 2Xing triangle,
  //follow non-visited, regular, 2Xing triangles in both sens to
  //create a Ridge line.  Each time an edge is added the strength and
  //sharpness(if Ridge_order_4) are updated.
  template <class OutputIterator>
  OutputIterator compute_ridges(Ridge_interrogation_type r_type, 
			  OutputIterator ridge_lines_it,
			  Ridge_order ord = Ridge_order_3);

 protected:
  const TriangulatedSurfaceMesh& P;
  FT squared_model_size;//squared radius of the smallest enclosing sphere of the TriangulatedSurfaceMesh
		//used to make the sharpness scale independant and iso indep
  Ridge_order tag_order;

  //tag to visit faces
  struct Facet_cmp{ //comparison is wrt facet addresses
    bool operator()(Facet_const_handle a,  Facet_const_handle b) const{
      return &*a < &*b;
    }
  };
  typedef std::map<Facet_const_handle, bool, Facet_cmp> Facet2bool_map_type;
  Facet2bool_map_type is_visited_map;

  //Property maps
  const Vertex2FTPropertyMap &k1, &k2, &b0, &b3, &P1, &P2;
  const Vertex2VectorPropertyMap &d1, &d2;

  //is a facet crossed by a BLUE, RED or CREST_RIDGE ridge? if so, return
  //the crossed edges and more precise type from MAX_ELLIPTIC_RIDGE,
  //MAX_HYPERBOLIC_RIDGE, MAX_CREST_RIDGE, MIN_ELLIPTIC_RIDGE,
  //MIN_HYPERBOLIC_RIDGE, MIN_CREST_RIDGE or NO_RIDGE
  Ridge_type facet_ridge_type(const Facet_const_handle f, 
			      Halfedge_const_handle& he1, 
			      Halfedge_const_handle& he2,
			      Ridge_interrogation_type r_type);
  
  //is an edge crossed by a BLUE/RED ridge? (color is MAX_RIDGE or
  //MIN_RIDGE ).  As we only test edges of regular triangles, the ppal
  //direction at endpoints d_p and d_q cannot be orthogonal. If both
  //extremalities vanish, we consider no crossing occurs. If only one
  //of them vanishes, we consider it as an positive infinitesimal and
  //apply the general rule. The general rule is that for both
  //non-vanishing extremalities, a crossing occurs if their sign
  //differ; Assuming the accute rule to orient the ppal directions,
  //there is a crossing iff d_p.d_q * b_p*b_q < 0
  void xing_on_edge(const Halfedge_const_handle he, 
		    bool& is_crossed, 
		    Ridge_interrogation_type color);
 
  //given a ridge segment of a given color, in a triangle crossing he1
  //(v_p1 -> v_q1) and he2 (v_p2 -> v_q2) return true if it is
  //elliptic, false if it is hyperbolic.
  bool tag_as_elliptic_hyperbolic(const Ridge_interrogation_type color,
				  const Halfedge_const_handle he1, 
				  const Halfedge_const_handle he2);

  //for the computation with tag_order == 3 only
  //for a ridge segment [r1,r2] in a triangle (v1,v2,v3), let r = r2 -
  //r1 and normalize, the projection of a point p on the line (r1,r2)
  //is pp=r1+tr, with t=(p-r1)*r then the vector v starting at p is
  //pointing to the ridge line (r1,r2) if (pp-p)*v >0. Return the sign
  //of b, for a ppal direction pointing to the ridge segment,
  //appearing at least at two vertices of the facet.
  //
  // for color = MAX_RIDGE, sign = 1 if MAX_ELLIPTIC_RIDGE, -1 if
  // MAX_HYPERBOLIC_RIDGE 
  //
  // for color = MIN_RIDGE, sign = -1 if MIN_ELLIPTIC_RIDGE, 1 if
  // MIN_HYPERBOLIC_RIDGE
  int b_sign_pointing_to_ridge(const Vertex_const_handle v1, 
			       const Vertex_const_handle v2,
			       const Vertex_const_handle v3,
			       const Vector_3 r1, const Vector_3 r2, 
			       const Ridge_interrogation_type color);

  //a ridge line begins with a segment in a triangle given by the 2 he
  //crossed
  void init_ridge_line(Ridge_line* ridge_line, 
		       const Halfedge_const_handle h1, 
		       const Halfedge_const_handle h2, 
		       const Ridge_type r_type);
  //When the line is extended with a he, the bary coord of the
  //crossing point is computed, the pair (he,coord) is added and the
  //weights are updated 
  void addback(Ridge_line* ridge_line, 
	       const Halfedge_const_handle he, 
	       const Ridge_type r_type);
  void addfront(Ridge_line* ridge_line, 
		const Halfedge_const_handle he,
		const Ridge_type r_type);

  //compute the barycentric coordinate of the xing point (blue or red)
  //for he: p->q (wrt the extremality values b0/3).  coord is st
  //xing_point = coord*p + (1-coord)*q
  FT bary_coord(const Halfedge_const_handle he, 
		const Ridge_type r_type);
};


// IMPLEMENTATION OF Ridge_approximation members
/////////////////////////////////////////////////////////////////////////////
 //contructor
template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap,
  class Vertex2VectorPropertyMap > 
  Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap >::
  Ridge_approximation(const TriangulatedSurfaceMesh &p,
		      const Vertex2FTPropertyMap& vertex2k1_pm, 
		      const Vertex2FTPropertyMap& vertex2k2_pm,
		      const Vertex2FTPropertyMap& vertex2b0_pm, 
		      const Vertex2FTPropertyMap& vertex2b3_pm,
		      const Vertex2VectorPropertyMap& vertex2d1_pm, 
		      const Vertex2VectorPropertyMap& vertex2d2_pm,
		      const Vertex2FTPropertyMap& vertex2P1_pm, 
		      const Vertex2FTPropertyMap& vertex2P2_pm)
    : P(p), k1(vertex2k1_pm), k2(vertex2k2_pm), b0(vertex2b0_pm), b3(vertex2b3_pm), 
      P1(vertex2P1_pm), P2(vertex2P2_pm), d1(vertex2d1_pm), d2(vertex2d2_pm)
{
  //init the is_visited_map and check that the mesh is a triangular one.
  Facet_const_iterator itb = P.facets_begin(), ite = P.facets_end();
  for(;itb!=ite;itb++) {
    is_visited_map[itb] = false;
    CGAL_precondition( itb->is_triangle() );
  }

  CGAL::Min_sphere_d<CGAL::Optimisation_d_traits_3<typename TriangulatedSurfaceMesh::Traits> > 
    min_sphere(P.points_begin(), P.points_end());
  squared_model_size = min_sphere.squared_radius();
  //maybe better to use CGAL::Min_sphere_of_spheres_d ?? but need to create spheres?

  tag_order = Ridge_order_3;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
  template <class OutputIterator>
  OutputIterator Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  compute_max_ridges(OutputIterator it, Ridge_order ord)
{
  compute_ridges(MAX_RIDGE, it, ord);
  return it;
}
template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
  template <class OutputIterator>
  OutputIterator Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  compute_min_ridges(OutputIterator it, Ridge_order ord)
{
  compute_ridges(MIN_RIDGE, it, ord);
  return it;
}
template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
  template <class OutputIterator>
  OutputIterator Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  compute_crest_ridges(OutputIterator it, Ridge_order ord)
{
  compute_ridges(CREST_RIDGE, it, ord);
  return it;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
  template <class OutputIterator>
  OutputIterator Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  compute_ridges(Ridge_interrogation_type r_type, OutputIterator ridge_lines_it, Ridge_order ord)
{
  tag_order = ord;

  //reinit the is_visited_map
  Facet_const_iterator itb = P.facets_begin(), ite = P.facets_end();
  for(;itb!=ite;itb++) is_visited_map[itb] = false;
  
  itb = P.facets_begin();
  for(;itb!=ite;itb++)
    {
      Facet_const_handle f = itb;
      if (is_visited_map.find(f)->second) continue;
      is_visited_map.find(f)->second = true;
      Halfedge_const_handle h1, h2, curhe1, curhe2, curhe;
      
      //h1 h2 are the hedges crossed if any, r_type should be
      //MAX_RIDGE, MIN_RIDGE or CREST_RIDGE ; cur_ridge_type should be
      //MAX_ELLIPTIC_RIDGE, MAX_HYPERBOLIC_RIDGE, MAX_CREST_RIDGE,
      //MIN_ELLIPTIC_RIDGE, MIN_HYPERBOLIC_RIDGE, MIN_CREST_RIDGE or NO_RIDGE
      Ridge_type cur_ridge_type = facet_ridge_type(f,h1,h2,r_type);
      if ( cur_ridge_type == NO_RIDGE ) continue;
      
      //a ridge_line is begining and stored
      Ridge_line* cur_ridge_line = new Ridge_line();
      init_ridge_line(cur_ridge_line, h1, h2, cur_ridge_type);
      *ridge_lines_it++ = cur_ridge_line;
    
      //next triangle adjacent to h1 (push_front)
      if ( !(h1->is_border_edge()) ) 
	{
	  f = h1->opposite()->facet();
	  curhe = h1;
	  while (cur_ridge_type == facet_ridge_type(f,curhe1,curhe2,r_type))
	    {
	      //follow the ridge from curhe
	      if (is_visited_map.find(f)->second) break;
	      is_visited_map.find(f)->second = true;
	      if (curhe->opposite() == curhe1) curhe = curhe2;
	      else curhe = curhe1;//curhe stays at the ridge extremity
	      addfront(cur_ridge_line, curhe, cur_ridge_type);
	      if ( !(curhe->is_border_edge()) ) f =
						  curhe->opposite()->facet();
	      else break;
	    }
	  //exit from the while if
	  //1. border or already visited (this is a ridge loop)
	  //2. not same type, then do not set visisted cause a MAX_ELLIPTIC_RIDGE
	  //	  follows a MAX_HYPERBOLIC_RIDGE
	}

      //next triangle adjacent to h2 (push_back)
      if ( !(h2->is_border_edge()) ) 
	{
	  f = h2->opposite()->facet();
	  curhe = h2;
	  while (cur_ridge_type ==
		 facet_ridge_type(f,curhe1,curhe2,r_type))
	    {
	      //follow the ridge from curhe
	      if (is_visited_map.find(f)->second) break;
	      is_visited_map.find(f)->second = true;
	      if (curhe->opposite() == curhe1) curhe = curhe2;
	      else curhe = curhe1;
	      addback(cur_ridge_line, curhe, cur_ridge_type);
	      if ( !(curhe->is_border_edge()) ) f =
						  curhe->opposite()->facet();
	      else break;
	    }
	} 
    }
  return ridge_lines_it;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
Ridge_type Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
facet_ridge_type(const Facet_const_handle f, Halfedge_const_handle& he1, Halfedge_const_handle&
		 he2, Ridge_interrogation_type r_type)
{
  //polyhedral data
  //we have v1--h1-->v2--h2-->v3--h3-->v1
  const Halfedge_const_handle h1 = f->halfedge();
  const Vertex_const_handle v2 = h1->vertex();
  const Halfedge_const_handle h2 = h1->next();
  const Vertex_const_handle v3 = h2->vertex();
  const Halfedge_const_handle h3 = h2->next();
  const Vertex_const_handle v1 = h3->vertex();

  //check for regular facet
  //i.e. if there is a coherent orientation of ppal dir at the facet vertices
  if ( d1[v1]*d1[v2] * d1[v1]*d1[v3] * d1[v2]*d1[v3] < 0 ) 
    return NO_RIDGE;
   
  //determine potential crest color
  //MAX_CREST_RIDGE if |sum(k1)|>|sum(k2)| sum over facet vertices vi
  //MIN_CREST_RIDGE if |sum(k1)|<|sum(k2)|
  Ridge_type crest_color = NO_RIDGE;
  if (r_type == CREST_RIDGE) 
    {
      if ( CGAL::abs(k1[v1]+k1[v2]+k1[v3]) > CGAL::abs(k2[v1]+k2[v2]+k2[v3]) ) 
	crest_color = MAX_CREST_RIDGE; 
      if ( CGAL::abs(k1[v1]+k1[v2]+k1[v3]) < CGAL::abs(k2[v1]+k2[v2]+k2[v3]) ) 
	crest_color = MIN_CREST_RIDGE;
      if ( CGAL::abs(k1[v1]+k1[v2]+k1[v3]) == CGAL::abs(k2[v1]+k2[v2]+k2[v3]) ) 
	return NO_RIDGE;
    }
  
  //compute Xing on the 3 edges
  bool h1_is_crossed = false, h2_is_crossed = false, h3_is_crossed = false;
  if ( r_type == MAX_RIDGE || crest_color == MAX_CREST_RIDGE ) 
    {
      xing_on_edge(h1, h1_is_crossed, MAX_RIDGE);
      xing_on_edge(h2, h2_is_crossed, MAX_RIDGE);
      xing_on_edge(h3, h3_is_crossed, MAX_RIDGE);
    }
  if ( r_type == MIN_RIDGE || crest_color == MIN_CREST_RIDGE ) 
    {
      xing_on_edge(h1, h1_is_crossed, MIN_RIDGE);
      xing_on_edge(h2, h2_is_crossed, MIN_RIDGE);
      xing_on_edge(h3, h3_is_crossed, MIN_RIDGE);
    }

  //there are either 0 or 2 crossed edges
  if ( !h1_is_crossed && !h2_is_crossed && !h3_is_crossed ) 
    return NO_RIDGE; 
  if (h1_is_crossed && h2_is_crossed && !h3_is_crossed)
    {
      he1 = h1; 
      he2 = h2;
    }
  if (h1_is_crossed && !h2_is_crossed && h3_is_crossed)
    {
      he1 = h1; 
      he2 = h3;
    }
  if (!h1_is_crossed && h2_is_crossed && h3_is_crossed)
    {
      he1 = h2; 
      he2 = h3;
    }
  //check there is no other case (just one edge crossed)
  CGAL_postcondition ( !( (h1_is_crossed && !h2_is_crossed && !h3_is_crossed)
			  || (!h1_is_crossed && h2_is_crossed && !h3_is_crossed)
			  || (!h1_is_crossed && !h2_is_crossed && h3_is_crossed)) );

  //There is a ridge segment in the triangle, determine its type elliptic/hyperbolic
  bool is_elliptic;  
  if ( r_type == MAX_RIDGE || crest_color == MAX_CREST_RIDGE ) 
    is_elliptic = tag_as_elliptic_hyperbolic(MAX_RIDGE, he1, he2);
  else is_elliptic = tag_as_elliptic_hyperbolic(MIN_RIDGE, he1, he2);
  
  if (r_type == MAX_RIDGE) 
    {if (is_elliptic) return MAX_ELLIPTIC_RIDGE;
    else return MAX_HYPERBOLIC_RIDGE; }
  if (crest_color == MAX_CREST_RIDGE && is_elliptic) return MAX_CREST_RIDGE;

  if (r_type == MIN_RIDGE) 
    {if (is_elliptic) return MIN_ELLIPTIC_RIDGE;
    else return MIN_HYPERBOLIC_RIDGE; }
  if (crest_color == MIN_CREST_RIDGE && is_elliptic) return MIN_CREST_RIDGE;
  
  return NO_RIDGE;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
void Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
xing_on_edge(const Halfedge_const_handle he, bool& is_crossed, Ridge_interrogation_type color)
{
  is_crossed = false;
  FT sign = 0;
  FT b_p, b_q; // extremalities at p and q for he: p->q
  Vector_3  d_p = d1[he->opposite()->vertex()],
    d_q = d1[he->vertex()]; //ppal dir
  if ( color == MAX_RIDGE ) {
    b_p = b0[he->opposite()->vertex()];
    b_q = b0[he->vertex()];
  }
  else {     
    b_p = b3[he->opposite()->vertex()];
    b_q = b3[he->vertex()];
  }
  if ( b_p == 0 && b_q == 0 ) return;
  if ( b_p == 0 && b_q !=0 ) sign = d_p*d_q * b_q;
  if ( b_p != 0 && b_q ==0 ) sign = d_p*d_q * b_p;
  if ( b_p != 0 && b_q !=0 ) sign = d_p*d_q * b_p * b_q;
  if ( sign < 0 ) is_crossed = true;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
bool Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  tag_as_elliptic_hyperbolic(const Ridge_interrogation_type color,
			     const Halfedge_const_handle he1, 
			     const Halfedge_const_handle he2)
{
  const Vertex_const_handle v_p1 = he1->opposite()->vertex(), v_q1 = he1->vertex(),
    v_p2 = he2->opposite()->vertex(), v_q2 = he2->vertex(); // hei: pi->qi

  FT coord1, coord2;
  if (color == MAX_RIDGE) 
    {
      coord1 = CGAL::abs(b0[v_q1]) / ( CGAL::abs(b0[v_p1]) + CGAL::abs(b0[v_q1]) );
      coord2 = CGAL::abs(b0[v_q2]) / ( CGAL::abs(b0[v_p2]) + CGAL::abs(b0[v_q2]) ); 
    }
  else 
    {
      coord1 = CGAL::abs(b3[v_q1]) / ( CGAL::abs(b3[v_p1]) + CGAL::abs(b3[v_q1]) );
      coord2 = CGAL::abs(b3[v_q2]) / ( CGAL::abs(b3[v_p2]) + CGAL::abs(b3[v_q2]) ); 
    }

  if ( tag_order == Ridge_order_3 ) {
    Vector_3 r1 = CGAL::barycenter(v_p1->point(), coord1, v_q1->point()) - ORIGIN,
             r2 = CGAL::barycenter(v_p2->point(), coord2, v_q2->point()) - ORIGIN; 
    //identify the 3 different vertices v_p1, v_q1 and v3 = v_p2 or v_q2
    Vertex_const_handle v3;
    if (v_p2 == v_p1 || v_p2 == v_q1) v3 = v_q2;
    else v3 = v_p2;

    int b_sign = b_sign_pointing_to_ridge(v_p1, v_q1, v3, r1, r2, color); 

    if (color == MAX_RIDGE) 
      if (b_sign == 1) return true; 
      else return false;
    else if (b_sign == -1) return true; 
      else return false;
  }
  else {//tag_order == Ridge_order_4, check the sign of the meanvalue of the signs
    //      of Pi at the two crossing points
    FT sign_P;
    if (color == MAX_RIDGE) 
      sign_P =  P1[v_p1]*coord1 + P1[v_q1]*(1-coord1) 
	+ P1[v_p2]*coord2 + P1[v_q2]*(1-coord2);
    else sign_P =  P2[v_p1]*coord1 + P2[v_q1]*(1-coord1) 
	+ P2[v_p2]*coord2 + P2[v_q2]*(1-coord2);

    if ( sign_P < 0 ) return true; else return false;
  }
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
int Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
  b_sign_pointing_to_ridge(const Vertex_const_handle v1, 
			       const Vertex_const_handle v2,
			       const Vertex_const_handle v3,
			       const Vector_3 r1, const Vector_3 r2, 
			       const Ridge_interrogation_type color)
{
  Vector_3 r = r2 - r1, dv1, dv2, dv3;
  FT bv1, bv2, bv3;
  if ( color == MAX_RIDGE ) {
    bv1 = b0[v1];
    bv2 = b0[v2];
    bv3 = b0[v3];
    dv1 = d1[v1];
    dv2 = d1[v2];
    dv3 = d1[v3];
  }
  else {
    bv1 = b3[v1];
    bv2 = b3[v2];
    bv3 = b3[v3];
    dv1 = d2[v1];
    dv2 = d2[v2];
    dv3 = d2[v3];    
  }
  if ( r != CGAL::NULL_VECTOR ) r = r/CGAL::sqrt(r*r);
  FT sign1, sign2, sign3;
  sign1 = bv1*(r1 - (v1->point()-ORIGIN) + (((v1->point()-ORIGIN)-r1)*r)*r )*dv1;
  sign2 = bv2*(r1 - (v2->point()-ORIGIN) + (((v2->point()-ORIGIN)-r1)*r)*r )*dv2;
  sign3 = bv3*(r1 - (v3->point()-ORIGIN) + (((v3->point()-ORIGIN)-r1)*r)*r )*dv3;
  
  int compt = 0;
  if ( sign1 > 0 ) compt++; else if (sign1 < 0) compt--;
  if ( sign2 > 0 ) compt++; else if (sign2 < 0) compt--;
  if ( sign3 > 0 ) compt++; else if (sign3 < 0) compt--;
  
  if (compt > 0) return 1; else return -1;
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
void Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
init_ridge_line(Ridge_line* ridge_line, 
		const Halfedge_const_handle h1, 
		const Halfedge_const_handle h2, 
		const Ridge_type r_type)
{
  ridge_line->line_type() = r_type;
  ridge_line->line()->push_back(Ridge_halfhedge(h1, bary_coord(h1,r_type)));
  addback(ridge_line, h2, r_type);
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
void Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
addback(Ridge_line* ridge_line, const Halfedge_const_handle he,
	const Ridge_type r_type)
{
  Halfedge_const_handle he_cur = ( --(ridge_line->line()->end()) )->first;
  FT coord_cur = ( --(ridge_line->line()->end()) )->second;//bary_coord(he_cur);
  FT coord = bary_coord(he,r_type);
  Vertex_const_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he_cur->vertex(); // he: p->q
  Vector_3 segment = CGAL::barycenter(v_p->point(), coord, v_q->point()) -
                     CGAL::barycenter(v_p_cur->point(), coord_cur, v_q_cur->point());

  FT k1x, k2x; //abs value of the ppal curvatures at the Xing point on he.
  FT k_second = 0; // abs value of the second derivative of the curvature
               // along the line of curvature
  k1x = CGAL::abs(k1[v_p]) * coord + CGAL::abs(k1[v_q]) * (1-coord) ;   
  k2x = CGAL::abs(k2[v_p]) * coord + CGAL::abs(k2[v_q]) * (1-coord) ;   

  if ( (ridge_line->line_type() == MAX_ELLIPTIC_RIDGE) 
       || (ridge_line->line_type() == MAX_HYPERBOLIC_RIDGE) 
       || (ridge_line->line_type() == MAX_CREST_RIDGE) ) {
    ridge_line->strength() += k1x * CGAL::sqrt(segment * segment); 
    if (tag_order == Ridge_order_4) { 
      if (k1x != k2x) 
	k_second =CGAL::abs(( CGAL::abs(P1[v_p]) * coord + CGAL::abs(P1[v_q]) * (1-coord) )/(k1x-k2x));
      ridge_line->sharpness() += k_second * CGAL::sqrt(segment * segment) * squared_model_size; }
  }
  if ( (ridge_line->line_type() == MIN_ELLIPTIC_RIDGE) 
       || (ridge_line->line_type() == MIN_HYPERBOLIC_RIDGE) 
       || (ridge_line->line_type() == MIN_CREST_RIDGE) ) {
   ridge_line->strength() += k2x * CGAL::sqrt(segment * segment); 
   if (tag_order == Ridge_order_4) {
     if (k1x != k2x) 
       k_second =CGAL::abs(( CGAL::abs(P2[v_p]) * coord + CGAL::abs(P2[v_q]) * (1-coord) )/(k1x-k2x));
     ridge_line->sharpness() += k_second * CGAL::sqrt(segment * segment) * squared_model_size; }
   } 
  ridge_line->line()->push_back( Ridge_halfhedge(he, coord));
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
void Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
addfront(Ridge_line* ridge_line, 
	 const Halfedge_const_handle he, 
	 const Ridge_type r_type)
{
  Halfedge_const_handle he_cur = ( ridge_line->line()->begin() )->first;
  FT coord_cur = ( ridge_line->line()->begin() )->second;
  FT coord = bary_coord(he,r_type);
  Vertex_const_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he_cur->vertex(); // he: p->q
  Vector_3 segment = CGAL::barycenter(v_p->point(), coord, v_q->point()) -
                     CGAL::barycenter(v_p_cur->point(), coord_cur, v_q_cur->point());

  FT k1x, k2x; //abs value of the ppal curvatures at the Xing point on he.
  FT k_second = 0.; // abs value of the second derivative of the curvature
               // along the line of curvature
  k1x = CGAL::abs(k1[v_p]) * coord + CGAL::abs(k1[v_q]) * (1-coord) ;   
  k2x = CGAL::abs(k2[v_p]) * coord + CGAL::abs(k2[v_q]) * (1-coord) ;   

  if ( (ridge_line->line_type() == MAX_ELLIPTIC_RIDGE) 
       || (ridge_line->line_type() == MAX_HYPERBOLIC_RIDGE) 
       || (ridge_line->line_type() == MAX_CREST_RIDGE) ) {
    ridge_line->strength() += k1x * CGAL::sqrt(segment * segment); 
   if (tag_order == Ridge_order_4) {
     if (k1x != k2x) 
       k_second =CGAL::abs(( CGAL::abs(P1[v_p]) * coord + CGAL::abs(P1[v_q]) * (1-coord) )/(k1x-k2x));
     ridge_line->sharpness() += k_second * CGAL::sqrt(segment * segment) * squared_model_size; }
  }
  if ( (ridge_line->line_type() == MIN_ELLIPTIC_RIDGE) 
       || (ridge_line->line_type() == MIN_HYPERBOLIC_RIDGE) 
       || (ridge_line->line_type() == MIN_CREST_RIDGE) ) {
   ridge_line->strength() += k2x * CGAL::sqrt(segment * segment); 
   if (tag_order == Ridge_order_4) {
     if (k1x != k2x) 
       k_second =CGAL::abs(( CGAL::abs(P2[v_p]) * coord + CGAL::abs(P2[v_q]) * (1-coord) )/(k1x-k2x));
     ridge_line->sharpness() += k_second * CGAL::sqrt(segment * segment) * squared_model_size; }
   } 
  ridge_line->line()->push_front( Ridge_halfhedge(he, coord));
}

template < class TriangulatedSurfaceMesh,  
           class Vertex2FTPropertyMap,
           class Vertex2VectorPropertyMap > 
typename TriangulatedSurfaceMesh::Traits::FT 
Ridge_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap , Vertex2VectorPropertyMap  >::
bary_coord(const Halfedge_const_handle he, const Ridge_type r_type)
{
  FT b_p = 0., b_q = 0.; // extremalities at p and q for he: p->q
  if ( (r_type == MAX_ELLIPTIC_RIDGE) 
       || (r_type == MAX_HYPERBOLIC_RIDGE) 
       || (r_type == MAX_CREST_RIDGE) ) {
    b_p = b0[he->opposite()->vertex()];
    b_q = b0[he->vertex()];    
  }
  if ( (r_type == MIN_ELLIPTIC_RIDGE) 
       || (r_type == MIN_HYPERBOLIC_RIDGE) 
       || (r_type == MIN_CREST_RIDGE) ) {
    b_p = b3[he->opposite()->vertex()];
    b_q = b3[he->vertex()];    
  }
  //denominator cannot be 0 since there is no crossing when both extremalities are 0
  return CGAL::abs(b_q) / ( CGAL::abs(b_q) + CGAL::abs(b_p) );
}
  

//---------------------------------------------------------------------------
//Global functions
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap,
  class Vertex2VectorPropertyMap,
  class OutputIterator>
  OutputIterator compute_max_ridges(const TriangulatedSurfaceMesh &P,
				    const Vertex2FTPropertyMap& vertex2k1_pm, 
				    const Vertex2FTPropertyMap& vertex2k2_pm,
				    const Vertex2FTPropertyMap& vertex2b0_pm, 
				    const Vertex2FTPropertyMap& vertex2b3_pm,
				    const Vertex2VectorPropertyMap& vertex2d1_pm, 
				    const Vertex2VectorPropertyMap& vertex2d2_pm,
				    const Vertex2FTPropertyMap& vertex2P1_pm, 
				    const Vertex2FTPropertyMap& vertex2P2_pm,
				    OutputIterator it, 
				    Ridge_order order = Ridge_order_3)
{
  typedef Ridge_approximation < TriangulatedSurfaceMesh, 
    Vertex2FTPropertyMap, Vertex2VectorPropertyMap > Ridge_approximation;
  
  Ridge_approximation ridge_approximation(P, 
					  vertex2k1_pm, vertex2k2_pm,
					  vertex2b0_pm, vertex2b3_pm,
					  vertex2d1_pm, vertex2d2_pm,
					  vertex2P1_pm, vertex2P2_pm );
  return ridge_approximation.compute_max_ridges(it, order);  
}

template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap,
  class Vertex2VectorPropertyMap,
  class OutputIterator>
  OutputIterator compute_min_ridges(const TriangulatedSurfaceMesh &P,
				    const Vertex2FTPropertyMap& vertex2k1_pm, 
				    const Vertex2FTPropertyMap& vertex2k2_pm,
				    const Vertex2FTPropertyMap& vertex2b0_pm, 
				    const Vertex2FTPropertyMap& vertex2b3_pm,
				    const Vertex2VectorPropertyMap& vertex2d1_pm, 
				    const Vertex2VectorPropertyMap& vertex2d2_pm,
				    const Vertex2FTPropertyMap& vertex2P1_pm, 
				    const Vertex2FTPropertyMap& vertex2P2_pm,
				    OutputIterator it, 
				    Ridge_order order = Ridge_order_3)
{
  typedef Ridge_approximation < TriangulatedSurfaceMesh, 
    Vertex2FTPropertyMap, Vertex2VectorPropertyMap > Ridge_approximation;
  
  Ridge_approximation ridge_approximation(P, 
					  vertex2k1_pm, vertex2k2_pm,
					  vertex2b0_pm, vertex2b3_pm,
					  vertex2d1_pm, vertex2d2_pm,
					  vertex2P1_pm, vertex2P2_pm );
  return ridge_approximation.compute_min_ridges(it, order);  
}

template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap,
  class Vertex2VectorPropertyMap,
  class OutputIterator>
  OutputIterator compute_crest_ridges(const TriangulatedSurfaceMesh &P,
				    const Vertex2FTPropertyMap& vertex2k1_pm, 
				    const Vertex2FTPropertyMap& vertex2k2_pm,
				    const Vertex2FTPropertyMap& vertex2b0_pm, 
				    const Vertex2FTPropertyMap& vertex2b3_pm,
				    const Vertex2VectorPropertyMap& vertex2d1_pm, 
				    const Vertex2VectorPropertyMap& vertex2d2_pm,
				    const Vertex2FTPropertyMap& vertex2P1_pm, 
				    const Vertex2FTPropertyMap& vertex2P2_pm,
				    OutputIterator it, 
				    Ridge_order order = Ridge_order_3)
{
  typedef Ridge_approximation < TriangulatedSurfaceMesh, 
    Vertex2FTPropertyMap, Vertex2VectorPropertyMap > Ridge_approximation;
  
  Ridge_approximation ridge_approximation(P, 
					  vertex2k1_pm, vertex2k2_pm,
					  vertex2b0_pm, vertex2b3_pm,
					  vertex2d1_pm, vertex2d2_pm,
					  vertex2P1_pm, vertex2P2_pm );
  return ridge_approximation.compute_crest_ridges(it, order);  
}


} //namespace CGAL

#endif
