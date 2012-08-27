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
#ifndef CGAL_UMBILIC_H_
#define CGAL_UMBILIC_H_

#include <list>
#include <vector>
#include <math.h>
#include <CGAL/basic.h>
#include <CGAL/PolyhedralSurf_neighbors.h>

namespace CGAL {

enum Umbilic_type { NON_GENERIC_UMBILIC = 0, ELLIPTIC_UMBILIC, HYPERBOLIC_UMBILIC};

//-------------------------------------------------------------------
//Umbilic : stores umbilic data, its location given by a vertex, its
//type and a circle of edges bording a disk containing the vertex
//------------------------------------------------------------------
template < class TriangulatedSurfaceMesh >
class Umbilic
{
 public:
  typedef typename TriangulatedSurfaceMesh::Vertex_const_handle    Vertex_const_handle;
  typedef typename TriangulatedSurfaceMesh::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3 Vector_3;
  
  //contructor
  Umbilic(const Vertex_const_handle v_init,
	  const std::list<Halfedge_const_handle> contour_init); 
  //access fct
  Vertex_const_handle vertex() const { return v;}
  Umbilic_type umbilic_type() const { return umb_type;}
  Umbilic_type& umbilic_type() { return umb_type;}
  const std::list<Halfedge_const_handle>& contour_list() const { return contour;}

 protected:
  const Vertex_const_handle v;
  Umbilic_type umb_type;
  const std::list<Halfedge_const_handle> contour;
};

//constructor
template <class TriangulatedSurfaceMesh>
Umbilic<TriangulatedSurfaceMesh>::
Umbilic(const Vertex_const_handle v_init,
	const std::list<Halfedge_const_handle> contour_init) 
  : v(v_init), contour(contour_init) {} 


template <class TriangulatedSurfaceMesh>
std::ostream& 
operator<<(std::ostream& out_stream, const Umbilic<TriangulatedSurfaceMesh>& umbilic)
{
  out_stream << "Umbilic at location (" << umbilic.vertex()->point() << ") of type ";
  switch (umbilic.umbilic_type())
    {
    case CGAL::NON_GENERIC_UMBILIC: out_stream << "non generic" << std::endl; break;
    case CGAL::ELLIPTIC_UMBILIC: out_stream << "elliptic" << std::endl; break;
    case CGAL::HYPERBOLIC_UMBILIC: out_stream << "hyperbolic" << std::endl; break;
    default : out_stream << "Something wrong occured for sure..." << std::endl; break;
    }
  return out_stream;
}
//---------------------------------------------------------------------------
//Umbilic_approximation : enable computation of umbilics of a
//TriangulatedSurfaceMesh. It uses the class
//T_PolyhedralSurf_neighbors to compute topological disk patches
//around vertices
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  class Umbilic_approximation
{
 public:
  typedef typename TriangulatedSurfaceMesh::Traits::FT       FT;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3 Vector_3;
  typedef typename TriangulatedSurfaceMesh::Vertex_const_handle    Vertex_const_handle;
  typedef typename TriangulatedSurfaceMesh::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename TriangulatedSurfaceMesh::Facet_const_iterator   Facet_const_iterator;
  typedef typename TriangulatedSurfaceMesh::Vertex_const_iterator  Vertex_const_iterator;

  //requirements for the templates TriangulatedSurfaceMesh and Vertex2FTPropertyMap or Vertex2VectorPropertyMap
  CGAL_static_assertion((boost::is_same<Vertex_const_handle, typename Vertex2FTPropertyMap::key_type>::value));
  CGAL_static_assertion((boost::is_same<Vertex_const_handle, typename Vertex2VectorPropertyMap::key_type>::value));
  CGAL_static_assertion((boost::is_same<FT, typename Vertex2FTPropertyMap::value_type>::value));
  CGAL_static_assertion((boost::is_same<Vector_3, typename Vertex2VectorPropertyMap::value_type>::value));

  typedef CGAL::Umbilic<TriangulatedSurfaceMesh> Umbilic;

  //constructor : sets propertymaps and the poly_neighbors
  Umbilic_approximation(const TriangulatedSurfaceMesh& P, 
			const Vertex2FTPropertyMap& vertex2k1_pm, 
			const Vertex2FTPropertyMap& vertex2k2_pm,
			const Vertex2VectorPropertyMap& vertex2d1_pm, 
			const Vertex2VectorPropertyMap& vertex2d2_pm);
  //identify umbilics as vertices minimizing the function k1-k2 on
  //their patch and for which the index is not 0. We avoid
  //potential umbilics whose contours touch the border.
  template <class OutputIterator>
  OutputIterator compute(OutputIterator it, FT size);

 protected:
  const TriangulatedSurfaceMesh& P;
  
  typedef T_PolyhedralSurf_neighbors<TriangulatedSurfaceMesh> Poly_neighbors;
  Poly_neighbors* poly_neighbors;

  CGAL::Abs<FT> cgal_abs;
  CGAL::To_double<FT> To_double;

  //Property maps
  const Vertex2FTPropertyMap &k1, &k2;
  const Vertex2VectorPropertyMap &d1, &d2;

  // index: following CW the contour, we choose an orientation for the
  // max dir of an arbitrary starting point, the max dir field is
  // oriented on the next point so that the scalar product of the
  // consecutive vectors is positive.  Adding all the angles between
  // consecutive vectors around the contour gives ~ -/+180 for a
  // wedge/trisector, ~ 0 gives a false umbilic, everything else gives
  // a non_generic umbilic.
  int compute_type(Umbilic& umb);
};

template < class TriangulatedSurfaceMesh,  class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  Umbilic_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
Umbilic_approximation(const TriangulatedSurfaceMesh& p, 
		      const Vertex2FTPropertyMap& vertex2k1_pm, 
		      const Vertex2FTPropertyMap& vertex2k2_pm,
		      const Vertex2VectorPropertyMap& vertex2d1_pm, 
		      const Vertex2VectorPropertyMap& vertex2d2_pm)
  : P(p), k1(vertex2k1_pm), k2(vertex2k2_pm), 
    d1(vertex2d1_pm), d2(vertex2d2_pm)
{
  //check that the mesh is a triangular one.
  Facet_const_iterator itb = P.facets_begin(), ite = P.facets_end();
  for(;itb!=ite;itb++) CGAL_precondition( itb->is_triangle() );

  poly_neighbors = new Poly_neighbors(P);
}

template < class TriangulatedSurfaceMesh,  class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  template <class OutputIterator>
  OutputIterator Umbilic_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
  compute(OutputIterator umbilics_it, FT size)
{
  CGAL_precondition( size >= 1 );
  
  std::vector<Vertex_const_handle> vces;
  std::list<Halfedge_const_handle> contour;
  FT umbilicEstimatorVertex, umbilicEstimatorNeigh;
  
  bool is_umbilic = true;

  //MAIN loop on P vertices
  Vertex_const_iterator itb = P.vertices_begin(), ite = P.vertices_end();
  for (;itb != ite; itb++) {
    Vertex_const_handle vh = itb;
    umbilicEstimatorVertex = cgal_abs(k1[vh]-k2[vh]);
    //reset vector, list and bool
    vces.clear();
    contour.clear();
    is_umbilic = true;
    //the size of neighbourhood is (size * OneRingSize)
    poly_neighbors->compute_neighbors(vh, vces, contour, size);
    
    
    // avoid umbilics whose contours touch the border (Note may not be
    // desirable?)
    typename std::list<Halfedge_const_handle>::const_iterator itb_cont = contour.begin(),
      ite_cont = contour.end();
    for (; itb_cont != ite_cont; itb_cont++)
      if ( (*itb_cont)->is_border() ) {is_umbilic = false; continue;}
    if (is_umbilic == false) continue;
    
    //is v an umbilic?
    //a priori is_umbilic = true, and it switches to false as soon as a 
    //  neigh vertex has a lower umbilicEstimator value
    typename std::vector<Vertex_const_handle>::const_iterator itbv = vces.begin(),
      itev = vces.end();
    itbv++;
    for (; itbv != itev; itbv++)
      {	umbilicEstimatorNeigh = cgal_abs( k1[*itbv] - k2[*itbv] );
	if ( umbilicEstimatorNeigh < umbilicEstimatorVertex ) 
	  {is_umbilic = false; break;}
      }
    if (is_umbilic == false) continue;
    
    //v is an umbilic (wrt the min of k1-k2), compute the index. If
    //the index is not 0 then we have actually an umbilic which is output
    Umbilic*  cur_umbilic = new Umbilic(vh, contour);
    if (compute_type(*cur_umbilic) != 0)  *umbilics_it++ = cur_umbilic;
  }
  return umbilics_it;
}

template < class TriangulatedSurfaceMesh,  class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  int Umbilic_approximation< TriangulatedSurfaceMesh, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
  compute_type(Umbilic& umb)
{
  Vector_3 dir, dirnext, normal;
  double cosinus, angle=0, angleSum=0;
  const double  pi=3.141592653589793;
  Vertex_const_handle v;
  typename std::list<Halfedge_const_handle>::const_iterator itb = umb.contour_list().begin(),
    itlast = --umb.contour_list().end();
  v = (*itb)->vertex();

  dir = d1[v];
  normal = CGAL::cross_product(d1[v], d2[v]);

  //sum angles along the contour
  do{
    itb++;
    v=(*itb)->vertex();
    dirnext = d1[v];
    cosinus = To_double(dir*dirnext);
    if (cosinus < 0) {dirnext = dirnext*(-1); cosinus *= -1;}
    if (cosinus>1) cosinus = 1;
    //orientation of (dir, dirnext, normal)
    if ( (dir * CGAL::cross_product(dirnext, normal)) > 0) angle = acos(cosinus);
    else angle = -acos(cosinus);
    angleSum += angle;
    dir = dirnext;
    normal = CGAL::cross_product(d1[v], d2[v]);
  }
  while (itb != (itlast));
  
  //angle (v_last, v_0)
  v=(*umb.contour_list().begin())->vertex();
   dirnext = d1[v];
   cosinus = To_double(dir*dirnext);
  if (cosinus < 0) {dirnext = dirnext*(-1); cosinus *= -1;}
  if (cosinus>1) cosinus = 1;
  if ( (dir * CGAL::cross_product(dirnext, normal)) > 0) angle = acos(cosinus);
  else angle = -acos(cosinus);
  angleSum += angle;

  if ((angleSum > (pi/2)) && (angleSum < (3*pi/2))) umb.umbilic_type() = HYPERBOLIC_UMBILIC ;
  else if ((angleSum < (-pi/2)) && (angleSum > (-3*pi/2))) umb.umbilic_type() = ELLIPTIC_UMBILIC;
  else if ((angleSum <= (pi/2)) && (angleSum >= (-pi/2))) return 0;//is not considered as an umbilic
  else umb.umbilic_type() = NON_GENERIC_UMBILIC;
  return 1;
}

//Global function

template < class TriangulatedSurfaceMesh,  
  class Vertex2FTPropertyMap,
  class Vertex2VectorPropertyMap,
  class OutputIterator>
  OutputIterator compute_umbilics(const TriangulatedSurfaceMesh &P,
				  const Vertex2FTPropertyMap& vertex2k1_pm, 
				  const Vertex2FTPropertyMap& vertex2k2_pm,
				  const Vertex2VectorPropertyMap& vertex2d1_pm, 
				  const Vertex2VectorPropertyMap& vertex2d2_pm,
				  OutputIterator it, 
				  double size)
{
  typedef Umbilic_approximation < TriangulatedSurfaceMesh, 
    Vertex2FTPropertyMap, Vertex2VectorPropertyMap > Umbilic_approximation;
  
  Umbilic_approximation umbilic_approximation(P, 
					      vertex2k1_pm, vertex2k2_pm,
					      vertex2d1_pm, vertex2d2_pm);
  return umbilic_approximation.compute(it, size);  
}




} //namespace CGAL

#endif
