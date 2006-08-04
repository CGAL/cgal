#ifndef _UMBILIC_H_
#define _UMBILIC_H_

#include <math.h>
#include <CGAL/basic.h>
#include <CGAL/PolyhedralSurf_neighbors.h>

CGAL_BEGIN_NAMESPACE

enum Umbilic_type { NON_GENERIC = 0, WEDGE, TRISECTOR};

//-------------------------------------------------------------------
// Umbilic : stores umbilic data
//------------------------------------------------------------------
template < class Poly >
class Umbilic
{
public:
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef typename Poly::Traits::Vector_3 Vector_3;

public:
  Vertex_handle v;
  Umbilic_type umb_type;
  std::list<Halfedge_handle> contour;
  //contructor
  Umbilic(Vertex_handle v_init,
	  std::list<Halfedge_handle> contour_init) 
    : v(v_init), contour(contour_init) {}; 
};


//---------------------------------------------------------------------------
//Umbilic_approximation
//--------------------------------------------------------------------------
template < class Poly, class OutputIt, class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  class Umbilic_approximation
{
public:
  typedef typename Poly::Traits::FT FT;
  typedef typename Poly::Traits::Vector_3 Vector_3;
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef typename Poly::Facet_handle Facet_handle;
  typedef typename Poly::Facet_iterator Facet_iterator;
  typedef typename Poly::Vertex_iterator Vertex_iterator;
  typedef Umbilic<Poly> Umbilic;
  static FT neigh_size;//the size of neighbourhood for umbilic
  //  computation is (neigh_size * OneRingSize)
 
  Umbilic_approximation(Poly &P, 
			Vertex2FTPropertyMap vertex2k1_pm, Vertex2FTPropertyMap vertex2k2_pm,
			Vertex2VectorPropertyMap vertex2d1_pm, Vertex2VectorPropertyMap vertex2d2_pm);

  OutputIt compute(Poly &P, OutputIt it, FT size);

 protected:
  typedef T_PolyhedralSurf_neighbors<Poly> Poly_neighbors;
  Poly_neighbors* poly_neighbors;

  CGAL::Abs<FT> cgal_abs;

  //Property maps
  Vertex2FTPropertyMap k1, k2;
  Vertex2VectorPropertyMap d1, d2;

  // index: following CW the contour, we choose an orientation for the
  // max dir of an arbitrary starting point, the max dir field is
  // oriented on the next point so that the scalar product of the
  // consecutive vectors is positive.  Adding all the angles between
  // consecutive vectors around the contour gives -/+180 for a
  // wedge/trisector
  void compute_type(Umbilic& umb);
};

template < class Poly, class OutputIt, class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  Umbilic_approximation< Poly, OutputIt, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
Umbilic_approximation(Poly &P, 
		      Vertex2FTPropertyMap vertex2k1_pm, Vertex2FTPropertyMap vertex2k2_pm,
		      Vertex2VectorPropertyMap vertex2d1_pm, Vertex2VectorPropertyMap vertex2d2_pm)
  : k1(vertex2k1_pm), k2(vertex2k2_pm), 
    d1(vertex2d1_pm), d2(vertex2d2_pm)
{
  poly_neighbors = new Poly_neighbors(P);
}

template < class Poly, class OutputIt, class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  OutputIt Umbilic_approximation< Poly, OutputIt, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
compute(Poly &P, OutputIt umbilics_it, FT size)
{
  std::vector<Vertex_handle> vces;
  std::list<Halfedge_handle> contour;
  double umbilicEstimatorVertex, umbilicEstimatorNeigh;
  
  bool is_umbilic = true;

  //MAIN loop on P vertices
  Vertex_iterator itb = P.vertices_begin(), ite = P.vertices_end();
  for (;itb != ite; itb++) {
    Vertex_handle vh = itb;
    umbilicEstimatorVertex = cgal_abs(k1[vh]-k2[vh]);
    //reset vector, list and bool
    vces.clear();
    contour.clear();
    is_umbilic = true;
    poly_neighbors->compute_neighbors(vh, vces, contour, size);
    
    
    // OPTIONAL: avoid umbilics whose contours touch the border
    typename std::list<Halfedge_handle>::iterator itb_cont = contour.begin(),
      ite_cont = contour.end();
    for (; itb_cont != ite_cont; itb_cont++)
      if ( (*itb_cont)->is_border() ) {is_umbilic = false; continue;}
    if (is_umbilic == false) continue;
    
    //is v an umbilic?
    //a priori is_umbilic = true, and it switches to false as soon as a 
    //  neigh vertex has a lower umbilicEstimator value
    typename std::vector<Vertex_handle>::iterator itbv = vces.begin(),
      itev = vces.end();
    assert(*itbv ==  vh);
    itbv++;
    for (; itbv != itev; itbv++)
      {	umbilicEstimatorNeigh = cgal_abs( k1[*itbv] - k2[*itbv] );
	if ( umbilicEstimatorNeigh < umbilicEstimatorVertex ) 
	  {is_umbilic = false; break;}
      }
    if (is_umbilic == false) continue;
    
    //v is an umbilic, compute the index
    Umbilic*  cur_umbilic = new Umbilic(vh, contour);
    compute_type(*cur_umbilic);
    *umbilics_it++ = cur_umbilic;
  }
  return umbilics_it;
}

template < class Poly, class OutputIt, class Vertex2FTPropertyMap, class Vertex2VectorPropertyMap >
  void Umbilic_approximation< Poly, OutputIt, Vertex2FTPropertyMap, Vertex2VectorPropertyMap >::
compute_type(Umbilic& umb)
{
  Vector_3 dir, dirnext, normal;
  double cosinus, angle=0, angleSum=0;
  const double  pi=3.141592653589793;
  Vertex_handle v;
  typename std::list<Halfedge_handle>::iterator itb = umb.contour.begin(),
    itlast = --umb.contour.end();
  v = (*itb)->vertex();

  dir = d1[v];
  normal = CGAL::cross_product(d1[v], d2[v]);

  //sum angles along the contour
  do{
    itb++;
    v=(*itb)->vertex();
    dirnext = d1[v];
    cosinus = dir*dirnext;
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
  v=(*umb.contour.begin())->vertex();
   dirnext = d1[v];
  cosinus = dir*dirnext;
  if (cosinus < 0) {dirnext = dirnext*(-1); cosinus *= -1;}
  if (cosinus>1) cosinus = 1;
  if ( (dir * CGAL::cross_product(dirnext, normal)) > 0) angle = acos(cosinus);
  else angle = -acos(cosinus);
  angleSum += angle;

  if ((angleSum > (pi/2)) && (angleSum < (3*pi/2))) umb.umb_type = TRISECTOR ;
  else if ((angleSum < (-pi/2)) && (angleSum > (-3*pi/2))) umb.umb_type = WEDGE;
  else umb.umb_type = NON_GENERIC;
}

CGAL_END_NAMESPACE

#endif
