#ifndef _UMBILIC_H_
#define _UMBILIC_H_

#include <math.h>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

enum Umbilic_type { NON_GENERIC = 0, WEDGE, TRISECTOR};

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
	  std::list<Halfedge_handle> contour_init); 
  // index: following CW the contour, we choose an orientation for the
  // max dir of an arbitrary starting point, the max dir field is
  // oriented on the next point so that the scalar product of the
  // consecutive vectors is positive.  Adding all the angles between
  // consecutive vectors around the contour gives -/+180 for a
  // wedge/trisector
  void compute_type();
};

//-------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------

template < class Poly >
Umbilic<Poly >::
Umbilic(Vertex_handle v_init,
	std::list<Halfedge_handle> contour_init)
  : v(v_init), contour(contour_init)
{}

template < class Poly >
void Umbilic<Poly>::
compute_type()
{
  Vector_3 dir, dirnext, normal;
  double cosinus, angle=0, pi=3.141592653589793, angleSum=0;
  Vertex_handle v;
  typename std::list<Halfedge_handle>::iterator itb = contour.begin(),
    itlast = --contour.end();
  v = (*itb)->vertex();

  dir = v->d1();
  normal = CGAL::cross_product(v->d1(), v->d2());

  //sum angles along the contour
  do{
    itb++;
    v=(*itb)->vertex();
    dirnext = v->d1();
    cosinus = dir*dirnext;
    if (cosinus < 0) {dirnext = dirnext*(-1); cosinus *= -1;}
    if (cosinus>1) cosinus = 1;
    //orientation of (dir, dirnext, normal)
    if ( (dir * CGAL::cross_product(dirnext, normal)) > 0) angle = acos(cosinus);
    else angle = -acos(cosinus);
    angleSum += angle;
    dir = dirnext;
   normal = CGAL::cross_product(v->d1(), v->d2());
  }
  while (itb != (itlast));
  
  //angle (v_last, v_0)
  v=(*contour.begin())->vertex();
   dirnext = v->d1();
  cosinus = dir*dirnext;
  if (cosinus < 0) {dirnext = dirnext*(-1); cosinus *= -1;}
  if (cosinus>1) cosinus = 1;
  if ( (dir * CGAL::cross_product(dirnext, normal)) > 0) angle = acos(cosinus);
  else angle = -acos(cosinus);
  angleSum += angle;

  if ((angleSum > (pi/2)) && (angleSum < (3*pi/2))) umb_type = TRISECTOR ;
  else if ((angleSum < (-pi/2)) && (angleSum > (-3*pi/2))) umb_type = WEDGE;
  else umb_type = NON_GENERIC;
}

CGAL_END_NAMESPACE

#endif
