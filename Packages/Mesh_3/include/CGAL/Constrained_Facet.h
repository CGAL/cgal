#ifndef FACET
#define FACET

#include <string>
#include <CGAL/basic.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <utility>
#include <map> 
#include <set>
#include <list> 

#include "Constrained_Element.h"
#include "Constrained_Vertex.h"
#include "Constrained_Edge.h"

#include "CGAL/squared_distance_3.h"


CGAL_BEGIN_NAMESPACE

template <class Gt>
class Constrained_Facet : public  virtual  Constrained_Element<Gt>
{

public :
	typedef Gt                                       Geom_traits;
	typedef typename Gt::FT                          Nb;

	typedef typename Gt::Point_3                     Point; 
	typedef typename Gt::Vector_3                    Vector;
	typedef typename Gt::Segment_3                   Segment;
	typedef typename Gt::Iso_cuboid_3                Iso_Cuboid;
	typedef typename Gt::Tetrahedron_3               Tetrahedron;
	typedef typename Gt::Triangle_3                  Triangle;
	typedef typename Gt::Plane_3                     Plane;
	typedef typename Gt::Line_3                      Line;
	typedef typename Gt::Object_3                    Object;
	typedef typename Gt::Ray_3                       Ray;

	// typedef typename Gt::Bbox_3                   Bbox_3;
	typedef CGAL::Bbox_3                             Bbox_3;

	typedef Constrained_Element<Gt>              C_elt;
	typedef Constrained_Vertex<Gt>               C_vertex;

	typedef Constrained_Edge<Gt>                 C_edge;


private :
	Geom_traits       gt;

public :
	Triangle t;

public :  
	//constructeur
	Constrained_Facet(): C_elt("facet"){ }

	Constrained_Facet( Triangle tr) 
		:	C_elt("face")
	{
		t = tr;
	} 

	Constrained_Facet( Point p, Point q, Point r) : C_elt("facet")
	{
		typename Geom_traits::Construct_triangle_3  construct_triangle =
			gt.construct_triangle_3_object();
		t = construct_triangle(p,q,r) ;
	} 


	//~Constrained_Facet();


	Triangle getTriangle() { return t; }


	//Declaration Des Fonctions 
	Bbox_3           getBBox();
	bool             does_intersect					( Iso_Cuboid c);
	bool             does_intersect					( Tetrahedron tet)   ;
	bool             does_intersect					( Segment seg );
	bool             does_intersect					( Ray ray );
	bool             does_intersect					( C_elt* ce);

	bool             does_intersect_vertex	(	C_vertex* s);
	bool             does_intersect_edge		(	C_edge* a);
	bool             does_intersect_facet		(	Constrained_Facet* f);

	bool             do_intersect2					(	Iso_Cuboid c, Triangle t);
	bool             do_intersect2					( Tetrahedron tet, Triangle t);


	Nb               squared_distance2			(	Point p,Triangle t); 

	Nb               sqared_distance_to_sheet(	Point p);

	Object           intersection						 (	Segment seg); 
	Object           intersection						 (	Ray ray);

	bool             do_intersect_on_a_point (	Triangle t,Segment s);
	bool             do_intersect_on_a_point (	Triangle t,Ray s);
};


//-------------------------------------------------------------------
// fonction qui retourne la bounding box de la contrainte 
template <class Gt>
typename Constrained_Facet<Gt> :: Bbox_3 
Constrained_Facet<Gt> :: getBBox()  { // ????????
	return t.bbox();
}
//-------------------------------------------------------------------
// fonction qui retourne true si le triangle t intersecte
// l'isocuboide c
template <class Gt>
bool
Constrained_Facet<Gt> :: do_intersect2(Iso_Cuboid c, Triangle t) 
{ 

	typename Geom_traits::Construct_vertex_3  construct_vertex =
		gt.construct_vertex_3_object();

	typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
		=gt.has_on_bounded_side_3_object();

	typename Geom_traits::Construct_triangle_3  construct_triangle =
		gt.construct_triangle_3_object();

	typename Geom_traits:: Has_on_boundary_3  has_on_boundary
		=gt.has_on_boundary_3_object();

	//std :: cout<< "do_intersect2 "<<endl;

	// on commence par tester si un des sommets du triangle  est
	// dans c
	Point st[3];
	for (int i=0; i<3; i++){
		st[i] = construct_vertex(t,i);
		if (has_on_bounded_side(c,st[i]) ||  has_on_boundary(c,st[i]) )
		{
			//std::cerr<<"(_)ooo"<<endl;
			return true;
		}
	}


	//auccun des sommets n'est dans c.
	// on construit les triangles qui forment les faces de c

	Point sc[8];
	for (int i=0; i<8;i++)
		sc[i] = construct_vertex(c,i);

	Triangle tf;

	//on va regarder pour toutes les facets du cube.
	//on coupe une facet en deux triangle
	tf = construct_triangle(sc[1],sc[3],sc[0]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[1],sc[3],sc[2]);
	if (do_intersect(tf,t)) return true;

	tf = construct_triangle(sc[1],sc[6],sc[0]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[0],sc[6],sc[5]);
	if (do_intersect(tf,t)) return true;

	tf = construct_triangle(sc[1],sc[7],sc[6]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[1],sc[7],sc[2]);
	if (do_intersect(tf,t)) return true;

	tf = construct_triangle(sc[2],sc[4],sc[3]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[2],sc[4],sc[7]);
	if (do_intersect(tf,t)) return true;     

	tf = construct_triangle(sc[3],sc[5],sc[0]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[3],sc[5],sc[4]);
	if (do_intersect(tf,t)) return true;

	tf = construct_triangle(sc[5],sc[4],sc[6]);
	if (do_intersect(tf,t)) return true;
	tf = construct_triangle(sc[6],sc[4],sc[7]);
	if (do_intersect(tf,t)) return true;

	return false;

}


//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte un
//isocuboide donne C
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect(Iso_Cuboid c) 
{ 

	typename Geom_traits::Construct_triangle_3
		construct_triangle =  gt.construct_triangle_3_object();

	return do_intersect2(c,t);

}
//-------------------------------------------------------------------
// fonction qui retourne true si le triangle t intersecte
// l'isocuboide c
template <class Gt>
bool
Constrained_Facet<Gt> :: do_intersect2(Tetrahedron tet, Triangle t) 
{ 

	typename Geom_traits::Construct_vertex_3  construct_vertex =
		gt.construct_vertex_3_object();

	typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
		=gt.has_on_bounded_side_3_object();

	typename Geom_traits::Construct_triangle_3  construct_triangle =
		gt.construct_triangle_3_object();

	typename Geom_traits:: Has_on_boundary_3  has_on_boundary
		=gt.has_on_boundary_3_object();


	// on commence par tester si un des sommets du triangle  est
	// dans c
	Point st[3];
	for (int i=0; i<3; i++){
		st[i] = construct_vertex(t,i);
		if (has_on_bounded_side(tet,st[i]) ||  has_on_boundary(tet,st[i]) )
		{
			return true;
		}
	}

	//auccune des extremites n'est dans tet
	// on construit les triangles qui forment les facets de tet
	Point stet[4];
	for (int i=0; i<4;i++)
		stet[i] = construct_vertex(tet,i);

	Triangle tf;

	// face 3
	tf = construct_triangle(stet[0],stet[1],stet[2]);
	if (do_intersect(tf,t)) return true;

	// face 2
	tf = construct_triangle(stet[0],stet[1],stet[3]);
	if (do_intersect(tf,t)) return true;

	// face 1
	tf = construct_triangle(stet[0],stet[2],stet[3]);
	if (do_intersect(tf,t)) return true;

	// face 0
	tf = construct_triangle(stet[1],stet[2],stet[3]);
	if (do_intersect(tf,t)) return true;

	return false;
}

//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//cellule donnee c
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect( Tetrahedron tet)  
{ 


	return (do_intersect2(tet,t));

}

//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//cellule donnee c
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect(Segment seg)  
{ 

	//std::cout<<t<<endl;

	return (do_intersect(t,seg));


}
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//cellule donnee c
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect(Ray ray)  
{ 

	//std::cout<<t<<endl;

	return (do_intersect_on_a_point(t,ray));


}
//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect_vertex(C_vertex* s) {

	return s->does_intersect_facet(this);
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect_edge(C_edge* a) {
	//std::cerr<<"****"<<endl;
	return a->does_intersect_facet(this);
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect_facet(Constrained_Facet*  facet) {

	std::cout<<" dans intersection facet-facet"<<std::endl;

	std::list<int> l1,l2;
	Triangle t_facet = facet->getTriangle();

	typename Geom_traits:: Equal_3  equal
		=gt.equal_3_object();
	return equal(t_facet,t);

}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Facet<Gt> :: does_intersect( C_elt* ce) 
{

	char* type = ce->getType();
	/*std::cout<<"intersection facet-"<<type <<std::endl;
	std::cout<<"C_elt before " <<type<<" ? "<<typeid(C_elt).before(typeid(Constrained_Facet))<<std::endl;;
	std::cout<<"type de ce "<<typeid(ce).name()<<std::endl;*/

	if ((strcmp(type,"vertex") == 0))
		return does_intersect_vertex(dynamic_cast<C_vertex*> (ce) );
	else  if ((strcmp(type,"edge") == 0))
		return does_intersect_edge(dynamic_cast<C_edge*> (ce) );
	else 
	{
		Constrained_Facet* cf = reinterpret_cast<Constrained_Facet*> (ce);
		if (cf)
		{
			//std::cout<<"type de cf "<<typeid(cf).name()<<std::endl;
			return does_intersect_facet(cf);
		}
		else
		{
			//std::cout<<"pb !!! "<<std::endl;
			return false;
		}
	}


}




//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Facet<Gt> :: Nb 
Constrained_Facet<Gt> :: squared_distance2(Point p,
																					 Triangle t)
{
	typename Geom_traits::Construct_vertex_3  construct_vertex =
		gt.construct_vertex_3_object(); 

	typename Geom_traits::Construct_segment_3 
		construct_segment =  gt.construct_segment_3_object();

	typename Geom_traits:: Compute_squared_distance_3
		compute_squared_distance = gt.compute_squared_distance_3_object();

	Nb res,restemp;
	Point 
		a0 = construct_vertex(t,0),
		a1 = construct_vertex(t,1),
		a2 = construct_vertex(t,2);

	Vector v0,v1,v2;
	v0 = a2 - a1;
	v1 = a0 - a2; 
	v2 = a1 - a0; 

	//normal of the plane a0a1a2
	Vector N = cross_product(v0,v1);

	//let p' the projected point of p on the plane a0a1a2

	bool b0=0,b1=0,b2=0;
	//test if p' is anticlockwise of v0
	if (cross_product(v0,p - a0)*N >0)
		b0 = 1;
	//test if p' is anticlockwise of v1
	if (cross_product(v1,p - a1)*N >0)
		b0 = 1;
	//test if p' is anticlockwise of v2
	if (cross_product(v2,p - a2)*N >0)
		b0 = 1;

	if (b0*b1*b2) {
		res = (N*(p - a1)*
			(N*(p - a1)/(N*N)));
	}
	else{
		Segment
			e1 = construct_segment (a0,a1),
			e2 = construct_segment (a0,a2),
			e3 = construct_segment (a1,a2);

		res = squared_distance(p,e1);
		restemp = squared_distance(p,e2);

		if (res > restemp) res = restemp;

		restemp =  squared_distance(p,e3);
		if (res > restemp) res = restemp; 
	}
	return res;

}



//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Facet<Gt> :: Nb 
Constrained_Facet<Gt> :: sqared_distance_to_sheet(Point p)
{
	return  squared_distance2(p,t);
}

////----------------------------------------------------------------------------
//fonction qui construit le point d'intersection entre le segment et
//la face
template <class Gt>
typename Constrained_Facet<Gt> :: Object
Constrained_Facet<Gt> :: intersection (Segment seg) {

	//ATTENTION !! je sais que le segment intersecte la facet

	typename Geom_traits::Construct_plane_3  construct_plane =
		gt.construct_plane_3_object();
	typename Geom_traits::Construct_vertex_3  construct_vertex =
		gt.construct_vertex_3_object();
	typename Geom_traits::Intersect_3  intersect =
		gt.intersect_3_object();


	Point p,q,r;
	p = construct_vertex(t,0);
	q = construct_vertex(t,1);
	r = construct_vertex(t,2);

	Plane plan = construct_plane(p,q,r);
	Object o = intersect(plan,seg);

	return  o; 

}
////----------------------------------------------------------------------------
//fonction qui construit le point d'intersection entre le rayon et
//la facet
template <class Gt>
typename Constrained_Facet<Gt> :: Object
Constrained_Facet<Gt> :: intersection (Ray ray) {

	//ATTENTION !! je sait que le segment intersecte la facet

	typename Geom_traits::Construct_plane_3  construct_plane =
		gt.construct_plane_3_object();
	typename Geom_traits::Construct_vertex_3  construct_vertex =
		gt.construct_vertex_3_object();
	typename Geom_traits::Intersect_3  intersect =
		gt.intersect_3_object();


	Point p,q,r;
	p = construct_vertex(t,0);
	q = construct_vertex(t,1);
	r = construct_vertex(t,2);

	Plane plan = construct_plane(p,q,r);
	Object o = intersect(plan,ray);

	return  o; 

}
////----------------------------------------------------------------------------
//fonction qui construit le point d'intersection entre le segment et
//la facet
template <class Gt>
bool
Constrained_Facet<Gt> ::  do_intersect_on_a_point (Triangle t,Segment s)
{
	int i,j,k;
	Geom_traits ker;

	// On regarde si le segment intersecte le plan support du triangle

	i=ker.orientation_3_object()(t.vertex(0),t.vertex(1),t.vertex(2),s.vertex(0));
	j=ker.orientation_3_object()(t.vertex(0),t.vertex(1),t.vertex(2),s.vertex(1));

	if (i*j>=0)
		return(false);


	// On regarde si l'intersection est dans le triangle

	i=ker.orientation_3_object()(s.vertex(0),s.vertex(1),t.vertex(0),t.vertex(2));

	j=ker.orientation_3_object()(s.vertex(0),s.vertex(1),t.vertex(2),t.vertex(1));

	k=ker.orientation_3_object()(s.vertex(0),s.vertex(1),t.vertex(1),t.vertex(0));

	return(i*j>=0 && i*k>=0 && j*k>=0 && (i!=0 || j!=0 || k!=0));
}

////----------------------------------------------------------------------------
//fonction qui construit le point d'intersection entre le rayon et
//la facet
template <class Gt>
bool
Constrained_Facet<Gt> ::  do_intersect_on_a_point (Triangle t,Ray s)
{
	int i,j,k;
	Point v1,v2;
	Geom_traits ker;

	//      cout << "rayon " << s << "-> ";

	// On prend deux points arbitraires sur le rayon
	v1=s.point(0);
	v2=s.point(1);

	// On regarde si le rayon intersecte le plan support du triangle
	if (ker.do_intersect_3_object()(t.supporting_plane(),s)==false)
	{
		//	  cout << "le rayon n'intersecte pas le plan " << t << endl;
		return(false);
	}

	//      cout << "le rayon intersecte le plan " << t << "-> ";

	// On regarde si l'intersection est dans le triangle
	i=ker.orientation_3_object()(v1,v2,t.vertex(0),t.vertex(2));
	j=ker.orientation_3_object()(v1,v2,t.vertex(2),t.vertex(1));
	k=ker.orientation_3_object()(v1,v2,t.vertex(1),t.vertex(0));

	//      cout << "(i=" << i << ",j=" << j << ",k=" << k << ")" << endl;

	//      if (i*j>=0 && i*k>=0 && j*k>=0)
	//	cout << "youpi" << endl;
	//      else
	//	cout << "perdu" << endl;

	return(i*j>=0 && i*k>=0 && j*k>=0 && (i!=0 || j!=0 || k!=0));
}


CGAL_END_NAMESPACE

#endif
