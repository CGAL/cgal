#ifndef CONTRAINTE
#define CONTRAINTE

#include <CGAL/basic.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath> 
#include <utility>
#include <map> 
#include <set>
#include <list>

CGAL_BEGIN_NAMESPACE

template <class Gt>
class Constrained_Element
{
  
public:  
  typedef Gt                                       Geom_traits_;
  typedef typename Gt::FT                          Nb;

  typedef typename Gt::Point_3                     Point; 
  typedef typename Gt::Vector_3                    Vector;
  typedef typename Gt::Segment_3                   Segment;
  typedef typename Gt::Iso_cuboid_3                Iso_Cuboid;
  typedef typename Gt::Tetrahedron_3               Tetrahedron;
  // structure de donnees pour la liste d'incidence
  typedef std::map < Constrained_Element ,int >             Contrainte_map;
  typedef typename Contrainte_map::iterator        Contrainte_it;

private: 
  typedef typename Gt::Triangle_3                             Triangle;
  typedef typename Gt::Plane_3                                Plane;
  typedef typename Gt::Ray_3                                  Ray;
 
  typedef CGAL::Bbox_3                              Bbox_3;
  
  typedef char* Mot;

public :
  char* type;
  
public :

  Constrained_Element(){ type = new char; }
  
  Constrained_Element(char* nom) { type = new char; strcpy(type,nom);}
  
  virtual ~Constrained_Element(){}
 
  //les methodes virtuelles
  virtual Bbox_3  getBBox() = 0;
  virtual bool    does_intersect( Iso_Cuboid c) = 0;
  virtual bool    does_intersect( Tetrahedron tet)   = 0;
  virtual bool    does_intersect( Constrained_Element* ce) = 0;
  virtual Nb      sqared_distance_to_sheet(Point p) = 0;
  virtual bool    does_intersect( Segment seg )   = 0;
  virtual Object  intersection (Segment seg) = 0;
  virtual bool    does_intersect(Ray ray )   = 0;
  virtual Object  intersection (Ray ray) = 0;
  
  void setType(char *nom){type = nom;}
  Mot getType() { return type;}
};

CGAL_END_NAMESPACE

#endif
 
