// we have to provide CGAL-compatible circles ...

#ifndef CEP_LEDA_RAT_CIRCLE_H
#define CEP_LEDA_RAT_CIRCLE_H

// use CGALs orientation enum ...
#include <CGAL/enum.h>

#include <LEDA/rat_point.h>


#if !defined(LEDA_BEGIN_NAMESPACE)
#define LEDA_BEGIN_NAMESPACE
#endif

#if !defined(LEDA_END_NAMESPACE)
#define LEDA_END_NAMESPACE
#endif


LEDA_BEGIN_NAMESPACE

class __exportC cgal_rat_circle_rep  : public handle_rep {

   friend class __exportC cgal_rat_circle;

   static leda_mutex mutex_id_counter;
   static unsigned long id_counter;

   leda_rat_point center;       // center point
   leda_rational  sq_rad;       // squared radius
   CGAL::Orientation       orientation;  // orientation of the circle ...

   unsigned long id;

public:

   cgal_rat_circle_rep() : center(0,0,1), sq_rad(0), orientation(CGAL::COUNTERCLOCKWISE)
   { 
     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }
   
   // p1 is the center ...
   cgal_rat_circle_rep(const leda_rat_point& p1, 
                       const leda_rat_point& p2,
		       CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE)
   {
     center = p1;
     sq_rad = p1.sqr_dist(p2);
     orientation = ori;
          
     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();               
   }
   
   // p1 is the center ...
   cgal_rat_circle_rep(const leda_rat_point& p1, 
                       const leda_rational&  sq,
		       CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE)
   {
     center = p1;
     sq_rad = sq;
     orientation = ori;
          
     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();               
   }   

  ~cgal_rat_circle_rep() {}

   friend inline unsigned long ID_Number(const cgal_rat_circle& c);

};

#if !defined(CEP_LEDA_RAT_CIRCLE_NO_STATIC_INIT)

leda_mutex cgal_rat_circle_rep::mutex_id_counter;

unsigned long cgal_rat_circle_rep::id_counter = 0;

#endif


#if (__LEDA__ >= 440)
#define __HANDLE_BASE_TYPE  HANDLE_BASE(cgal_rat_circle_rep)
#else
#define __HANDLE_BASE_TYPE  handle_base
#endif 


/*{\Manpage {cgal_rat_circle} {} {Rational circles} {C}}*/

class __exportC cgal_rat_circle : public __HANDLE_BASE_TYPE {

/*{\Mdefinition
An instance of the data type |\Mname| is an oriented circle in the plane. 
The circle is defined by a center with rational coordinates 
(|rat_point|) and a squared radius. 
}*/

//#if (__LEDA__ < 440)
  cgal_rat_circle_rep* ptr() const { return (cgal_rat_circle_rep*)PTR; } 
//#endif  
  
public:

/*{\Mcreation}*/

  // constructors ...
  cgal_rat_circle() { PTR = new cgal_rat_circle_rep; }
/*{\Mcreate introduces a variable |C| of type |\Mname|.}*/   
  
  cgal_rat_circle(const leda_rat_point& p1, const leda_rat_point& p2, 
                  CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE)
/*{\Mcreate introduces a variable |C| of type |\Mname| with center |p1|, a point
|p2| on the circle and orientation |ori|. }*/	
  { PTR = new cgal_rat_circle_rep(p1,p2,ori); }
  
  cgal_rat_circle(const leda_rat_point& p1, 
                  CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE)
/*{\Mcreate introduces a variable |C| of type |\Mname| with center |p1|, orientation |ori|
and radius = 0.
}*/	
  { PTR = new cgal_rat_circle_rep(p1,p1,ori); } 
  
  cgal_rat_circle(const leda_rat_point& p1, const leda_rational& sqrad, 
                  CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE)
/*{\Mcreate introduces a variable |C| of type |\Mname| with center |p1|, orientation |ori|
and squared radius |sqrad|.
}*/
  { PTR = new cgal_rat_circle_rep(p1,sqrad,ori); }   

  // copy constructor ...
  cgal_rat_circle(const cgal_rat_circle& c) : __HANDLE_BASE_TYPE(c) {}
  ~cgal_rat_circle() {}

  cgal_rat_circle& operator=(const cgal_rat_circle& c)
  { __HANDLE_BASE_TYPE::operator=(c); return *this; }


/*{\Moperations 4 5}*/  
  
  leda_rat_point  center() const 
/*{\Mop returns the center of the circle.}*/    
  { return ptr()->center; }
  
  leda_rational   squared_radius() const 
/*{\Mop returns the square of the radius.}*/   
  { return ptr()->sq_rad; }

  leda_rational   sqr_radius() const 
/*{\Mop returns the square of the radius.}*/   
  { return ptr()->sq_rad; }
  
  CGAL::Orientation orientation() const
/*{\Mop returns the orientation of the circle.}*/    
  { return ptr()->orientation; }
  
  bool is_degenerate() const
/*{\Mop returns |true| if the radius is 0, |false| otherwise.}*/    
  {
    return (squared_radius() == 0);
  }  
  
  // comparison operators
  bool operator==(const cgal_rat_circle& c2) const
  {
    if (this->center() != c2.center()) return false;
    if (this->sqr_radius() != c2.sqr_radius()) return false;  
    if (this->orientation() != c2.orientation()) return false;  
    return true;
  }
  
  bool operator!=(const cgal_rat_circle& c2) const
  { return (! (*this == c2)); }  

  friend unsigned long ID_Number(const cgal_rat_circle&);
};

unsigned long ID_Number(const cgal_rat_circle& c) { return c.ptr()->id; }

#if defined(LEDA_STD_IO_HEADERS)
std::ostream& operator<<(std::ostream& O, const cgal_rat_circle& C)
{
  O << "(" << C.center() << "," << C.squared_radius() << ")";
  return O;
}

std::istream& operator>>(std::istream& I, cgal_rat_circle& C)
{
  leda_rat_point p1,p2;
  I >> p1;
  I >> p2;
  C = cgal_rat_circle(p1,p2);
  return I;
}
#else
ostream& operator<<(ostream& O, const cgal_rat_circle& C)
{
  O << "(" << C.center() << "," << C.squared_radius() << ")";
  return O;
}

istream& operator>>(istream& I, cgal_rat_circle& C)
{
  leda_rat_point p1,p2;
  I >> p1;
  I >> p2;
  C = cgal_rat_circle(p1,p2);
  return I;
}
#endif

#undef __HANDLE_BASE_TYPE

LEDA_END_NAMESPACE

#endif
