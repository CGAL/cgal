#ifndef CEP_LEDA_D3_RAT_ISO_CUBOID_H
#define CEP_LEDA_D3_RAT_ISO_CUBOID_H

#include <LEDA/d3_rat_point.h>

#if !defined(LEDA_BEGIN_NAMESPACE)
#define LEDA_BEGIN_NAMESPACE
#endif

#if !defined(LEDA_END_NAMESPACE)
#define LEDA_END_NAMESPACE
#endif


LEDA_BEGIN_NAMESPACE

class __exportC d3_rat_iso_cuboid_rep  : public handle_rep {

   friend class __exportC d3_rat_iso_cuboid;

   static leda_mutex mutex_id_counter;
   static unsigned long id_counter;

   leda_integer xmin;
   leda_integer xmax;
   leda_integer ymin;
   leda_integer ymax;
   leda_integer zmin;
   leda_integer zmax;
   leda_integer w; // common w value

   double xmin_d;
   double xmax_d;
   double ymin_d;
   double ymax_d;
   double zmin_d;
   double zmax_d;
   double w_d;

   unsigned long id;

public:

   d3_rat_iso_cuboid_rep() : xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0), w(1)
   { xmin_d = 0; xmax_d = 0;
     ymin_d = 0; ymax_d = 0;
     zmin_d = 0; zmax_d = 0;
     w_d = 1;
     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }

   // w coord 1:
   d3_rat_iso_cuboid_rep(leda_integer xmi, leda_integer xma,
                         leda_integer ymi, leda_integer yma,
			 leda_integer zmi, leda_integer zma) : xmin(xmi), xmax(xma), ymin(ymi), ymax(yma), zmin(zmi), zmax(zma), w(1)
   {
     xmin_d = xmin.to_double();
     xmax_d = xmax.to_double();
     ymin_d = ymin.to_double();
     ymax_d = ymax.to_double();
     zmin_d = zmin.to_double();
     zmax_d = zmax.to_double();
     w_d    = w.to_double();

     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }

   d3_rat_iso_cuboid_rep(leda_integer xmi, leda_integer xma,
                         leda_integer ymi, leda_integer yma,
			 leda_integer zmi, leda_integer zma,
			 leda_integer wv) : xmin(xmi), xmax(xma), ymin(ymi), ymax(yma), zmin(zmi), zmax(zma), w(wv)
   {
     xmin_d = xmin.to_double();
     xmax_d = xmax.to_double();
     ymin_d = ymin.to_double();
     ymax_d = ymax.to_double();
     zmin_d = zmin.to_double();
     zmax_d = zmax.to_double();
     w_d    = w.to_double();

     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }

   d3_rat_iso_cuboid_rep(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2)
   {
     leda_integer w1 = p1.W();
     leda_integer w2 = p2.W();
     w = w1 * w2;
     leda_integer ix1 = p1.X() * w2;
     leda_integer iy1 = p1.Y() * w2;
     leda_integer iz1 = p1.Z() * w2;
     leda_integer ix2 = p2.X() * w1;
     leda_integer iy2 = p2.Y() * w1;
     leda_integer iz2 = p2.Z() * w1; 
     
     if (ix1 < ix2) { xmin=ix1; xmax=ix2; } else { xmax=ix1; xmin=ix2; }   
     if (iy1 < iy2) { ymin=iy1; ymax=iy2; } else { ymax=iy1; ymin=iy2; }  
     if (iz1 < iz2) { zmin=iz1; zmax=iz2; } else { zmax=iz1; zmin=iz2; }      
     
     // compute double approximations ...
     xmin_d = xmin.to_double(); 
     xmax_d = xmax.to_double();
     ymin_d = ymin.to_double(); 
     ymax_d = ymax.to_double();
     zmin_d = zmin.to_double(); 
     zmax_d = zmax.to_double();
     w_d    = w.to_double();
          
     mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();               
   }

  ~d3_rat_iso_cuboid_rep() {}


   friend inline unsigned long ID_Number(const d3_rat_iso_cuboid& c);

};

#if !defined(CEP_LEDA_D3_RAT_ISO_CUBOID_NO_STATIC_INIT)
leda_mutex d3_rat_iso_cuboid_rep::mutex_id_counter;

unsigned long d3_rat_iso_cuboid_rep::id_counter = 0;
#endif

#if (__LEDA__ >= 440)
#define __HANDLE_BASE_TYPE  HANDLE_BASE(d3_rat_iso_cuboid_rep)
#else
#define __HANDLE_BASE_TYPE  handle_base
#endif 

/*{\Manpage {d3_rat_iso_cuboid} {} {Iso-oriented boxes in space} {IC} }*/

class __exportC d3_rat_iso_cuboid : public __HANDLE_BASE_TYPE {

/*{\Mdefinition
An instance |\Mvar| of the data type |\Mname| is an iso-oriented cuboid (box)
in 3d space with rational coordinates. }*/

//#if (__LEDA__ < 440)
  d3_rat_iso_cuboid_rep* ptr() const { return (d3_rat_iso_cuboid_rep*)PTR; } 
//#endif  
  
public:

/*{\Mcreation}*/

  d3_rat_iso_cuboid() { PTR = new d3_rat_iso_cuboid_rep; }
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|.}*/   
  
  d3_rat_iso_cuboid(leda_integer xmin, leda_integer xmax, leda_integer ymin, leda_integer ymax, leda_integer zmin, leda_integer zmax)
/*{\Mcreate introduces a variable |\Mvar| with minimal/maximal x/y/z - coordinates |xmin,xmax,ymin,ymax,zmin,zmax|
of type |\Mname|. }*/   
  { PTR = new d3_rat_iso_cuboid_rep(xmin,xmax,ymin,ymax,zmin,zmax); }
  
  d3_rat_iso_cuboid(leda_integer xmin, leda_integer xmax, leda_integer ymin, leda_integer ymax, leda_integer zmin, leda_integer zmax, leda_integer w)
/*{\Mcreate introduces a variable |\Mvar| with diagonal opposite vertices |(xmin,ymin,zmin,w)| and |(xmax,ymax,zmax,w)|
of type |\Mname|. }*/
  { PTR = new d3_rat_iso_cuboid_rep(xmin,xmax,ymin,ymax,zmin,zmax,w); }
  
  d3_rat_iso_cuboid(const leda_d3_rat_point& p,const leda_d3_rat_point& q)
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with diagonal opposite vertices |p| and |q|. }*/  
  { PTR = new d3_rat_iso_cuboid_rep(p,q); }
  
  d3_rat_iso_cuboid(const d3_rat_iso_cuboid& c) : __HANDLE_BASE_TYPE(c) {}
  ~d3_rat_iso_cuboid() {}

 d3_rat_iso_cuboid& operator=(const d3_rat_iso_cuboid& c)
 { __HANDLE_BASE_TYPE::operator=(c); return *this; }


  // the ordering is compatible to the CGAL iso cuboid ...

/*{\Moperations 2.5 4}*/  

  leda_d3_rat_point vertex(int i) const
/*{\Mop returns a vertex (corner) of |\Mvar|. \precond |0<=i<=7| .
The ordering of the vertices is compatible to the CGAL |Iso_cuboid_3|.
The following vertices are returned:
\begin{itemize}
\item i==0: lower left front 
\item i==1: lower right front
\item i==2: lower right back
\item i==3: lower left back
\item i==4: upper left front
\item i==5: upper left back
\item i==6: upper right front
\item i==7: upper right back
\end{itemize}
}*/  
  {
    // precondition: 0 <= i <= 7
    switch(i){
      case 0: { return leda_d3_rat_point(ptr()->xmin,ptr()->ymin,ptr()->zmin,ptr()->w); }
      case 1: { return leda_d3_rat_point(ptr()->xmax,ptr()->ymin,ptr()->zmin,ptr()->w); }
      case 2: { return leda_d3_rat_point(ptr()->xmax,ptr()->ymax,ptr()->zmin,ptr()->w); }
      case 3: { return leda_d3_rat_point(ptr()->xmin,ptr()->ymax,ptr()->zmin,ptr()->w); }
      case 4: { return leda_d3_rat_point(ptr()->xmin,ptr()->ymin,ptr()->zmax,ptr()->w); }
      case 5: { return leda_d3_rat_point(ptr()->xmin,ptr()->ymax,ptr()->zmax,ptr()->w); }
      case 6: { return leda_d3_rat_point(ptr()->xmax,ptr()->ymin,ptr()->zmax,ptr()->w); }
      case 7: { return leda_d3_rat_point(ptr()->xmax,ptr()->ymax,ptr()->zmax,ptr()->w); }
    }

    // precondition violated ...
    return leda_d3_rat_point(ptr()->xmax,ptr()->ymax,ptr()->zmax,ptr()->w);
  }
  
  leda_rational xmin() const { return leda_rational(ptr()->xmin,ptr()->w); }
/*{\Mop returns the smallest x-coordinate of |\Mvar|.}*/  
  
  leda_rational xmax() const { return leda_rational(ptr()->xmax,ptr()->w); }
/*{\Mop returns the largest x-coordinate of |\Mvar|.}*/   
  
  leda_rational ymin() const { return leda_rational(ptr()->ymin,ptr()->w); }
/*{\Mop returns the smallest y-coordinate of |\Mvar|.}*/   
  
  leda_rational ymax() const { return leda_rational(ptr()->ymax,ptr()->w); }
/*{\Mop returns the largest y-coordinate of |\Mvar|.}*/   
  
  leda_rational zmin() const { return leda_rational(ptr()->zmin,ptr()->w); }
/*{\Mop returns the smallest z-coordinate of |\Mvar|.}*/   
  
  leda_rational zmax() const { return leda_rational(ptr()->zmax,ptr()->w); } 
/*{\Mop returns the largest z-coordinate of |\Mvar|.}*/ 
  
  // homogeneous coords ...
  leda_integer  Xmin() const { return ptr()->xmin; }
  leda_integer  Xmax() const { return ptr()->xmax; }
  leda_integer  Ymin() const { return ptr()->ymin; }
  leda_integer  Ymax() const { return ptr()->ymax; }
  leda_integer  Zmin() const { return ptr()->zmin; }
  leda_integer  Zmax() const { return ptr()->zmax; } 
  leda_integer  W()    const { return ptr()->w; }
  
  // ... and double approximations
  double  XminD() const { return ptr()->xmin_d; }
  double  XmaxD() const { return ptr()->xmax_d; }
  double  YminD() const { return ptr()->ymin_d; }
  double  YmaxD() const { return ptr()->ymax_d; }
  double  ZminD() const { return ptr()->zmin_d; }
  double  ZmaxD() const { return ptr()->zmax_d; } 
  double  WD()    const { return ptr()->w_d; }     

  friend unsigned long ID_Number(const d3_rat_iso_cuboid&);
};

unsigned long ID_Number(const d3_rat_iso_cuboid& c) { return c.ptr()->id; }

#if defined(LEDA_STD_IO_HEADERS)
std::ostream& operator<<(std::ostream& O, const d3_rat_iso_cuboid& C)
{
  O << "(" << C.vertex(0) << "," << C.vertex(7) << ")";
  return O;
}

std::istream& operator>>(std::istream& I, d3_rat_iso_cuboid& C)
{
  leda_d3_rat_point p1,p2;
  I >> p1;
  I >> p2;
  C = d3_rat_iso_cuboid(p1,p2);
  return I;
}
#else
ostream& operator<<(ostream& O, const d3_rat_iso_cuboid& C)
{
  O << "(" << C.vertex(0) << "," << C.vertex(7) << ")";
  return O;
}

istream& operator>>(istream& I, d3_rat_iso_cuboid& C)
{
  leda_d3_rat_point p1,p2;
  I >> p1;
  I >> p2;
  C = d3_rat_iso_cuboid(p1,p2);
  return I;
}
#endif

#undef __HANDLE_BASE_TYPE

LEDA_END_NAMESPACE

#endif
