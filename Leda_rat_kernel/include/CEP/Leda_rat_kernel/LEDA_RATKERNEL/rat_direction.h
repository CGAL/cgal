#ifndef CEP_LEDA_RAT_DIRECTION_H
#define CEP_LEDA_RAT_DIRECTION_H

#include <LEDA/rat_vector.h>
#include <LEDA/rat_segment.h>

LEDA_BEGIN_NAMESPACE

class __exportC rat_direction_rep  : public handle_rep {

   friend class __exportC rat_direction;

   static leda_mutex mutex_id_counter;
   static unsigned long id_counter;
   
   leda_rat_vector V;
   unsigned long id;

public:

   rat_direction_rep() : V(2)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }

   rat_direction_rep(const leda_rat_vector& vec) : V(vec)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }

   rat_direction_rep(int dim) : V(dim)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }
   
   rat_direction_rep(leda_integer x, leda_integer y, leda_integer w) : V(x,y,w)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }   

   rat_direction_rep(leda_rational x, leda_rational y) : V(x,y)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }
   
   rat_direction_rep(leda_integer x, leda_integer y, leda_integer z, leda_integer w) : V(x,y,z,w)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }   

   rat_direction_rep(leda_rational x, leda_rational y, leda_rational z) : V(x,y,z)
   { mutex_id_counter.lock();
     id = id_counter++;
     mutex_id_counter.unlock();
   }   

  ~rat_direction_rep() {}

   friend inline unsigned long ID_Number(const rat_direction& c);
};

#if !defined(CEP_LEDA_RAT_DIRECTION_NO_STATIC_INIT)

leda_mutex rat_direction_rep::mutex_id_counter;

unsigned long rat_direction_rep::id_counter = 0;

#endif

#if (__LEDA__ >= 440)
#define __HANDLE_BASE_TYPE  HANDLE_BASE(rat_direction_rep)
#else
#define __HANDLE_BASE_TYPE  handle_base
#endif 


/*{\Manpage {rat_direction} {} {Rational directions} {D} }*/

class __exportC rat_direction : public __HANDLE_BASE_TYPE  {

/*{\Mdefinition
An instance |\Mvar| of the data type |\Mname| is a direction (in the 2d plane or in 3d space).
It stores internally a vector with rational coordiantes. }*/

//#if (__LEDA__ < 440)
  rat_direction_rep* ptr() const { return (rat_direction_rep*)PTR; } 
//#endif

public:

/*{\Mcreation}*/

   rat_direction() 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|.}*/      
   { PTR = new rat_direction_rep; }
   
   rat_direction(const leda_rat_vector& vec)
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction of vector vec.}*/   
   { PTR = new rat_direction_rep(vec); }

   rat_direction(const leda_rat_segment& s)
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction of |s|.}*/   
   { PTR = new rat_direction_rep(s.to_vector()); }

   rat_direction(int dim) 
   { PTR = new rat_direction_rep(dim); }
   
   rat_direction(leda_integer x, leda_integer y, leda_integer w) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction passing
through the origin and point |p=(x/w,y/w)|.}*/   
   { PTR = new rat_direction_rep(x,y,w); }

   rat_direction(leda_rational x, leda_rational y) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction passing
through the origin and point |p=(x,y)|.}*/    
   { PTR = new rat_direction_rep(x,y); }
   
   rat_direction(leda_integer x, leda_integer y, leda_integer z, leda_integer w) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction passing
through the origin and point |p=(x/w,y/w,z/w)|.}*/    
   { PTR = new rat_direction_rep(x,y,z,w); }   

   rat_direction(leda_rational x, leda_rational y, leda_rational z) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| with the direction passing
through the origin and point |p=(x,y,z)|.}*/   
   { PTR = new rat_direction_rep(x,y,z); }
  
   rat_direction(const rat_direction& d) : __HANDLE_BASE_TYPE(d) {}
   ~rat_direction() {}

   rat_direction& operator=(const rat_direction& d) 
   { __HANDLE_BASE_TYPE::operator=(d); return *this; }     

/*{\Moperations 3 4}*/  

   int dim() const { return (ptr()->V).dim(); }
/*{\Mop returns the dimension of the underlying vector. }*/   
   
   const leda_rat_vector& get_vector() const { return ptr()->V; }
/*{\Mop returns the the underlying vector. }*/  

   friend unsigned long ID_Number(const rat_direction&); 
};

unsigned long ID_Number(const rat_direction& dir) { return dir.ptr()->id; }


#if defined(LEDA_STD_IO_HEADERS)
std::ostream& operator<<(std::ostream& O, const rat_direction& D)
{
  O << D.get_vector();
  return O;
}

std::istream& operator>>(std::istream& I, rat_direction& D)
{
  leda_rat_vector vec;
  I >> vec;
  D = rat_direction(vec);
  return I;
}
#else
ostream& operator<<(ostream& O, const rat_direction& D)
{
  O << D.get_vector();
  return O;
}

istream& operator>>(istream& I, rat_direction& D)
{
  leda_rat_vector vec;
  I >> vec;
  D = rat_direction(vec);
  return I;
}
#endif


#undef __HANDLE_BASE_TYPE

LEDA_END_NAMESPACE

#endif
