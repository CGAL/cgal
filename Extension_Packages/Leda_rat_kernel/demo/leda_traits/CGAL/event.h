// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/event.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================

/*
the CGAL events are based on the LEDA event implementation by Oliver Zlotowski 
*/


#ifndef CGAL_EVENT_H
#define CGAL_EVENT_H

#include <CGAL/basic.h>
#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

// ---------------------------------------------------------------------------
// event_item : BASE_RECEIVER*;
// ---------------------------------------------------------------------------

struct BASE_RECEIVER;
typedef BASE_RECEIVER* event_item;

struct BASE_EVENT
{ virtual ~BASE_EVENT() {}
  virtual void detach(event_item) = 0;
};

struct BASE_RECEIVER
{ BASE_EVENT* e;
  
  void* data;
  bool  enabled;

  BASE_RECEIVER(BASE_EVENT* x) : e(x), data(0), enabled(true) {}
  virtual ~BASE_RECEIVER() {}
  
  void enable()  { enabled = true; }
  void disable() { enabled = false; }
  bool is_enabled() const { return enabled; }
};


template <class PT>
struct base_receiver : public BASE_RECEIVER
{ base_receiver(BASE_EVENT* e) : BASE_RECEIVER(e) {}
  virtual void notify(const PT&) const = 0;
  virtual ~base_receiver() {}
};

// ---------------------------------------------------------------------------
// parameter structs
// ---------------------------------------------------------------------------

struct param0
{ enum { number = 0 }; 
  explicit param0() {}
};

template <class T> 
struct param1
{ enum { number = 1 }; 
  T val1;    
  explicit param1(T x1) : val1(x1) {}
};

template <class T1, class T2> 
struct param2 : public param1<T1>
{ enum { number = 2 };
  T2 val2;
  explicit param2(T1 x1, T2 x2) : param1<T1>(x1), val2(x2) {}
};

template <class T1, class T2, class T3> 
struct param3 : public param2<T1,T2>
{ enum { number = 3 };
  T3 val3;
  explicit param3(T1 x1, T2 x2, T3 x3) : param2<T1,T2>(x1,x2), val3(x3) {}
};

template <class T1, class T2, class T3, class T4> 
struct param4 : public param3<T1,T2,T3>
{ enum { number = 4 };
  T4 val4;
  explicit param4(T1 x1, T2 x2, T3 x3, T4 x4) 
    : param3<T1,T2,T3>(x1,x2,x3), val4(x4) {}
};

template <class T1, class T2, class T3, class T4, class T5> 
struct param5 : public param4<T1,T2,T3,T4>
{ enum { number = 5 };
  T5 val5;
  explicit param5(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5) 
    : param4<T1,T2,T3,T4>(x1,x2,x3,x4), val5(x5) {}
};

template <class T1, class T2, class T3, class T4, class T5, class T6> 
struct param6 : public param5<T1,T2,T3,T4,T5>
{ enum { number = 6 };
  T6 val6;
  explicit param6(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6) 
    : param5<T1,T2,T3,T4,T5>(x1,x2,x3,x4,x5), val6(x6) {}
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7> 
struct param7 : public param6<T1,T2,T3,T4,T5,T6>
{ enum { number = 7 };
  T7 val7;
  explicit param7(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7) 
    : param6<T1,T2,T3,T4,T5,T6>(x1,x2,x3,x4,x5,x6), val7(x7) {}
};

template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8> 
struct param8 : public param7<T1,T2,T3,T4,T5,T6,T7>
{ enum { number = 8 };
  T8 val8;
  explicit param8(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8) 
    : param7<T1,T2,T3,T4,T5,T6,T7>(x1,x2,x3,x4,x5,x6,x7), val8(x8) {}
};

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9> 
struct param9 : public param8<T1,T2,T3,T4,T5,T6,T7,T8>
{ enum { number = 9 };
  T9 val9;
  explicit param9(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9) 
    : param8<T1,T2,T3,T4,T5,T6,T7,T8>(x1,x2,x3,x4,x5,x6,x7,x8), val9(x9) {}
};

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10> 
struct param10 : public param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>
{ enum { number = 10 };
  T10 val10;
  explicit param10(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10) 
    : param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>(x1,x2,x3,x4,x5,x6,x7,x8,x9), val10(x10) {}
};

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class T11> 
struct param11 : public param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>
{ enum { number = 11 };
  T11 val11;
  explicit param11(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11) 
    : param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10), val11(x11) {}
};

template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7, class T8, class T9, class T10, class T11, class T12> 
struct param12 : public param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>
{ enum { number = 12 };
  T12 val12;
  explicit param12(T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11, T12 x12) 
    : param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11), val12(x12) {}
};

// ---------------------------------------------------------------------------
// caller classes
// ---------------------------------------------------------------------------

template <unsigned>
struct caller
{ template <class P, class R, class F>
  caller(const P& p, R r, F f);

  template <class P, class F>
  caller(const P& p, F f);
};  

template <>
struct caller<0>
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) { (r->*f)(); }

  template <class P, class F>
  caller(const P& p, F f) { f(); } 
};

template <>  
struct caller<1> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) { (r->*f)(p.val1); } 

  template <class P, class F>
  caller(const P& p, F f) { f(p.val1); } 
};

template <>
struct caller<2> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) { (r->*f)(p.val1,p.val2); } 
  
  template <class P, class F>
  caller(const P& p, F f) { f(p.val1,p.val2); } 
};

template <>
struct caller<3> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) { (r->*f)(p.val1,p.val2,p.val3); } 

  template <class P, class F>
  caller(const P& p, F f) { f(p.val1,p.val2,p.val3); } 
};

template <>
struct caller<4> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4); } 

  template <class P, class F>
  caller(const P& p, F f) { f(p.val1,p.val2,p.val3,p.val4); } 
};

template <>
struct caller<5> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5); } 

  template <class P, class F>
  caller(const P& p, F f) { f(p.val1,p.val2,p.val3,p.val4,p.val5); } 
};

template <>
struct caller<6> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6); } 
};

template <>
struct caller<7> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7); } 
};

template <>
struct caller<8> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8); } 
};

template <>
struct caller<9> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9); } 
};

template <>
struct caller<10> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10); } 
};

template <>
struct caller<11> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10,p.val11); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10,p.val11); } 
};

template <>
struct caller<12> 
{ template <class P, class R, class F>
  caller(const P& p, R r, F f) 
  { (r->*f)(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10,p.val11,p.val12); } 

  template <class P, class F>
  caller(const P& p, F f) 
  { f(p.val1,p.val2,p.val3,p.val4,p.val5,p.val6,p.val7,p.val8,p.val9,p.val10,p.val11,p.val12); } 
};

// ---------------------------------------------------------------------------
// receiver class
// ---------------------------------------------------------------------------

struct _NoClassReceiver {};

template <class PT, class F, class R = _NoClassReceiver>
struct receiver : public base_receiver<PT>
{
  enum { num = PT::number };
  
  typedef PT param_type;

  typedef base_receiver<PT> base;
  typedef receiver<PT,F,R>  self;
  
  R* r; // receiver object
  F  f; // handler function    

  void (self::*call)(const PT& p);

  void c_call(const PT& p) { caller<num>(p,r,f); }
  void f_call(const PT& p) { caller<num>(p,f); }

  receiver(BASE_EVENT& e, const self& x) : base(&e), r(x.r), f(x.f), call(x.call) {}

  receiver(R& _r, F _f) : base(0), r(&_r), f(_f), call(&self::c_call) {}
  receiver(F _f) : base(0), f(_f), call(&self::f_call) {}  

  void notify(const PT& p) const 
  { if (!is_enabled()) return;
    (const_cast<self*>(this)->*call)(p); 
  } 
};


// ---------------------------------------------------------------------------
// event
// ---------------------------------------------------------------------------

class event;

namespace detail {

class base_singleton;

struct data_struct
{ base_singleton* instance;

  typedef std::list< std::list<event_item> >::iterator ev_list_list_it;
  typedef std::list<event_item>::iterator              ev_list_it;

  ev_list_list_it       global;
  ev_list_it            local;  
  data_struct(base_singleton* x) : instance(x) {}
};
  
class base_singleton
{ protected:
  std::list< std::list<event_item> > L;
  
  typedef std::list< std::list<event_item> >::iterator ev_list_list_it;
  typedef std::list<event_item>::iterator              ev_list_it;
    
  typedef std::list< std::list<event_item> >::const_iterator ev_list_list_cit;
  typedef std::list<event_item>::const_iterator              ev_list_cit;
      
  public:
      
  template <class RT>
  event_item attach(event& e, event_item& x, const RT& r) 
  { data_struct& T = *(new data_struct(this));
    
    RT* R = new RT(e,r); R->data = &T;    
    if (x != 0)
    { data_struct& S = *(data_struct*) x->data;    
      T.global = S.global;        
    } 
    else
    { L.push_back(std::list<event_item>());   
      ev_list_list_it iter = --L.end();
      T.global = iter;     
      x = R;
    } 
    
    (*(T.global)).push_back(R);
    
    T.local = --(*(T.global)).end();    
    return R;
  }  
  
  template <class PT>     
  void occur(event_item x, const PT& p) const
  { data_struct& D = *(data_struct*) x->data;    
    event_item it;
    
    ev_list_list_it iter = D.global;
    const std::list<event_item>& act_list = *iter;
    
    // iteration on act_list
    ev_list_cit act_it = act_list.begin();
    
    for(;act_it != act_list.end();act_it++){
     it = *act_it;
     ((base_receiver<PT>*)it)->notify(p);  
    }
  }

  virtual event_item detach(event_item x) = 0;
  virtual void clear(event_item x) = 0;
  virtual ~base_singleton() {}
};


template <class T>
class singleton : public base_singleton
{ static void* instance;
  singleton() {}  
  
  typedef std::list< std::list<event_item> >::iterator ev_list_list_it;
  typedef std::list<event_item>::iterator              ev_list_it;    
  
  public:

  event_item detach(event_item x) 
  { data_struct& D = *(data_struct*) x->data;    
    ev_list_list_it item = D.global;    
    
    delete *(D.local);
    (*item).erase(D.local); 
    delete &D;    
    
    if (! (*item).empty()) return (*item).front();
    
    L.erase(item);    
      
    if (L.empty()) 
    { instance = 0;
      delete this;
    }      
    return 0;      
  }  

  void clear(event_item x)
  { data_struct& D = *(data_struct*) x->data;    
    ev_list_list_it item = D.global;    
    
    while (! (*item).empty())
    { event_item x = (*item).front();
      (*item).pop_front(); 
      delete (data_struct*) x->data;
      delete x;
    }
    
    L.erase(item);
    
    if (L.empty()) 
    { instance = 0;
      delete this;
    }      
  }

  static singleton<T>& get()
  { if (instance == 0) instance = new singleton<T>;
    return *(singleton<T>*) instance;
  }
};

template <class T> void* singleton<T>::instance = 0;

}

// ---------------------------------------------------------------------------
// event
// ---------------------------------------------------------------------------

/*{\Manpage {event} {} {Events}}*/

/*{\Mdefinition
Events are designed to support all kinds of Observer Pattern. An event |E|
is able to store arbitrary receivers. Each receiver contains a reference to 
a global function or to a member function of an arbitrary data type. 
Functions can be added/removed by calling the non-member functions |attach|/|detach|. 
If an event occurs, all attached functions with a suitable parameter list
are called. 
}*/

class event : public BASE_EVENT
{ std::map<detail::base_singleton*,event_item> M;
  bool enabled;
  
  typedef std::map<detail::base_singleton*,event_item>::const_iterator  const_iterator;
  
  event(const event& e);
  event& operator=(const event& e);
 
  public:

/*{\Mcreation E }*/
  event() : enabled(true) {}
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and initializes it
            to the empty receiver list. }*/
   
 ~event() 
  { event_item x;
  
    std::map<detail::base_singleton*,event_item>::iterator it = M.begin();
  
    for(;it != M.end();it++)
    { 
      x = (*it).second;
      if (!x) continue;
      detail::data_struct& d = *(detail::data_struct*) x->data;    
      detail::base_singleton& s = *d.instance;
      s.clear(x);      
    }
  }
  
  template <class RT>
  event_item attach(const RT& r)
  { typedef typename RT::param_type PT;
    detail::singleton<PT>& s = detail::singleton<PT>::get();    
    return s.attach(*this,M[&s],r);
  }

  template <class PT>
  void occur(const PT& x) const 
  { if (!enabled) return;
    detail::singleton<PT>& s = detail::singleton<PT>::get();   
    
    const_iterator cit = M.find(&s);
    
    if (cit != M.end() ) {
      event_item ev_it = (*cit).second;
      if (ev_it != NULL) s.occur(ev_it,x);
    }   
  }
  
  void detach(event_item x)
  { detail::data_struct& d = *(detail::data_struct*) x->data;    
    detail::base_singleton& s = *d.instance;
    M[&s] = s.detach(x);
  }
    
  void enable()  { enabled = true; }   
  void disable() { enabled = false; }    
  bool is_enabled() const { return enabled; } 
};

/*{\Moperations 1.3 2.7 }*/


/*{\Moptions nextwarning=no}*/
/*
void enable();
*/
/*{\Mop enables |E|, i. e., all attached functions will be notified
        when |E| occurs.}*/

/*{\Moptions nextwarning=no}*/
/*
void disable();
*/
/*{\Mop disables |E|, no function will be notified when |E| occurs.}*/

/*{\Moptions nextwarning=no}*/
/*
bool is_enabled();
*/
/*{\Mop returns |true| if |E| is enabled, otherwise |false|.}*/


/*{\Mtext
\bigskip
{\bf Non-Member Functions} 
}*/

/*{\Moptions nextwarning=no}*/
/*
event_item attach(event& E, observer& o, function f);
*/
/*{\Mfunc    appends a new member function |f| of object |o| to the  
             receiver list of |E| and returns an |event_item| for it. 
             If |E| occurs with suitable parameters, |f| is called 
             for object |o|. }*/ 

/*{\Moptions nextwarning=no}*/
/*
event_item attach(event& E, function f);
*/
/*{\Mfunc    appends a new function |f| to the list of receivers 
             of |E| and returns an |event_item| for it. }*/ 

/*{\Moptions nextwarning=no}*/
/*
template<class T>
void occur(event& E, T p);
*/
/*{\Mfunc  calls all attached functions with |p| as parameter
           of the associated event |E|. There are also |occur| functions 
           for more than one parameter.
}*/ 

/*{\Moptions nextwarning=no}*/
/*
void detach(event_item x)
*/
/*{\Mfunc   removes |x| from the list of receivers of the associated event.}*/

/*{\Moptions nextwarning=no}*/
/*
void enable(event_item x)
*/
/*{\Mfunc   enables |x| in the list of receivers of the associated event.}*/

/*{\Moptions nextwarning=no}*/
/*
void disable(event_item x)
*/
/*{\Mfunc   disables |x| in the list of receivers of the associated event.}*/

/*{\Moptions nextwarning=no}*/
/*
bool is_enabled(event_item x)
*/
/*{\Mfunc   returns |true| if |x| is an enabled receiver in the list of 
            the associated event, otherwise |false|.}*/


/*{\Moptions nextwarning=no}*/
/*
void detach(event_item x, int c) 
*/
/*{\Mfunc  removes |c| receivers handled by vector |x| from the list 
           of receivers of the associated event.}*/

/*{\Moptions nextwarning=no}*/
/*
void enable(event_item x, int c) 
*/
/*{\Mfunc  enables |c| receivers handled by vector |x| of the 
           associated event.}*/

/*{\Moptions nextwarning=no}*/
/*
void disable(event_item x, int c) 
*/
/*{\Mfunc  disables |c| receivers handled by vector |x| of the 
           associated event.}*/

/*{\Mexample  
  \begin{verbatim}
  #include <CGAL/event.h>
  #include <iostream>

  using namespace CGAL;
  using std::cout;

  void by_reference(int& val)
  { cout << "Call by reference:  " << val <<  " Change value! \n";
    val = 3;
  }

  void by_value(int val)
  { cout << "Call by value:      " << val << "\n"; }

  int main()
  {
    event e;  
    attach(e,by_value);
    attach(e,by_reference); 
    int i = 10;  
    occur(e,i);
    occur<int&>(e,i);  
    occur(e,i); 
    return 0;
  }
  \end{verbatim}
}*/

// ---------------------------------------------------------------------------
// non-member functions
// ---------------------------------------------------------------------------

inline void detach(event_item x)
{ BASE_EVENT* e = x->e;
  if (e) e->detach(x);
}

inline void enable(event_item x)     { x->enable(); }

inline void disable(event_item x)    { x->disable(); }

inline bool is_enabled(event_item x) { return x->is_enabled(); }

inline void detach(event_item* x, int c) 
{ while (c > 0) 
  { BASE_EVENT* e = (x[--c]->e);
    if (e) e->detach(x[c]);
  }
}

inline void enable(event_item* x, int c) 
{ while (c > 0) x[--c]->enable(); }

inline void disable(event_item* x, int c) 
{ while (c > 0) x[--c]->disable(); }

//-----------------------------------------------------------------------------
// handler without arguments
//-----------------------------------------------------------------------------

template <class R, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)()) 
{ typedef void (R::*F)(); 
  return e.attach(receiver<param0,F,R>(r,f)); 
}

template <class R, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)() const) 
{ typedef void (R::*F)() const; 
  return e.attach(receiver<param0,F,R>(r,f)); 
}

template <class ET>
inline event_item attach(ET& e, void (*f)()) 
{ typedef void (*F)();
  return e.attach(receiver<param0,F>(f));
}

inline void occur(event& e) { e.occur(param0()); }

//-----------------------------------------------------------------------------
// handler with one argument
//-----------------------------------------------------------------------------

template <class R, class T, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T)) 
{ typedef void (R::*F)(T);
  return e.attach(receiver<param1<T>,F,R>(r,f)); 
}

template <class R, class T, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T) const) 
{ typedef void (R::*F)(T) const;
  return e.attach(receiver<param1<T>,F,R>(r,f)); 
}

template <class T, class ET>
inline event_item attach(ET& e, void (*f)(T)) 
{ typedef void (*F)(T);
  return e.attach(receiver<param1<T>,F>(f)); 
}

template <class T>
inline void occur(event& e, T x) { e.occur(param1<T>(x)); }

//-----------------------------------------------------------------------------
// handler with two arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2)) 
{ typedef void (R::*F)(T1,T2);
  return e.attach(receiver<param2<T1,T2>,F,R>(r,f)); 
}

template <class R, class T1, class T2, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2) const) 
{ typedef void (R::*F)(T1,T2) const;
  return e.attach(receiver<param2<T1,T2>,F,R>(r,f)); 
}

template <class T1, class T2, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2)) 
{ typedef void (*F)(T1,T2);
  return e.attach(receiver<param2<T1,T2>,F>(f)); 
}

template <class T1, class T2>
inline void occur(event& e, T1 x1, T2 x2) { e.occur(param2<T1,T2>(x1,x2)); }

//-----------------------------------------------------------------------------
// handler with three arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3)) 
{ typedef void (R::*F)(T1,T2,T3);
  return e.attach(receiver<param3<T1,T2,T3>,F,R>(r,f)); 
}

template <class R, class T1, class T2, class T3, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3) const) 
{ typedef void (R::*F)(T1,T2,T3) const;
  return e.attach(receiver<param3<T1,T2,T3>,F,R>(r,f)); 
}

template <class T1, class T2, class T3, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3)) 
{ typedef void (*F)(T1,T2,T3);
  return e.attach(receiver<param3<T1,T2,T3>,F>(f)); 
}

template <class T1, class T2, class T3>
inline void occur(event& e, T1 x1, T2 x2, T3 x3) 
{ e.occur(param3<T1,T2,T3>(x1,x2,x3)); }

//-----------------------------------------------------------------------------
// handler with four arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, class T4, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4)) 
{ typedef void (R::*F)(T1,T2,T3,T4);
  return e.attach(receiver<param4<T1,T2,T3,T4>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, class T4, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4) const) 
{ typedef void (R::*F)(T1,T2,T3,T4) const;
  return e.attach(receiver<param4<T1,T2,T3,T4>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4)) 
{ typedef void (*F)(T1,T2,T3,T4);
  return e.attach(receiver<param4<T1,T2,T3,T4>,F>(f));
}

template <class T1, class T2, class T3, class T4>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4) 
{ e.occur(param4<T1,T2,T3,T4>(x1,x2,x3,x4)); }

//-----------------------------------------------------------------------------
// handler with five arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, class T4, class  T5, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5);
  return e.attach(receiver<param5<T1,T2,T3,T4,T5>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, class T4, class  T5, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5) const;
  return e.attach(receiver<param5<T1,T2,T3,T4,T5>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, class T5, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5)) 
{ typedef void (*F)(T1,T2,T3,T4,T5);
  return e.attach(receiver<param5<T1,T2,T3,T4,T5>,F>(f));
}

template <class T1, class T2, class T3, class T4, class T5>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5) 
{ e.occur(param5<T1,T2,T3,T4,T5>(x1,x2,x3,x4,x5)); }

//-----------------------------------------------------------------------------
// handler with six arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6);
  return e.attach(receiver<param6<T1,T2,T3,T4,T5,T6>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6) const;
  return e.attach(receiver<param6<T1,T2,T3,T4,T5,T6>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, class T5, class T6, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6);
  return e.attach(receiver<param6<T1,T2,T3,T4,T5,T6>,F>(f));
}

template <class T1, class T2, class T3, class T4, class T5, class T6>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6) 
{ e.occur(param6<T1,T2,T3,T4,T5,T6>(x1,x2,x3,x4,x5,x6)); }


//-----------------------------------------------------------------------------
// handler with seven arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7);
  return e.attach(receiver<param7<T1,T2,T3,T4,T5,T6,T7>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7) const;
  return e.attach(receiver<param7<T1,T2,T3,T4,T5,T6,T7>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7);
  return e.attach(receiver<param7<T1,T2,T3,T4,T5,T6,T7>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7) 
{ e.occur(param7<T1,T2,T3,T4,T5,T6,T7>(x1,x2,x3,x4,x5,x6,x7)); }


//-----------------------------------------------------------------------------
// handler with eight arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8);
  return e.attach(receiver<param8<T1,T2,T3,T4,T5,T6,T7,T8>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8) const;
  return e.attach(receiver<param8<T1,T2,T3,T4,T5,T6,T7,T8>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7,T8)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7,T8);
  return e.attach(receiver<param8<T1,T2,T3,T4,T5,T6,T7,T8>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8) 
{ e.occur(param8<T1,T2,T3,T4,T5,T6,T7,T8>(x1,x2,x3,x4,x5,x6,x7,x8)); }


//-----------------------------------------------------------------------------
// handler with nine arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9);
  return e.attach(receiver<param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9) const;
  return e.attach(receiver<param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9);
  return e.attach(receiver<param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9) 
{ e.occur(param9<T1,T2,T3,T4,T5,T6,T7,T8,T9>(x1,x2,x3,x4,x5,x6,x7,x8,x9)); }


//-----------------------------------------------------------------------------
// handler with ten arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10);
  return e.attach(receiver<param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10) const;
  return e.attach(receiver<param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10);
  return e.attach(receiver<param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10) 
{ e.occur(param10<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)); }


//-----------------------------------------------------------------------------
// handler with eleven arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11);
  return e.attach(receiver<param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11) const;
  return e.attach(receiver<param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class T11, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11);
  return e.attach(receiver<param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class T11>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11) 
{ e.occur(param11<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)); }


//-----------------------------------------------------------------------------
// handler with twelve arguments
//-----------------------------------------------------------------------------

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12)) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12);
  return e.attach(receiver<param12<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>,F,R>(r,f));
}

template <class R, class T1, class T2, class T3, 
          class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12, class ET>
inline event_item attach(ET& e, R& r, void (R::*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12) const) 
{ typedef void (R::*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12) const;
  return e.attach(receiver<param12<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>,F,R>(r,f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12, class ET>
inline event_item attach(ET& e, void (*f)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12)) 
{ typedef void (*F)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12);
  return e.attach(receiver<param12<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>,F>(f));
}

template <class T1, class T2, class T3, class T4, 
          class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12>
inline void occur(event& e, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11, T12 x12) 
{ e.occur(param12<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)); }


CGAL_END_NAMESPACE

#endif


