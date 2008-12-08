#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_integer.h>
#include <CGAL/CORE_BigRat.h>

template <typename Rational> struct Rational_traits; 
template <> struct Rational_traits<leda::rational>{
  typedef leda::integer Integer; 
  static inline Integer numerator(const leda_rational& r)
  { return r.numerator(); }
  // ...
};
template <> struct Rational_traits<CORE::BigRat>{
  // ...
};

struct Tag_false {};
struct Tag_true  {}; 

template <typename T> struct Traits{
  typedef Tag_true Has_property;
};

template<typename T> inline 
int some_function(const T& x, const Tag_false&){
  // default code for function 
  return int(0);
}
template<typename T> inline 
int some_function(const T& x, const Tag_true&){
  // fast code using property
  return int(1);
}

template<typename T> inline
int some_function(const T& x){
  typedef typename Traits<T>::Has_property Tag;
  return function(x, Tag());
}

struct My_order{
  bool state;
  My_order(bool state_ = false): state(state_){};
  bool operator()(int a, int b){
    return  state ? a < b : a > b ; 
  }
};

struct Object_base{
  virtual ~Object_base(){}
};

template <typename T>
struct Wrapper: public Object_base{
  Wrapper(const T& t): _value(t){}
  // Wrapper(){}
  T get_value(){return _value;}
  virtual ~Wrapper(){};
private:
  T _value;
};

struct Object{
  template <typename T> 
  Object(const T& t):_base(new Wrapper<T>(t)){}
  
  ~Object(){ delete _base; }

  template <typename T>
  bool assign(T& t) const{
    Wrapper<T>* w_ptr = dynamic_cast<Wrapper<T>*>(_base);
    if( w_ptr == NULL ) return false; 
    t = w_ptr->get_value();
    return true;
  }
private:
  Object_base* _base;
};

template <typename T>
Object make_object(const T& t){
  return Object(t);
}

template <typename T>
bool assign_object(T& t, const Object& obj){
  return obj.assign(t);
}

int main() { 
  double d = -1; 
  int j =-1; 
  Object obj = make_object(int(1)); 
  assert(!assign_object(d,obj) && d == -1);
  assert( assign_object(j,obj) && j ==  1);

  std::cout << "assign double "<< assign_object(d,obj)<<std::endl;
  std::cout << "assign int    "<< assign_object(j,obj)<<std::endl;
  std::cout << "new j :       "<< j << std::endl;

  return 0; 


  //std::cout << some_function(4) << std::endl;
  //leda::rational r(2);
  //std::cout << Rational_traits<leda::rational>::numerator(r) << std::endl;

  int array[4]={10,20,5,8}; // could be vector<int>,list<int> 
  std::sort(&array[0],&array[4]);
  std::cout << array[0] <<","<< array[1] <<","<< array[2] <<","<< array[3] <<std::endl;
  My_order less(true); // less(3,4) == true 
  std::sort(&array[0],&array[4], less);
  My_order greater(false); 
  std::cout << array[0] <<","<< array[1] <<","<< array[2] <<","<< array[3] <<std::endl;
  std::sort(&array[0],&array[4], My_order(false));
  std::cout << array[0] <<","<< array[1] <<","<< array[2] <<","<< array[3] <<std::endl;
  
  return 0;
}
