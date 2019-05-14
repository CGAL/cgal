#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/function_objects.h>
#include <boost/bind.hpp>

// functor int x int -> Quotient<int>, (a,b) -> a/b
// ------------------------------------------------
typedef CGAL::Creator_2<int, int, CGAL::Quotient<int> > 
Quotient_creator;

// functor Quotient<int> ->  Quotient<int>, a/b -> b/a
// ---------------------------------------------------
struct Quotient_inverter
{
  typedef CGAL::Quotient<int> result_type;
  CGAL::Quotient<int> operator() (const CGAL::Quotient<int>& q) const
  {
    return CGAL::Quotient<int> (q.denominator(), q.numerator());
  }
}; 

int main() 
{
  // create composed functor (a,b) -> b/a...
  // ---------------------------------------
  int three = 3;
  int two = 2;
  std::cout << boost::bind
    (Quotient_inverter(), boost::bind
     (Quotient_creator(), _1, _2))
  // ...and apply it to (3, 2)
  // -------------------------
  (three, two);

  return 0;
}
