namespace A {

template <class T>
const T&
mix(const T& a1, const T& a2)
{ return a1 < a2 ? a1 : a2; } 

template <class T>
const T&
use(const T& a1, const T& a2)
{ return mix( a1, a2); }

} // namespace A

namespace B {

template <class T>
const T&
mix(const T& a1, const T& a2)
{ return a2 < a1 ? a1 : a2; } 

double
use_use( const double& t1, const double& t2)
{ return A::use(t1,t2); }

} // namespace B

int
main()
{
  B::use_use( 0.0, 1.0);
  return 0;
}
