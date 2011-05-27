#include <CGAL/tuple.h>
#include <iostream>

template<typename T>
struct Convert_void
{ typedef T type; };

template<>
struct Convert_void<void>
{ typedef int type; };


#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template<typename ... Items>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<Items...> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<Items>::type... > type;
};
#else //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class T1>
struct Convert_tuple_with_void;

template <class T1>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type > type;
};
template <class T1, class T2>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type> type;
};
template <class T1, class T2, class T3>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type> type;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES




typedef CGAL::cpp0x::tuple<void,int> testtuple;
typedef Convert_tuple_with_void<testtuple>::type converttuple;

int main()
{
  // The following line does not compile on windows:
  // std::cout<<CGAL::cpp0x::tuple_size<testtuple>::value<<std::endl;

  // Idea: convert the tuple with possible void into a tuple without void.
  converttuple t;
  std::cout<<CGAL::cpp0x::tuple_size<converttuple>::value<<std::endl;
}
