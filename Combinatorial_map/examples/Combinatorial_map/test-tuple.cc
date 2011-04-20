#include <CGAL/tuple.h>
#include <iostream>
typedef CGAL::cpp0x::tuple<void,int> testtuple;

int main()
{
  std::cout<<CGAL::cpp0x::tuple_size<testtuple>::value<<std::endl;
}
