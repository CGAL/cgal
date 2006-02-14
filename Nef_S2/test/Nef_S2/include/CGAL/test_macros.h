#ifndef CGAL_TEST_MACROS_H
#define CGAL_TEST_MACROS_H

#include <CGAL/basic.h>
#include <iostream>
#include <sstream>
#include <vector>

#define CGAL_TEST_START int cgal_test_res=0

#define CGAL_TEST(b) if (!(b)) { ++cgal_test_res; \
std::cerr<<"ERROR: ("<<__LINE__ <<") test "<<#b<<" failed."<<std::endl; } \
else

#define CGAL_IO_TEST(datao,datai) { \
std::ostrstream OS; OS<<datao<<'\n'<<'\0'; \
std::istrstream IS(OS.str()); OS.freeze(0); IS>>datai; \
if (datao != datai) { ++cgal_test_res; \
std::cerr<<"ERROR in IO of "<<#datao<<" "<<#datai<<" : "\
         <<OS.str()<<" failed."<<std::endl;\
OS.freeze(0); }}
   
#define CGAL_TEST_END return cgal_test_res

#undef CGAL_NEF_TRACE
#undef CGAL_NEF_TRACEN
#undef CGAL_NEF_TRACEV
#define CGAL_NEF_TRACE(t)  std::cerr << t
#define CGAL_NEF_TRACEN(t) std::cerr << t << std::endl
#define CGAL_NEF_TRACEV(t) std::cerr << #t << " = " << (t) << std::endl

template <class T>
std::vector<T> make_vector(const T& t1, const T& t2)
{ std::vector<T> V(2); V[0]=t1; V[1]=t2; return V; }

template <class T>
std::vector<T> make_vector(const T& t1, const T& t2, const T& t3)
{ std::vector<T> V(3); V[0]=t1; V[1]=t2; V[2]=t3; return V; }

template <class T>
std::vector<T> make_vector(const T& t1, const T& t2, 
	                   const T& t3, const T& t4)
{ std::vector<T> V(4); V[0]=t1; V[1]=t2; V[2]=t3; V[3]=t4; return V; }


#endif //CGAL_TEST_MACROS_H
