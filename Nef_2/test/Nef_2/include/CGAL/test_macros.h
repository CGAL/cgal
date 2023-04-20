#ifndef CGAL_TEST_MACROS_H
#define CGAL_TEST_MACROS_H

#include <CGAL/IO/io.h>
#include <iostream>
#include <sstream>
#include <vector>

#define CGAL_TEST_START int cgal_test_res=0

#define CGAL_TEST(b) if (!(b)) { ++cgal_test_res; \
std::cerr<<"ERROR: ("<<__LINE__ <<") test "<<#b<<" failed."<<std::endl; }


#define CGAL_IO_TEST(datao,datai,iomode) {                                \
    std::stringstream S;                                                  \
    CGAL::IO::set_mode(S,iomode);                                             \
    S << datao;                                                           \
    if ( iomode != CGAL::IO::BINARY)                                      \
        S << '\n';                                                        \
    S << datao;                                                           \
    S >> datai;                                                           \
    if (datao != datai) {                                                 \
       ++cgal_test_res;                                                   \
       std::cerr << "ERROR in 1.IO " << #iomode << " of " << #datao << " "\
                 << #datai << " : " << S.str() << " failed." <<std::endl; \
    }                                                                     \
    S >> datai;                                                           \
    if (datao != datai) {                                                 \
       ++cgal_test_res;                                                   \
       std::cerr << "ERROR in 2.IO " << #iomode << " of " << #datao << " "\
                 << #datai << " : " << S.str() << " failed." <<std::endl; \
    }                                                                     \
}


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
