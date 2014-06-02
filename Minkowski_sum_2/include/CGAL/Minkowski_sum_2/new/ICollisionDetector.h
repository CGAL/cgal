#ifndef ICOLLISIONDETECTOR_HEADER
#define ICOLLISIONDETECTOR_HEADER

#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>

template <class Kernel_, class Container_> class ICollisionDetector {
public:
    typedef Kernel_                                        Kernel;
    typedef CGAL::Polygon_2<Kernel, Container_>            Polygon_2;

    virtual bool checkCollision(const Polygon_2 &p, const Polygon_2 &q) = 0;
};

#endif