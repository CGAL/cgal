#ifndef CGAL_BBOX_3_ADAPTER_H
#define CGAL_BBOX_3_ADAPTER_H

#include <algorithm>
#include <functional>

/*
class Bbox_3
{
    double x_lo, y_lo, z_lo;
    double x_hi, y_hi, z_hi;
public:
    Bbox_3(double x_lo, double y_lo, double z_lo,
           double x_hi, double y_hi, double z_hi)
        : x_lo( x_lo ), y_lo( y_lo ), z_lo( z_lo ),
          x_hi( x_hi ), y_hi( y_hi ), z_hi( z_hi )
    {}

    double xlo() const { return x_lo; }
    double ylo() const { return y_lo; }
    double zlo() const { return z_lo; }
    double xhi() const { return x_hi; }
    double yhi() const { return y_hi; }
    double zhi() const { return z_hi; }

    Bbox_3 operator+(const Bbox_3& b) const {
        return Bbox_3(std::lo( xlo(), b.xlo() ),
                      std::lo( ylo(), b.ylo() ),
                      std::lo( zlo(), b.zlo() ),
                      std::hi( xhi(), b.xhi() ),
                      std::hi( yhi(), b.yhi() ),
                      std::hi( zhi(), b.zhi() ));
    }
};

class Bbox_3_Adapter : public Bbox_3 {
public:
    Bbox_3_Adapter( double x_lo, double y_lo, double z_lo,
                    double x_hi, double y_hi, double z_hi  )
        : Bbox_3( x_lo, y_lo, z_lo, x_hi, y_hi, z_hi )
    {}

    double get_lo( unsigned int n ) const {
        switch( n ) {
        case 0:  return xlo();
        case 1:  return ylo();
        case 2:  return zlo();
        default: return 0;
        }
    }

    double get_hi( unsigned int n ) const {
        switch( n ) {
        case 0:  return xhi();
        case 1:  return yhi();
        case 2:  return zhi();
        default: return 0;
        }
    }


};
*/

template< class T >
class Bbox_3
{
    T __lo[3], __hi[3];
    unsigned int __num;
public:
    typedef T NumberType;

    Bbox_3() {}
    Bbox_3( T x_lo, T y_lo, T z_lo,
            T x_hi, T y_hi, T z_hi )
    {
        __lo[0]= x_lo;   __lo[1]= y_lo;   __lo[2]= z_lo;
        __hi[0]= x_hi;   __hi[1]= y_hi;   __hi[2]= z_hi;
        __num = getCounter();
    }

    static unsigned int getCounter( bool reset = false ) {
       static unsigned int counter = 0;
       if( reset )
         counter = 0;
       else
         ++counter;
       return counter;
    }

    T lo( unsigned int dim ) const { return __lo[dim]; }
    T hi( unsigned int dim ) const { return __hi[dim]; }
    unsigned int num()       const { return __num;     }
    //unsigned int num()       const { return (unsigned int )this;  }
};

template< class _Box >
struct Bbox_3_Adapter {
    typedef _Box Box;
    typedef typename _Box::NumberType NumberType;

    static NumberType get_lo( const Box& b, unsigned int dim )
    { return b.lo( dim ); }

    static NumberType get_hi( const Box& b, unsigned int dim )
    { return b.hi( dim ); }

    static unsigned int get_num( const Box& b )
    { return b.num();     }
};

template< class _Box >
struct Bbox_3_Pointer_Adapter {
    typedef _Box* Box;
    typedef typename _Box::NumberType NumberType;

    static NumberType get_lo( const Box b, unsigned int dim )
    { return b->lo( dim ); }

    static NumberType get_hi( const Box b, unsigned int dim )
    { return b->hi( dim ); }

    static unsigned int get_num( const Box b )
    { return b->num();     }
};

#endif
