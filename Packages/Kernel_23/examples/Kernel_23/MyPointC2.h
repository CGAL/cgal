#ifndef MY_POINTC2_H
#define MY_POINTC2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_2.h>




template < class R_ >
class MyPointC2

{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  FT x_, y_;
  int c_;

public:
  typedef CGAL::Cartesian_coordinate_iterator_2<R_> Cartesian_const_iterator;

  typedef R_   R;

  MyPointC2()
    : x_(FT(0)),y_(FT(0)), c_(0) 
  {}

  
  MyPointC2(const FT &x, const FT &y)
    : x_(x),y_(y), c_(0) {}

  
  const FT& x() const
  {
    return x_;
  }
  const FT& y() const
  {
    return y_;
  }

  FT& x()
  {
    return x_;
  }

  FT& y()
  {
    return y_;
  }

  int c() const
  {
    return c_;
  }

  int& c()
  {
    return c_;
  }
  
  
  const FT& cartesian(int i) const;
  FT homogeneous(int i) const;
  const FT& operator[](int i) const
  {
      return cartesian(i);
  }


  Cartesian_const_iterator cartesian_begin() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* >(this),0);
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* >(this), 2);
  }

  int dimension() const
  {
      return 2;
  }

  bool operator==(const MyPointC2 &p) const
  {
    return ( x_ == p.x_ )  && ( y_ == p.y_ );
  }

  bool operator!=(const MyPointC2 &p) const
  {
      return !(*this == p);
  }

  CGAL::Bbox_2 bbox() const;

  Point_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
CGAL_KERNEL_INLINE
const typename MyPointC2<R>::FT &
MyPointC2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
CGAL_KERNEL_INLINE
typename MyPointC2<R>::FT
MyPointC2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2)
    return cartesian(i);
  return FT(1);
}

template < class R >
CGAL_KERNEL_INLINE
CGAL::Bbox_2
MyPointC2<R>::bbox() const
{
  std::pair<double,double> xp = CGAL::to_interval(x());
  std::pair<double,double> yp = CGAL::to_interval(y());
  return CGAL::Bbox_2(xp.first, yp.first,  xp.second, yp.second);
}




template < class R >
std::ostream &
operator<<(std::ostream &os, const MyPointC2<R> &p)
{
    switch(os.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
        return os << p.x() << ' ' << p.y() << ' ' << p.c();
    case CGAL::IO::BINARY :
        CGAL::write(os, p.x());
        CGAL::write(os, p.y());
        CGAL::write(os, p.c());
        return os;
    default:
        return os << "MyPointC2(" << p.x() << ", " << p.y() << ", " << p.c() << ')';
    }
}



template < class R >
std::istream &
operator>>(std::istream &is, MyPointC2<R> &p)
{
    typename R::FT x, y;
    int c;
    switch(is.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
      is >> x >> y >> c;
        break;
    case CGAL::IO::BINARY :
        CGAL::read(is, x);
        CGAL::read(is, y);
	CGAL::read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is) {
      p = MyPointC2<R>(x, y);
      p.c() = c;
    }
    return is;
}


#endif // MY_POINTC2_H


