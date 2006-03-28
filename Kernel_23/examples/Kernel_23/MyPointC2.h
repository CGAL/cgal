#ifndef MY_POINTC2_H
#define MY_POINTC2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>


class MyPointC2 {

private:
  double vec[2];
  int col;

public:

  MyPointC2()
    : col(0)
  {
    *vec = 0;
    *(vec+1) = 0;
  }

  
  MyPointC2(const double x, const double y, int c)
    : col(c)
  {
    *vec = x;
    *(vec+1) = y;
  }

  const double& x() const  { return *vec; }

  const double& y() const { return *(vec+1); }

  double & x() { return *vec; }

  double& y() { return *(vec+1); }

  int color() const { return col; }

  int& color() { return col; }
  
  
  bool operator==(const MyPointC2 &p) const
  {
    return ( *vec == *(p.vec) )  && ( *(vec+1) == *(p.vec + 1) && ( col == p.col) );
  }

  bool operator!=(const MyPointC2 &p) const
  {
      return !(*this == p);
  }

};

template <class ConstructBbox_2>
class MyConstruct_bbox_2 : public ConstructBbox_2 {
public:
  CGAL::Bbox_2 operator()(const MyPointC2& p) const {
    return CGAL::Bbox_2(p.x(), p.y(), p.x(), p.y());
  }
};


class MyConstruct_coord_iterator {
public:
  const double* operator()(const MyPointC2& p)
  {
    return &p.x();
  }

  const double* operator()(const MyPointC2& p, int)
  {
    const double* pyptr = &p.y();
    pyptr++;
    return pyptr;
  }
};

  template <typename K, typename OldK>
  class MyConstruct_point_2
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2    Line_2;
  public:
    typedef Point_2          result_type;
    typedef CGAL::Arity_tag< 1 >   Arity;

    Point_2
    operator()(CGAL::Origin o) const
    { return MyPointC2(0, 0, 0); }

    Point_2
    operator()(const RT& x, const RT& y) const
    { 
      return MyPointC2(x, y, 0); 
    }

    Point_2
    operator()(const Line_2& l) const
    {
      typename OldK::Construct_point_2 base_operator;
      Point_2 p = base_operator(l);
      return p;
    }
    
    Point_2
    operator()(const Line_2& l, int i) const
    {
      typename OldK::Construct_point_2 base_operator;
      return base_operator(l, i);
    }

    // We need this one, as such a functor is in the Filtered_kernel
    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { 
      if(w != 1){
	return MyPointC2(x/w, y/w, 0); 
      } else {
	return MyPointC2(x,y, 0);
      }
    }
  };

std::ostream &
operator<<(std::ostream &os, const MyPointC2 &p)
{
    switch(os.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
        return os << p.x() << ' ' << p.y() << ' ' << p.color();
    case CGAL::IO::BINARY :
        CGAL::write(os, p.x());
        CGAL::write(os, p.y());
        CGAL::write(os, p.color());
        return os;
    default:
        return os << "MyPointC2(" << p.x() << ", " << p.y() << ", " << p.color() << ')';
    }
}



std::istream &
operator>>(std::istream &is, MyPointC2 &p)
{
    double x, y;
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
      p = MyPointC2(x, y, c);
    }
    return is;
}


#endif // MY_POINTC2_H


