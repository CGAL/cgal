#ifndef CGAL_BOUNDING_BOX_3_H
#define CGAL_BOUNDING_BOX_3_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>

CGAL_BEGIN_NAMESPACE

template <typename Tag, typename Traits> class Bounding_box_3;

template <typename Traits> 
class Bounding_box_rep_3 {

  typedef typename Traits::Point_3   Point_3;
  
 public:
  Bounding_box_rep_3() {
    min = Point_3(0,0,0);
    max = Point_3(0,0,0);
  }
  Bounding_box_rep_3(Point_3 init_min, Point_3 init_max) {
    min = init_min;
    max = init_max;
  }

  Point_3 min;
  Point_3 max;
};

template <typename Traits>
class Bounding_box_3<Cartesian_tag,Traits> : 
public Handle_for< Bounding_box_rep_3<Traits> >
{
  typedef Handle_for< Bounding_box_rep_3<Traits> > BBox_handle_3;
  typedef typename BBox_handle_3::element_type     BBox_ref_3;

  typedef typename Traits::Point_3             Point_3;
  typedef typename Traits::Vector_3            Vector_3;

public:
        Bounding_box_3()
	  : BBox_handle_3(BBox_ref_3()) {}

        Bounding_box_3(Point_3 init_min, Point_3 init_max) 
	  : BBox_handle_3(BBox_ref_3(init_min, init_max)) {}
       

	const Point_3& get_min() const { return this->Ptr()->min; }
	const Point_3& get_max() const { return this->Ptr()->max; }

	Bounding_box_3  operator+(const Bounding_box_3& b) const {
	  Point_3 res_min(b.get_min()), res_max(b.get_max());
	  
	  if(res_min.x() > get_min().x())
	    res_min = Point_3(res_min.x(),get_min().y(),get_min().z());
	  if(res_min.y() > get_min().y()) 
	    res_min = Point_3(get_min().x(),res_min.y(),get_min().z());
	  if(res_min.z() > get_min().z())
	    res_min = Point_3(get_min().x(),get_min().y(),res_min.z());
	  
	  if(res_max.x() < get_max().x())
	    res_max = Point_3(res_max.x(),get_max().y(),get_max().z());
	  if(res_max.y() < get_max().y())
	    res_max = Point_3(get_max().x(),res_max.y(),get_max().z());
	  if(res_max.z() < get_max().z())
	    res_max = Point_3(get_max().x(),get_max().y(),res_max.z());
	  
	  return Bounding_box_3(normalized(res_min),normalized(res_max));
	}
};

template <typename Traits>
class Bounding_box_3<Homogeneous_tag,Traits> : 
public Handle_for< Bounding_box_rep_3<Traits> >
{
  typedef Handle_for< Bounding_box_rep_3<Traits> > BBox_handle_3;
  typedef typename BBox_handle_3::element_type     BBox_ref_3;

  typedef typename Traits::Point_3             Point_3;
  typedef typename Traits::Vector_3            Vector_3;

public:
        Bounding_box_3()
	  : BBox_handle_3(BBox_ref_3()) {}

        Bounding_box_3(Point_3 init_min, Point_3 init_max) 
	  : BBox_handle_3(BBox_ref_3(init_min, init_max)) {}
       

	const Point_3& get_min() const { return this->Ptr()->min; }
	const Point_3& get_max() const { return this->Ptr()->max; }

	Bounding_box_3  operator+(const Bounding_box_3& b) const {

	  Point_3 res_min(b.get_min()), res_max(b.get_max());
	  
	  if(res_min.x() > get_min().x()) {
	    Point_3 src(res_min.hx(),0,0,res_min.hw());
	    Point_3 tgt(get_min().hx(),0,0,get_min().hw());
	    Vector_3 delta(tgt - src);
	    res_min = res_min + delta;
	  }
	  if(res_min.y() > get_min().y()) {
	    Point_3 src(0,res_min.hy(),0,res_min.hw());
	    Point_3 tgt(0,get_min().hy(),0,get_min().hw());
	    Vector_3 delta(tgt - src);
	    res_min = res_min + delta;
	  }
	  if(res_min.z() > get_min().z()) {
	    Point_3 src(0,0,res_min.hz(),res_min.hw());
	    Point_3 tgt(0,0,get_min().hz(),get_min().hw());
	    Vector_3 delta(tgt - src);
	    res_min = res_min + delta;
	  }
	  
	  if(res_max.x() < get_max().x()) {
	    Point_3 src(res_max.hx(),0,0,res_max.hw());
	    Point_3 tgt(get_max().hx(),0,0,get_max().hw());
	    Vector_3 delta(tgt - src);
	    res_max = res_max + delta;
	  }
	  if(res_max.y() < get_max().y()) {
	    Point_3 src(0,res_max.hy(),0,res_max.hw());
	    Point_3 tgt(0,get_max().hy(),0,get_max().hw());
	    Vector_3 delta(tgt - src);
	    res_max = res_max + delta;
	  }
	  if(res_max.z() < get_max().z()) {
	    Point_3 src(0,0,res_max.hz(),res_max.hw());
	    Point_3 tgt(0,0,get_max().hz(),get_max().hw());
	    Vector_3 delta(tgt - src);
	    res_max = res_max + delta;
	  }
	  
	  return Bounding_box_3(normalized(res_min),normalized(res_max));
	}
};

/*
template <typename NT>
inline
Bounding_box_3<Cartesian_tag,NT>
Bounding_box_3<Cartesian_tag,NT>::
operator+(const Bounding_box_3<Cartesian_tag,NT>& b) const {

  Point_3 res_min(b.get_min()), res_max(b.get_max());

  if(res_min.x() > get_min().x())
    res_min = Point_3(res_min.x(),get_min().y(),get_min().z());
  if(res_min.y() > get_min().y()) 
    res_min = Point_3(get_min().x(),res_min.y(),get_min().z());
  if(res_min.z() > get_min().z())
    res_min = Point_3(get_min().x(),get_min().y(),res_min.z());

  if(res_max.x() < get_max().x())
    res_max = Point_3(res_max.x(),get_max().y(),get_max().z());
  if(res_max.y() < get_max().y())
    res_max = Point_3(get_max().x(),res_max.y(),get_max().z());
  if(res_max.z() < get_max().z())
    res_max = Point_3(get_max().x(),get_max().y(),res_max.z());

  return Bounding_box_3(normalized(res_min),normalized(res_max));
}

template <typename NT>
inline
Bounding_box_3<Homogeneous_tag,NT>
Bounding_box_3<Tag_true,NT>::operator+(const Bounding_box_3<Tag_true,NT>& b) const {

  Point_3 res_min(b.get_min()), res_max(b.get_max());

  if(res_min.x() > get_min().x()) {
    Point_3 src(res_min.hx(),0,0,res_min.hw());
    Point_3 tgt(get_min().hx(),0,0,get_min().hw());
    Vector_3 delta(tgt - src);
    res_min = res_min + delta;
  }
  if(res_min.y() > get_min().y()) {
    Point_3 src(0,res_min.hy(),0,res_min.hw());
    Point_3 tgt(0,get_min().hy(),0,get_min().hw());
    Vector_3 delta(tgt - src);
    res_min = res_min + delta;
  }
  if(res_min.z() > get_min().z()) {
    Point_3 src(0,0,res_min.hz(),res_min.hw());
    Point_3 tgt(0,0,get_min().hz(),get_min().hw());
    Vector_3 delta(tgt - src);
    res_min = res_min + delta;
  }

  if(res_max.x() < get_max().x()) {
    Point_3 src(res_max.hx(),0,0,res_max.hw());
    Point_3 tgt(get_max().hx(),0,0,get_max().hw());
    Vector_3 delta(tgt - src);
    res_max = res_max + delta;
  }
  if(res_max.y() < get_max().y()) {
    Point_3 src(0,res_max.hy(),0,res_max.hw());
    Point_3 tgt(0,get_max().hy(),0,get_max().hw());
    Vector_3 delta(tgt - src);
    res_max = res_max + delta;
  }
  if(res_max.z() < get_max().z()) {
    Point_3 src(0,0,res_max.hz(),res_max.hw());
    Point_3 tgt(0,0,get_max().hz(),get_max().hw());
    Vector_3 delta(tgt - src);
    res_max = res_max + delta;
  }

  return Bounding_box_3(normalized(res_min),normalized(res_max));
}
*/

template <typename Tag,typename NT>
inline
bool
do_overlap(const Bounding_box_3<Tag,NT>& bb1, const Bounding_box_3<Tag,NT>& bb2)
{
    if (bb1.get_max().x() < bb2.get_min().x() || 
	bb2.get_max().x() < bb1.get_min().x())
        return false;
    if (bb1.get_max().y() < bb2.get_min().y() || 
	bb2.get_max().y() < bb1.get_min().y())
        return false;
    if (bb1.get_max().z() < bb2.get_min().z() || 
	bb2.get_max().z() < bb1.get_min().z())
        return false;
    return true;
}

#define CGAL_NO_ISTREAM_EXTRACT_BOUNDING_BOX_3

#ifndef CGAL_NO_OSTREAM_INSERT_BOUNDING_BOX_3
template <typename Tag,typename NT>
inline
std::ostream&
operator<<(std::ostream &os, const Bounding_box_3<Tag,NT>& b)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
      return os << b.get_min().hx() << ' ' << b.get_min().hy() << ' ' 
		<< b.get_min().hz() << ' ' << b.get_min().hw() << std::endl
		<< b.get_max().hx() << ' ' << b.get_max().hy() << ' ' 
		<< b.get_max().hz() << ' ' << b.get_max().hz();
    case IO::BINARY :
        write(os, b.get_min().hx());
        write(os, b.get_min().hy());
        write(os, b.get_min().hz());
        write(os, b.get_min().hw());
        write(os, b.get_max().hx());
        write(os, b.get_max().hy());
        write(os, b.get_max().hz());
        write(os, b.get_max().hw());
        return os;
    default:
      os << "Bounding_box_3<Tag,NT>((" << b.get_min().hx()
	 << ", "       << b.get_min().hy()
	 << ", "       << b.get_min().hz()
	 << ", "       << b.get_min().hw() << "), (";
      os <<               b.get_max().hx()
	 << ", "       << b.get_max().hy()
	 << ", "       << b.get_max().hz()
	 << ", "       << b.get_max().hw() << "))";
        return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_BOUNDING_BOX_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_BOUNDING_BOX_3
template <typename Tag, typename NT>
inline
std::istream&
operator>>(std::istream &is, Bounding_box_3<Tag,NT>& b)
{
  NT xmin, ymin, zmin, xmax, ymax, zmax;

  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> xmin >> ymin >> xmax >> ymax;
        break;
    case IO::BINARY :
        read(is, xmin);
        read(is, ymin);
        read(is, zmin);
        read(is, xmax);
        read(is, ymax);
        read(is, zmax);
        break;
  }
  b = Bounding_box_3<Tag,NT>(xmin, ymin, zmin, xmax, ymax, zmax);
  return is;
}

#endif // CGAL_NO_ISTREAM_EXTRACT_BOUNDING_BOX_3

CGAL_END_NAMESPACE

#endif // CGAL_BOUNDING_BOX_3_H
