#ifndef CGAL_BOUNDING_BOX_3_H
#define CGAL_BOUNDING_BOX_3_H

#include <CGAL/basic.h>
#include <CGAL/Sixtuple.h>
#include <CGAL/Simple_Handle_for.h>

CGAL_BEGIN_NAMESPACE

template <typename NT>
class Bounding_box_3 : public Simple_Handle_for< Sixtuple<NT> >
{
  typedef Simple_Handle_for< Sixtuple<NT> > BBox_handle_3;
  typedef typename BBox_handle_3::element_type BBox_ref_3;

public:
        Bounding_box_3()
	  : BBox_handle_3(BBox_ref_3()) {}

        Bounding_box_3(const NT& x_min, const NT& y_min, const NT& z_min,
               const NT& x_max, const NT& y_max, const NT& z_max)
	  : BBox_handle_3(BBox_ref_3(x_min, y_min, z_min,
                                     x_max, y_max, z_max)) {}

  NT  xmin() const;
  NT  ymin() const;
  NT  zmin() const;
  NT  xmax() const;
  NT  ymax() const;
  NT  zmax() const;

  Bounding_box_3  operator+(const Bounding_box_3& b) const;
};

template <typename NT>
inline
NT
Bounding_box_3<NT>::xmin() const
{ return Ptr()->e0; }

template <typename NT>
inline
NT
Bounding_box_3<NT>::ymin() const
{ return Ptr()->e1; }

template <typename NT>
inline
NT
Bounding_box_3<NT>::zmin() const
{ return Ptr()->e2; }

template <typename NT>
inline
NT
Bounding_box_3<NT>::xmax() const
{ return Ptr()->e3; }

template <typename NT>
inline
NT
Bounding_box_3<NT>::ymax() const
{ return Ptr()->e4; }

template <typename NT>
inline
NT
Bounding_box_3<NT>::zmax() const
{ return Ptr()->e5; }

template <typename NT>
inline
Bounding_box_3<NT>
Bounding_box_3<NT>::operator+(const Bounding_box_3<NT>& b) const
{
  return Bounding_box_3<NT>(std::min(xmin(), b.xmin()),
                std::min(ymin(), b.ymin()),
                std::min(zmin(), b.zmin()),
                std::max(xmax(), b.xmax()),
                std::max(ymax(), b.ymax()),
                std::max(zmax(), b.zmax()));
}

template <typename NT>
inline
bool
do_overlap(const Bounding_box_3<NT>& bb1, const Bounding_box_3<NT>& bb2)
{
    // check for emptiness ??
    if (bb1.xmax() < bb2.xmin() || bb2.xmax() < bb1.xmin())
        return false;
    if (bb1.ymax() < bb2.ymin() || bb2.ymax() < bb1.ymin())
        return false;
    if (bb1.zmax() < bb2.zmin() || bb2.zmax() < bb1.zmin())
        return false;
    return true;
}

#ifndef CGAL_NO_OSTREAM_INSERT_BOUNDING_BOX_3
template <typename NT>
inline
std::ostream&
operator<<(std::ostream &os, const Bounding_box_3<NT>& b)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
      return os << b.xmin() << ' ' << b.ymin() << ' ' << b.zmin() << std::endl
		  << b.xmax() << ' ' << b.ymax() << ' ' << b.zmax();
    case IO::BINARY :
        write(os, b.xmin());
        write(os, b.ymin());
        write(os, b.zmin());
        write(os, b.xmax());
        write(os, b.ymax());
        write(os, b.zmax());
        return os;
    default:
        os << "Bounding_box_3<NT>((" << b.xmin()
           << ", "       << b.ymin()
           << ", "       << b.zmin() << "), (";
        os <<               b.xmax()
           << ", "       << b.ymax()
           << ", "       << b.zmax() << "))";
        return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_BOUNDING_BOX_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_BOUNDING_BOX_3
template <typename NT>
inline
std::istream&
operator>>(std::istream &is, Bounding_box_3<NT>& b)
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
  b = Bounding_box_3<NT>(xmin, ymin, zmin, xmax, ymax, zmax);
  return is;
}

#endif // CGAL_NO_ISTREAM_EXTRACT_BOUNDING_BOX_3

CGAL_END_NAMESPACE

#endif // CGAL_BOUNDING_BOX_3_H
