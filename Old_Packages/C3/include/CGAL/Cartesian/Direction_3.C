// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_DIRECTION_3_C
#define CGAL_CARTESIAN_DIRECTION_3_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
DirectionC3<R CGAL_CTAG>::
DirectionC3()
{
   new ( static_cast< void*>(ptr)) Threetuple<FT>();
}

template < class R >
DirectionC3<R CGAL_CTAG>::
DirectionC3(const DirectionC3<R CGAL_CTAG> &d)
  : Handle_for<Threetuple<typename R::FT> >(d)
{}

template < class R >
DirectionC3<R CGAL_CTAG>::
DirectionC3(const typename DirectionC3<R CGAL_CTAG>::Vector_3 &v)
  : Handle_for<Threetuple<typename R::FT> >(v)
{}

template < class R >
DirectionC3<R CGAL_CTAG>::
DirectionC3(const typename DirectionC3<R CGAL_CTAG>::FT &x,
            const typename DirectionC3<R CGAL_CTAG>::FT &y,
            const typename DirectionC3<R CGAL_CTAG>::FT &z)
{
  new ( static_cast< void*>(ptr)) Threetuple<FT>(x, y, z);
}

template < class R >
DirectionC3<R CGAL_CTAG>::~DirectionC3()
{}

template < class R >
bool
DirectionC3<R CGAL_CTAG>::operator==(const DirectionC3<R CGAL_CTAG> &d) const
{
  if (ptr == d.ptr) return true;
  return equal_directionC3(dx(), dy(), dz(), d.dx(), d.dy(), d.dz());
}

template < class R >
inline
bool
DirectionC3<R CGAL_CTAG>::operator!=(const DirectionC3<R CGAL_CTAG> &d) const
{
  return !(*this == d);
}


template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::Vector_3
DirectionC3<R CGAL_CTAG>::to_vector() const
{
  return Vector_3(*this);
}

template < class R >
inline
DirectionC3<R CGAL_CTAG>
DirectionC3<R CGAL_CTAG>::
transform
  (const typename DirectionC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return t.transform(*this);
}

template < class R >
inline
DirectionC3<R CGAL_CTAG> 
DirectionC3<R CGAL_CTAG>::operator-() const
{
  return DirectionC3<R>(-dx(), -dy(), -dz());
}

template < class R >
typename DirectionC3<R CGAL_CTAG>::FT 
DirectionC3<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( i >= 0 && i <= 2 );
  // return (i==0) ? dx() :
//          (i==1) ? dy() : dz() ;
  if (i==0) return dx();
  else if (i==1) return dy();
  return dz() ;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::dx() const
{
  return ptr->e0;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::dy() const
{
  return ptr->e1;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::dz() const
{
  return ptr->e2;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::hdx() const
{
  return ptr->e0;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::hdy() const
{
  return ptr->e1;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::hdz() const
{
  return ptr->e2;
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::FT
DirectionC3<R CGAL_CTAG>::hw() const
{
  return FT(1);
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const DirectionC3<R CGAL_CTAG> &d)
{
  typename DirectionC3<R CGAL_CTAG>::Vector_3 v = d.to_vector();
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "DirectionC3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
      return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3
template < class R >
std::istream &
operator>>(std::istream &is, DirectionC3<R CGAL_CTAG> &p)
{
  typename R::FT x, y, z;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
      is >> x >> y >> z;
      break;
    case IO::BINARY :
      read(is, x);
      read(is, y);
      read(is, z);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  p = DirectionC3<R CGAL_CTAG>(x, y, z);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_DIRECTION_3_C
