#ifndef CGAL_POINTHD_C
#define CGAL_POINTHD_C
CGAL_BEGIN_NAMESPACE

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Vector res = t(vector_rep());
  return PointHd<RT,LA>(dimension(),res.begin(),res.end()); }

template <class RT, class LA>
VectorHd<RT,LA> PointHd<RT,LA>::operator-(const Origin& o) const 
{ return VectorHd<RT,LA>(Base(*this)); }

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator+(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr->homogeneous_add(ptr, v.ptr); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator-(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr->homogeneous_sub(ptr, v.ptr); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator+= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr->homogeneous_add(old.ptr, v.ptr); 
  return *this; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator-= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr->homogeneous_sub(old.ptr, v.ptr); 
  return *this; 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, PointHd<RT,LA>& p)
{ p.copy_on_write(); p.ptr->read(I); 
  CGAL_assertion_msg((p.homogeneous(p.dimension()) > 0), 
  "operator>>: denominator of point must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const PointHd<RT,LA>& p)
{ p.ptr->print(O,"PointHd"); return O; } 

template <class RT, class LA>
inline CGAL::io_Operator io_tag(const PointHd<RT,LA>&) 
{ return CGAL::io_Operator(); }



CGAL_END_NAMESPACE
#endif // CGAL_POINTHD_C

