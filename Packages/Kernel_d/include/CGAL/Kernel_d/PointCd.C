#ifndef CGAL_POINTCDXXX_C
#define CGAL_POINTCDXXX_C
CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ typename LA::Vector res = t(vector_rep());
  return PointCd<FT,LA>(dimension(),res.begin(),res.end()); }

template <class FT, class LA>
VectorCd<FT,LA> PointCd<FT,LA>::operator-(const Origin& o) const 
{ return VectorCd<FT,LA>(Base(*this)); }

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::operator+(const VectorCd<FT,LA> &v) const
{ PointCd<FT,LA> res(dimension()); 
  res.ptr->cartesian_add(ptr, v.ptr);
  return res; 
}

template <class FT, class LA>
PointCd<FT,LA> PointCd<FT,LA>::operator-(const VectorCd<FT,LA> &v) const
{ PointCd<FT,LA> res(dimension()); 
  res.ptr->cartesian_sub(ptr, v.ptr); 
  return res; 
}

template <class FT, class LA>
PointCd<FT,LA>& PointCd<FT,LA>::operator+= (const VectorCd<FT,LA>& v)
{ int d = dimension(); 
  PointCd<FT,LA> old(*this); 
  *this = PointCd<FT,LA>(d); 
  ptr->cartesian_add(old.ptr, v.ptr); 
  return *this; 
}

template <class FT, class LA>
PointCd<FT,LA>& PointCd<FT,LA>::operator-= (const VectorCd<FT,LA>& v)
{ int d = dimension(); 
  PointCd<FT,LA> old(*this); 
  *this = PointCd<FT,LA>(d); 
  ptr->cartesian_sub(old.ptr, v.ptr); 
  return *this; 
}

template <class FT, class LA>
std::istream& operator>>(std::istream& I, PointCd<FT,LA>& p)
{ p.copy_on_write(); p.ptr->read(I);
  return I; 
}

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const PointCd<FT,LA>& p)
{ p.ptr->print(O,"PointCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const PointCd<FT,LA>&) 
{ return CGAL::io_Operator(); }

#undef PointCd
CGAL_END_NAMESPACE
#endif // CGAL_POINTCDXXX_C

