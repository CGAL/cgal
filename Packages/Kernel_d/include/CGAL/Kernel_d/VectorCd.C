#ifndef CGAL_VECTORCD_C
#define CGAL_VECTORCD_C
CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <class FT,class LA>
PointCd<FT,LA> VectorCd<FT,LA>::to_point() const
{ return PointCd<FT,LA>(Base(*this)); }

template <class FT,class LA>
PointCd<FT,LA> 
operator+ (const Origin& o, const VectorCd<FT,LA>& v)
{ return v.to_point(); }

template <class FT, class LA>
DirectionCd<FT,LA>  VectorCd<FT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorCd::direction: \
  zero vector cannot be a direction."); 
  return DirectionCd<FT,LA>(*this);
}

template <class FT, class LA>
VectorCd<FT,LA> VectorCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ typename LA::Matrix m_at = t.matrix();
  int d = t.dimension();
  for (int i = 0; i < d; i++) m_at(i,d) = 0;
  typename LA::Vector res(m_at*vector_rep());
  return VectorCd<FT,LA>(dimension(),res.begin(),res.end());
}

template <class FT, class LA>
VectorCd<FT,LA> operator*(int n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
VectorCd<FT,LA> operator*(const FT& n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
std::istream& operator>>(std::istream& I, VectorCd<FT,LA>& v)
{ v.copy_on_write(); v.ptr->read(I);
  return I; 
}

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorCd<FT,LA>& v)
{ v.ptr->print(O,"VectorCd"); return O; } 

template <class FT, class LA>
inline CGAL::io_Operator io_tag(const VectorCd<FT,LA>&) 
{ return CGAL::io_Operator(); }

#undef PointCd
CGAL_END_NAMESPACE
#endif // CGAL_VECTORCD_C

