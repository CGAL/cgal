#ifndef CGAL_VECTORHD_C
#define CGAL_VECTORHD_C
CGAL_BEGIN_NAMESPACE

template <class RT,class LA>
PointHd<RT,LA> VectorHd<RT,LA>::to_point() const
{ return PointHd<RT,LA>(Base(*this)); }

template <class RT,class LA>
PointHd<RT,LA> 
operator+ (const Origin& o, const VectorHd<RT,LA>& v)
{ return v.to_point(); }

template <class RT, class LA>
DirectionHd<RT,LA>  VectorHd<RT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorHd::direction: \
  zero vector cannot be a direction."); 
  return DirectionHd<RT,LA>(*this); 
}

template <class RT, class LA>
VectorHd<RT,LA> VectorHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Matrix m_at = t.matrix(); 
  int d = t.dimension(); 
  for (int i = 0; i < d; i++) m_at(i,d) = 0;
  typename LA::Vector res(m_at*vector_rep());
  return VectorHd<RT,LA>(dimension(),res.begin(),res.end()); 
}

template <class RT, class LA>
VectorHd<RT,LA> VectorHd<RT,LA>::scale(const RT& m, const RT& n) const
{ int d = dimension(); 
  VectorHd<RT,LA> result(d); 
  result.entry(d) = entry(d) * n; 
  for (int i = 0; i < d; i++) 
    result.entry(i) = entry(i) * m; 
  return result; 
}

template <class RT, class LA>
void VectorHd<RT,LA>::self_scale(const RT& m, const RT& n)
{ int d = dimension(); 
  copy_on_write();
  entry(d) *= n; 
  for (int i = 0; i < d; i++) entry(i) *= m; 
}

template <class RT, class LA>
VectorHd<RT,LA> operator*(int n, const VectorHd<RT,LA>& v) 
{ return v.scale(n,1); }

template <class RT, class LA>
VectorHd<RT,LA> operator*(const RT& n, const VectorHd<RT,LA>& v) 
/*{\Mbinopfunc returns the vector with Cartesian coordinates $n v_i$.}*/
{ return v.scale(n,1); }

template <class RT, class LA>
VectorHd<RT,LA> operator*(const Quotient<RT>& r, const VectorHd<RT,LA>& v)
/*{\Mbinopfunc returns the vector with Cartesian coordinates 
$r v_i, 0 \leq i < d$.}*/
{ return v.scale(r.numerator(),r.denominator()); }

template <class RT, class LA>
std::istream& operator>>(std::istream& I, VectorHd<RT,LA>& v)
{ v.copy_on_write(); v.ptr->read(I); 
  CGAL_assertion_msg((v.homogeneous(v.dimension()) > 0),
  "operator>>: denominator of vector must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorHd<RT,LA>& v)
{ v.ptr->print(O,"VectorHd"); return O; } 

template <class RT, class LA>
inline CGAL::io_Operator io_tag(const VectorHd<RT,LA>&) 
{ return CGAL::io_Operator(); }


CGAL_END_NAMESPACE
#endif // CGAL_VECTORHD_C

