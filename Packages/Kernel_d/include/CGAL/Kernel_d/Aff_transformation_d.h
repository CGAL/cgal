#ifndef CGAL_AFF_TRANSFORMATION_D_H
#define CGAL_AFF_TRANSFORMATION_D_H

CGAL_BEGIN_NAMESPACE

template <class pR>
class Aff_transformation_d : public pR::Aff_transformation_d_base
{ public:
  typedef typename pR::Aff_transformation_d_base Base;
  typedef Aff_transformation_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Aff_transformation_d(int d = 0, bool identity=false) 
    : Base(d,identity) {}
  Aff_transformation_d(const typename LA::Matrix& M) 
    : Base(M) {}
  template <typename Forward_iterator>
  Aff_transformation_d(
    Forward_iterator start, Forward_iterator end) : Base(start,end) {}
  Aff_transformation_d(const Vector_d<R>& v) : Base(v) {}
  Aff_transformation_d(int d, const RT& num, const RT& den) 
    : Base(d,num,den) {}
  Aff_transformation_d(int d, const RT& sin_num, const RT& cos_num, 
                       const RT& den, int e1 = 0, int e2 = 1)
    : Base(d,sin_num,cos_num,den,e1,e2) {}
  Aff_transformation_d(int d, const Direction_d<R>& dir, 
                       const RT& num, const RT& den, 
                       int e1 = 0, int e2 = 1)
    : Base(d,dir,num,den,e1,e2) {}
  Aff_transformation_d(const Base& a) : Base(a) {}
  Aff_transformation_d(const Self& a) : Base(a) {}

  Self operator*(const Self& a)
  { return Base::operator*(a); }
  Self inverse() const { return Base::inverse(); }

  bool operator==(const Self& a) const
  { return Base::operator==(a); }
  bool operator!=(const Self& a) const
  { return Base::operator!=(a); } 

};

CGAL_END_NAMESPACE
#endif //CGAL_AFF_TRANSFORMATION_D_H
