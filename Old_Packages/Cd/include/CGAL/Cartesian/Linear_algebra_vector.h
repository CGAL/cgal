// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class _LA >
class LA_vectorCd : public Handle
{
public:
  typedef _LA                              LA;
  typedef typename LA::FT                  FT;
  typedef typename LA::RT                  RT;
  typedef LA_vectorCd<LA>                  Self;
  typedef FT*                              iterator;
  typedef const FT*                        const_iterator;
  
  LA_vectorCd(int dim = 0);
  LA_vectorCd(const FT &a);
  LA_vectorCd(const FT &a, const FT &b);
  LA_vectorCd(const FT &a, const FT &b, const FT &c);
  LA_vectorCd(const FT &a, const FT &b, const FT &c, const FT &d);
  LA_vectorCd(const Self &v);
  ~LA_vectorCd();

  template < class InputIterator >
  LA_vectorCd(InputIterator first, InputIterator last)
  {
    PTR = new _d_tuple<FT>(last-first);
    std::copy(first,last,begin());
  }

  Self&          operator=(const Self &v);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;

  FT             operator[](int i) const;

  Self           operator+(const Self &w) const;
  Self           operator-(const Self &w) const;
  Self           operator-() const;
  FT             operator*(const Self &w) const;
  Self           operator*(const FT &c) const;
  Self           operator/(const FT &c) const;

  int            dimension() const { return ptr()->d; }
  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e + dimension(); }
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e + dimension(); }

private:
  const _d_tuple<FT>* ptr()  const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*       ptr()        { return (_d_tuple<FT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_vector_d.C>
#endif 

#endif // CGAL_CARTESIAN_VECTOR_D_H
