// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>
#include <algorithm>
#include <functional>
#include <utility>

CGAL_BEGIN_NAMESPACE

template < class _LA >
class LA_matrixCd : public Handle
{
public:
  typedef _LA                              LA;
  typedef typename LA::FT                  FT;
  typedef typename LA::RT                  RT;
  typedef typename LA::Vector              Vector;
  typedef LA_matrixCd<LA>                  Self;
  typedef FT*                              iterator;
  typedef const FT*                        const_iterator;
  
  LA_matrixCd(int n = 0);
  LA_matrixCd(int n, int m = n);
  LA_matrixCd(const std::pair<int,int> &d);
  LA_matrixCd(const Self &v);
  template < class InputIterator >
  LA_matrixCd(int n, const InputIterator &first, const InputIterator &last)
    {
      CGAL_kernel_precondition( last-first == n);
      int m = first->dimension();
      new_rep(n,m);
      InputIterator i;
      ptrdiff_t j;
      for (i=first,j=0; i!=last; ++i, j+=m)
        std::copy_n(i->begin(),m,begin()+j);
    }
  template < class InputIterator >
  LA_matrixCd(int n, int m,
              const InputIterator &first, const InputIterator &last)
    {
      CGAL_kernel_precondition( last-first == n*m);
      new_rep(n,m);
      std::copy(first,last,begin());
    }
  ~LA_matrixCd();

  template < class InputIterator >
  LA_matrixCd(InputIterator first, InputIterator last)
  {
    PTR = new _d_tuple<FT>(last-first);
    std::copy(first,last,begin());
  }

  Self&          operator=(const Self &v);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;

  Self           operator+(const Self &w) const;
  Self           operator-(const Self &w) const;
  Self           operator+() const { return *this; }
  Self           operator-() const;
  Self           operator*(const FT &c) const;
  Self           operator/(const FT &c) const;

  Self           operator*(const Self &w) const;
  Vector         operator*(const Vector &w) const;

  std::pair<int,int> dimension(); // (row,column)

  int            row_dimension()    const { return _rdim; }
  const Vector  &row(int i);        
  const Vector  &operator[](int i); 

  int            column_dimension() const { return _cdim; }
  Vector         column(int j) const;

  const FT      &operator[](int i, int j) const;
  const FT      &operator()(int i, int j) const { return operator[](i,j); }

  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e + ptr()->d; }
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e + ptr()->d; }

private:
  int _rdim, _cdim;
  void new_rep(int n, int m);
  const _d_tuple<FT>* ptr() const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*       ptr()       { return (_d_tuple<FT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_matrix_d.C>
#endif 

#endif // CGAL_CARTESIAN_VECTOR_D_H
