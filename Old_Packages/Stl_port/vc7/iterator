#ifndef CGAL_ITER_VC7
#define CGAL_ITER_VC7
#include <../include/iterator>
#include <stl_iterator_base.h>
namespace std {
template<class C__> inline
typename iterator_traits<C__>::iterator_category
_Iter_cat(const C__&)
  {
    typedef typename iterator_traits<C__>::iterator_category c;
    return c();
  }

template <class _Iter> inline 
typename iterator_traits<_Iter>::difference_type*
  _Dist_type(const _Iter&)
  {
    typedef typename iterator_traits<_Iter>::difference_type _diff_type;
    return static_cast<_diff_type*>(0);
  }

template <class _Iter> inline 
typename iterator_traits<_Iter>::value_type*
  _Val_type(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::value_type _value_type;
  return static_cast<_value_type*>(0);
}


}
#endif // CGAL_ITER_VC7
