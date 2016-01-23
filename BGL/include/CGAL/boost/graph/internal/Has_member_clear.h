#ifndef CGAL_HAS_MEMBER_CLEAR_H
#define CGAL_HAS_MEMBER_CLEAR_H

namespace CGAL {
namespace internal {

template<class T>
class Has_member_clear
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static char f(check<void(C::*)(void), &C::clear>*);

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

}  // internal
}  // cgal

#endif /* CGAL_HAS_MEMBER_CLEAR_H */
