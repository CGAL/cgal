#include <cassert>
#include <CGAL/Handle_with_policy.h>

struct Int_rep {
    int val;
    Int_rep( int i = 0) : val(i) {}
};

template < class Unify>
struct Int_t : public ::CGAL::Handle_with_policy< Int_rep, Unify > {
    typedef ::CGAL::Handle_with_policy< Int_rep, Unify > Base;
    Int_t( int i = 0) : Base( i) {}

    // This is needed to prevent VC7.1 and VC8 to call
    // the explicit templated constructor in Base instead of its copy-ctor.
    Int_t( Int_t const& rhs ) : Base( static_cast<Base const&>(rhs) ) {}

    int  value() const { return this->ptr()->val; }

    bool operator==( const Int_t<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            Base::unify(i);
        return equal;
    }
};


int main() {
  typedef Int_t< ::CGAL::Handle_policy_union> Int;

  Int j(5);
  Int k(5);

  // pump up the union_size counter for j and k

  Int j1(j);
  Int k1(k);

  assert( j == k);
  return 0;
}
