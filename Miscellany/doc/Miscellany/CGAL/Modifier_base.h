
namespace CGAL {

/*!
\ingroup MiscellanyRef

`Modifier_base` is an abstract base class providing the
interface for any modifier. A modifier is a function object derived
from `Modifier_base` that implements the pure virtual member
function `operator()`, which accepts a single reference parameter
`R&` on which the modifier is allowed to work. `R` is the
type of the internal representation that is to be modified.

\cgalHeading{Example}

The following fragment defines a class <TT>A</TT> with an internal
representation <TT>i</TT> of type <TT>int</TT>. It provides a member
function <TT>delegate()</TT>, which gives a modifier access to the
internal variable and checks validity thereafter. The
example modifier sets the internal variable to 42. The example
function applies the modifier to an instance of class <TT>A</TT>.

\code
class A {
  int i; // protected internal representation
public:
  void delegate( CGAL::Modifier_base<int>& modifier) {
    modifier(i);
    CGAL_postcondition( i > 0); // check validity
  }
};

struct Modifier : public CGAL::Modifier_base<int> {
  void operator()( int& rep) { rep = 42;}
};

void use_it() {
  A a;
  Modifier m;
  a.delegate(m); // a.i == 42 and A has checked that A::i > 0.
}
\endcode

*/
template< typename R >
class Modifier_base {
public:

/// \name Types
/// @{

/*!
the internal representation type.
*/
typedef R Representation;

/// @}

/// \name Operations
/// @{

/*!
\post `rep` is a valid representation.
*/
virtual void operator()( R& rep);

/// @}

}; /* end Modifier_base */
} /* end namespace CGAL */
