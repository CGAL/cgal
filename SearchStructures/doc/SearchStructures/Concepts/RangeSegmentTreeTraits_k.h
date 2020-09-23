
/*!
\ingroup PkgSearchStructuresConcepts
\cgalConcept

A tree traits class gives the range tree and segment tree class the necessary
type information of the keys and intervals. Further more, they define function objects that allow to access
the keys and intervals, and provide comparison functions that
are needed for window queries.


\cgalHasModel `CGAL::Range_segment_tree_set_traits_2`
\cgalHasModel `CGAL::Range_segment_tree_set_traits_3`
\cgalHasModel `CGAL::Range_tree_map_traits_2`
\cgalHasModel `CGAL::Range_tree_map_traits_3`
\cgalHasModel `CGAL::Segment_tree_map_traits_2`
\cgalHasModel `CGAL::Segment_tree_map_traits_3`

\cgalHeading{Example}

The following piece of code gives an example of how a traits class
might look like, if you have keys that are of the type `int`
in the first and that are of the type `double` in the second
dimension.

\code{.cpp}
class Int_double_tree_traits_2{
public:
  typedef std::pair<int, double> Key;
  typedef int Key_1;
  typedef double Key_2;
  typedef std::pair<Key,Key> Interval;

  class C_Key_1{
  public:
    Key_1 operator()(const Key& k)
    { return k.first;}
  };
  class C_Key_2{
  public:
    Key_2 operator()(const Key& k)
    { return k.second;}
  };
  class C_Low_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.first.first;}
  };
  class C_High_1{
  public:
    Key_1 operator()(const Interval& i)
    { return i.second.first;}
  };
  class C_Low_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.first.second;}
  };
  class C_High_2{
  public:
    Key_2 operator()(const Interval& i)
    { return i.second.second;}
  };
  class C_Compare_1{
  public:
    bool operator()(Key_1 k1, Key_1 k2)
    { return less<int>()(k1,k2);}
  };
  class C_Compare_2{
  public:
    bool operator()(Key_2 k1, Key_2 k2)
    { return less<double>()(k1,k2);}
  };
  typedef C_Compare_1 compare_1;
  typedef C_Compare_2 compare_2;
  typedef C_Low_1 low_1;
  typedef C_High_1 high_1;
  typedef C_Key_1 key_1;
  typedef C_Low_2 low_2;
  typedef C_High_2 high_2;
  typedef C_Key_2 key_2;
};
\endcode

*/

class RangeSegmentTreeTraits_k {
public:

/// \name Types
/// @{

/*!
The k-dimensional key type.
*/
typedef unspecified_type Key;

/*!
The k-dimensional interval type.
*/
typedef unspecified_type Interval;

/*!
The type in dimension \f$ i\f$, with \f$ 1\leq i
\leq k\f$.
*/
typedef unspecified_type Key_i;

/*!
function object providing an
`operator()` that takes an argument of type `Key` and returns
a component of type `Key_i`.
*/
typedef unspecified_type key_i;

/*!
function object providing an
`operator()` that takes an argument of type `Interval` and returns
a component of type `Key_i`.
*/
typedef unspecified_type low_i;

/*!
function object providing an
`operator()` that takes an argument of type `Interval` and returns
a component of type `Key_i`.
*/
typedef unspecified_type high_i;

/*!
function object providing an
`operator()` that takes two arguments argument \f$ a\f$, \f$ b\f$ of type `Key_i` and returns
`true` if \f$ a<b\f$, `false` otherwise.
*/
typedef unspecified_type compare_i;

/// @}

}; /* end RangeSegmentTreeTraits_k */

