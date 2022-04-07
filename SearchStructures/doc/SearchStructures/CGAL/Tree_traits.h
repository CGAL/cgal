
namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

`tree_interval_traits` is a template class
that provides an interface to data items. It is similar to
`tree_point_traits`, except that it provides
access to two data slots of the same type of each container class
(`Data, Window`) instead of providing access to one
data slot of container class `Data` and two data slots
of class `Window`.
*/
template<typename Data, typename Window, typename Key,
         typename Data_left_func, typename Data_right_func,
         typename Window_left_func, typename Window_right_func,
         typename Compare>
class tree_interval_traits {
public:



/// \name Types
/// @{

/*!
the container `Data` -
the data type. It may consist of
several data slots. Two of these data slots have to be of
type `Key`.
*/
typedef unspecified_type Data;

/*!
the container
`Window` - the query window type. It may consist of
several data slots. Two of these data slots have to be of
type `Key`.
*/
typedef unspecified_type Window;

/*!
the type
`Key` of the data
slot this traits class provides access to.
*/
typedef unspecified_type Key;

/*!

`Data_left_func` is a
function object providing an
`operator()` that takes an argument of type `Data`
and returns
a (the left) component of type `Key`.
*/
typedef unspecified_type Data_left_func;

/*!

`Data_right_func` is a
function object providing an
`operator()` that takes an argument of type `Data`
and returns
a (the right) component of type `Key`.
*/
typedef unspecified_type Data_right_func;

/*!

`Window_left_func` is a function objects that
allow to access the
left data slot of container
`Window` which has type `Key`
*/
typedef unspecified_type Window_left_func;

/*!

`Window_right_func` is a function objects that
allow to access the
right data slot of container
`Window` which has type `Key`
*/
typedef unspecified_type Window_right_func;

/*!
defines a comparison relation which must
define a strict ordering of the objects of type
`Key`. If defined, `less<Key>`
is sufficient.
*/
typedef unspecified_type Compare;

/// @}

/// \name Operations
/// @{

/*!
The data slot of
the data item of `d` of type `Key` is
accessed by function object
`Data_left_func`.
*/
Key get_left(Data d);

/*!
The data slot of
the data item of `d` of type `Key` is
accessed by function object
`Data_right_func`.
*/
Key get_right(Data d);

/*!
The data slot of
the data item of `w` of type `Key` is
accessed by function object
`Window_left_func`.
*/
Key get_left_win(Window w);

/*!
The data slot of
the data item of `w` of type `Key` is
accessed by function object `Window_right_func`.
*/
Key get_right_win(Window w);

/*!
returns Compare(key1, key2).
*/
static bool comp(Key& key1, Key& key2);

/// @}

}; /* end tree_interval_traits */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSearchStructuresTraitsClasses

`tree_point_traits` is a template class
that provides an interface to data items.

*/
template<typename Data, typename Window, typename Key, typename Data_func,
         typename Window_left_func, typename Window_right_func,
         typename Compare>
class tree_point_traits {
public:

/// \name Types
/// @{

/*!
the container
`Data` - defines the Data type. It may consist of
several data slots. One of these data slots has to be of
type `Key`.
*/
typedef unspecified_type Data;

/*!
the container
`Window` - defines the type of the query rectangle. It
may consist of
several data slots. Two of these data slots has to be of
type `Key`
*/
typedef unspecified_type Window;

/*!
the type
`Key` of the data
slot this traits class provides access to.
*/
typedef unspecified_type Key;

/*!

`Data_func` is a
function object providing an
`operator()` that takes an argument of type `Data` and returns
a component of type `Key`.
*/
typedef unspecified_type Data_func;

/*!

`Window_left_func` is a function objects that
allow to access the
left data slot of container
`Window` which has type `Key`
*/
typedef unspecified_type Window_left_func;

/*!

`Window_right_func` is a function objects that
allow to access the
right data slot of container
`Window` which has type `Key`
*/
typedef unspecified_type Window_right_func;

/*!
defines a comparison relation which must
define a strict ordering of the objects of type
`Key`. If defined, `less<Key>`
is sufficient.
*/
typedef unspecified_type Compare;

/// @}

/// \name Operations
/// @{

/*!
The data slot of
the data item of `d` of type `Key` is
accessed by function object `Data_func`.
*/
Key get_key(Data d);

/*!
The data slot of
the data item of `w` of type `Key` is
accessed by function object
`Window_left_func`.
*/
Key get_left(Window w);

/*!
The data slot of
the data item of `w` of type `Key` is
accessed by function object `Window_right_func`.
*/
Key get_right(Window w);

/*!
returns `Compare(key1, key2)`.
*/
static bool comp(Key& key1, Key&
key2);

/*!
returns `Compare(get_key(data1), get_key(data2))`.
*/
static bool key_comp(Data& data1, Data&
data2);

/// @}

}; /* end Point_traits */
} /* end namespace CGAL */
