/// 
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  Not all values of a type need to be valid. The function object 
///  class `Is_valid` checks this. 
///
///  For example, an expression like
///  `NT(0)/NT(0)` can result in an invalid number. 
///  Routines may have as a precondition that all values are valid.
///  
///  \models ::AdaptableFunctor
template<class T>
class Is_valid {
/*!
  returns if the argument is valid.
*/
bool operator()(const T& x);
};

///  
///  The function object class `Max` returns the larger of two values.
///  The order is induced by the second template argument <TT>Less</TT>.
///  The default value for `Less` is `std::less`.
///  Note that `T` must be a model of `LessThanComparable`
///  in case `std::less` is used.
///  
///  \models ::AdaptableFunctor
template<typename T, typename Less>
class Max {
/*!
 default constructor. 
*/

Max();
/*!
 The constructed object will use `c` to compare the arguments. 
*/
Max(Less c);

/*!
 returns the larger of `x` and `y`,  with respect to the order induced by `Less`. 
*/
T operator()(const T& x, const T& y);
};
template 


/// @}
} // namespace CGAL

namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function object class `Min` returns the smaller of two values.
///  The order is induced by the second template argument <TT>Less</TT>.
///  The default value for `Less` is `std::less`.
///  Note that `T` must be a model of `LessThanComparable`
///  in case `std::less` is used.
///  
template<typename T, typename Less>
class Min {
/*!
 default constructor. 
*/
Min();
  
/*!
 The constructed object will use `c` to compare the arguments. 
*/
Min(Less c);

/*!
 returns the larger of `x` and `y`,  with respect to the order induced by `Less`. 
*/
T operator()(const T& x, const T& y);

};

///  \models ::AdaptableFunctor
/// @}
} // namespace CGAL
