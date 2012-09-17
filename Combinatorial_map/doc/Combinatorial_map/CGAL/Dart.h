
namespace CGAL {

/*!
\ingroup PkgCombinatorialMaps

The class `Dart` represents a <I>d</I>D dart. 

\f$ \beta_i\f$ pointers are coded in a array of <I>d+1</I> `Dart_handle` 
(because we describe also the \f$ \beta_0\f$ link). Attributes are 
associated to each dart by `Attribute_handle<i>`, one for each 
non void <I>i</I>-attribute. 

\models ::Dart 

\tparam d an integer for the dimension of the dart. 

\tparam CMap must be a model of the `CombinatorialMap` concept. 

###Complexity ###

Each \f$ \beta_i\f$ link is initialized to `CMap::null_dart_handle`, and each 
attribute handle of non void <I>i</I>-attribute is initialized to `NULL` 
at the creation of the dart, thus the complexity of the creation is in 
<I>O</I>(<I>d+1</I>). 

The complexity of `opposite` and `other_extremity` methods is in 
<I>O</I>(<I>d+1</I>). 

Other methods have all a constant time complexity. 

\sa `CombinatorialMap` 

*/
template< typename d, typename CMap >
class Dart {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef CMap::Dart_handle Dart_handle; 

/*! 

*/ 
typedef CMap::Dart_const_handle Dart_const_handle; 

/*! 

*/ 
typedef CMap::Attribute_handle<i>::type Attribute_handle<i>::type; 

/*! 

*/ 
typedef CMap::Attribute_const_handle<i>::type Attribute_const_handle<i>::type; 

/// @}

}; /* end Dart */
} /* end namespace CGAL */
