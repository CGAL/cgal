
namespace CGAL {

/*!
\ingroup PkgProperty_map

Property map that converts a `T*` pointer (or in general an iterator
over `T` elements) to the `T` object.

\models `LvaluePropertyMap`

*/
template< typename T >
class Dereference_property_map : 
    public boost::put_get_helper< T&, 
                                  Dereference_property_map<T> > {
public:

  /// \name Types 
  /// @{

  /*! 

  typedef to `T*` 

  */ 
  typedef Hidden_type key_type; 

  /*! 

  typedef to `T` 

  */ 
  typedef Hidden_type value_type; 

  /*! 

  typedef to `T&` 

  */ 
  typedef Hidden_type reference; 

  /*! 

  `boost::lvalue_property_map_tag` 

  */ 
  typedef Hidden_type category; 

  /// @} 

  /// \name Creation 
  /// @{

  /*! 

  %Default constructor. 

  */ 
  Dereference_property_map(); 

  /// @} 

  /// \name Operations 
  /// @{

  /*! 

  Access a property map element. 

  \tparam Iter Type convertible to `key_type`. 
  */ 
  template<class Iter> reference operator[](Iter it) const; 


  /// @}

}; /* end Dereference_property_map */

/*! 

Free function to create a `Dereference_property_map` property map. 

\relates Dereference_property_map 
*/ 
Dereference_property_map<typename value_type_traits<Iter>::type> make_dereference_property_map(Iter); 


/*!
\ingroup PkgProperty_map

Property map that accesses the first item of a `std::pair`. 

\tparam Pair Instance of `std::pair`. 

\models `LvaluePropertyMap`

\sa `CGAL::Second_of_pair_property_map<Pair>`


*/
template< typename Pair >
class First_of_pair_property_map : 
    public boost::put_get_helper< Pair::first_type&, First_of_pair_property_map<Pair> >
{
public:

  /// \name Types 
  /// @{

  /*! 

  typedef to `Pair*` 

  */ 
  typedef Hidden_type key_type; 

  /*! 

  typedef to `Pair::first_type` 

  */ 
  typedef Hidden_type value_type; 

  /*! 

  typedef to `value_type`& 

  */ 
  typedef Hidden_type reference; 

  /*! 

  `boost::lvalue_property_map_tag` 

  */ 
  typedef Hidden_type category; 

  /// @} 

  /// \name Creation 
  /// @{

  /*! 

  Default constructor. 

  */ 
  First_of_pair_property_map(); 

  /// @} 

  /// \name Operations 
  /// @{

  /*! 

  Access a property map element. 

  \tparam Iter Type convertible to `key_type`. 

  */ 
  template<class Iter> reference operator[](Iter pair) const; 

  /// @}

}; /* end First_of_pair_property_map */

/*! 
Free function to create a `First_of_pair_property_map` property map. 

\relates First_of_pair_property_map 
*/ 
First_of_pair_property_map<typename value_type_traits<Iter>::type> make_first_of_pair_property_map(Iter); 


/*!
\ingroup PkgProperty_map

Property map that accesses the Nth item of a `boost::tuple`. 

\tparam N Index of the item to access.
\tparam Tuple Instance of `boost::tuple`.

\models `LvaluePropertyMap`

*/
template< typename N, typename Tuple >
class Nth_of_tuple_property_map : 
    public boost::put_get_helper< boost::tuples::element<N, Tuple>::type& , Nth_of_tuple_property_map<N, Tuple>> {
public:

  /// \name Types 
  /// @{

  /*! 

  typedef to `Tuple*` 

  */ 
  typedef Hidden_type key_type; 

  /*! 

  typedef to `boost::tuples::element<N, Tuple>::type` 

  */ 
  typedef Hidden_type value_type; 

  /*! 

  typedef to `value_type`& 

  */ 
  typedef Hidden_type reference; 

  /*! 

  `boost::lvalue_property_map_tag` 

  */ 
  typedef Hidden_type category; 

  /// @} 

  /// \name Creation 
  /// @{

  /*! 

  %Default constructor. 

  */ 
  Nth_of_tuple_property_map(); 

  /// @} 

  /// \name Operations 
  /// @{

  /*! 

  Access a property map element. 

  \tparam Iter Type convertible to `key_type`. 

  */ 
  template<class Iter> reference operator[](Iter tuple) const; 


  /// @}

}; /* end Nth_of_tuple_property_map */

/*! 

Free function to create a `Nth_of_tuple_property_map` property map. 

\relates Nth_of_tuple_property_map 
*/ 
Nth_of_tuple_property_map<N, typename value_type_traits<Iter>::type> make_nth_of_tuple_property_map(Iter); 


/*!
\ingroup PkgProperty_map

Property map that accesses the second item of a `std::pair`. 

\tparam Pair Instance of `std::pair`. 

\models `LvaluePropertyMap`

\sa `CGAL::First_of_pair_property_map<Pair>`

*/
template< typename Pair >
class Second_of_pair_property_map : public 
boost::put_get_helper< Pair::second_type& , Second_of_pair_property_map<Pair>> 
{
public:

  /// \name Types 
  /// @{

  /*! 

  typedef to `Pair*` 

  */ 
  typedef Hidden_type key_type; 

  /*! 

  typedef to `Pair::second_type` 

  */ 
  typedef Hidden_type value_type; 

  /*! 

  typedef to `value_type`& 

  */ 
  typedef Hidden_type reference; 

  /*! 

  `boost::lvalue_property_map_tag` 

  */ 
  typedef Hidden_type category; 

  /// @} 

  /// \name Creation 
  /// @{

  /*! 

  %Default constructor. 

  */ 
  Second_of_pair_property_map(); 

  /// @} 

  /// \name Operations 
  /// @{

  /*! 
  Access a property map element. 

  \tparam Iter Type convertible to `key_type`. 

  */ 
  template<class Iter> reference operator[](Iter pair) const; 


  /// @}

}; /* end Second_of_pair_property_map */

/*! 
Free function to create a `Second_of_pair_property_map` property map. 

\relates Second_of_pair_property_map 
*/ 
Second_of_pair_property_map<typename value_type_traits<Iter>::type> make_second_of_pair_property_map(Iter); 


} /* end namespace CGAL */
