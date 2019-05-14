namespace CGAL {

/*!
\ingroup PkgEnvelope3

The class `Env_surface_data_traits_3` is a model of the `EnvelopeTraits_3` concept 
and serves as a decorator class that allows the extension of the surfaces 
defined by the base traits-class (the `Traits` parameter), which serves 
as a geometric traits-class (a model of the `EnvelopeTraits_3` concept), 
with extraneous (non-geometric) data fields. 

The traits-class decorator extends the `Surface_3` and the 
`Xy_monotone_surface_3` types as detailed below. 
Each `Surface_3` object is associated with a single data field of type 
`SData`, and each `Xy_monotone_surface_3` object is associated with 
a single data field of type `XyData`. When a surface is 
subdivided into \f$ xy\f$-monotone surfaces, its data field is converted using 
the conversion functor, which is specified by the `Cnv` template-parameter, 
and the resulting objects is copied to all `Xy_monotone_surface_3` objects 
induced by this surface. The conversion functor should provide an operator with 
the following prototype: 

`XyData operator() (const SData& d) const;` 

By default, the two data types are the same, so the conversion operator 
is trivial and just copies the data object: 

<TABLE><TR><TD>
`SData` = 
</TD>
<TD>
`XyData` 
</TD></TR>
<TR><TD>
`Cnv` = 
</TD>
<TD>
`_Default_convert_functor<SData,XyData>` 
</TD></TR>
</TABLE> 

\cgalModels `EnvelopeTraits_3`

*/
template< typename Traits, typename XyData, typename SData, typename Cnv >
class Env_surface_data_traits_3 : public Traits {
public:

/// \name Types 
/// @{

/*!
the base traits-class. 
*/ 
typedef Traits Base_traits_3; 

/*!
the base surface. 
*/ 
typedef typename Base_traits_3::Surface_3 Base_surface_3; 

/*!
the base \f$ xy\f$-monotone surface surface. 
*/ 
typedef typename Base_traits_3::Xy_monotone_surface_3 Base_xy_monotone_surface_3; 

/*!
the conversion functor. 
*/ 
typedef Cnv Convert; 

/*!
the type of data associated with surfaces. 
*/ 
typedef SData Surface_data; 

/*!
the type of data associated with \f$ xy\f$-monotone surfaces. 
*/ 
typedef XyData Xy_monotone_surface_data; 
  
  
  /// The `Surface_3` class nested within the surface-data traits 
  /// extends the `Base_traits_3::Surface_3` type with an extra data field. 
  class Surface_3 : public Base_surface_3 {
  public:
    /// \name Creation 
    /// @{
    
    /*!
    default constructor. 
    */ 
    Surface_3 (); 
    
    /*!
    constructs surface from the given `base` surface with uninitialized 
    data field. 
    */ 
    Surface_3 (const Base_surface_3& base); 
    
    /*!
    constructs surface from the given `base` surface with an attached 
    `data` field. 
    */ 
    Surface_3 (const Base_surface_3& base, const Surface_data& data); 
    
    /// @} 
    
    /// \name Access Functions 
    /// @{
    
    /*!
    returns the data field (a non-const version, which returns a reference 
    to the data object, is also available). 
    */ 
    const Surface_data& data () const; 
    
    /*!
    sets the data field. 
    */ 
    void set_data (const Surface_data& data); 
    
    /// @} 
  };

  /// The `Xy_monotone_surface_3` class nested within the surface-data traits
  /// extends the `Base_traits_3::Xy_monotone_surface_3` type with an extra data
  /// field.
  class  Xy_monotone_surface_3 : public Base_xy_monotone_surface_3 {
  public:
    /// \name Creation 
    /// @{
    
    /*!
    default constructor. 
    */ 
    Xy_monotone_surface_3 (); 
    
    /*!
    constructs an \f$ xy\f$-monotone surface from the given `base` surface with 
    uninitialized data field. 
    */ 
    Xy_monotone_surface_3 (const Base_xy_monotone_surface_3& base); 
    
    /*!
    constructs an \f$ x\f$-monotone surface from the given `base` \f$ x\f$-monotone 
    surface with an attached `data` field. 
    */ 
    Xy_monotone_surface_3 (const Base_xy_monotone_surface_3& base, 
    const Xy_monotone_surface_data& data); 
    
    /// @} 
    
    /// \name Access Functions 
    /// @{
    
    /*!
    returns the field (a non-const version, which returns a reference 
    to the data object, is also available). 
    */ 
    const Xy_monotone_surface_data& data () const; 
    
    /*!
    sets the data field. 
    */ 
    void set_data (const Xy_monotone_surface_data& data); 
    
    /// @}
  };

/// @} 

}; /* end */
} /* end namespace CGAL */
