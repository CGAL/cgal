namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * The class `Arr_curve_data_traits_2` is a model of the `ArrangementTraits_2`
 * concept and serves as a decorator class that allows the extension of the
 * curves defined by the base traits-class (the `Tr` parameter), which serves as
 * a geometric traits-class (a model of the `ArrangementTraits_2` concept), with
 * extraneous (non-geometric) data fields.
 *
 * The traits class inherits its point type from `Traits::Point_2`, and defines
 * an extended `Curve_2` and `X_monotone_curve_2` types, as detailed below.
 *
 * Each `Curve_2` object is associated with a single data field of type `CData`,
 * and each `X_monotone_curve_2` object is associated with a single data field
 * of type `XData`. When a curve is subdivided into \f$ x\f$-monotone subcurves,
 * its data field is converted using the conversion functor, which is specified
 * by the `Cnv` template-parameter, and the resulting objects is copied to all
 * `X_monotone_curve_2` objects induced by this curve. The conversion functor
 * should provide an operator with the following prototype:
 *
 * `XData operator() (const CData& d) const;`
 *
 * By default, the two data types are the same, so the conversion operator
 * is trivial:
 *
 * <TABLE><TR><TD>
 * `CData` =
 * </TD>
 * <TD>
 * `XData`
 * </TD></TR>
 * <TR><TD>
 * `Cnv` =
 * </TD>
 * <TD>
 * `_Default_convert_functor<CData,XData>`
 * </TD></TR>
 * </TABLE>
 *
 * In case two (or more) \f$ x\f$-monotone curves overlap, their data fields are
 * merged to a single field, using the merge functor functor, which is specified
 * by the `Mrg` template-parameter. This functor should provide an operator with
 * the following prototype:
 *
 * `XData operator() (const XData& d1, const XData& d2) const;`
 *
 * which returns a single data object that represents the merged data field of
 * `d1` and `d2`. The \f$ x\f$-monotone curve that represents the overlap is
 * associated with the output of this functor.
 *
 * \cgalModels{ArrangementTraits_2}
 */
template <typename Tr, typename XData, typename Mrg, typename CData, typename Cnv>
class Arr_curve_data_traits_2 : public Tr {
public:

  /// \name Types
  /// @{

  /*! the base traits-class.
   */
  typedef Tr                                            Base_traits_2;

  /*! the base curve.
   */
  typedef typename Base_traits_2::Curve_2               Base_curve_2;

  /*! the base \f$ x\f$-monotone curve curve.
   */
  typedef typename Base_traits_2::X_monotone_curve_2    Base_x_monotone_curve_2;

  /*! the point type.
   */
  typedef typename Base_traits_2::Point_2               Point_2;

  /*! the merge functor.
   */
  typedef Mrg                                           Merge;

  /*! the conversion functor.
   */
  typedef Cnv                                           Convert;

  /*! the type of data associated with curves.
   */
  typedef CData                                         Curve_data;

  /*! the type of data associated with \f$ x\f$-monotone curves.
   */
  typedef XData                                         X_monotone_curve_data;

  /// @}

  /*! The `Curve_2` class nested within the curve-data traits extends the
   * `Base_traits_2::Curve_2` type with an extra data field of type `Data`.
   */
  class Curve_2 : public Base_curve_2 {
  public:

    /// \name Creation
    /// @{

    /*! default constructor.
     */
    Curve_2();

    /*! constructs curve from the given `base` curve with uninitialized
     * data field.
     */
    Curve_2(const Base_curve_2& base);

    /*! constructs curve from the given `base` curve with an attached
     * `data` field.
     */
    Curve_2(const Base_curve_2& base, const Data& data);

    /// @}

    /// \name Access Functions
    /// @{

    /*! returns the data field (a non-const version, which returns a reference
     * to the data object, is also available).
     */
    const Curve_data& data() const;

    /*! sets the data field.
     */
    void set_data(const Curve_data& data);

    /// @}

}; /* end Arr_curve_data_traits_2::Curve_2 */

/*! The `X_monotone_curve_2` class nested within the curve-data traits extends
 * the `Base_traits_2::X_monotone_curve_2` type with an extra data field.
 */
class X_monotone_curve_2 : public Base_x_monotone_curve_2 {
public:

  /// \name Creation
  /// @{

  /*! default constructor.
   */
  X_monotone_curve_2();

  /*! constructs an \f$ x\f$-monotone curve from the given `base` curve with
   * uninitialized data field.
   */
  X_monotone_curve_2(const Base_x_monotone_curve_2& base);

  /*! constructs an \f$ x\f$-monotone curve from the given `base` \f$
   * x\f$-monotone curve with an attached `data` field.
   */
  X_monotone_curve_2(const Base_x_monotone_curve_2& base,
                     const X_monotone_curve_data& data);

  /// @}

  /// \name Access Functions
  /// @{

  /*! returns the field (a non-const version, which returns a reference
   * to the data object, is also available).
   */
  const X_monotone_curve_data& data() const;

  /*! sets the data field.
   */
  void set_data(const X_monotone_curve_data& data);

  /// @}

}; /* end Arr_curve_data_traits_2::X_monotone_curve_2 */

}; /* end Arr_curve_data_traits_2 */
} /* end namespace CGAL */
