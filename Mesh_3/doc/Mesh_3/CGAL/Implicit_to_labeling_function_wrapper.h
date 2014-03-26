namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Implicit_multi_domain_to_labeling_function_wrapper` is an helping class to get a function with integer values
labeling the components of a multi-domain. This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.
Each component corresponds to a sign vector [s1, s2, ..., sn] where si is the sign of the function fi(p) at a point p of the component.

\par Example
<pre>
\e Creation \e of \e the \e wrapper\n
    [f1, f2, f3]
    [
      [ -,  -,  +]
      [ +,  -,  +]
    ]
\e Input\n
    %Point_3 p
\e Output\n
    int N = 1  if f1(p)<0 and f2(p)<0 and f3(p)>0
            2  if f1(p)>0 and f2(p)<0 and f3(p)>0
            0  else
</pre>

\tparam Function is the type of the input implicit functions.

\tparam BGT is a geometric traits class that provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\sa `Implicit_mesh_domain_3`.

*/
template<class Function, class BGT>
class Implicit_multi_domain_to_labeling_function_wrapper
{
public:
  /// \name Types
  /// @{
  //!
  typedef std::vector<Function*>   Function_vector;
  //!
  typedef typename BGT::Point_3     Point_3;
  /// @}

  /// \name Creation
  /// @{
  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions);

  /*!
   * \brief Construction from a vector of implicit functions and a vector of vector of signs.
   * \param implicit_functions the vector of implicit functions.
   * \param positions_vectors the vector of vector of signs.
   * \sa `Sign`
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions, const std::vector<std::vector<Sign> >& positions_vectors);

  /*!
   * \brief Construction from a vector of implicit functions and a vector of strings.
   * \param implicit_functions the vector of implicit functions.
   * \param positions_strings the vector of string. The strings contained in this vector must contain '+' or '-' only.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions, const std::vector<std::string>& positions_strings);
/// @}

}; /* end Implicit_multi_domain_to_labeling_function_wrapper */
/// @}
} /* end namespace CGAL */
