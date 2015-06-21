namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Implicit_multi_domain_to_labeling_function_wrapper` is an helping class to get a function with integer values
labeling the components of a multi-domain. The multidomain is described through a set of functions {fi(p), i=1, ...n}.
Each component corresponds to a sign vector [s1, s2, ..., sn] where si is the sign of the function fi(p) at a point p of the component.
This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.

\par Example
For example, the multidomain described by the three functions [f1,f2,f3] and the two sign vectors [-,-,+] and [+,-,+]
 includes two components.<br />
The first one matches the locus of points satisfying f1(p)<0 and f2(p)<0 and f3(p)>0.<br />
The second one matches the locus of points satisfying f1(p)>0 and f2(p)<0 and f3(p)>0.<br />

\tparam Function provides the definition of the function.
This parameter stands for a model of the concept ImplicitFunction described in the surface mesh generation package.
The number types Function::FT and BGT::FT are required to match.

\sa `Implicit_mesh_domain_3`.

*/
template<class Function>
class Implicit_multi_domain_to_labeling_function_wrapper
{
public:
  /// \name Types
  /// @{
  //!
  typedef std::vector<Function*>   Function_vector;
  //!
  typedef typename Function::Point     Point_3;
  /// @}

  /// \name Creation
  /// @{
  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   *
   * Position vectors are built automatically so that the union of components equals the union of the functions.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions);

  /*!
   * \brief Construction from a vector of implicit functions and a vector of vector of signs.
   * \param implicit_functions the vector of implicit functions.
   * \param position_vectors the vector of vector of signs. Each vector of positions describes a component.
   * \sa `Sign`
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions, const std::vector<std::vector<Sign> >& position_vectors);

  /*!
   * \brief Construction from a vector of implicit functions and a vector of strings.
   * \param implicit_functions the vector of implicit functions.
   * \param position_strings the vector of strings. The strings contained in this vector must contain '+' or '-' only. Each string (vector of positions) describes a component.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& implicit_functions, const std::vector<std::string>& position_strings);
/// @}

}; /* end Implicit_multi_domain_to_labeling_function_wrapper */
/// @}
} /* end namespace CGAL */
