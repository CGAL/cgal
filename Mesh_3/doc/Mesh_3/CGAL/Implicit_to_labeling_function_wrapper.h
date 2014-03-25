namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Implicit_multi_domain_to_labeling_function_wrapper` wraps several implicit functions [f1,f2,...] to one labeling function F
that takes its values into N. This wrapper allows the user to define several domains from implicit functions and a vector of sign.
Implicit functions must take a point as input parameter and return an arithmetic value.

We associate a power of two to each implicit function of the set.
For each implicit function f, we look at the sign \e s of f(p). If \e s is the same as
the sign given in the vector of signs, we add the power of two associated to f to F(p).

\par Example
<pre>
\e Creation \e of \e the \e wrapper\n
    [f1, f2, f3]
    [ -,  -,  +]
\e Input\n
    Point_3 p
\e Output\n
    int N = 0b00...00( f1(p)<0 )( f2(p)<0 )( f3(p)>0 )
</pre>

This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.

\tparam Function is the type of the input implicit functions.

\tparam BGT is a geometric traits class that provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\sa `Implicit_mesh_domain_3`.
\sa `Sign`

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
   *
   * The array of signs is built automatically with `NEGATIVE` signs.
   * \sa Sign
   */
  Implicit_multi_domain_to_labeling_function_wrapper(const std::vector<Function>& implicit_functions);

  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   * \param mask the vector of signs expected by the user.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const std::vector<Function>& implicit_functions, const std::vector<Sign>& mask);

  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   * \param mask_str a string from which the vector of signs will be built. It must contains '+' or '-' only.
   */
  Implicit_multi_domain_to_labeling_function_wrapper (const std::vector<Function>& implicit_functions, const std::string& mask_str);
/// @}

}; /* end Implicit_multi_domain_to_labeling_function_wrapper */
/// @}
} /* end namespace CGAL */
