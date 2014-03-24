namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Implicit_multi_domain_to_labeled_function_wrapper` wraps a set of implicit function [f1,f2,...] to one function F which
takes its values into N. This wrapper allows the user to define several domain from implicit functions and a vector of sign.

We associate a power of two to each function of the set. For a point p, F(p) will be equal to the sum
of powers of two of the functions where f(p) has the sign expected by the user.

\b Creation \b of \b the \b wrapper<br />
  [f1, f2, f3]<br />
  [ -,  -,  +]

\b Input<br />
  Point_3 p

\b Output<br />
  int N = 0b00..00( f1(p)<0 )( f2(p)<0 )( f3(p)>0 )

(Functions (f1, f2, f3) must take a point as input parameter and return an arithmetic value.)


This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.

\tparam Function is the type of the input function.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\sa `Implicit_mesh_domain_3`.
\sa `Sign`

*/
template<class Function, class BGT>
class Implicit_multi_domain_to_labeled_function_wrapper
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
   * The array of signs is built automatically with NEGATIVE signs.
   */
  Implicit_multi_domain_to_labeled_function_wrapper(const std::vector<Function>& implicit_functions);

  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   * \param mask the vector of signs expected by the user.
   */
  Implicit_multi_domain_to_labeled_function_wrapper (const std::vector<Function>& implicit_functions, const std::vector<Sign>& mask);

  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   * \param mask_str a string from which the vector of sign will be built. It must contains '+' or '-' only.
   */
  Implicit_multi_domain_to_labeled_function_wrapper (const std::vector<Function>& implicit_functions, const std::string& mask_str);
/// @}

}; /* end Implicit_multi_domain_to_labeled_function_wrapper */

/*!
\ingroup PkgMesh_3Domains

\deprecated
Use `Implicit_multi_domain_to_labeled_function_wrapper` instead.

The class `Implicit_vector_to_labeled_function_wrapper` wraps a set of implicit function [f1,f2,...] to one function F which
takes its values into N.

Let p be a point, f1 and f2 two functions:
F(p) = 0b000000(f2(p)<0)(f1(p)<0)
(Functions must take a point as input parameter and return an arithmetic value.)

It can handle at most 8 functions.

This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.

\tparam Function is the type of the input function.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\sa `Implicit_mesh_domain_3`.

*/
template<class Function, class BGT>
class Implicit_vector_to_labeled_function_wrapper
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
  Implicit_vector_to_labeled_function_wrapper(const std::vector<Function*>& implicit_functions);
}; /* end Implicit_vector_to_labeled_function_wrapper */
/// @}
} /* end namespace CGAL */
