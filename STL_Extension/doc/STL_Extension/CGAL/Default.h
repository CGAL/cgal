namespace CGAL {

/*!
\ingroup PkgSTLExtensionRef



`Default` is a tag class. It can be used to state that one wants
to use the default argument of a template parameter of a class template.

This can be useful in several cases: (a) when one needs a non-default value for
another template parameter coming next (since \cpp only supports defaults at
the end of lists), (b) when the default is actually a complex expression, e.g.
referring to previous template parameters (in this case, it shortens compiler
error messages and mangled symbol names), (c) when defining the default
involves circular dependencies of type instantiations (there, it breaks the
cycle in a nice way).

Using the mechanism is easy: just plug `Default` as template argument
in the place where you would like to use the default. You should refer
to the documentation of the template class you are using in order to know
whether this functionality is offered.

Also beware that the type of the instantiated template class will not be the
same when instantiating it using `Default` instead of the type of the default
argument, even though their interfaces will otherwise be the same. This may
have consequences in some cases.

\cgalModels{DefaultConstructible,CopyConstructible}

In order to help the template class writer, `Default` provides a convenient
way to extract the desired type for a template parameter which may be defaulted
using `Default`. It is enough to fetch the type as
`Default::Get<Parameter, Type>::%type`, as in the example program below.

\cgalHeading{Example}

\cgalExample{STL_Extension/Default.cpp}

*/
struct Default {


/// \name Types
/// @{
/*!
A nested template providing a typedef `type` which equals `Type` if
`Parameter` is `Default`, and `Parameter` otherwise.
*/
template <typename Parameter, typename Type>
struct Get {
  typedef unspecified_type type;
};
/*!
A nested template providing a typedef `type` which equals `MetaFct::type` if
`Parameter` is `Default`, and `Parameter` otherwise. Use this template
instead of `CGAL::Default::Get` when the default parameter should not be
instantiated by default.
*/
template <typename Parameter, typename MetaFct>
struct Lazy_get {
  typedef unspecified_type type;
};/// @}


}; /* end Default */
} /* end namespace CGAL */
