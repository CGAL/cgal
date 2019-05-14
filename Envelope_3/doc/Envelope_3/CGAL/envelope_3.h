namespace CGAL {

/*!
\ingroup PkgEnvelope3

The class-template `Envelope_diagram_2` represents the minimization 
diagram that corresponds to the lower envelope of a set of curves, or the 
maximization diagram that corresponds to their upper envelope. It is 
parameterized by a traits class that must be a model of the 
`EnvelopeTraits_3` concept, and is basically a planar arrangement of 
\f$ x\f$-monotone curves, as defined by this traits class. These \f$ x\f$-monotone 
curves are the projections of boundary curves of \f$ xy\f$-monotone surfaces, 
or the intersection curves between such surfaces, onto the \f$ xy\f$-plane. 
Thus, it is possible to traverse the envelope diagram using the 
methods inherited from the `Arrangement_2` class. 

The envelope diagram extends the arrangement features (namely the vertices, 
halfedges, and faces), such that each feature stores a container of 
originators - namely, the \f$ xy\f$-monotone surfaces (instances of the type 
`EnvTraits::Xy_monotone_surface_3`) that induce the lower envelope 
(or the upper envelope, in case of a maximization diagram) over this 
feature. The envelope diagram provides access methods to these originators. 

*/
template< typename EnvTraits >
class Envelope_diagram_2 : public Arrangement_2<EnvTraits> {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef Envelope_diagram_2<EnvTraits> Self; 

/*!

*/ 
typedef Arrangement_2<EnvTraits> Base; 

/*!
an iterator for the \f$ xy\f$-monotone surfaces that induce a diagram feature. 
Its value-type is `EnvTraits::Xy_monotone_surface_3`. 
*/ 
typedef unspecified_type Surface_const_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
constructs an empty diagram containing one unbounded face, 
which corresponds to the entire plane and has no originators. 
*/ 
Envelope_diagram_2(); 

/*!
copy constructor. 
*/ 
Envelope_diagram_2 (const Self& other); 

/*!
constructs an empty diagram that uses the given `traits` 
instance for performing the geometric predicates. 
*/ 
Envelope_diagram_2 (EnvTraits *traits); 

/// @} 

  class Vertex : public Arrangement_2<EnvTraits>::Vertex {
  public:
    /// \name Access Functions 
    /// @{
    
    /*!
      returns the number of \f$ xy\f$-monotone surfaces that induce `v`. 
    */ 
    size_t number_of_surfaces () const; 

    /*!
      returns an iterator for the first \f$ xy\f$-monotone surface that induces `v`. 
    */ 
    Surface_const_iterator surfaces_begin () const; 

    /*!
      returns a past-the-end iterator for the \f$ xy\f$-monotone surfaces that induce 
      `v`. 
    */ 
    Surface_const_iterator surfaces_end () const; 

    /*!
      returns the first \f$ xy\f$-monotone surface that induce `v`. 
      \pre The number of surfaces is not 0. 
    */ 
    Xy_monotone_surface_3 surface () const; 

    /// @}
  };

  class Halfedge : public Arrangement_2<EnvTraits>::Halfedge {
  public:
    /// \name Access Functions 
    /// @{
    
    /*!
    returns the number of \f$ xy\f$-monotone surfaces that induce the halfedge. 
    */ 
    size_t number_of_surfaces () const; 
    
    /*!
    returns an iterator for the first \f$ xy\f$-monotone surface that induces the halfedge. 
    */ 
    Surface_const_iterator surfaces_begin () const; 
    
    /*!
    returns a past-the-end iterator for the \f$ xy\f$-monotone surfaces that induce 
    the halfedge. 
    */ 
    Surface_const_iterator surfaces_end () const; 
    
    /*!
    returns the first \f$ xy\f$-monotone surface that induce `e`. 
    \pre The number of surfaces is not 0. 
    */ 
    Xy_monotone_surface_3 surface () const; 
    /// @}
  };

  class Face : public Arrangement_2<EnvTraits>::Face {
  public:
    /// \name Access Functions 
    /// @{
    
    /*!
    returns the number of \f$ xy\f$-monotone surfaces that induce the face. 
    */ 
    size_t number_of_surfaces () const; 
    
    /*!
    returns an iterator for the first \f$ xy\f$-monotone surface that induces the face. 
    */ 
    Surface_const_iterator surfaces_begin () const; 
    
    /*!
    returns a past-the-end iterator for the \f$ xy\f$-monotone surfaces that induce 
    the face. 
    */ 
    Surface_const_iterator surfaces_end () const; 
    
    /*!
    returns the first \f$ xy\f$-monotone surface that induce the face. 
    \pre The number of surfaces is not 0. 
    */ 
    Xy_monotone_surface_3 surface () const; 
    
    /// @}
  };

/// @}

}; /* end  */
} /* end namespace CGAL */
namespace CGAL {

/*!
\ingroup PkgEnvelope3

Computes the lower envelope of a set of surfaces in \f$ \mathbb{R}^3\f$,
as given by the range `[begin, end)`. The lower envelope is
represented using the output minimization diagram `diag`.
\tparam Traits must be a model of `EnvelopeTraits_3`.
InputIterator must be an input iterator with value type `Traits::Surface_3`.
*/
template<class InputIterator, class Traits>
void lower_envelope_3 (InputIterator begin, InputIterator end,
Envelope_diagram_2<Traits>& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope3

Computes the lower envelope of a set of \f$ xy\f$-monotone surfaces in
\f$ \mathbb{R}^3\f$, as given by the range `[begin, end)`. The lower 
envelope is represented using the output minimization diagram `diag`.
\tparam Traits must be a model of `EnvelopeTraits_3`.
InputIterator must be an input iterator with value type `Traits::Xy_monotone_surface_3`.
*/
template<class InputIterator, class Traits>
void lower_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
Envelope_diagram_2<Traits>& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope3

Computes the upper envelope of a set of surfaces in \f$ \mathbb{R}^3\f$,
as given by the range `[begin, end)`. The upper envelope is
represented using the output maximization diagram `diag`.
\tparam Traits must be a model of `EnvelopeTraits_3`.
InputIterator must be an input iterator with value type `Traits::Surface_3`.
*/
template<class InputIterator, class Traits>
void upper_envelope_3 (InputIterator begin, InputIterator end,
Envelope_diagram_2<Traits>& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope3

Computes the upper envelope of a set of \f$ xy\f$-monotone surfaces in 
\f$ \mathbb{R}^3\f$, as given by the range `[begin, end)`. The lower 
envelope is represented using the output maximization diagram `diag`.
\tparam Traits must be a model of `EnvelopeTraits_3`.
InputIterator must be an input iterator with value type `Traits::Xy_monotone_surface_3`.
*/
template<class InputIterator, class Traits>
void upper_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
Envelope_diagram_2<Traits>& diag);

} /* namespace CGAL */
