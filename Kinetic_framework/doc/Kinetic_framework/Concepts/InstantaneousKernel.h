namespace Kinetic {
/*!
\ingroup PkgKdsFrameworkConcepts
\cgalConcept

The concept `Kinetic::InstantaneousKernel` covers models that act as adaptors allowing 
%CGAL static data structures to act on snapshots of kinetic 
data. Different methods for evaluating predicates are used depending 
on whether time is set using an `NT` or a `Time` 
object. Evaluating predicates when time is the former is much cheaper. 

\cgalHasModel `CGAL::Kinetic::Default_instantaneous_kernel`

*/

class InstantaneousKernel {
public:

/// \name Types 
/// @{

/*!
A number type which can be used to represent the current time. This must be a ring or field type. 
*/ 
typedef unspecified_type NT; 

/*!
The type used to represent the current time. This type must be comparable. 
*/ 
typedef unspecified_type Time; 

/// @} 

/// \name Operations 
/// @{

/*!
Return the current time. 
*/ 
Time time(); 

/*!
Return the current time as an `NT`. As a precondition, `time_is_nt` must be true. 
*/ 
NT time_as_nt(); 

/*!
Return true if the last time time was set, it was using an object of type `NT`. 
*/ 
bool time_is_nt(); 

/*!
Set the current time to have a certain value. All existing predicates are updated automatically. 
*/ 
void set_time(Time); 

/*!
The the current time to be an instance of `NT`. With this more efficient techniques can be used. `time_is_nt()` must be true. 
*/ 
void set_time(NT); 

/*!
Return a static object corresponding to the kinetic object at this instant in time. `time_is_nt()` must be true. 
*/ 
Static_object static_object(Key); 

/// @}

}; /* end Kinetic::InstantaneousKernel */

} /* end namespace KineticConcepts */
