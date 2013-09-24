namespace CGAL {

/*!
\ingroup PkgStlExtensionAssertions

*/

/// @{

/*!
 */
enum Failure_behaviour { ABORT, EXIT, EXIT_WITH_SUCCESS, CONTINUE,
                         THROW_EXCEPTION };


/*!

\param type is a string that contains one of the words precondition, postcondition, assertion or warning. 
\param expression contains the expression that was violated. 
\param file, line contain the place where the check was made. 
\param explanation contains an explanation of what was checked. It can be `NULL`, in which case the expression is thought to be descriptive enough.


 */
typedef
    void
    (*Failure_function)(
        const char* type, const char* expression, const char* file, int line, const char* explanation);


/*!
 */
Failure_function
set_error_handler( Failure_function handler);


/*!
 */
Failure_function
set_warning_handler( Failure_function handler);


/*!
 */
Failure_behaviour
set_error_behaviour(Failure_behaviour eb);


/*!
 */
Failure_behaviour
set_warning_behaviour(Failure_behaviour eb);

/// @}

} //namespace CGAL
