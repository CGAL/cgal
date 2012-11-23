
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkClasses

MovingObjects are stored in a vector. This means that access is constant 
time, but storage is not generally freed. The only way to be sure is 
to remove all reference counts for the table or to call `clear()`. 

\cgalModels `Kinetic::ActiveObjectsTable`

*/
template< typename MovingObject >
class Active_objects_vector {

}; /* end Kinetic::Active_objects_vector */
} /* end namespace Kinetic */
} /* end namespace CGAL */
