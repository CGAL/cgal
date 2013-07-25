
namespace CGAL {

/*!
\ingroup PkgNef3

This class helps to perform the union of a set of 3D Nef polyhedra 
efficiently. It succesively applies the binary union operation of 
`Nef_polyhedron_3`, but schedules these union operations in an 
opportune way. The class is most efficient, if the polyhedra are added 
in sorted order. Any order that reflects proximity in the 
three-dimensional space is helpful. To allow saving memory space, the 
sorting is left to the user. This way the user can generate the 
polyhedra in a sorted way and add them one by one to `Nef_nary_union_3`. 



\tparam NefPolyhedron_3 must be  an instantiation of the class template `Nef_polyhedron_3`. 



*/
template< typename Nef_polyhedron_3 >
class Nef_nary_union_3 {
public:

/// \name Creation 
/// @{

/*!
initialization only. 
*/ 
Nef_nary_union_3<NefPolyhedron_3>(); 

/// @} 

/// \name Member Functions 
/// @{

/*!
returns the union of the polyhedra previously added to the class. 
*/ 
NefPolyhedron_3 get_union() const; 

/*!
adds a polyhedron. 
*/ 
void add_polyhedron(const NefPolyhedron_3& N); 

/// @}

}; /* end Nef_nary_union_3 */
} /* end namespace CGAL */
