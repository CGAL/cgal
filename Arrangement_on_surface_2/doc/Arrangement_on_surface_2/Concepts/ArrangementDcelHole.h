
/*!
\ingroup PkgArrangement2ConceptsDCEL
\cgalConcept

A hole record in a <span class="textsc">Dcel</span> data structure, which stores the face that contains 
the hole in its interior, along with an iterator for the hole in the holes' 
container of this face. 

\sa `ArrangementDcel` 
\sa `ArrangementDcelFace` 

*/

class ArrangementDcelHole {
public:

/// \name Types 
/// @{

/*!
the corresponding <span class="textsc">Dcel</span> face type. 
*/ 
typedef unspecified_type Face; 

/*!

*/ 
typedef Face::Hole_iterator Hole_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Arr_dcel_hole(); 

/// @} 

/// \name Access Functions 
/// All functions below also have `const` counterparts, returning non-mutable pointers or iterators:
/// @{

/*!
returns the incident face, which contains `ho` in its interior. 
*/ 
Face* face (); 

/*!
returns an iterator for the hole. 
*/ 
Hole_iterator iterator(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
sets the incident face. 
*/ 
void set_face (Face* f); 

/*!
sets the hole iterator. 
*/ 
void set_iterator (Hole_iterator it); 

/// @}

}; /* end ArrangementDcelHole */

