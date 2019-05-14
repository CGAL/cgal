namespace CGAL {

/*!
\ingroup PkgHDS_Decorators

 The class 
`CGAL::HalfedgeDS_items_decorator<HDS>` provides additional functions 
for vertices, halfedges, and faces of a halfedge data structure 
without knowing the containing halfedge data structure. The class 
`CGAL::HalfedgeDS_decorator<HDS>` stores a reference to the halfedge 
data structure and provides functions that modify the halfedge data 
structure, for example Euler-operators. The class 
`CGAL::HalfedgeDS_const_decorator<HDS>` stores a const reference to 
the halfedge data structure. It contains non-modifying functions, for 
example the test for validness of the data structure. 

All these additional functions take care of the different capabilities 
a halfedge data structure may have or may not have. The functions 
evaluate the type tags of the halfedge data structure to decide on the 
actions. If a particular feature is not supported nothing is done. 
Note that for example the creation of new halfedges is mandatory for 
all halfedge data structures and will not appear here again. 

\tparam HDS must be a model of `HalfedgeDS`


\cgalHeading{Example}

The following program fragment illustrates the implementation of a 
`is_valid()` member function for a simplified polyhedron class. 
We assume here that the level three check is the appropriate default 
for polyhedral surfaces. 

\code{.cpp} 

namespace CGAL { 
  template <class Traits> 
  class Polyhedron { 
    typedef HalfedgeDS_default<Traits> HDS; 
    HDS hds; 
  public: 
    // ... 
    bool 
    is_valid( bool verb = false, int level = 0) const { 
      Verbose_ostream verr(verb); 
      verr << "begin Polyhedron::is_valid( verb=true, level = " << level << "):" 
           << std::endl; 
      HalfedgeDS_const_decorator<HDS> decorator(hds); 
      bool valid = decorator.is_valid( verb, level + 3); 
      // further checks ... 
    } 
  }; 
} 

\endcode 

\sa `CGAL::HalfedgeDS_items_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_decorator<HDS>` 

*/
template< typename HDS >
class HalfedgeDS_const_decorator : public CGAL::HalfedgeDS_items_decorator<HDS> {
public:

/// \name Creation 
/// @{

/*!
keeps internally a const reference to `hds`. 
*/ 
HalfedgeDS_const_decorator( const HDS& hds); 

/// @} 

/*! \name Validness Checks 
A halfedge data structure has no definition of validness of its own, 
but a useful set of tests is defined with the following levels: 
<DL> 
<DT><B>Level 0</B><DD> 
The number of halfedges is even. All pointers except 
the vertex pointer and the face pointer for border halfedges are 
unequal to their respective default construction value. For all 
halfedges `h`: The opposite halfedge is different from `h` and the 
opposite of the opposite is equal to `h`. The next of the previous 
halfedge is equal to `h`. For all vertices `v`: the incident vertex 
of the incident halfedge of `v` is equal to `v`. The halfedges 
around `v` starting with the incident halfedge of `v` form a cycle. 
For all faces `f`: the incident face of the incident halfedge of `f` 
is equal to `f`. The halfedges around `f` starting with the incident 
halfedge of `f` form a cycle. Redundancies among internal variables 
are tested, e.g., that iterators enumerate as many items as the 
related size value indicates. 
<DT><B>Level 1</B><DD> 
All tests of level 0. For all halfedges `h`: The 
incident vertex of `h` exists and is equal to the incident vertex of 
the opposite of the next halfedge. The incident face (or hole) of 
`h` is equal to the incident face (or hole) of the next halfedge. 
<DT><B>Level 2</B><DD> 
All tests of level 1. The sum of all halfedges that can 
be reached through the vertices must be equal to the number of all 
halfedges, i.e., all halfedges incident to a vertex must form a single 
cycle. 
<DT><B>Level 3</B><DD> 
All tests of level 2. The sum of all halfedges that can 
be reached through the faces must be equal to the number of all 
halfedges, i.e., all halfedges surrounding a face must form a single 
cycle (no holes in faces). 
<DT><B>Level 4</B><DD> 
All tests of level 3 and `normalized_border_is_valid()`. 
</DL> 
*/
/// @{

/*!
returns `true` if the halfedge data structure `hds` is 
valid with respect to the `level` value as defined above. 
If `verbose` is `true`, statistics are written to `cerr`. 
*/ 
bool is_valid( bool verbose = false, int level = 0) const; 

/*!
returns `true` if the border halfedges are in normalized 
representation, which is when enumerating all halfedges with the 
halfedge iterator the following holds: The non-border edges precede the 
border edges. For border edges, the second halfedge is a border halfedge. 
(The first halfedge may or may not be a border halfedge.) The halfedge 
iterator `HalfedgeDS::border_halfedges_begin()` denotes the first border 
edge. If `verbose` is `true`, statistics are written to `cerr`. 

*/ 
bool normalized_border_is_valid( bool verbose = false) const; 

/// @}

}; /* end HalfedgeDS_const_decorator */
} /* end namespace CGAL */
