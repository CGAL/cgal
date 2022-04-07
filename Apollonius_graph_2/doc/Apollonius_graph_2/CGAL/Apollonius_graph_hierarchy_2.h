
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

We provide an alternative to the class
`Apollonius_graph_2<Gt,Agds>` for the dynamic
construction of the Apollonius graph. The `Apollonius_graph_hierarchy_2` class maintains
a hierarchy of Apollonius graphs. The bottom-most level of the
hierarchy contains the full Apollonius diagram. A site that
is in level \f$ i\f$, is in level \f$ i+1\f$ with probability \f$ 1/\alpha\f$
where \f$ \alpha > 1\f$ is some constant. The difference between the
`Apollonius_graph_2<Gt,Agds>` class and the
`Apollonius_graph_hierarchy_2` is on how the nearest neighbor location is done. Given a
point \f$ p\f$ the location is done as follows: at the top most level we
find the nearest neighbor of \f$ p\f$ as in the
`Apollonius_graph_2<Gt,Agds>` class. At every subsequent level \f$ i\f$
we use the nearest neighbor found at level \f$ i+1\f$ to find the nearest
neighbor at level \f$ i\f$. This is a variant of the corresponding
hierarchy for points found in \cgalCite{d-iirdt-98}.

The class has two template parameters which have essentially the same
meaning as in the `Apollonius_graph_2<Gt,Agds>` class.

\tparam Gt is the geometric traits class and must be a model of `ApolloniusGraphTraits_2`.

\tparam Agds is the Apollonius graph data structure and must be a model of `ApolloniusGraphDataStructure_2`
whose vertex and face must be models of `ApolloniusGraphHierarchyVertexBase_2` and `TriangulationFaceBase_2`, respectively.
It defaults to:
\code
  CGAL::Triangulation_data_structure_2<
    CGAL::Apollonius_graph_hierarchy_vertex_base_2<CGAL::Apollonius_graph_vertex_base_2<Gt,true> >,
    CGAL::Triangulation_face_base_2<Gt> >
\endcode

\cgalHeading{Heritage}

The `Apollonius_graph_hierarchy_2` class derives publicly from the
`Apollonius_graph_2<Gt,Agds>` class. The interface is
the same with its base class. In the sequel only the methods
overridden are documented.

`Apollonius_graph_hierarchy_2` does not introduce other types than those introduced by
its base class `Apollonius_graph_2<Gt,Agds>`.

\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
*/
template< typename Gt, typename Agds >
class Apollonius_graph_hierarchy_2 : public CGAL::Apollonius_graph_2<Gt,Agds> {
public:

/// \name Creation
/// @{

/*!
Creates an hierarchy of Apollonius graphs using `gt` as
geometric traits.
*/
Apollonius_graph_hierarchy_2(Gt gt=Gt());

/*!
Creates an Apollonius graph hierarchy using
`gt` as geometric traits and inserts all sites in the
range [`first`, `beyond`).
*/
template< class Input_iterator >
Apollonius_graph_hierarchy_2(Input_iterator first, Input_iterator beyond, Gt gt=Gt());

/*!
Copy constructor. All faces, vertices, and inter-level pointers
are duplicated. After the construction, `agh` and `other` refer
to two different Apollonius graph hierarchies: if
`other` is modified, `agh` is not.
*/
Apollonius_graph_hierarchy_2(const Apollonius_graph_hierarchy_2<Gt,Agds>& other);

/*!
Assignment. All faces, vertices and inter-level pointers
are duplicated. After the construction, `agh` and `other` refer
to two different Apollonius graph hierarchies: if
`other` is modified, `agh` is not.
*/
Apollonius_graph_hierarchy_2<Gt,Agds>
operator=(Apollonius_graph_hierarchy_2<Gt,Agds>
other);

/// @}

/// \name Insertion
/// @{

/*!
Inserts the sites in the range
[`first`,`beyond`). The number of sites in the range
[`first`, `beyond`) is returned.
\pre `Input_iterator` must be a model of `InputIterator` and its value type must be `Site_2`.
*/
template< class Input_iterator >
unsigned int insert(Input_iterator first, Input_iterator beyond);

/*!
Inserts the
site `s` in the Apollonius graph hierarchy. If `s`
is visible then the vertex handle of `s` is returned, otherwise
`Vertex_handle(nullptr)` is returned.
*/
Vertex_handle insert(const Site_2& s);

/*!
Inserts `s` in the Apollonius graph hierarchy using the
site associated with `vnear` as
an estimate for the nearest neighbor of the center of `s`.
If `s` is visible then the vertex handle of `s` is
returned, otherwise `Vertex_handle(nullptr)` is returned.
A call to this method is equivalent to `agh.insert(s);` and it has
been added for the sake of conformity with the interface of the
`Apollonius_graph_2<Gt,Agds>` class.
*/
Vertex_handle insert(const Site_2& s, Vertex_handle vnear);

/// @}

/// \name Removal
/// @{

/*!
Removes the site
associated to the vertex handle `v` from the Apollonius
graph hierarchy.
\pre `v` must correspond to a valid finite vertex of the Apollonius graph hierarchy.
*/
void remove(Vertex_handle v);

/// @}

/// \name Nearest Neighbor Location
/// @{

/*!
Finds the nearest neighbor of the point `p`. In other words it
finds the site whose Apollonius cell contains `p`. Ties are broken
arbitrarily and one of the nearest neighbors of `p` is
returned. If there are no visible sites in the Apollonius diagram
`Vertex_handle(nullptr)` is returned.
*/
Vertex_handle nearest_neighbor(const Point_2& p) const;

/*!
Finds the nearest neighbor of the point
`p`. If there are no visible sites in the Apollonius diagram
`Vertex_handle(nullptr)` is returned.
A call to this method is equivalent to
`agh.nearest_neighbor(p);` and it has been added for the sake of
conformity with the interface of the
`Apollonius_graph_2<Gt,Agds>` class.
*/
Vertex_handle nearest_neighbor(const Point_2& p, Vertex_handle vnear) const;

/// @}

/// \name I/O
/// @{

/*!
Writes the current
state of the Apollonius graph hierarchy to an output stream. In particular,
all visible and hidden sites are written as well as the
underlying combinatorial hierarchical data structure.
*/
void file_output(std::ostream& os) const;

/*!
Reads the state of the
Apollonius graph hierarchy from an input stream.
*/
void file_input(std::istream& is);

/*!
Writes the current state of the Apollonius graph hierarchy to an
output stream.
*/
std::ostream& operator<<(std::ostream& os, Apollonius_graph_hierarchy_2<Gt,Agds> agh) const;

/*!
Reads the state of the Apollonius graph hierarchy from an input stream.
*/
std::istream& operator>>(std::istream& is, Apollonius_graph_hierarchy_2<Gt,Agds> agh);

/// @}

/// \name Validity Check
/// @{

/*!
Checks the validity of the Apollonius graph hierarchy. If
`verbose` is `true` a short message is sent to
`std::cerr`. If `level` is 0, the data structure at all levels
is validated, as well as the inter-level pointers. If `level` is
1, then the data structure at all levels is validated, the inter-level
pointers are validated and all levels of the Apollonius graph
hierarchy are also validated. Negative values of `level` always
return `true`, and values greater than 1 are equivalent to
`level` being 1.
*/
bool is_valid(bool verbose = false, int level = 1) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
Clears all contents of the Apollonius graph
hierarchy.
*/
void clear();

/*!
The Apollonius graph hierarchies `other` and `agh` are
swapped. `agh.swap(other)` should be preferred to `agh = other`
or to `agh(other)` if `other` is deleted afterwards.
*/
void swap(Apollonius_graph_hierarchy_2<Gt,Agds>& other);

/// @}

}; /* end Apollonius_graph_hierarchy_2 */
} /* end namespace CGAL */
