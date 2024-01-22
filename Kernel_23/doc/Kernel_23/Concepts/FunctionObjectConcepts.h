/// \cgalConceptNamespace
namespace Kernel {

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `angle_grp`

*/
class Angle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
  on the angle formed by the two vectors `u` and `v`.
  */
  Angle operator()(const Kernel::Vector_2&u,
                   const Kernel::Vector_2&v);

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
    on the angle formed by the three points `p`, `q`, `r` (`q` being the vertex of
    the angle). The returned value is the same as `operator()(p - q, r - q)`.
  */
  Angle operator()(const Kernel::Point_2&p,
                   const Kernel::Point_2&q,
                   const Kernel::Point_2&r);

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
    on the angle formed by the two vectors `pq`, `rs`. The returned value is
    the same as `operator()(q - p, s - r)`.
  */
  Angle operator()(const Kernel::Point_2&p,
                   const Kernel::Point_2&q,
                   const Kernel::Point_2&r,
                   const Kernel::Point_2&s);


  /// @}

}; /* end Kernel::Angle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

\sa `angle_grp`

*/
class Angle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
  on the angle formed by the two vectors `u` and `v`.
  */
  Angle operator()(const Kernel::Vector_3&u,
                   const Kernel::Vector_3&v);

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
    on the angle formed by the three points `p`, `q`, `r` (`q` being the vertex of
    the angle). The returned value is the same as `operator()(p - q, r - q)`.
  */
  Angle operator()(const Kernel::Point_3&p,
                   const Kernel::Point_3&q,
                   const Kernel::Point_3&r);

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
    on the angle formed by the two vectors `pq`, `rs`. The returned value is
    the same as `operator()(q - p, s - r)`.
  */
  Angle operator()(const Kernel::Point_3&p,
                   const Kernel::Point_3&q,
                   const Kernel::Point_3&r,
                   const Kernel::Point_3&s);

  /*!
    returns \ref CGAL::OBTUSE, \ref CGAL::RIGHT or \ref CGAL::ACUTE depending
    on the angle formed by the normal of the triangle `pqr` and the vector `v`.
  */
  Angle operator()(const Kernel::Point_3&p,
                   const Kernel::Point_3&q,
                   const Kernel::Point_3&r,
                   const Kernel::Vector_3&v);
  /// @}

}; /* end Kernel::Angle_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `are_ordered_along_line_grp`

*/
class AreOrderedAlongLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff the three points are collinear and
    `q` lies between `p` and `r`.
    Note that `true` is returned, if `q==p` or
    `q==r`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);

  /// @}

}; /* end Kernel::AreOrderedAlongLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `are_ordered_along_line_grp`

*/
class AreOrderedAlongLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff the three points are collinear and
    `q` lies between `p` and `r`.
    Note that `true` is returned, if `q==p` or
    `q==r`.
    */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);

  /// @}

}; /* end Kernel::AreOrderedAlongLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `parallel_grp`

*/
class AreParallel_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, if `l1` and `l2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Line_2&l1,
                  const Kernel::Line_2&l2);

  /*!
    returns `true`, if `r1` and `r2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Ray_2&r1,
                  const Kernel::Ray_2&r2);

  /*!
    returns `true`, if `s1` and `s2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Segment_2&s1,
                  const Kernel::Segment_2&s2);

  /// @}

}; /* end Kernel::AreParallel_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `parallel_grp`

*/
class AreParallel_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, if `l1` and `l2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Line_3&l1,
                  const Kernel::Line_3&l2);

  /*!
    returns `true`, if `h1` and `h2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Plane_3&h1,
                  const Kernel::Plane_3&h2);

  /*!
    returns `true`, if `r1` and `r2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Ray_3&r1,
                  const Kernel::Ray_3&r2);

  /*!
    returns `true`, if `s1` and `s2` are parallel or if one
    of those (or both) is degenerate.
  */
  bool operator()(const Kernel::Segment_3&s1,
                  const Kernel::Segment_3&s2);


  /// @}

}; /* end Kernel::AreParallel_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `are_strictly_ordered_along_line_grp`

*/
class AreStrictlyOrderedAlongLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns `true`, iff the three points are collinear and
    `q` lies strictly between `p` and `r`.
    Note that `false` is returned, if `q==p` or
    `q==r`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);


  /// @}

}; /* end Kernel::AreStrictlyOrderedAlongLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

\sa `are_strictly_ordered_along_line_grp`

*/
class AreStrictlyOrderedAlongLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff the three points are collinear and
    `q` lies strictly between `p` and `r`.
    Note that `false` is returned, if `q==p` or
    `q==r`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);

  /// @}

}; /* end Kernel::AreStrictlyOrderedAlongLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

 \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred.

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Object`
  \sa `Kernel::Object_2`
  \sa `Kernel::Intersect_2`

*/
class Assign_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    assigns `o` to `t` if `o`
    was constructed from an object of type `T`.
    Returns `true`, if the assignment was possible.
  */
  template <class T>
  bool operator()(T& t, const Kernel::Object_2&o);


  /// @}

}; /* end Kernel::Assign_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

 \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred.
  \sa `CGAL::Object`
  \sa `Kernel::Object_3`
  \sa `Kernel::Intersect_3`

*/
class Assign_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    assigns `o` to `t` if `o`
    was constructed from an object of type `T`.
    Returns `true`, if the assignment was possible.
  */
  template <class T>
  bool operator()(T& t, const Kernel::Object_3&o);

  /// @}

}; /* end Kernel::Assign_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`

*/
class BoundedSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is relative to
    circle `c`.
  */
  Bounded_side operator()(const Kernel::Circle_2&c,
                          const Kernel::Point_2&p);

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is relative to
    triangle `t`.
  */
  Bounded_side operator()(const Kernel::Triangle_2& t,
                          const Kernel::Point_2&p);

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is relative to
    rectangle `r`.
  */
  Bounded_side operator()(const Kernel::Iso_rectangle_2& r,
                          const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::BoundedSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class BoundedSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is with
    respect to sphere `s`.
  */
  Bounded_side operator()(const Kernel::Sphere_3& s,
                          const Kernel::Point_3&p);

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is with
    respect to tetrahedron `t`.
  */
  Bounded_side operator()(const Kernel::Tetrahedron_3& t,
                          const Kernel::Point_3&p);

  /*!
    returns either \ref CGAL::ON_UNBOUNDED_SIDE,
    \ref CGAL::ON_BOUNDED_SIDE, or the constant
    \ref CGAL::ON_BOUNDARY, depending on where point `p` is with
  respect to iso-cuboid `c`.
  */
  Bounded_side operator()(const Kernel::Iso_cuboid_3& c,
                          const Kernel::Point_3&p);


  /// @}

}; /* end Kernel::BoundedSide_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  A type representing an iterator to the %Cartesian coordinates of a point
  in two dimensions.

  \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

  \sa `Kernel::ConstructCartesianConstIterator_2`

*/
class CartesianConstIterator_2 {
public:

}; /* end Kernel::CartesianConstIterator_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  A type representing an iterator to the %Cartesian coordinates of a point
  in three dimensions.

  \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

  \sa `Kernel::ConstructCartesianConstIterator_3`

*/
class CartesianConstIterator_3 {
public:

}; /* end Kernel::CartesianConstIterator_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_are_ordered_along_line_grp`

*/
class CollinearAreOrderedAlongLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff `q` lies between `p`
    and `r`. \pre `p, q` and `r` are collinear.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);


  /// @}

}; /* end Kernel::CollinearAreOrderedAlongLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_are_ordered_along_line_grp`

*/
class CollinearAreOrderedAlongLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff `q` lies between `p`
    and `r`. \pre `p, q` and `r` are collinear.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);

  /// @}

}; /* end Kernel::CollinearAreOrderedAlongLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_are_strictly_ordered_along_line_grp`

*/
class CollinearAreStrictlyOrderedAlongLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff `q` lies strictly between
    `p` and `r`. \pre `p, q` and `r` are collinear.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);


  /// @}

}; /* end Kernel::CollinearAreStrictlyOrderedAlongLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_are_strictly_ordered_along_line_grp`

*/
class CollinearAreStrictlyOrderedAlongLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff `q` lies strictly between
    `p` and `r`. \pre `p, q` and `r` are collinear.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);


  /// @}

}; /* end Kernel::CollinearAreStrictlyOrderedAlongLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class CollinearHasOn_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    checks if point `p` is on `r`.
    \pre `p` is on the supporting line of `r`.
  */
  bool operator()(const Kernel::Ray_2& r,
                  const Kernel::Point_2&p);

  /*!
    checks if point `p` is on `s`.
    \pre `p` is on the supporting line of `s`.
  */
  bool operator()(const Kernel::Segment_2& s,
                  const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::CollinearHasOn_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_grp`

*/
class Collinear_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, if `p`, `q`, and `r` are collinear.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);

  /// @}

}; /* end Kernel::Collinear_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `collinear_grp`

*/
class Collinear_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, if `p`, `q`, and `r` are collinear.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);

  /// @}

}; /* end Kernel::Collinear_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

*/
class CompareAngleWithXAxis_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares
  the angles between the positive \f$ x\f$-axis and the directions in
  counterclockwise order.
  */
  Comparison_result operator()(const
                               Kernel::Direction_2& d, const Kernel::Direction_2& e);


  /// @}

}; /* end Kernel::CompareAngleWithXAxis_2 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

*/
class CompareAngle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
  compares the angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
  \f$ \theta_1\f$ is the angle in \f$ [0, \pi]\f$ of the triangle
  \f$ (a, b, c)\f$ at the vertex `b`, and \f$ \theta_2\f$ is
  the angle in \f$ [0, \pi]\f$ such that \f$ cos(\theta_2) = cosine\f$.
  \pre `a!=b && c!=b`.
  */
  Comparison_result operator()(const K::Point_3& a,
                               const K::Point_3& b,
                               const K::Point_3& c,
                               const K::FT& cosine);
};

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

*/
class CompareDihedralAngle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
  \f$ \theta_1\f$ is the dihedral angle, in \f$ [0, \pi]\f$, of the tetrahedron
  \f$ (a_1, b_1, c_1, d_1)\f$ at the edge `(a_1, b_1)`, and \f$ \theta_2\f$ is
  the angle in \f$ [0, \pi]\f$ such that \f$ cos(\theta_2) = cosine\f$.
  The result is the same as `operator()(b1-a1, c1-a1, d1-a1, cosine)`.
  \pre `a_1`, `b_1`, `c_1` are not collinear, and `a_1`, `b_1`, `d_1` are not collinear.
  */
  Comparison_result operator()(const K::Point_3& a1,
                               const K::Point_3& b1,
                               const K::Point_3& c1,
                               const K::Point_3& d1,
                               const K::FT& cosine);

  /*!
    compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
    \f$ \theta_i\f$ is the dihedral angle in the tetrahedron `(a_i, b_i,
    c_i, d_i)` at the edge `(a_i, b_i)`. These two angles are computed
    in \f$ [0, \pi]\f$.
    The result is the same as `operator()(b1-a1, c1-a1, d1-a1, b2-a2, c2-a2, d2-a2)`.
    \pre For \f$ i \in\{1,2\}\f$, `a_i`, `b_i`, `c_i` are not collinear, and `a_i`, `b_i`, `d_i` are not collinear.
  */
  Comparison_result operator()(const K::Point_3& a1,
                               const K::Point_3& b1,
                               const K::Point_3& c1,
                               const K::Point_3& d1,
                               const K::Point_3& a2,
                               const K::Point_3& b2,
                               const K::Point_3& c2,
                               const K::Point_3& d2);

  /*!
    compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
    \f$ \theta_1\f$ is the dihedral angle, in \f$ [0, \pi]\f$, between the
    vectorial planes defined by `(u_1, v_1)` and `(u_1, w_1)`, and
    \f$ \theta_2\f$ is the angle in \f$ [0, \pi]\f$ such that \f$ cos(\theta_2) =
    cosine\f$.
    \pre `u_1` and `v_1` are not collinear, and `u_1` and `w_1` are not collinear.
  */
  Comparison_result operator()(const K::Vector_3& u1,
                               const K::Vector_3& v1,
                               const K::Vector_3& w1,
                               const K::FT& cosine);

  /*!
    compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
    \f$ \theta_i\f$ is the dihedral angle between the vectorial planes
    defined by `(u_i, v_i)` and `(u_i, w_i)`. These two angles are
    computed in \f$ [0, \pi]\f$.
    \pre For \f$ i \in\{1,2\}\f$, `u_i` and `v_i` are not collinear, and `u_i` and `w_i` are not collinear.
  */
  Comparison_result operator()(const K::Vector_3& u1,
                               const K::Vector_3& v1,
                               const K::Vector_3& w1,
                               const K::Vector_3& u2,
                               const K::Vector_3& v2,
                               const K::Vector_3& w2);

  /// @}

}; /* end Kernel::CompareDihedralAngle_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `Kernel::CompareSquaredDistance_2`
  \sa `compare_distance_to_point_grp`
  \sa `compare_squared_distance_grp`

  \cgalRefines{AdaptableTernaryFunction}

*/
class CompareDistance_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the squared distance between `obj1` and `obj2` to
    the squared distance between `obj1` and `obj3`,
    for all triples of types `Type1`, `Type2` and`Type3`
    in the following set of types:

    - `Kernel::Point_2`
    - `Kernel::Line_2`
    - `Kernel::Ray_2`
    - `Kernel::Segment_2`
    - `Kernel::Triangle_2`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Type3& obj3);

  /*!
    compares the squared distance between `obj1` and `obj2` to
    the squared distance between `obj3` and `obj4`,
    for all tuples of types `Type1`, `Type2`, `Type3`
    and `Type4` in the following set of types:

    - `Kernel::Point_2`
    - `Kernel::Line_2`
    - `Kernel::Ray_2`
    - `Kernel::Segment_2`
    - `Kernel::Triangle_2`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Type3& obj3,
                               const Type4& obj4);

  /// @}

}; /* end Kernel::CompareDistance_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::CompareSquaredDistance_3`
  \sa `compare_distance_to_point_grp`
  \sa `compare_squared_distance_grp`

*/
class CompareDistance_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the squared distance between `obj1` and `obj2` to
    the squared distance between `obj1` and `obj3`, for all triples of types `Type1`, `Type2` and `Type3`
    in the following set of types:
    - `Kernel::Point_3`
    - `Kernel::Line_3`
    - `Kernel::Ray_3`
    - `Kernel::Segment_3`
    - `Kernel::Plane_3`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Type3& obj3);

  /*!
    compares the squared distance between `obj1` and `obj2` to
    the squared distance between `obj3` and `obj4`, for all tuples of types `Type1`, `Type2`, `Type3`
    and `Type4` in the following set of types:
    - `Kernel::Point_3`
    - `Kernel::Line_3`
    - `Kernel::Ray_3`
    - `Kernel::Segment_3`
    - `Kernel::Plane_3`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Type3& obj3,
                               const Type4& obj4);

  /// @}

}; /* end Kernel::CompareDistance_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Weighted_point_2<Kernel>`
  \sa `ComputePowerProduct_2` for the definition of power distance.

*/
class ComparePowerDistance_2
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the power distance between `p` and `q` to the power distance
    between `p` and `r`.
  */
  Comparison_result operator()(const Kernel::Point_2& p,
                               const Kernel::Weighted_point_2& q,
                               const Kernel::Weighted_point_2& r);
  /// @}
};

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Weighted_point_3<Kernel>`
  \sa `ComputePowerProduct_3` for the definition of power distance.
*/
class ComparePowerDistance_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the power distance between `p` and `q` to the power distance
    between `p` and `r`.
  */
  Comparison_result operator()(const Kernel::Point_3& p,
                               const Kernel::Weighted_point_3& q,
                               const Kernel::Weighted_point_3& r);
  /// @}
}; /* end Kernel::ComparePowerDistance_3 */




/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}
*/
class CompareSignedDistanceToLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the signed distance of `r` and `s` to the directed line through `p` and `q`.
  */
  Comparison_result operator()(const Kernel::Point_2& p,
                               const Kernel::Point_2& q,
                               const Kernel::Point_2& r,
                               const Kernel::Point_2& s);

  /*!
    compares the signed distance of `r` and `s` to the directed line `l`.
  */
  Comparison_result operator()(const Kernel::Line_2& l,
                               const Kernel::Point_2& r,
                               const Kernel::Point_2& s);
  /// @}
}; /* end Kernel::CompareSignedDistanceToLine_2 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_slopes_grp`

*/
class CompareSlope_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    compares the slopes of the lines `l1` and `l2`
  */
  Comparison_result operator()(const Kernel::Line_2& l1,
                               const Kernel::Line_2& l2);

  /*!
    compares the slopes of the segments `s1` and `s2`,
    where the slope is the variation of the `y`-coordinate
    from the left to the right endpoint of the segments.
  */
  Comparison_result operator()(const Kernel::Segment_2& s1,
                               const Kernel::Segment_2& s2);

  /*!
    compares the slopes of the segments `(s1s,s1t)` and `(s2s,s2t)`,
    where the slope is the variation of the `y`-coordinate
    from the left to the right endpoint of the segments.
  */
  Comparison_result operator()(const Kernel::Point_2& s1s,
                               const Kernel::Point_2& s1t,
                               const Kernel::Point_2& s2s,
                               const Kernel::Point_2& s2t));

  /// @}

}; /* end Kernel::CompareSlope_2 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_slopes_grp`

*/
class CompareSlope_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    compares the slopes of the segments `(p,q)` and `(r,s)`,
    where the slope is the variation of the `z`-coordinate
    from the first to the second point of the segment divided
    by the length of the segment.
  */
  Comparison_result operator()(const Kernel::Point_3& p,
                               const Kernel::Point_3& q,
                               const Kernel::Point_3& r,
                               const Kernel::Point_3& s);


  /// @}

}; /* end Kernel::CompareSlope_3 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `compare_distance_to_point_grp`
  \sa `compare_squared_distance_grp`

*/
class CompareSquaredDistance_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the squared distance between the two geometrical objects
    `obj1` and `obj2` to the value `d2`, where the types `Type1` and `Type2` can be any of the
    following:

    - `Kernel::Point_2`
    - `Kernel::Line_2`
    - `Kernel::Ray_2`
    - `Kernel::Segment_2`
    - `Kernel::Triangle_2`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Kernel::FT&d2);

  /// @}

}; /* end Kernel::CompareSquaredDistance_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::CompareDistance_3`
  \sa `compare_distance_to_point_grp`
  \sa `compare_squared_distance_grp`

*/
class CompareSquaredDistance_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the squared distance between the two geometrical objects
    `obj1` and `obj2` to the value `d2`, for all pairs `Type1` and `Type2`, where
    the types `Type1` and `Type2` can be any of the following:
    - `Kernel::Point_3`
    - `Kernel::Line_3`
    - `Kernel::Ray_3`
    - `Kernel::Segment_3`
    - `Kernel::Plane_3`
  */
  Comparison_result operator()(const Type1& obj1,
                               const Type2& obj2,
                               const Kernel::FT&d2);

  /// @}

}; /* end Kernel::CompareSquaredDistance_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `compare_squared_radius_grp`

*/
class CompareSquaredRadius_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the squared radius of the sphere of radius 0 centered
    at `p` to `sr`.
    This returns the opposite sign of `sr`.
  */
  Comparison_result operator()(const Kernel::Point_3& p,
                               const Kernel::FT& sr);

  /*!
    compares the squared radius of the sphere defined by the
    points `p` and `q` to `sr`.
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q,
                               const Kernel::FT& sr);

  /*!
    compares the squared radius of the sphere defined by the
    points `p`, `q`, and `r` to `sr`.
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q,
                               const Kernel::Point_3&r,
                               const Kernel::FT& sr);

  /*!
    compares the squared radius of the sphere defined by the
    points `p`, `q`, `r`, and `s` to `sr`.
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q,
                               const Kernel::Point_3&r,
                               const Kernel::Point_3&s,
                               const Kernel::FT& sr);


  /// @}

}; /* end Kernel::CompareSquaredRadius_3 */

/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableFunctor}

\sa `ComputePowerProduct_3` for the definition of of orthogonality for power distances.

*/
class CompareWeightedSquaredRadius_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the weight of the smallest sphere orthogonal to the
    input weighted point(s) with the input weight.
    */
    Comparison_result operator()(const Kernel::Weighted_point_3& pw,
                                 const Kernel::Weighted_point_3& qw,
                                 const Kernel::Weighted_point_3& rw,
                                 const Kernel::Weighted_point_3& sw,
                                 const Kernel::FT& w);

    Comparison_result operator()(const Kernel::Weighted_point_3& pw,
                                 const Kernel::Weighted_point_3& qw,
                                 const Kernel::Weighted_point_3& rw,
                                 const Kernel::FT& w);

    Comparison_result operator()(const Kernel::Weighted_point_3& pw,
                                 const Kernel::Weighted_point_3& qw,
                                 const Kernel::FT& w);

    Comparison_result operator()(const Kernel::Weighted_point_3& pw,
                                 const Kernel::FT& w);

    /// @}
}; /* end Kernel::CompareWeightedSquaredRadius_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \anchor fig-compare_x_at_y_2
  \image html compare_x_at_y.png
  \image latex compare_x_at_y.png

  \cgalRefines{AdaptableTernaryFunction}

  \sa `compare_x_at_y_grp`

*/
class CompareXAtY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the \f$ x\f$-coordinates of `p` and the horizontal projection
    of `p` on `h`. See Figure \ref fig-compare_x_at_y_2 (a).

  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &h);

  /*!
    compares the \f$ x\f$-coordinates of the horizontal projection
    of `p` on `h1` and on `h2`.
    See Figure \ref fig-compare_x_at_y_2 (b).

  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /*!
    Let `p` be the intersection of lines `l1` and `l2`.
    This function compares the \f$ x\f$-coordinates of `p` and
    the horizontal projection of `p` on `h`.

    See Figure \ref fig-compare_x_at_y_2 (c).
  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h);

  /*!
    Let `p` be the intersection of lines `l1` and `l2`. This
    function compares the \f$ x\f$-coordinates of the horizontal projection of
    `p` on `h1` and on `h2`.

    See Figure \ref fig-compare_x_at_y_2 (d)).
  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /// @}

}; /* end Kernel::CompareXAtY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_xyz_grp`

*/
class CompareXYZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian coordinates of points `p` and
    `q` lexicographically in \f$ xyz\f$ order: first
    \f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
    are compared. If they are equal, \f$ z\f$-coordinates are compared.
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::CompareXYZ_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_xy_grp`
  \sa `Kernel::CompareYX_2`

*/
class CompareXY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian coordinates of points `p` and
    `q` lexicographically in \f$ xy\f$ order: first
    \f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
    are compared.
  */
  Comparison_result operator()(const Kernel::Point_2&p,
                               const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::CompareXY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_xy_grp`

*/
class CompareXY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    Compares the %Cartesian coordinates of points `p` and
    `q` lexicographically in \f$ xy\f$ order: first
    \f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
    are compared.
    */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::CompareXY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \anchor fig-compare12
  \image html compare1.png
  \image latex compare1.png

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_x_grp`

*/
class CompareX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the %Cartesian \f$ x\f$-coordinates of points `p` and `q`
  */
  Comparison_result operator()(const Kernel::Point_2&p,
                               const Kernel::Point_2&q);

  /*!
    compares the \f$ x\f$-coordinates of `p` and the intersection
    of lines `l1` and `l2`.
    See Figure \ref fig-compare12 (a).
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2);

  /*!
    compares the `x`-coordinates of the intersection of line `l`
    with line `h1` and with line `h2`.

    See Figure \ref fig-compare12 (b).
  */
  Comparison_result operator()(const Kernel::Line_2 &l,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /*!
    compares the \f$ x\f$-coordinates of the intersection of lines `l1`
    and `l2` and the intersection of lines `h1` and `h2`.

    See Figure \ref fig-compare12 (c).
  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /// @}

}; /* end Kernel::CompareX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_x_grp`

*/
class CompareX_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian \f$ x\f$-coordinates of points `p` and
    `q`
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q);


  /// @}

}; /* end Kernel::CompareX_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \image html compare2.png
  \image latex compare2.png

  \cgalRefines{AdaptableTernaryFunction}

  \anchor fig-compare2
  \sa `compare_y_at_x_grp`

*/
class CompareYAtX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compares the \f$ y\f$-coordinates of `p` and the vertical projection
    of `p` on `h`. See Figure \ref fig-compare2 (e).
    \pre `h` is not vertical.
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &h);

  /*!
    This function compares the `y`-coordinates of the vertical projection
    of `p` on `h1` and on `h2`. See Figure \ref fig-compare2 (e).
    \pre `h1` and `h2` are not vertical.
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /*!
    Let `p` be the intersection of lines `l1` and `l2`.
    This function compares the \f$ y\f$-coordinates of `p` and
    the vertical projection of `p` on `h`.
    See (Figure \ref fig-compare2 (f)).
    \pre `l1`, `l2` intersect and `h` is not vertical.

  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h);

  /*!
    Let `p` be the intersection of lines `l1` and `l2`. This function
    compares the `y`-coordinates of the vertical projection of `p` on
    `h1` and on `h2`.
    See (Figure \ref fig-compare2 (g)).
    \pre `l1` and `l2` intersect; `h1` and `h2` are not vertical.

  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /*!
  compares the `y`-coordinates of `p` and the vertical projection
  of `p` on `s`. If `s` is vertical, then return
  \ref CGAL::EQUAL when `p` lies on `s`, \ref CGAL::SMALLER when `p` lies
  under s, and \ref CGAL::LARGER otherwise.

  \pre `p` is within the x range of `s`.
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Segment_2 &s);

  /*!
    This function compares the \f$ y\f$-coordinates of the vertical projection
    of `p` on `s1` and on `s2`. If `s1` or `s2`
    is vertical, then return \ref CGAL::EQUAL if they intersect, otherwise return
    \ref CGAL::SMALLER if `s1` lies below `s2`, and return \ref CGAL::LARGER
    otherwise.

    \pre `p` is within the x range of `s1` and `s2`.
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Segment_2 &s1,
                               const Kernel::Segment_2 &s2);

  /// @}

}; /* end Kernel::CompareYAtX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_yx_grp`
  \sa `Kernel::CompareXY_2`

*/
class CompareYX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian coordinates of points `p` and
    `q` lexicographically in \f$ yx\f$ order: first
    \f$ y\f$-coordinates are compared, if they are equal, \f$ x\f$-coordinates
    are compared.
  */
  Comparison_result operator()(const Kernel::Point_2&p,
                               const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::CompareYX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \anchor fig-compare14
  \image html compare1.png
  \image latex compare1.png

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_y_grp`

*/
class CompareY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian \f$ y\f$-coordinates of points `p` and
    `q`
  */
  Comparison_result operator()(const Kernel::Point_2&p,
                               const Kernel::Point_2&q);

  /*!
    compares the \f$ y\f$-coordinates of `p` and the
    intersection of lines `l1` and `l2`.

    See Figure \ref fig-compare14 (a).
  */
  Comparison_result operator()(const Kernel::Point_2 &p,
                               const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2);

  /*!
    compares the \f$ y\f$-coordinates of the intersection of line `l`
    with line `h1` and with line `h2`.

    See Figure \ref fig-compare14 (b).
  */
  Comparison_result operator()(const Kernel::Line_2 &l,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);

  /*!
    compares the \f$ y\f$-coordinates of the intersection of lines `l1`
    and `l2` and the intersection of lines `h1` and `h2`.

    See Figure \ref fig-compare14 (c).
  */
  Comparison_result operator()(const Kernel::Line_2 &l1,
                               const Kernel::Line_2 &l2,
                               const Kernel::Line_2 &h1,
                               const Kernel::Line_2 &h2);


  /// @}

}; /* end Kernel::CompareY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_y_grp`

*/
class CompareY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian \f$ y\f$-coordinates of points `p` and
    `q`
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::CompareY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_z_grp`

*/
class CompareZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Compares the %Cartesian \f$ z\f$-coordinates of points `p` and
    `q`
  */
  Comparison_result operator()(const Kernel::Point_3&p,
                               const Kernel::Point_3&q);


  /// @}

}; /* end Kernel::CompareZ_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeA_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `a` of the line with equation `ax +by + c = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_2& l) const;

  /// @}

}; /* end Kernel::ComputeA_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeA_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `a` of the line with equation `ax +by + cz+ d = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_3& l) const;

  /// @}

}; /* end Kernel::ComputeA_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `CGAL::Circle_3<Kernel>`

  \cgalRefines{AdaptableFunctor}

*/
class ComputeApproximateArea_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an approximation of the area of `c`.
  */
  double operator()(const Kernel::Circle_3& c);

  /// @}

}; /* end Kernel::ComputeApproximateArea_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeApproximateAngle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an approximation of the angle between `u` and `v`.
    The angle is given in degrees.
    \pre `u` and `v` are not equal to the null vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& u,
                        const Kernel::Vector_3& v) const;

  /*!
    returns an approximation of the angle between `p-q` and `r-q`.
    The angle is given in degrees.
    \pre `p` and `r` are not equal to `q`.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r) const;

  /// @}

}; /* end Kernel::ComputeApproximateAngle_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeApproximateDihedralAngle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an approximation of the signed dihedral angle in the tetrahedron `pqrs` of edge `pq`.
    The sign is negative if `orientation(p,q,r,s)` is `CGAL::NEGATIVE` and positive otherwise.
    The angle is given in degrees.
    \pre `p,q,r` and `p,q,s` are not collinear.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r,
                        const Kernel::Point_3& s) const;

  /// @}

}; /* end Kernel::ComputeApproximateDihedralAngle_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_3<Kernel>`

*/
class ComputeApproximateSquaredLength_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an approximation of the squared length (i.e.\ perimeter) of `c`.
  */
  double operator()(const Kernel::Circle_3& c);

  /// @}

}; /* end Kernel::ComputeApproximateSquaredLength_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_3<Kernel>`

*/
class ComputeAreaDividedByPi_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the area of `c`, divided by \f$ \pi\f$.
  */
  Kernel::FT operator()(const Kernel::Circle_3& c);

  /// @}

}; /* end Kernel::ComputeAreaDividedByPi_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_rectangle_2<Kernel>`
\sa `CGAL::Triangle_2<Kernel>`

*/
class ComputeArea_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the signed area of the triangle defined by the points `p`,
    `q` and `r`.
  */
  Kernel::FT operator()(const Kernel::Point_2& p,
                        const Kernel::Point_2& q,
                        const Kernel::Point_2& r);

  /*!
    returns the area of `r`.
  */
  Kernel::FT operator()(const Kernel::Iso_rectangle_2& r);

  /*!
    returns the signed area of `t`.
  */
  Kernel::FT operator()(const Kernel::Triangle_2& t);

  /// @}

}; /* end Kernel::ComputeArea_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Triangle_3<Kernel>`

*/
class ComputeArea_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the area of  `t`. This requires that `Kernel::FT`
    supports the `sqrt` operation.
  */
  Kernel::FT operator()(const Kernel::Triangle_3& t);

  /*!
    returns the area of the triangle `p`, `q`, `r`.
    This requires that `Kernel::FT` supports the `sqrt` operation.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r);

  /// @}

}; /* end Kernel::ComputeArea_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeB_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `b` of the line with equation `ax +by + c = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_2& l) const;

  /// @}

}; /* end Kernel::ComputeB_2 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeB_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `b` of the line with equation `ax +by + cz+ d = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_3& l) const;

  /// @}

}; /* end Kernel::ComputeB_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeC_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `c` of the line with equation `ax +by + c = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_2& l) const;

  /// @}

}; /* end Kernel::ComputeC_2 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeC_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `c` of the line with equation `ax +by + cz+ d = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_3& l) const;

  /// @}

}; /* end Kernel::ComputeC_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeD_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the coefficient `d` of the line with equation `ax +by + cz+ d = 0`.
  */
  Kernel::FT operator()(const Kernel::Line_3& l) const;

  /// @}

}; /* end Kernel::ComputeD_3 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`
  \sa `determinant_grp`

*/
class ComputeDeterminant_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the determinant of the two vectors `v` and `w`.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v,
                        const Kernel::Vector_2& w);

  /// @}

}; /* end Kernel::ComputeDeterminant_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Vector_3<Kernel>`
  \sa `determinant_grp`

*/
class ComputeDeterminant_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the determinant of the three vectors `u`, `v` and `w`.
  */
  Kernel::FT operator()(const Kernel::Vector_3& u,
                        const Kernel::Vector_3& v,
                        const Kernel::Vector_3& w);


  /// @}

}; /* end Kernel::ComputeDeterminant_3 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeDx_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns an \f$ x\f$-coordinate of the direction.
  */
  Kernel::FT operator()(const Kernel::Direction_2& d) const;

  /// @}

}; /* end Kernel::ComputeDx_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeDx_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns an \f$ x\f$-coordinate of the direction.
  */
  Kernel::FT operator()(const Kernel::Direction_3& d) const;

  /// @}

}; /* end Kernel::ComputeDx_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeDy_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an \f$ y\f$-coordinate of the direction.
  */
  Kernel::FT operator()(const Kernel::Direction_2& d) const;

  /// @}

}; /* end Kernel::ComputeDy_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeDy_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an \f$ y\f$-coordinate of the direction.
  */
  Kernel::FT operator()(const Kernel::Direction_3& d) const;

  /// @}

}; /* end Kernel::ComputeDy_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeDz_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an \f$ z\f$-coordinate of the direction.
  */
  Kernel::FT operator()(const Kernel::Direction_3& d) const;

  /// @}

}; /* end Kernel::ComputeDz_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHx_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogeneous \f$ x\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_2& p) const;

  /*!
    returns the homogeneous \f$ x\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v) const;

  /// @}

}; /* end Kernel::ComputeHx_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHx_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogeneous \f$ x\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the homogeneous \f$ x\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeHx_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHy_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogeneous \f$ y\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_2& p) const;

  /*!
    returns the homogeneous \f$ y\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v) const;

  /// @}

}; /* end Kernel::ComputeHy_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHy_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogeneous \f$ y\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the homogeneous \f$ y\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeHy_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHw_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogenizing coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_2& p) const;

  /*!
    returns the homogenizing coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v) const;

  /// @}

}; /* end Kernel::ComputeHw_2 */


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHw_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogenizing coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the homogenizing coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeHw_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeHz_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the homogeneous \f$ z\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the homogeneous \f$ z\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeHz_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuinaryFunction}

  \sa `CGAL::Weighted_point_3<Kernel>`
  \sa `ComputePowerProduct_3` for the definitions of power distance and orthogonality.
*/
class ComputePowerDistanceToPowerSphere_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared radius of the sphere centered in `t`
    and orthogonal to the sphere orthogonal to `p`, `q`, `r` ,and `s`.
  */
  Kernel::FT operator()(const Kernel::Weighted_point_3& p,
                        const Kernel::Weighted_point_3& q,
                        const Kernel::Weighted_point_3& r,
                        const Kernel::Weighted_point_3& s,
                        const Kernel::Weighted_point_3& t
                        ) const;
  /// @}
};

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Weighted_point_2<Kernel>`
  \sa `ComputePowerProduct_3`
*/
class ComputePowerProduct_2
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the power product of `pw` and `qw`.
    Let\f$ {p}^{(w)} = (p,w_p), p\in\mathbb{R}^2, w_p\in\mathbb{R}\f$ and
    \f$ {q}^{(w)}=(q,w_q), q\in\mathbb{R}^2, w_q\in\mathbb{R}\f$ be two weighted points.

    The <I>power product</I>, also called <i>power distance</i>
    between \f$ {p}^{(w)}\f$ and \f$ {q}^{(w)}\f$ is defined as
    \f[ \Pi({p}^{(w)},{q}^{(w)}) = {\|{p-q}\|^2-w_p-w_q} \f]
    where \f$ \|{p-q}\|\f$ is the Euclidean distance between \f$ p\f$ and \f$ q\f$.

    The weighted points \f$ {p}^{(w)}\f$ and \f$ {q}^{(w)}\f$
    are said to be <I>orthogonal</I> iff \f$ \Pi{({p}^{(w)},{q}^{(w)})}
    = 0\f$.

    Three weighted points have, in 2D, a unique common orthogonal weighted point
    called the <I>power circle</I>. The <I>power segment</I> will denote the
    weighted point orthogonal to two weighted points on the line defined by
    these two points.
  */
  Kernel::FT operator()(const Kernel::Weighted_point_2& pw,
                        const Kernel::Weighted_point_2& qw) const;

  /// @}
}; /* end Kernel::ComputePowerProduct_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Weighted_point_3<Kernel>`
  \sa `ComputePowerProduct_2`
*/
class ComputePowerProduct_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the power product of `pw` and `qw`.
    Let\f$ {p}^{(w)}=(p,w_p), p\in\mathbb{R}^3, w_p\in\mathbb{R}\f$ and
    \f$ {q}^{(w)}=(q,w_q), q\in\mathbb{R}^3, w_q\in\mathbb{R}\f$ be two weighted points.

    The <I>power product</I>, also called <i>power distance</i>
    between \f$ {p}^{(w)}\f$ and \f$ {q}^{(w)}\f$ is defined as
    \f[ \Pi({p}^{(w)},{q}^{(w)}) = {\|{p-q}\|^2-w_p-w_q} \f]
    where \f$ \|{p-q}\|\f$ is the Euclidean distance between \f$ p\f$ and \f$ q\f$.


    The weighted points \f$ {p}^{(w)}\f$ and \f$ {q}^{(w)}\f$
    are said to be <I>orthogonal</I> iff \f$ \Pi{({p}^{(w)},{q}^{(w)})}
    = 0\f$.

    Four weighted points have, in 3D, a unique common orthogonal weighted point
    called the <I>power sphere</I>. The weighted point orthogonal to
    three weighted points in the plane defined by these three points is
    called the <I>power circle</I>. The
    <I>power segment</I> will denote the weighted point orthogonal to
    two weighted points on the line defined by these two points.
  */
  Kernel::FT operator()(const Kernel::Weighted_point_3& pw,
                        const Kernel::Weighted_point_3& qw) const;

  /// @}
}; /* end Kernel::ComputePowerProduct_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

*/
class ComputeLInfinityDistance_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the distance between the two points in the L-infinity metric.
  */
  Kernel::FT operator()(const Kernel::Point_2& p,
                        const Kernel::Point_2& q) const;


  /// @}

}; /* end Kernel::ComputeLInfinityDistance_2 */
/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeLInfinityDistance_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the distance between the two points in the L-infinity metric.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q) const;


  /// @}

}; /* end Kernel::ComputeLInfinityDistance_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `CGAL::Vector_2<Kernel>`
  \sa scalar_product_grp

  \cgalRefines{AdaptableBinaryFunction}

*/
class ComputeScalarProduct_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the scalar (inner) product of the two vectors `v` and `w`.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v,
                        const Kernel::Vector_2& w);

  /// @}

}; /* end Kernel::ComputeScalarProduct_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`
  \sa scalar_product_grp

*/
class ComputeScalarProduct_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the scalar (inner) product of the two vectors `v` and `w`.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v,
                        const Kernel::Vector_3& w);


  /// @}

}; /* end Kernel::ComputeScalarProduct_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Triangle_3<Kernel>`

*/
class ComputeSquaredArea_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the square of the area of `t`.
  */
  Kernel::FT operator()(const Kernel::Triangle_3& t);

  /*!
    returns the square of the area of the triangle `p`, `q`, `r`.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r);

  /// @}

}; /* end Kernel::ComputeSquaredArea_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `squared_distance_grp`

*/
class ComputeSquaredDistance_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared distance between two geometrical objects of type
    for all pairs `Type1` and `Type2`, where  the types `Type1` and `Type2` can be any of the
    following:

    - `Kernel::Point_2`
    - `Kernel::Line_2`
    - `Kernel::Ray_2`
    - `Kernel::Segment_2`
    - `Kernel::Triangle_2`

    as well as any combination of `Kernel::Point_2` and `Kernel::Weighted_point_2`
  */
  Kernel::FT operator()(Type1 obj1, Type2 obj2);

  /// @}

}; /* end Kernel::ComputeSquaredDistance_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `squared_distance_grp`

*/
class ComputeSquaredDistance_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared distance between two geometrical objects of type
    `Type1` and `Type2`, for all pairs `Type1` and `Type2`, where
    the types `Type1` and `Type2` can be any of the
    following:

    - `Kernel::Point_3`
    - `Kernel::Line_3`
    - `Kernel::Ray_3`
    - `Kernel::Segment_3`
    - `Kernel::Plane_3`

    as well as any combination of `Kernel::Point_3` and `Kernel::Weighted_point_3`
  */
  Kernel::FT operator()(Type1 obj1, Type2 obj2);

  /// @}

}; /* end Kernel::ComputeSquaredDistance_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_3<Kernel>`

*/
class ComputeSquaredLengthDividedByPiSquare_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared length of `c`, divided by \f$ \pi^2\f$.
  */
  Kernel::FT operator()(const Kernel::Circle_3& c);


  /// @}

}; /* end Kernel::ComputeSquaredLengthDividedByPiSquare_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Vector_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class ComputeSquaredLength_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared length of `v`.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v);

  /*!
    returns the squared length of `s`.
  */
  Kernel::FT operator()(const Kernel::Segment_2& s);

  /// @}

}; /* end Kernel::ComputeSquaredLength_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Vector_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`

*/
class ComputeSquaredLength_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared length of `v`.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v);

  /*!
    returns the squared length of `s`.
  */
  Kernel::FT operator()(const Kernel::Segment_3& s);

  /// @}

}; /* end Kernel::ComputeSquaredLength_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `squared_radius_grp`

*/
class ComputeSquaredRadius_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the squared radius of `c`.
  */
  Kernel::FT operator()(const Kernel::Circle_2& c);

  /*!
    returns the squared radius of the circle passing through `p`, `q`
    and `r`. \pre `p, q` and `r` are not collinear.
  */
  Kernel::FT operator()(const Kernel::Point_2& p,
                        const Kernel::Point_2& q,
                        const Kernel::Point_2& r);

  /*!
    returns the squared radius of the smallest circle passing through `p`,
    and `q`, i.e.\ one fourth of the squared distance between `p` and `q`.
  */
  Kernel::FT operator()(const Kernel::Point_2& p,
                        const Kernel::Point_2& q);

  /*!
    returns the squared radius of the smallest circle passing through `p`, i.e.\ \f$ 0\f$.
  */
  Kernel::FT operator()(const Kernel::Point_2& p);


  /// @}

}; /* end Kernel::ComputeSquaredRadius_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Circle_3<Kernel>`
  \sa `squared_radius_grp`

*/
class ComputeSquaredRadius_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the squared radius of `s`.
  */
  Kernel::FT operator()(const Kernel::Sphere_3& s);

  /*!
    returns the squared radius of `c`.
  */
  Kernel::FT operator()(const Kernel::Circle_3& c);

  /*!
    returns the squared radius of the sphere passing through `p`, `q`, `r`
    and `s`. \pre `p, q, r` and `s` are not coplanar.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r,
                        const Kernel::Point_3& s);

  /*!
    returns the squared radius of the sphere passing through `p`, `q` and
  `r`, and whose center is in the plane defined by these three points.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q,
                        const Kernel::Point_3& r);

  /*!
    returns the squared radius of the smallest circle passing through `p`,
    and `q`, i.e.\ one fourth of the squared distance between `p` and `q`.
  */
  Kernel::FT operator()(const Kernel::Point_3& p,
                        const Kernel::Point_3& q);

  /*!
    returns the squared radius of the smallest circle passing through `p`, i.e.\ \f$ 0\f$.
  */
  Kernel::FT operator()(const Kernel::Point_3& p);


  /// @}

}; /* end Kernel::ComputeSquaredRadius_3 */

/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\sa `CGAL::Weighted_point_2<Kernel>`
\sa `ComputePowerProduct_2` for the definition of orthogonality for power distances.

\cgalRefines{AdaptableFunctor}

*/
class ComputeSquaredRadiusSmallestOrthogonalCircle_2
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
  returns the squared radius of the
  smallest sphere circle to the argument(s).
  */
  Kernel::FT operator() (const Kernel::Weighted_point_2& pw,
                         const Kernel::Weighted_point_2& qw,
                         const Kernel::Weighted_point_2& rw) const;

  Kernel::FT operator() (const Kernel::Weighted_point_2& pw,
                         const Kernel::Weighted_point_2& qw) const;

  Kernel::FT operator() (const Kernel::Weighted_point_2& pw) const;

  /// @}
}; /* end Kernel::ComputeSquaredRadiusSmallestOrthogonalCircle_2 */

/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\sa `CGAL::Weighted_point_3<Kernel>`
\sa `ComputePowerProduct_3` for the definition of of orthogonality for power distances.

\cgalRefines{AdaptableFunctor}

*/
class ComputeSquaredRadiusSmallestOrthogonalSphere_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
  returns the squared radius of the
  smallest sphere orthogonal to the argument(s).
  */
  Kernel::FT operator() (const Kernel::Weighted_point_3& pw,
                         const Kernel::Weighted_point_3& qw,
                         const Kernel::Weighted_point_3& rw,
                         const Kernel::Weighted_point_3& sw) const;

  Kernel::FT operator() (const Kernel::Weighted_point_3& pw,
                         const Kernel::Weighted_point_3& qw,
                         const Kernel::Weighted_point_3& rw) const;

  Kernel::FT operator() (const Kernel::Weighted_point_3& pw,
                         const Kernel::Weighted_point_3& qw) const;

  Kernel::FT operator() (const Kernel::Weighted_point_3& pw) const;

  /// @}
}; /* end Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class ComputeVolume_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the volume of `c`.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& c);

  /*!
    returns the signed volume of `t`.
  */
  Kernel::FT operator()(const Kernel::Tetrahedron_3& t);

  /*!
    returns the signed volume of the tetrahedron defined by the four
    points `p0`, `p1`, `p2`, `p3`.
  */
  Kernel::FT operator()(const Kernel::Point_3& p0,
                        const Kernel::Point_3& p1,
                        const Kernel::Point_3& p2,
                        const Kernel::Point_3& p3);


  /// @}

}; /* end Kernel::ComputeVolume_3 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeWeight_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the weight of the weighted point.
  */
  Kernel::FT operator()(const Kernel::WeightedPoint_2& p) const;
  /// @}

}; /* end Kernel::ComputeWeight_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeWeight_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the weight of the weighted point.
  */
  Kernel::FT operator()(const Kernel::WeightedPoint_3& p) const;

  /// @}

}; /* end Kernel::ComputeWeight_3 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ x\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_2& p) const;

  /*!
    returns the \f$ x\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v) const;


  /// @}

}; /* end Kernel::ComputeX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeX_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ x\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the \f$ x\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeX_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeXmax_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the largest \f$ x\f$-coordinate of the iso-rectangle.
  */
  Kernel::FT operator()(const Kernel::Iso_rectangle_2& r) const;

  /// @}

}; /* end Kernel::ComputeXmax_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeXmax_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the largest \f$ x\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;

  /// @}

}; /* end Kernel::ComputeXmax_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeXmin_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the smallest \f$ x\f$-coordinate of the iso-rectangle.
  */
  Kernel::FT operator()(const Kernel::Iso_rectangle_2& r) const;

  /// @}

}; /* end Kernel::ComputeXmin_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeXmin_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the smallest \f$ x\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;

  /// @}

}; /* end Kernel::ComputeXmin_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `compare_y_at_x_grp`

  \cgalRefines{AdaptableFunctor}

*/
class ComputeYAtX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ y\f$-coordinate of the point at `l` with
    given \f$ x\f$-coordinate.
    \pre `l` is not vertical.
  */
  Kernel::FT operator()(const Kernel::Line_2& l,
                        const Kernel::FT &x) const;

  /// @}

}; /* end Kernel::ComputeYAtX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ y\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_2& p) const;

  /*!
    returns the \f$ y\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_2& v) const;

  /// @}

}; /* end Kernel::ComputeY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ y\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the \f$ y\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeYmax_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the largest \f$ y\f$-coordinate of the iso-rectangle.
  */
  Kernel::FT operator()(const Kernel::Iso_rectangle_2& r) const;

  /// @}

}; /* end Kernel::ComputeYmax_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeYmax_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the largest \f$ y\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;

  /// @}

}; /* end Kernel::ComputeYmax_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeYmin_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the smallest \f$ y\f$-coordinate of the iso-rectangle.
  */
  Kernel::FT operator()(const Kernel::Iso_rectangle_2& r) const;

  /// @}

}; /* end Kernel::ComputeYmin_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeYmin_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the smallest \f$ y\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;

  /// @}

}; /* end Kernel::ComputeYmin_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the \f$ z\f$-coordinate of the point.
  */
  Kernel::FT operator()(const Kernel::Point_3& p) const;

  /*!
    returns the \f$ z\f$-coordinate of the vector.
  */
  Kernel::FT operator()(const Kernel::Vector_3& v) const;

  /// @}

}; /* end Kernel::ComputeZ_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeZmax_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the largest \f$ z\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;
  /// @}
}; /* end Kernel::ComputeZmax_3 */
/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

*/
class ComputeZmin_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the smallest \f$ z\f$-coordinate of the iso-cuboid.
  */
  Kernel::FT operator()(const Kernel::Iso_cuboid_3& r) const;
  /// @}
}; /* end Kernel::ComputeZmin_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::ConstructCentroid_2`
  \sa `centroid_grp`
  \sa `barycenter_grp`

*/
class ConstructBarycenter_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    compute the barycenter of the points `p1` and `p2` with corresponding
    weights `w1` and `1-w1`.
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2);

  /*!
    compute the barycenter of the points `p1` and `p2` with corresponding
    weights `w1` and `w2`. \pre `w1+w2 != 0`.
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2, const Kernel::FT&w2);

  /*!
    compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
    weights `w1`, `w2` and `1-w1-w2`.
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2, const Kernel::FT&w2,
              const Kernel::Point_2& p3);

  /*!
    compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
    weights `w1`, `w2` and `w3`. \pre `w1+w2+w3 != 0`.
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2, const Kernel::FT&w2,
              const Kernel::Point_2& p3, const Kernel::FT&w3);

  /*!
    compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
    weights `w1`, `w2`, `w3` and `1-w1-w2-w3`.
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2, const Kernel::FT&w2,
              const Kernel::Point_2& p3, const Kernel::FT&w3,
              const Kernel::Point_2& p4);

  /*!
    compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
    weights `w1`, `w2`, `w3` and `w4`. \pre `1+w2+w3+w4 != 0.`
  */
  Kernel::Point_2
  operator()( const Kernel::Point_2& p1, const Kernel::FT&w1,
              const Kernel::Point_2& p2, const Kernel::FT&w2,
              const Kernel::Point_2& p3, const Kernel::FT&w3,
              const Kernel::Point_2& p4, const Kernel::FT&w4);

  /// @}

}; /* end Kernel::ConstructBarycenter_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::ConstructCentroid_3`
  \sa `centroid_grp`
  \sa `barycenter_grp`

*/
class ConstructBarycenter_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    compute the barycenter of the points `p1` and `p2` with corresponding
    weights `w1` and `1-w1`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2);

  /*!
    compute the barycenter of the points `p1` and `p2` with corresponding
    weights `w1` and `w2`. \pre `w1+w2 != 0`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2, const Kernel::FT&w2);

  /*!
    compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
    weights `w1`, `w2` and `1-w1-w2`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2, const Kernel::FT&w2,
              const Kernel::Point_3& p3);

  /*!
    compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
    weights `w1`, `w2` and `w3`. \pre `w1+w2+w3 != 0`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2, const Kernel::FT&w2,
              const Kernel::Point_3& p3, const Kernel::FT&w3);

  /*!
    compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
    weights `w1`, `w2`, `w3` and `1-w1-w2-w3`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2, const Kernel::FT&w2,
              const Kernel::Point_3& p3, const Kernel::FT&w3,
              const Kernel::Point_3& p4);

  /*!
    compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
    weights `w1`, `w2`, `w3` and `w4`. \pre `w1+w2+w3+w4 != 0`.
  */
  Kernel::Point_3
  operator()( const Kernel::Point_3& p1, const Kernel::FT&w1,
              const Kernel::Point_3& p2, const Kernel::FT&w2,
              const Kernel::Point_3& p3, const Kernel::FT&w3,
              const Kernel::Point_3& p4, const Kernel::FT&w4);


  /// @}

}; /* end Kernel::ConstructBarycenter_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructBaseVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    when `index` == 1, returns a vector `b1` that is orthogonal to the
    normal `n` to plane `h`; when `index` == 2, returns a vector
    `b2` that is orthogonal to `n` and `b1` and such that
    for an arbitrary point `p` on the plane `h`, the orientation of
    `p`, `p + b1`, `p + b2`, and `p + n` is positive.
  */
  Kernel::Vector_3 operator()(const Kernel::Plane_3& h,
                              int index);

  /// @}

}; /* end Kernel::ConstructBaseVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

*/
class ConstructBbox_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns a bounding box of `p`.
  */
  CGAL::Bbox_2 operator()(const Kernel::Point_2
                          &p);

  /*!
    returns a bounding box of `s`.
  */
  CGAL::Bbox_2 operator()(const Kernel::Segment_2
                          &s);

  /*!
    returns a bounding box of `t`.
  */
  CGAL::Bbox_2 operator()(const Kernel::Triangle_2
                          &t);

  /*!
    returns a bounding box of `i`.
  */
  CGAL::Bbox_2 operator()(const Kernel::Iso_rectangle_2
                          &i);

  /*!
    returns a bounding box of `c`.
  */
  CGAL::Bbox_2 operator()(const Kernel::Circle_2
                          &c);

  /// @}

}; /* end Kernel::ConstructBbox_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

*/
class ConstructBbox_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns a bounding box of `c`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Circle_3
                          &c);

  /*!
    returns a bounding box of `p`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Point_3
                          &p);

  /*!
    returns a bounding box of `s`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Segment_3
                          &s);

  /*!
    returns a bounding box of `t`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Triangle_3
                          &t);

  /*!
    returns a bounding box of `t`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Tetrahedron_3
                          &t);

  /*!
    returns a bounding box of `i`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Iso_Cuboid_3
                          &i);

  /*!
    returns a bounding box of `s`.
  */
  CGAL::Bbox_3 operator()(const Kernel::Sphere_3
                          &s);

  /// @}

}; /* end Kernel::ConstructBbox_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableFunctor}

  \sa `bisector_grp`

*/
class ConstructBisector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs the bisector of `p` and `q`.
    The bisector is oriented in such a way that `p` lies on its
    positive side. \pre `p != q`.
  */
  Kernel::Line_2 operator()(const Kernel::Point_2&p,
                            const Kernel::Point_2&q );

  /*!
    constructs the bisector of the two lines `l1` and `l2`.
    In the general case, the bisector has the direction of the vector which
    is the sum of the normalized directions of the two lines, and which passes
    through the intersection of `l1` and `l2`.
    If `l1` and `l2` are parallel, then the bisector is defined as the line
    which has the same direction as `l1`, and which is at the same distance
    from `l1` and `l2`.
  */
  Kernel::Line_2 operator()(const Kernel::Line_2&l1,
                            const Kernel::Line_2&l2);

  /// @}

}; /* end Kernel::ConstructBisector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `bisector_grp`

*/
class ConstructBisector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    constructs the bisector plane of `p` and `q`.
    The bisector is oriented in such a way that `p` lies on its
    positive side. \pre `p != q`.
  */
  Kernel::Plane_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q );

  /*!
    constructs the bisector of the two planes `h1` and `h2`.
    In the general case, the bisector has a normal vector which has the same
    direction as the sum of the normalized normal vectors of the two planes, and
    passes through the intersection of `h1` and `h2`.
    If `h1` and `h2` are parallel, then the bisector is defined as the
    plane which has the same oriented normal vector as `h1`, and which is at
    the same distance from `h1` and `h2`.
  */
  Kernel::Plane_3 operator()(const Kernel::Plane_3&h1,
                             const Kernel::Plane_3&h2);


  /// @}

}; /* end Kernel::ConstructBisector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `Kernel::CartesianConstIterator_2`

*/
class ConstructCartesianConstIterator_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns an iterator on the zeroth %Cartesian coordinate of `p`.
  */
  Kernel::Cartesian_const_iterator_2 operator()(const Kernel::Point_2
                                                &p);

  /*!
    returns the past the end iterator of the %Cartesian coordinates of `p`.
  */
  Kernel::Cartesian_const_iterator_2 operator()(const Kernel::Point_2
                                                &p, int);

  /*!
    returns an iterator on the zeroth %Cartesian coordinate of `v`.
  */
  Kernel::Cartesian_const_iterator_2 operator()(const Kernel::Vector_2
                                                &v);

  /*!
    returns the past the end iterator of the %Cartesian coordinates of `v`.
  */
  Kernel::Cartesian_const_iterator_2 operator()(const Kernel::Vector_2
                                                &v, int);


  /// @}

}; /* end Kernel::ConstructCartesianConstIterator_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `Kernel::CartesianConstIterator_3`

*/
class ConstructCartesianConstIterator_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an iterator on the zeroth %Cartesian coordinate of `p`.
  */
  Kernel::Cartesian_const_iterator_3 operator()(const Kernel::Point_3
                                                &p);

  /*!
    returns the past the end iterator of the %Cartesian coordinates of `p`.
  */
  Kernel::Cartesian_const_iterator_3 operator()(const Kernel::Point_3
                                                &p, int);

  /*!
    returns an iterator on the zeroth %Cartesian coordinate of `v`.
  */
  Kernel::Cartesian_const_iterator_3 operator()(const Kernel::Vector_3
                                                &v);

  /*!
    returns the past the end iterator of the %Cartesian coordinates of `v`.
  */
  Kernel::Cartesian_const_iterator_3 operator()(const Kernel::Vector_3
                                                &v, int);

  /// @}

}; /* end Kernel::ConstructCartesianConstIterator_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

\sa `CGAL::Circle_2<Kernel>`

*/
class ConstructCenter_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the center of the circle `c`.
  */
  Kernel::Point_2 operator()(const Kernel::Circle_2 & c);

  /// @}

}; /* end Kernel::ConstructCenter_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Circle_3<Kernel>`

*/
class ConstructCenter_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the center of the sphere `s`.
  */
  Kernel::Point_3 operator()(const Kernel::Sphere_3 & s);

  /*!
    compute the center of the circle `c`.
  */
  Kernel::Point_3 operator()(const Kernel::Circle_3 & c);


  /// @}

}; /* end Kernel::ConstructCenter_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::ConstructBarycenter_2`
  \sa `centroid_grp`
  \sa `barycenter_grp`

*/
class ConstructCentroid_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the centroid of the points `p`, `q`, and `r`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2&p,
                             const Kernel::Point_2&q,
                             const Kernel::Point_2&r);

  /*!
    compute the centroid of the points `p`, `q`, `r` and `s`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2&p,
                             const Kernel::Point_2&q,
                             const Kernel::Point_2&r,
                             const Kernel::Point_2&s);

  /*!
    compute the centroid of the triangle `t`.
  */
  Kernel::Point_2 operator()(const Kernel::Triangle_2&t);

  /// @}

}; /* end Kernel::ConstructCentroid_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `Kernel::ConstructBarycenter_3`
  \sa `centroid_grp`
  \sa `barycenter_grp`

*/
class ConstructCentroid_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the centroid of the points `p`, `q`, and `r`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q,
                             const Kernel::Point_3&r);

  /*!
    compute the centroid of the points `p`, `q`, `r` and `s`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q,
                             const Kernel::Point_3&r,
                             const Kernel::Point_3&s);

  /*!
    compute the centroid of the triangle `t`.
  */
  Kernel::Point_3 operator()(const Kernel::Triangle_3&t);

  /*!
    compute the centroid of the tetrahedron `t`.
  */
  Kernel::Point_3 operator()(const Kernel::Tetrahedron_3&t);

  /// @}

}; /* end Kernel::ConstructCentroid_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Circle_2<Kernel>`

*/
class ConstructCircle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!

    introduces a variable of type `Kernel::Circle_2`.
    It is initialized to the circle with center `center`,
    squared radius `squared_radius` and orientation
    `orientation`.
    \pre `orientation != CGAL::COLLINEAR` and `squared_radius >= 0`.
  */
  Kernel::Circle_2 operator()( Kernel::Point_2 const& center,
                               Kernel::FT const& squared_radius,
                               Orientation const& orientation
                               = COUNTERCLOCKWISE);

  /*!

    introduces a variable of type `Kernel::Circle_2`.
    It is initialized to the unique circle which passes through
    the points `p`, `q` and `r`. The orientation of
    the circle is the orientation of the point triple `p`,
    `q`, `r`.
    \pre `p`, `q`, and `r` are not collinear.
  */
  Kernel::Circle_2 operator()( Kernel::Point_2 const& p,
                               Kernel::Point_2 const& q,
                               Kernel::Point_2 const& r);

  /*!

    introduces a variable of type `Kernel::Circle_2`.
    It is initialized to the circle with diameter `pq`
    and orientation `orientation`.
    \pre `orientation != CGAL::COLLINEAR`.
  */
  Kernel::Circle_2 operator()( Kernel::Point_2 const& p,
                               Kernel::Point_2 const& q,
                               Orientation const& orientation
                               = COUNTERCLOCKWISE);

  /*!

    introduces a variable of type `Kernel::Circle_2`.
    It is initialized to the circle with center `center`, squared
    radius zero and orientation `orientation`.
    \pre `orientation != CGAL::COLLINEAR`.
    \post .`is_degenerate()` = `true`.
  */
  Kernel::Circle_2 operator()( Kernel::Point_2 const& center,
                               Orientation const& orientation
                               = COUNTERCLOCKWISE);


  /// @}

}; /* end Kernel::ConstructCircle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Circle_3<Kernel>`

*/
class ConstructCircle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a variable of type `Kernel::Circle_3`.
    It is initialized to the circle with center `center`,
    and squared radius `sq_r` in the plane `plane`.
    \pre `center` lies in `plane` and  `sq_r >= 0`.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Point_3 const& center,
    Kernel::FT const& sq_r,
    Kernel::Plane_3 const& plane);

  /*!
    introduces a variable of type `Kernel::Circle_3`.
    It is initialized to the circle with center `center`,
    and squared radius `sq_r` in the plane
    containing `center` and normal to `n`.
    \pre `sq_r >= 0`.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Point_3 const& center,
    Kernel::FT const& sq_r,
    Kernel::Vector_3 const& n);

  /*!
    introduces a variable of type `Kernel::Point_3`.
    It is initialized to the circle passing through the three points.
    \pre The three points are not collinear.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Point_3 const& p,
    Kernel::Point_3 const& q,
    Kernel::Point_3 const& r);

  /*!
    introduces a variable of type `Kernel::Circle_3`.
    It is initialized to the circle along which the two spheres intersect.
    \pre The two spheres intersect along a circle.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Sphere_3 const& sphere1,
    Kernel::Sphere_3 const& sphere2);

  /*!
    introduces a variable of type `Kernel::Circle_3`.
    It is initialized to the circle along which the sphere and the
    plane intersect.
    \pre The sphere and the plane intersect along a circle.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Sphere_3 const& sphere,
    Kernel::Plane_3 const& plane);

  /*!
    introduces a variable of type `Kernel::Circle_3`.
    It is initialized to the circle along which the sphere and the
    plane intersect.
    \pre The sphere and the plane intersect along a circle.
  */
  Kernel::Circle_3 operator()
  ( Kernel::Plane_3 const& plane,
    Kernel::Sphere_3 const& sphere);

  /// @}

}; /* end Kernel::ConstructCircle_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `circumcenter_grp`

*/
class ConstructCircumcenter_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the center of the smallest circle passing through the points `p` and `q`. Note : this is the same as `Kernel::ConstructMidpoint_2`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2&p,
                             const Kernel::Point_2&q);

  /*!
    compute the center of the circle passing through the points `p`, `q`, and `r`.
    \pre `p`, `q`, and `r` are not collinear.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2&p,
                             const Kernel::Point_2&q,
                             const Kernel::Point_2&r);

  /*!
    compute the center of the circle passing through the three vertices of `t`.
    \pre `t` is not degenerate.
  */
  Kernel::Point_2 operator()(const Kernel::Triangle_2&t);

  ///@}

}; /* end Kernel::ConstructCircumcenter_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `circumcenter_grp`

*/
class ConstructCircumcenter_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    compute the center of the smallest circle passing through the points `p` and
    `q`. Note : this is the same as `Kernel::ConstructMidpoint_3`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q);

  /*!
    compute the center of the sphere passing through the points `p`, `q`, `r`,
    and `s`. \pre `p`, `q`, `r`, and `s` are not coplanar.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q,
                             const Kernel::Point_3&r,
                             const Kernel::Point_3&s);

  /*!
    compute the center of the sphere passing through the vertices of `t`.
    \pre `t` is not degenerate.
  */
  Kernel::Point_3 operator()(const Kernel::Tetrahedron_3&t);

  /*!
    compute the center of the circle passing through the points `p`, `q` and `r`.
    \pre `p`, `q` and `r` are not collinear.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3&p,
                             const Kernel::Point_3&q,
                             const Kernel::Point_3&r);

  /*!
    compute the center of the circle passing through the vertices of `t`.
    \pre `t` is not degenerate.
  */
  Kernel::Point_3 operator()(const Kernel::Triangle_3&t);

  /// @}

}; /* end Kernel::ConstructCircumcenter_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `cross_product_grp`
  \sa `determinant_grp`

*/
class ConstructCrossProductVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the cross product of `v` and `w`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v,
                              const Kernel::Vector_3 &w);

  /// @}

}; /* end Kernel::ConstructCrossProductVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`

*/
class ConstructDifferenceOfVectors_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `v1 - v2`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2 &v1,
                              const Kernel::Vector_2 &v2);

  /// @}

}; /* end Kernel::ConstructDifferenceOfVectors_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`

*/
class ConstructDifferenceOfVectors_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `v1 - v2`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v1,
                              const Kernel::Vector_3 &v2);


  /// @}

}; /* end Kernel::ConstructDifferenceOfVectors_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Direction_2<Kernel>`

*/
class ConstructDirection_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the direction of vector `v`.
  */
  Kernel::Direction_2 operator()(const Kernel::Vector_2 &v);

  /*!
    introduces the direction of line `l`.
  */
  Kernel::Direction_2 operator()(const Kernel::Line_2 &l);

  /*!
    introduces the direction of ray `r`.
  */
  Kernel::Direction_2 operator()(const Kernel::Ray_2 &r);

  /*!
    introduces the direction of segment `s`.
  */
  Kernel::Direction_2 operator()(const Kernel::Segment_2 &s);


  /// @}

}; /* end Kernel::ConstructDirection_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Direction_3<Kernel>`

*/
class ConstructDirection_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a direction initialized with the
    direction of vector `v`.
  */
  Kernel::Direction_3 operator()(const Kernel::Vector_3 &v);

  /*!
    introduces the direction of line `l`.
  */
  Kernel::Direction_3 operator()(const Kernel::Line_3 &l);

  /*!
    introduces the direction of ray `r`.
  */
  Kernel::Direction_3 operator()(const Kernel::Ray_3 &r);

  /*!
    introduces the direction of segment `s`.
  */
  Kernel::Direction_3 operator()(const Kernel::Segment_3 &s);

  /// @}

}; /* end Kernel::ConstructDirection_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`

*/
class ConstructDividedVector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `v/s`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2 &v,
                              const Kernel::RT s);

  /// @}

}; /* end Kernel::ConstructDividedVector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`

*/
class ConstructDividedVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `v/s`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v,
                              const Kernel::RT s);

  /// @}

}; /* end Kernel::ConstructDividedVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `equidistant_line_grp`

*/
class ConstructEquidistantLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs the line which is at the same distance from the three points
    `p`, `q` and `r`.
    \pre `p`, `q` and `r` are not collinear.
  */
  Kernel::Line_3 operator()(const Kernel::Point_3&p,
                            const Kernel::Point_3&q,
                            const Kernel::Point_3&r );

  /// @}

}; /* end Kernel::ConstructEquidistantLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableBinaryFunction}

\sa `CGAL::Iso_cuboid_3<Kernel>`

*/
class ConstructIsoCuboid_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces an iso-oriented cuboid with diagonal
    opposite vertices `p` and `q` such that `p` is the
    lexicographically smallest point in the cuboid.
  */
  Kernel::Iso_cuboid_3 operator()(const Kernel::Point_3 &p,
                                  const Kernel::Point_3 &q);

  /*!
    introduces an iso-oriented cuboid with diagonal
    opposite vertices `p` and `q`. The `int` argument value is
    only used to distinguish the two overloaded functions.
    \pre `p.x()<=q.x()`, `p.y()<=q.y()` and `p.z()<=q.z()`.
  */
  Kernel::Iso_cuboid_3 operator()(const Kernel::Point_3 &p,
                                  const Kernel::Point_3 &q,
                                  int);

  /*!
    introduces an iso-oriented cuboid `fo` whose
    minimal \f$ x\f$ coordinate is the one of `left`, the
    maximal \f$ x\f$ coordinate is the one of `right`, the
    minimal \f$ y\f$ coordinate is the one of `bottom`, the
    maximal \f$ y\f$ coordinate is the one of `top`, the
    minimal \f$ z\f$ coordinate is the one of `far`, the
    maximal \f$ z\f$ coordinate is the one of `close`.
  */
  Kernel::Iso_cuboid_3
  operator()(const Kernel::Point_3 &left,
             const Kernel::Point_3 &right,
             const Kernel::Point_3 &bottom,
             const Kernel::Point_3 &top,
             const Kernel::Point_3 &far,
             const Kernel::Point_3 &close);

  /// @}

}; /* end Kernel::ConstructIsoCuboid_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_rectangle_2<Kernel>`

*/
class ConstructIsoRectangle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    introduces an iso-oriented rectangle with diagonal
    opposite vertices `p` and `q` such that `p` is the
    lexicographically smallest point in the rectangle.
  */
  Kernel::Iso_rectangle_2 operator()(const Kernel::Point_2 &p,
                                     const Kernel::Point_2 &q);

  /*!
    introduces an iso-oriented rectangle with diagonal
    opposite vertices `p` and `q`. The `int` argument value is
  only used to distinguish the two overloaded functions.
  \pre `p.x()<=q.x()` and `p.y()<=q.y()`.
  */
  Kernel::Iso_rectangle_2 operator()(const Kernel::Point_2 &p,
                                     const Kernel::Point_2 &q,
                                     int);

  /*!
    introduces an iso-oriented rectangle `fo` whose
    minimal \f$ x\f$ coordinate is the one of `left`, the
    maximal \f$ x\f$ coordinate is the one of `right`, the
    minimal \f$ y\f$ coordinate is the one of `bottom`, the
    maximal \f$ y\f$ coordinate is the one of `top`.
  */
  Kernel::Iso_rectangle_2
  operator()(const Kernel::Point_2 &left,
             const Kernel::Point_2 &right,
             const Kernel::Point_2 &bottom,
             const Kernel::Point_2 &top);


  /// @}

}; /* end Kernel::ConstructIsoRectangle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `ConstructProjectedXYPoint_2`
*/
class ConstructLiftedPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the image point of the projection of `p`
    under an affine transformation which maps the \f$ xy\f$-plane onto `h`.
    This affine transformation must be the inverse of the affine transformation used
    in `ConstructProjectedXYPoint_2`.
  */
  Kernel::Point_3 operator()(const Kernel::Plane_3& h,
                             const Kernel::Point_2& p);

  /// @}

}; /* end Kernel::ConstructLiftedPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_2<Kernel>`

*/
class ConstructLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a line passing through the points `p` and `q`.
    Line is directed from `p` to `q`.
  */
  Kernel::Line_2 operator()(const Kernel::Point_2 &p,
                            const Kernel::Point_2 &q);

  /*!
    introduces a line passing through point `p` with
    direction `d`.
  */
  Kernel::Line_2 operator()(const Kernel::Point_2 &p,
                            const Kernel::Direction_2&d);

  /*!
    introduces a line passing through point `p` and
    oriented by `v`.
  */
  Kernel::Line_2 operator()(const Kernel::Point_2 &p,
                            const Kernel::Vector_2&v);

  /*!
    introduces a line supporting the segment `s`,
    oriented from source to target.
  */
  Kernel::Line_2 operator()(const Kernel::Segment_2 &s);

  /*!
    introduces a line supporting the ray `r`,
    with same orientation.
  */
  Kernel::Line_2 operator()(const Kernel::Ray_2 &r);

  /// @}

}; /* end Kernel::ConstructLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_3<Kernel>`

*/
class ConstructLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a line passing through the points `p` and `q`.
    Line is directed from `p` to `q`.
  */
  Kernel::Line_3 operator()(const Kernel::Point_3 &p,
                            const Kernel::Point_3 &q);

  /*!
    introduces a line passing through point `p` and
    oriented by `v`.
  */
  Kernel::Line_3 operator()(const Kernel::Point_3 &p,
                            const Kernel::Vector_3&v);

  /*!
    introduces a line passing through point `p` with
    direction `d`.
  */
  Kernel::Line_3 operator()(const Kernel::Point_3 &p,
                            const Kernel::Direction_3&d);

  /*!
    returns the line supporting the segment `s`,
    oriented from source to target.
  */
  Kernel::Line_3 operator()(const Kernel::Segment_3 &s);

  /*!
    returns the line supporting the ray `r`, with the
    same orientation.
  */
  Kernel::Line_3 operator()(const Kernel::Ray_3 &r);


  /// @}

}; /* end Kernel::ConstructLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructMaxVertex_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the vertex of
    `r` with lexicographically largest coordinates.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Iso_rectangle_2 &r);

  /*!
    returns the vertex of
    `s` with lexicographically largest coordinates.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Segment_2 &s);

  /// @}

}; /* end Kernel::ConstructMaxVertex_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructMaxVertex_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the vertex of
    `c` with lexicographically largest coordinates.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Iso_cuboid_3 &c);

  /*!
    returns the vertex of
    `s` with lexicographically largest coordinates.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Segment_3 &s);


  /// @}

}; /* end Kernel::ConstructMaxVertex_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `midpoint_grp`

*/
class ConstructMidpoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the midpoint of the segment `pq`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2& p,
                             const Kernel::Point_2& q );
  /*!
    computes the midpoint of the segment `s`.
  */
  Kernel::Point_2 operator()(const Kernel::Segment_2& s);

  /// @}

}; /* end Kernel::ConstructMidpoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `midpoint_grp`

*/
class ConstructMidpoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the midpoint of the segment `pq`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3& p,
                             const Kernel::Point_3& q );

  /*!
    computes the midpoint of the segment `s`.
  */
  Kernel::Point_3 operator()(const Kernel::Segment_3& s);


  /// @}

}; /* end Kernel::ConstructMidpoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructMinVertex_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the vertex of
    `r` with lexicographically smallest coordinates.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Iso_rectangle_2 &r);

  /*!
    returns the vertex of
    `s` with lexicographically smallest coordinates.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Segment_2 &s);

  /// @}

}; /* end Kernel::ConstructMinVertex_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructMinVertex_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the vertex of
    `c` with lexicographically smallest coordinates.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Iso_cuboid_3 &c);

  /*!
    returns the vertex of
    `s` with lexicographically smallest coordinates.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Segment_3 &s);

  /// @}

}; /* end Kernel::ConstructMinVertex_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `normal_grp`
  \sa `unit_normal_grp`

*/
class ConstructNormal_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    computes the normal of the vectors `q-p` and `r-p`.
  */
  Kernel::Vector_3 operator()(const Kernel::Point_3& p,
                              const Kernel::Point_3& q,
                              const Kernel::Point_3& r );


  /// @}

}; /* end Kernel::ConstructNormal_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred.

  \sa `CGAL::Object`
  \sa `Kernel::Assign_2`
  \sa `Kernel::Assign_3`
  \sa `Kernel::Object_2`
  \sa `Kernel::Object_3`

*/
class ConstructObject_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs an object that contains `t` and returns it.
  */
  template <class T>
  Object_2 operator()(const T& t);

  /// @}

}; /* end Kernel::ConstructObject_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

 \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred.
  \sa `CGAL::Object`
  \sa `Kernel::Assign_2`
  \sa `Kernel::Assign_3`
  \sa `Kernel::Object_2`
  \sa `Kernel::Object_3`

*/
class ConstructObject_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs an object that contains `t` and returns it.
  */
  template <class T>
  Object_3 operator()(const T& t);


  /// @}

}; /* end Kernel::ConstructObject_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_2<Kernel>`

*/
class ConstructOppositeCircle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the circle with the same center and squared radius as
    `c`, but with opposite orientation.
  */
  Kernel::Circle_2 operator()(const Kernel::Circle_2& c);

  /// @}

}; /* end Kernel::ConstructOppositeCircle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

\sa `CGAL::Direction_2<Kernel>`

*/
class ConstructOppositeDirection_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the direction opposite to `d`.
  */
  Kernel::Direction_2 operator()(const
                                 Kernel::Direction_2& d);


  /// @}

}; /* end Kernel::ConstructOppositeDirection_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Direction_3<Kernel>`

*/
class ConstructOppositeDirection_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the direction opposite to `d`.
  */
  Kernel::Direction_3 operator()(const
                                 Kernel::Direction_3& d);


  /// @}

}; /* end Kernel::ConstructOppositeDirection_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Line_2<Kernel>`

*/
class ConstructOppositeLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the line representing the same set of points as `l`,
    but with opposite direction.
  */
  Kernel::Line_2 operator()(const Kernel::Line_2& l);


  /// @}

}; /* end Kernel::ConstructOppositeLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Line_3<Kernel>`

*/
class ConstructOppositeLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the line representing the same set of points as `l`,
    but with opposite direction.
  */
  Kernel::Line_3 operator()(const Kernel::Line_3& l);

  /// @}

}; /* end Kernel::ConstructOppositeLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructOppositePlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the plane representing the same set of points as `p`,
    but with opposite orientation.
  */
  Kernel::Plane_3 operator()(const Kernel::Plane_3& p);

  /// @}

}; /* end Kernel::ConstructOppositePlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Ray_2<Kernel>`

*/
class ConstructOppositeRay_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the ray with the same source as `r`, but in opposite direction.
  */
  Kernel::Ray_2 operator()(const Kernel::Ray_2& r);

  /// @}

}; /* end Kernel::ConstructOppositeRay_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Ray_3<Kernel>`

*/
class ConstructOppositeRay_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the ray with the same source as `r`, but in opposite direction.
  */
  Kernel::Ray_3 operator()(const Kernel::Ray_3& r);


  /// @}

}; /* end Kernel::ConstructOppositeRay_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructOppositeSegment_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the segment representing the same set of points as `s`,
    but with opposite orientation.
  */
  Kernel::Segment_2 operator()(const Kernel::Segment_2& s);

  /// @}

}; /* end Kernel::ConstructOppositeSegment_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructOppositeSegment_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the segment representing the same set of points as `s`,
    but with opposite orientation.
  */
  Kernel::Segment_3 operator()(const Kernel::Segment_3& s);


  /// @}

}; /* end Kernel::ConstructOppositeSegment_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Sphere_3<Kernel>`

*/
class ConstructOppositeSphere_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the sphere with the same center and squared radius as
    `s`, but with opposite orientation.
  */
  Kernel::Sphere_3 operator()(const Kernel::Sphere_3& s);

  /// @}

}; /* end Kernel::ConstructOppositeSphere_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Triangle_2<Kernel>`

*/
class ConstructOppositeTriangle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the triangle with opposite orientation to `t`
    (this flips the positive and the negative side, but
    not bounded and unbounded side).
  */
  Kernel::Triangle_2 operator()(const Kernel::Triangle_2& t);


  /// @}

}; /* end Kernel::ConstructOppositeTriangle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Vector_2<Kernel>`

*/
class ConstructOppositeVector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the vector `-v`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2& v);

  /// @}

}; /* end Kernel::ConstructOppositeVector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

\sa `CGAL::Vector_3<Kernel>`

*/
class ConstructOppositeVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the vector `-v`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3& v);

  /// @}

}; /* end Kernel::ConstructOppositeVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `Kernel::ConstructCrossProductVector_3`

*/
class ConstructOrthogonalVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns a vector that is orthogonal to the plane `p` and directed
    to the positive side of `p`.
  */
  Kernel::Vector_3 operator()(const Kernel::Plane_3& p);

  /*!
    returns a vector that is orthogonal to the plane defined by
    `Kernel::ConstructPlane_3()(p, q, r)` and directed
    to the positive side of this plane.
  */
  Kernel::Vector_3 operator()(const Kernel::Point_3& p,
                              const Kernel::Point_3& q, const Kernel::Point_3& r);


  /// @}

}; /* end Kernel::ConstructOrthogonalVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Direction_2<Kernel>`

*/
class ConstructPerpendicularDirection_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a direction orthogonal to `d`. If `o` is
    \ref CGAL::CLOCKWISE, `d` is rotated clockwise; if `o` is
    \ref CGAL::COUNTERCLOCKWISE, `d` is rotated counterclockwise.
    \pre `o != CGAL::COLLINEAR.`
  */
  Kernel::Direction_2 operator()(const Kernel::Direction_2& d,
                                 Orientation o);


  /// @}

}; /* end Kernel::ConstructPerpendicularDirection_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_2<Kernel>`

*/
class ConstructPerpendicularLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the line perpendicular to `l` and passing through `p`,
    where the direction is the direction of `l` rotated
    counterclockwise by 90 degrees.
  */
  Kernel::Line_2 operator()(const Kernel::Line_2& l,
                            const Kernel::Point_2& p);

  /// @}

}; /* end Kernel::ConstructPerpendicularLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructPerpendicularLine_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the line that is perpendicular to `pl` and that
    passes through point `p`. The line is oriented from
    the negative to the positive side of `pl`
  */
  Kernel::Line_3 operator()(const Kernel::Plane_3& pl,
                            const Kernel::Point_3& p);


  /// @}

}; /* end Kernel::ConstructPerpendicularLine_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructPerpendicularPlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the plane perpendicular to `l` passing through `p`,
    such that the normal direction of the plane coincides with the direction of
  the line.
  */
  Kernel::Plane_3 operator()(const Kernel::Line_3& l,
                             const Kernel::Point_3& p);

  /// @}

}; /* end Kernel::ConstructPerpendicularPlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`

*/
class ConstructPerpendicularVector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns `v` rotated clockwise by 90 degrees, if `o` is
    \ref CGAL::CLOCKWISE, and rotated counterclockwise otherwise.
    \pre `o != CGAL::COLLINEAR`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2& v,
                              Orientation o);


  /// @}

}; /* end Kernel::ConstructPerpendicularVector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructPlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    creates a plane defined by the equation
    \f$ a\, x +b\, y +c\, z + d = 0\f$.
    Notice that it is degenerate if \f$ a = b = c = 0\f$.
  */
  Kernel::Plane_3 operator()(const Kernel::RT &a,
                             const Kernel::RT &b,
                             const Kernel::RT &c,
                             const Kernel::RT &d);

  /*!
    creates a plane passing through the points `p`,
    `q` and `r`. The plane is oriented such that `p`,
    `q` and `r` are oriented in a positive sense
    (that is counterclockwise) when seen from the positive side of the plane.
    Notice that it is degenerate if the points are collinear.
  */
  Kernel::Plane_3 operator()(const Kernel::Point_3 &p,
                             const Kernel::Point_3 &q,
                             const Kernel::Point_3 &r);

  /*!
    introduces a plane that passes through point `p` and
    that has as an orthogonal direction equal to `d`.
  */
  Kernel::Plane_3 operator()(const Kernel::Point_3 &p,
                             const Kernel::Direction_3&d);

  /*!
    introduces a plane that passes through point `p` and
    that is orthogonal to `v`.
  */
  Kernel::Plane_3 operator()(const Kernel::Point_3 &p,
                             const Kernel::Vector_3 &v);

  /*!
    introduces a plane that is defined through the three points
    `l.point(0)`, `l.point(1)` and `p`.
  */
  Kernel::Plane_3 operator()(const Kernel::Line_3 &l,
                             const Kernel::Point_3 &p);

  /*!
    introduces a plane that is defined through the three points
    `r.point(0)`, `r.point(1)` and `p`.
  */
  Kernel::Plane_3 operator()(const Kernel::Ray_3 &r,
                             const Kernel::Point_3 &p);

  /*!
    introduces a plane that is defined through the three points
    `s.source()`, `s.target()` and `p`.
  */
  Kernel::Plane_3 operator()(const Kernel::Segment_3 &s,
                             const Kernel::Point_3 &p);

  /*!
    introduces a plane that is defined as the plane containing the circle.
  */
  Kernel::Plane_3 operator()(const Kernel::Circle_3 &c);

  /// @}

}; /* end Kernel::ConstructPlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructPointOn_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an arbitrary point on `l`. It holds
    `point(i) == point(j)`, iff `i==j`.
    Furthermore, is directed from `point(i)`
    to `point(j)`, for all `i < j`.
  */
  Kernel::Point_2 operator()(const Kernel::Line_2& l,
                             const Kernel::FT i);

  /*!
    returns a point on `r`. `point(0)` is the source,
    `point(i)`, with `i>0`, is different from the
    source. \pre `i>= 0`.
  */
  Kernel::Point_2 operator()(const Kernel::Ray_2& r,
                             const Kernel::FT i);

  /*!
    returns source or target of `s`: `point(0)` returns
    the source of `s`, `point(1)` returns the target of `s`.
    The parameter `i` is taken modulo 2, which gives
    easy access to the other end point.
  */
  Kernel::Point_2 operator()(const Kernel::Segment_2& s,
                             int i);


  /// @}

}; /* end Kernel::ConstructPointOn_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_3<Kernel>`
  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Ray_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructPointOn_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns an arbitrary point on `l`. It holds
    `point(i) == point(j)`, iff `i==j`.
    Furthermore, is directed from `point(i)`
    to `point(j)`, for all `i < j`.
  */
  Kernel::Point_3 operator()(const Kernel::Line_3& l,
                             const Kernel::FT i);

  /*!
    returns `point(0)` on `l`, identical to `operator()(l,0)`.
  */
  Kernel::Point_3 operator()(const Kernel::Line_3& l);

  /*!
    returns an arbitrary point on `h`.
  */
  Kernel::Point_3 operator()(const Kernel::Plane_3& h);

  /*!
    returns a point on `r`. `point(0)` is the source,
    `point(i)`, with `i>0`, is different from the
    source. \pre `i >= 0`.
  */
  Kernel::Point_3 operator()(const Kernel::Ray_3& r,
                             const Kernel::FT i);

  /*!
    returns source or target of `s`: `point(0)` returns
    the source of `s`, `point(1)` returns the target of `s`.
  The parameter `i` is taken modulo 2, which gives
  easy access to the other end point.
  */
  Kernel::Point_3 operator()(const Kernel::Segment_3& s,
                             int i);

  /// @}

}; /* end Kernel::ConstructPointOn_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Point_2<Kernel>`

*/
class ConstructPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a variable with %Cartesian coordinates
    \f$ (0,0)\f$.
  */
  Kernel::Point_2 operator()(const CGAL::Origin &CGAL::ORIGIN);

  /*!
    returns `p`.

    \note It is advised to return a const reference to `p` to avoid useless copies.

    \note This peculiar requirement is necessary because some \cgal structures such as triangulations
    internally manipulate points whose type might be `Point_2` or `Weighted_point_2`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2& p);

 /*!
    extracts the bare point from the weighted point.
  */
  Kernel::Point_2 operator()(const Kernel::Weighted_point_2& wp);

  ///@}

}; /* end Kernel::ConstructPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Point_3<Kernel>`

*/
class ConstructPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a point with %Cartesian coordinates\f$ (0,0,0)\f$.
  */
  Kernel::Point_3 operator()(const CGAL::Origin &CGAL::ORIGIN);

  /*!
    returns `p`.

    \note It is advised to return a const reference to `p` to avoid useless copies.

    \note This peculiar requirement is necessary because some \cgal structures such as triangulations
    internally manipulate points whose type might be `Point_3` or `Weighted_point_3`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3& p);

 /*!
    extracts the bare point from the weighted point.
  */
  Kernel::Point_3 operator()(const Kernel::Weighted_point_3& wp);

  ///@}

}; /* end Kernel::ConstructPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_2<Kernel>`

*/
class ConstructProjectedPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the orthogonal projection of `p` onto `l`.
  */
  Kernel::Point_2 operator()(const Kernel::Line_2& l,
                             const Kernel::Point_2& p);

  /// @}

}; /* end Kernel::ConstructProjectedPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_3<Kernel>`
  \sa `CGAL::Plane_3<Kernel>`

*/
class ConstructProjectedPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the orthogonal projection of `p` onto `l`.
  */
  Kernel::Point_3 operator()(const Kernel::Line_3& l,
                             const Kernel::Point_3& p);

  /*!
    returns the orthogonal projection of `p` onto `h`.
  */
  Kernel::Point_3 operator()(const Kernel::Plane_3& h,
                             const Kernel::Point_3& p);

  /*!
    returns the point of `s` that is the closest to `p`.
  */
  Kernel::Point_3 operator()(const Kernel::Segment_3& s,
                             const Kernel::Point_3& p);

  /*!
    returns the point of `r` that is the closest to `p`.
  */
  Kernel::Point_3 operator()(const Kernel::Ray_3& r,
                             const Kernel::Point_3& p);

  /*!
    returns the point of `t` that is the closest to `p`.
  */
  Kernel::Point_3 operator()(const Kernel::Triangle_3& t,
                             const Kernel::Point_3& p);

  /// @}

}; /* end Kernel::ConstructProjectedPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `ConstructLiftedPoint_3`
*/
class ConstructProjectedXYPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the image point of the projection of `p` under an affine
    transformation, which maps `h` onto the \f$ xy\f$-plane, with the
    \f$ z\f$-coordinate removed.
    This affine transformation must be the inverse of the affine transformation used
    in `ConstructLiftedPoint_3`.
  */
  Kernel::Point_2 operator()(const Kernel::Plane_3& h,
                             const Kernel::Point_3& p);

  /// @}

}; /* end Kernel::ConstructProjectedXYPoint_2 */



/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Weighted_point_2<Kernel>`

*/
class ConstructRadicalAxis_2
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the radical line of the weighted points.
  */
  Kernel::Line_2 operator()(const Kernel::Weighted_point_2& wp1,
                            const Kernel::Weighted_point_2& wp2);
  /// @}
};


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`

*/
class ConstructRadicalLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the radical line of the circles.
    \pre The two circles don't have the same center.
  */
  Kernel::Line_2 operator()
  (const Kernel::Circle_2& c1,
   const Kernel::Circle_2& c2);


  /// @}

}; /* end Kernel::ConstructRadicalLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Sphere_3<Kernel>`

*/
class ConstructRadicalPlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the radical plane of the spheres.
    \pre The two spheres don't have the same center.
  */
  Kernel::Plane_3 operator()
  (const Kernel::Sphere_3& sphere1,
   const Kernel::Sphere_3& sphere2);

  /// @}

}; /* end Kernel::ConstructRadicalPlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Ray_2<Kernel>`

*/
class ConstructRay_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a ray
    with source `p` and passing through point `q`.
  */
  Kernel::Ray_2 operator()(const Kernel::Point_2 &p,
                           const Kernel::Point_2 &q);

  /*!
    introduces a ray starting at source `p` with
    the direction given by `v`.
  */
  Kernel::Ray_2 operator()(const Kernel::Point_2 &p,
                           const Kernel::Vector_2 &v);

  /*!
    introduces a ray starting at source `p` with
    direction `d`.
  */
  Kernel::Ray_2 operator()(const Kernel::Point_2 &p,
                           const Kernel::Direction_2 &d);

  /*!
    introduces a ray starting at source `p` with
    the same direction as `l`.
  */
  Kernel::Ray_2 operator()(const Kernel::Point_2 &p,
                           const Kernel::Line_2 &l);

  ///@}

}; /* end Kernel::ConstructRay_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Ray_3<Kernel>`

*/
class ConstructRay_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a ray
    with source `p` and passing through point `q`.
  */
  Kernel::Ray_3 operator()(const Kernel::Point_3 &p,
                           const Kernel::Point_3 &q);

  /*!
    introduces a ray with source `p` and with
    the direction given by `v`.
  */
  Kernel::Ray_3 operator()(const Kernel::Point_3 &p,
                           const Kernel::Vector_3 &v);

  /*!
    introduces a ray with source `p` and with
    direction `d`.
  */
  Kernel::Ray_3 operator()(const Kernel::Point_3 &p,
                           const Kernel::Direction_3 &d);

  /*!
    introduces a ray with source `p` and with
    the same direction as `l`.
  */
  Kernel::Ray_3 operator()(const Kernel::Point_3 &p,
                           const Kernel::Line_3 &l);

  /// @}

}; /* end Kernel::ConstructRay_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`
*/
class ConstructScaledVector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    produces the vector `v` scaled by a factor `scale`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2 &v, const Kernel::RT& scale);

  /*!
    produces the vector `v` scaled by a factor `scale`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2 &v, const Kernel::FT& scale);

  /// @}

}; /* end Kernel::ConstructScaledVector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`

*/
class ConstructScaledVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    produces the vector `v` scaled by a factor `scale`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v,
                              const Kernel::RT& scale);

  /*!
    produces the vector `v` scaled by a factor `scale`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v,
                              const Kernel::FT& scale);


  /// @}

}; /* end Kernel::ConstructScaledVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Ray_2<Kernel>`

*/
class ConstructSecondPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns a point different from the source on
    the ray `r`.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Ray_2 &r);


  /// @}

}; /* end Kernel::ConstructSecondPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Ray_3<Kernel>`

*/
class ConstructSecondPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns a point different from the source on
    the ray `r`.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Ray_3 &r);


  /// @}

}; /* end Kernel::ConstructSecondPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructSegment_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a segment with source `p`
    and target `q`. The segment is directed from the source towards
  the target.
  */
  Kernel::Segment_2 operator()(const Kernel::Point_2 &p, const Kernel::Point_2 &q);

  /// @}

}; /* end Kernel::ConstructSegment_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructSegment_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a segment with source `p`
    and target `q`. It is directed from the source towards
    the target.
  */
  Kernel::Segment_3 operator()(const Kernel::Point_3 &p, const Kernel::Point_3 &q);

  /// @}

}; /* end Kernel::ConstructSegment_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`

*/
class ConstructSource_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the source of
    the segment `s`.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Segment_2 &s);

  /*!
    returns the source of
    the ray `r`.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Ray_2 &r);


  /// @}

}; /* end Kernel::ConstructSource_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_3<Kernel>`
  \sa `CGAL::Ray_3<Kernel>`

*/
class ConstructSource_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the source of
    the segment `s`.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Segment_3 &s);

  /*!
    returns the source of
    the ray `r`.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Ray_3 &r);


  /// @}

}; /* end Kernel::ConstructSource_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `CGAL::Sphere_3<Kernel>`

*/
class ConstructSphere_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a sphere initialized to the sphere with center `center`,
    squared radius `squared_radius` and orientation
    `orientation`.
    \pre `orientation != CGAL::COPLANAR` and `squared_radius >= 0`.
  */
  Kernel::Sphere_3 operator()(const Kernel::Point_3 & center,
                              const Kernel::FT & squared_radius,
                              const Orientation & orientation = COUNTERCLOCKWISE);

  /*!
    introduces a sphere initialized to the unique sphere which passes
    through the points `p`, `q`, `r` and `s`. The
    orientation of the sphere is the orientation of the point quadruple
    `p`, `q`, `r`, `s`.
    \pre `p`, `q`, `r`, and `s` are not coplanar.
  */
  Kernel::Sphere_3 operator()( const Kernel::Point_3 & p,
                               const Kernel::Point_3 & q,
                               const Kernel::Point_3 & r,
                               const Kernel::Point_3 & s);

  /*!
    introduces a sphere initialized to the smallest sphere which passes
    through the points `p`, `q`, and `r`. The orientation of
    the sphere is `o`. \pre `o != CGAL::COPLANAR`.
  */
  Kernel::Sphere_3 operator()(const Kernel::Point_3 & p,
                              const Kernel::Point_3 & q,
                              const Kernel::Point_3 & r,
                              const Orientation& o = COUNTERCLOCKWISE);

  /*!
    introduces a sphere initialized to the smallest sphere which passes
    through the points `p` and `q`. The orientation of
    the sphere is `o`. \pre `o != CGAL::COPLANAR`.
  */
  Kernel::Sphere_3 operator()(const Kernel::Point_3 & p,
                              const Kernel::Point_3 & q,
                              const Orientation & o = COUNTERCLOCKWISE);

  /*!
    introduces a sphere `s` initialized to the sphere with center
    `center`, squared radius zero and orientation `orientation`.
    \pre `orientation != CGAL::COPLANAR`.
    \post `s.is_degenerate()` = `true`.
  */
  Kernel::Sphere_3 operator()( const Kernel::Point_3 & center,
                               const Orientation & orientation = COUNTERCLOCKWISE);

  /*!
    introduces a sphere initialized to the diametral sphere of
    the circle.
  */
  Kernel::Sphere_3 operator()( const Kernel::Circle_3 &c);

  /// @}

}; /* end Kernel::ConstructSphere_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`

*/
class ConstructSumOfVectors_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    introduces the vector `v1 + v2`.
  */
  Kernel::Vector_2 operator()(const Kernel::Vector_2 &v1,
                              const Kernel::Vector_2 &v2);


  /// @}

}; /* end Kernel::ConstructSumOfVectors_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`

*/
class ConstructSumOfVectors_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `v1 + v2`.
  */
  Kernel::Vector_3 operator()(const Kernel::Vector_3 &v1,
                              const Kernel::Vector_3 &v2);

  /// @}

}; /* end Kernel::ConstructSumOfVectors_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Triangle_3<Kernel>`

*/
class ConstructSupportingPlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the supporting plane of `t`, with same orientation.
  */
  Kernel::Plane_3 operator()(const Kernel::Triangle_3& t);

  ///@}

}; /* end Kernel::ConstructSupportingPlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_2<Kernel>`

*/
class ConstructTarget_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the target of
    the segment `s`.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Segment_2 &s);

  /// @}

}; /* end Kernel::ConstructTarget_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Segment_3<Kernel>`

*/
class ConstructTarget_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the target of
    the segment `s`.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Segment_3 &s);

  /// @}

}; /* end Kernel::ConstructTarget_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class ConstructTetrahedron_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a tetrahedron with vertices `p_0`, `p_1`, `p_2` and `p_3`.
  */
  Kernel::Tetrahedron_3 operator()(const Kernel::Point_3 &p0,
                                   const Kernel::Point_3 &p1,
                                   const Kernel::Point_3 &p2,
                                   const Kernel::Point_3 &p3);

  /// @}

}; /* end Kernel::ConstructTetrahedron_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Point_2<Kernel>`

*/
class ConstructTranslatedPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the point obtained by translating `p` by the vector
    `v`.
  */
  Kernel::Point_2 operator()(const Kernel::Point_2& p,
                             const Kernel::Vector_2& v);

  /*!
    returns the point obtained by translating a point at the origin by the vector
    `v`.
  */
  Kernel::Point_2 operator()(const CGAL::Origin& o,
                             const Kernel::Vector_2& v);


  /// @}

}; /* end Kernel::ConstructTranslatedPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Point_3<Kernel>`

*/
class ConstructTranslatedPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the point obtained by translating `p` by the vector
    `v`.
  */
  Kernel::Point_3 operator()(const Kernel::Point_3& p,
                             const Kernel::Vector_3& v);

  /*!
    returns the point obtained by translating a point at the origin by the vector
    `v`.
  */
  Kernel::Point_3 operator()(const Origin& o,
                             const Kernel::Vector_3& v);


  /// @}

}; /* end Kernel::ConstructTranslatedPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Triangle_2<Kernel>`

*/
class ConstructTriangle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a triangle with vertices `p`, `q` and `r`.
  */
  Kernel::Triangle_2 operator()(const Kernel::Point_2 &p,
                                const Kernel::Point_2 &q,
                                const Kernel::Point_2 &r);

  /// @}

}; /* end Kernel::ConstructTriangle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Triangle_3<Kernel>`

*/
class ConstructTriangle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    introduces a triangle with vertices `p`, `q` and `r`.
  */
  Kernel::Triangle_3 operator()(const Kernel::Point_3 &p,
                                const Kernel::Point_3 &q,
                                const Kernel::Point_3 &r);


  /// @}

}; /* end Kernel::ConstructTriangle_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableTernaryFunction}

\sa `normal_grp`
\sa `unit_normal_grp`

*/
class ConstructUnitNormal_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the unit normal of the vectors `q-p` and `r-p`.
    This requires that `Kernel::FT` supports the `sqrt` operation.
  */
  Kernel::Vector_3 operator()(const Kernel::Point_3& p,
                              const Kernel::Point_3& q,
                              const Kernel::Point_3& r );


  /// @}

}; /* end Kernel::ConstructUnitNormal_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_2<Kernel>`
  \sa `Kernel::ConstructScaledVector_2`

*/
class ConstructVector_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    introduces the vector `b-a`.
  */
  Kernel::Vector_2 operator()(const Kernel::Point_2 &a,
                              const Kernel::Point_2 &b);

  /*!
    introduces the vector `b`.
  */
  Kernel::Vector_2 operator()(const CGAL::Origin &o,
                              const Kernel::Point_2 &b);

  /*!
    introduces the vector `-a`.
  */
  Kernel::Vector_2 operator()(const Kernel::Point_2 &a,
                              const CGAL::Origin &o);

  /*!
    introduces the vector `s.target()-s.source()`.
  */
  Kernel::Vector_2 operator()(const Kernel::Segment_2 &s);

  /*!
    introduces a vector having the same direction as `r`.
  */
  Kernel::Vector_2 operator()(const Kernel::Ray_2 &r);

  /*!
    introduces a vector having the same direction as `l`.
  */
  Kernel::Vector_2 operator()(const Kernel::Line_2 &l);

  /*!
    introduces a null vector.
  */
  Kernel::Vector_2 operator()(const CGAL::Null_vector &v);

  /// @}

}; /* end Kernel::ConstructVector_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Vector_3<Kernel>`
  \sa `Kernel::ConstructScaledVector_3`

*/
class ConstructVector_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces the vector `b-a`.
  */
  Kernel::Vector_3 operator()(const Kernel::Point_3 &a,
                              const Kernel::Point_3 &b);

  /*!
    introduces the vector `b`.
  */
  Kernel::Vector_3 operator()(const CGAL::Origin &o,
                              const Kernel::Point_3 &b);

  /*!
    introduces the vector `-a`.
  */
  Kernel::Vector_3 operator()(const Kernel::Point_3 &a,
                              const CGAL::Origin &o);

  /*!
    introduces the vector `s.target()-s.source()`.
  */
  Kernel::Vector_3 operator()(const Kernel::Segment_3 &s);

  /*!
    introduces a vector having the same direction as `r`.
  */
  Kernel::Vector_3 operator()(const Kernel::Ray_3 &r);

  /*!
    introduces a vector having the same direction as `l`.
  */
  Kernel::Vector_3 operator()(const Kernel::Line_3 &l);

  /*!
    introduces a null vector.
  */
  Kernel::Vector_3 operator()(const CGAL::Null_vector &v);


  /// @}

}; /* end Kernel::ConstructVector_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`

*/
class ConstructVertex_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns source or target of `s`: `fo``(s,0)`
    returns the source of `s`, `fo``(s,1)` returns the target
    of `s`. The parameter `i` is taken modulo 2.
  */
  Kernel::Point_2 operator()(const Kernel::Segment_2
                             &s, int i);

  /*!
    returns the i-th vertex of
    `r` in counterclockwise order, starting with the lower left
    vertex. The parameter `i` is taken modulo 4.
  */
  Kernel::Point_2 operator()(const
                             Kernel::Iso_rectangle_2 &r, int i);

  /*!
    returns the i-th vertex of `t`. The parameter
    `i` is taken modulo 3.
  */
  Kernel::Point_2 operator()(const Kernel::Triangle_2
                             &t, int i);

  /// @}

}; /* end Kernel::ConstructVertex_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\image html IsoCuboid.png
\image latex IsoCuboid.png

\cgalRefines{AdaptableBinaryFunction}

\sa `CGAL::Iso_cuboid_3<Kernel>`
\sa `CGAL::Segment_3<Kernel>`
\sa `CGAL::Tetrahedron_3<Kernel>`
\sa `CGAL::Triangle_3<Kernel>`

*/
class ConstructVertex_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns source or target of `s`: `fo``(s,0)`
    returns the source of `s`, `fo``(s,1)` returns the target
    of `s`. The parameter `i` is taken modulo 2.
  */
  Kernel::Point_3 operator()(const Kernel::Segment_3
                             &s, int i);

  /*!
    returns the i-th vertex of
    `c`, as indicated in the figure below. The parameter `i` is
    taken modulo 8.

  */
  Kernel::Point_3 operator()(const
                             Kernel::Iso_cuboid_3 &c, int i);

  /*!
    returns the i-th vertex of `t`. The parameter
    `i` is taken modulo 3.
  */
  Kernel::Point_3 operator()(const Kernel::Triangle_3
                             &t, int i);

  /*!
    returns the i-th vertex of
    `t`. The parameter `i` is taken modulo 4.
  */
  Kernel::Point_3 operator()(const
                             Kernel::Tetrahedron_3 &t, int i);

  /// @}

}; /* end Kernel::ConstructVertex_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableBinaryFunction}

\sa `CGAL::Weighted_point_2<Kernel>`

*/
class ConstructWeightedCircumcenter_2
{
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs the point which is the center of the smallest orthogonal circle to the input weighted points.
  */
  Kernel::Point_2 operator()(const Kernel::Weighted_point_2& p,
                             const Kernel::Weighted_point_2& q,
                             const Kernel::Weighted_point_2& s);
  /// @}
}; /* end Kernel::ConstructWeightedCircumcenter_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableBinaryFunction}

\sa `CGAL::Weighted_point_3<Kernel>`

*/
class ConstructWeightedCircumcenter_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    constructs the point which is the center of the smallest orthogonal sphere to the input weighted points.
  */
  Kernel::Point_3 operator()(const Kernel::Weighted_point_3& p,
                             const Kernel::Weighted_point_3& q,
                             const Kernel::Weighted_point_3& r,
                             const Kernel::Weighted_point_3& s);

  /*!
  constructs the point which is the center of the smallest orthogonal sphere to the input weighted points.
  */
  Kernel::Point_3 operator()(const Kernel::Weighted_point_3& p,
                             const Kernel::Weighted_point_3& q,
                             const Kernel::Weighted_point_3& r);

  /*!
  constructs the point which is the center of the smallest orthogonal sphere to the input weighted points.
  */
  Kernel::Point_3 operator()(const Kernel::Weighted_point_3& p,
                             const Kernel::Weighted_point_3& q);

  /// @}
}; /* end Kernel::ConstructWeightedCircumcenter_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Weighted_point_2<Kernel>`
*/
class ConstructWeightedPoint_2
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a weighted point with %Cartesian coordinates
    \f$ (0,0)\f$ and weight \f$ 0 \f$.
  */
  Kernel::Weighted_point_2 operator()(const CGAL::Origin &CGAL::ORIGIN);

 /*!
    introduces a weighted point with %Cartesian coordinates
    those of \f$ p \f$ and weight \f$ 0 \f$.
  */
  Kernel::Weighted_point_2 operator()(const Kernel::Point_2& p);

  /*!
     introduces a weighted point with %Cartesian coordinates
     those of \f$ p \f$ and weight \f$ w \f$.
   */
   Kernel::Weighted_point_2 operator()(const Kernel::Point_2& p, const Kernel::FT& w);
  ///@}

}; /* end Kernel::ConstructWeightedPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Weighted_point_3<Kernel>`

*/
class ConstructWeightedPoint_3
{
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    introduces a weighted point with %Cartesian coordinates
    \f$ (0,0,0)\f$ and weight \f$ 0 \f$.
  */
  Kernel::Weighted_point_3 operator()(const CGAL::Origin &CGAL::ORIGIN);

 /*!
    introduces a weighted point with %Cartesian coordinates
    those of \f$ p \f$ and weight \f$ 0 \f$.
  */
  Kernel::Weighted_point_3 operator()(const Kernel::Point_3& p);

  /*!
     introduces a weighted point with %Cartesian coordinates
     those of \f$ p \f$ and weight \f$ w \f$.
   */
   Kernel::Weighted_point_3 operator()(const Kernel::Point_3& p, const Kernel::FT& w);

  ///@}

}; /* end Kernel::ConstructWeightedPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `coplanar_orientation_grp`

*/
class CoplanarOrientation_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    If `p`, `q`, and `s` are collinear, then \ref CGAL::COLLINEAR is returned.

    If not, let `P` be the plane defined by the points `p`, `q`,
    and `r`. Note that the order defines the orientation of
    `P`. If `P` and the plane defined by `p`, `q`, and `s`
    have the same orientation, then \ref CGAL::POSITIVE is returned;
    otherwise \ref CGAL::NEGATIVE is returned.

    \pre `p`, `q`, `r`, and `s` are coplanar and `p`, `q`, and `r` are not collinear.
  */
  Orientation operator()(const Kernel::Point_3&p,
                         const Kernel::Point_3&q,
                         const Kernel::Point_3&r,
                         const Kernel::Point_3&s);

  /*!
    If `p`, `q`, and `r` are collinear, then \ref CGAL::COLLINEAR is returned.

    If not, let `P` be the plane defined by the points `p`, `q`, and `r`.
    The return value in this case is either \ref CGAL::POSITIVE or \ref CGAL::NEGATIVE,
    but we don't specify it explicitly.
    However, we guarantee that all calls to this predicate over 3 points in `P`
    will return a coherent orientation if considered as a 2D orientation in `P`.
  */
  Orientation operator()(const Kernel::Point_3&p,
                         const Kernel::Point_3&q,
                         const Kernel::Point_3&r);


  /// @}

}; /* end Kernel::CoplanarOrientation_3 */






/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `coplanar_side_of_bounded_circle_grp`

*/
class CoplanarSideOfBoundedCircle_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the bounded side of the circle defined
    by `p`, `q`, and `r` on which `s` lies.
    \pre `p`, `q`, `r`, and `s` are coplanar and `p`, `q`, and `r` are not collinear.
  */
  Bounded_side operator()(const Kernel::Point_3&p,
                          const Kernel::Point_3&q,
                          const Kernel::Point_3&r,
                          const Kernel::Point_3&s);

  /// @}

}; /* end Kernel::CoplanarSideOfBoundedCircle_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `coplanar_grp`

*/
class Coplanar_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, if `p`, `q`, `r`, and `s` are coplanar.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r,
                  const Kernel::Point_3&s);

  /// @}

}; /* end Kernel::Coplanar_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `CGAL::Direction_2<Kernel>`

*/
class CounterclockwiseInBetween_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns `true` iff `d` is not equal to `d1`, and
    while rotating counterclockwise starting at `d1`,
    `d` is reached strictly before `d2` is reached.
    Note that true is returned if `d1` == `d2`, unless
    also `d` == `d1`.
  */
  bool operator()(const Kernel::Direction_2&d,
                  const Kernel::Direction_2&d1,
                  const Kernel::Direction_2&d2);


  /// @}

}; /* end Kernel::CounterclockwiseInBetween_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `do_intersect_grp`

*/
class DoIntersect_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    determines if two geometrical objects of type
    `Type1` and `Type2` intersect or not, for all pairs `Type1` and `Type2`, where
    the types `Type1` and `Type2` can be any of the
    following:

    - `Kernel::Point_2`
    - `Kernel::Line_2`
    - `Kernel::Ray_2`
    - `Kernel::Segment_2`
    - `Kernel::Triangle_2`
    - `Kernel::Iso_rectangle_2`

  */
  bool operator()(Type1 obj1, Type2 obj2);

  /// @}

}; /* end Kernel::DoIntersect_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `do_intersect_grp`

*/
class DoIntersect_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    determines if two geometrical objects of type
    `Type1` and `Type2` intersect or not, for all pairs `Type1` and `Type2`, where
    the types `Type1` and
    `Type2` can be any of the following:

    - `Kernel::Point_2`
    - `Kernel::Plane_3`
    - `Kernel::Line_3`
    - `Kernel::Ray_3`
    - `Kernel::Segment_3`
    - `Kernel::Triangle_3`
    - `CGAL::Bbox_3`

    and also for `Type1` and `Type2` of respective types

    - `Kernel::Triangle_3` and `Kernel::Tetrahedron_3`
    - `Kernel::Plane_3` and `Kernel::Sphere_3` (or the contrary)
    - `Kernel::Sphere_3` and `Kernel::Sphere_3`.

  */
  bool operator()(Type1 obj1, Type2 obj2);

  /// @}

}; /* end Kernel::DoIntersect_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_xy_grp`

*/
class EqualXY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ x\f$-coordinate
    and the same %Cartesian \f$ y\f$-coordinate.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::EqualXY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `x_equal_grp`

*/
class EqualX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ x\f$-coordinate.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::EqualX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

\sa `x_equal_grp`

*/
class EqualX_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ x\f$-coordinate.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::EqualX_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `y_equal_grp`

*/
class EqualY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ y\f$-coordinate.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::EqualY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `y_equal_grp`

*/
class EqualY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ y\f$-coordinate.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::EqualY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `z_equal_grp`

*/
class EqualZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` and `q` have the same %Cartesian \f$ z\f$-coordinate.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::EqualZ_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Direction_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Point_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`
  \sa `CGAL::Vector_2<Kernel>`

*/
class Equal_2 {
public:

  /// \name Operations
  /// A model of this concept must provide the following operations. For each of
  /// them, it returns `true` iff `x` and `y` are equal.
  /// @{

  /*!

   */
  bool operator()(const Kernel::Point_2& x,
                  const Kernel::Point_2& y);

  /*!

   */
  bool operator()(const Kernel::Vector_2& x,
                  const Kernel::Vector_2& y);

  /*!

   */
  bool operator()(const Kernel::Direction_2& x,
                  const Kernel::Direction_2& y);

  /*!

   */
  bool operator()(const Kernel::Line_2& x,
                  const Kernel::Line_2& y);

  /*!

   */
  bool operator()(const Kernel::Ray_2& x,
                  const Kernel::Ray_2& y);

  /*!

   */
  bool operator()(const Kernel::Segment_2& x,
                  const Kernel::Segment_2& y);

  /*!

   */
  bool operator()(const Kernel::Circle_2& x,
                  const Kernel::Circle_2& y);

  /*!

   */
  bool operator()(const Kernel::Triangle_2& x,
                  const Kernel::Triangle_2& y);

  /*!

   */
  bool operator()(const Kernel::Iso_rectangle_2& x,
                  const Kernel::Iso_rectangle_2& y);

  /// @}

}; /* end Kernel::Equal_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Direction_3<Kernel>`
  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Line_3<Kernel>`
  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Point_3<Kernel>`
  \sa `CGAL::Ray_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`
  \sa `CGAL::Triangle_3<Kernel>`
  \sa `CGAL::Vector_3<Kernel>`

*/
class Equal_3 {
public:

  /// \name Operations
  /// A model of this concept must provide the following operations. For each of
  /// them, it returns `true` iff `x` and `y` are equal.
  /// @{


  /*!

   */
  bool operator()(const Kernel::Point_3& x,
                  const Kernel::Point_3& y);

  /*!

   */
  bool operator()(const Kernel::Vector_3& x,
                  const Kernel::Vector_3& y);

  /*!

   */
  bool operator()(const Kernel::Direction_3& x,
                  const Kernel::Direction_3& y);

  /*!

   */
  bool operator()(const Kernel::Line_3& x,
                  const Kernel::Line_3& y);

  /*!

   */
  bool operator()(const Kernel::Plane_3& x,
                  const Kernel::Plane_3& y);

  /*!

   */
  bool operator()(const Kernel::Ray_3& x,
                  const Kernel::Ray_3& y);

  /*!

   */
  bool operator()(const Kernel::Segment_3& x,
                  const Kernel::Segment_3& y);

  /*!

   */
  bool operator()(const Kernel::Circle_3& x,
                  const Kernel::Circle_3& y);

  /*!

   */
  bool operator()(const Kernel::Sphere_3& x,
                  const Kernel::Sphere_3& y);

  /*!

   */
  bool operator()(const Kernel::Triangle_3& x,
                  const Kernel::Triangle_3& y);

  /*!

   */
  bool operator()(const Kernel::Tetrahedron_3& x,
                  const Kernel::Tetrahedron_3& y);

  /*!

   */
  bool operator()(const Kernel::Iso_cuboid_3& x,
                  const Kernel::Iso_cuboid_3& y);

  /// @}

}; /* end Kernel::Equal_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`
  \sa `::Kernel::HasOnBoundedSide_2`
  \sa `::Kernel::HasOnUnboundedSide_2`
  \sa `::Kernel::BoundedSide_2`

*/
class HasOnBoundary_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the boundary of `c`.
  */
  bool operator()(const Kernel::Circle_2&c,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the boundary of `i`.
  */
  bool operator()(const Kernel::Iso_rectangle_2&i,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the boundary of `t`.
  */
  bool operator()(const Kernel::Triangle_2&t,
                  const Kernel::Point_2&p);


  /// @}

}; /* end Kernel::HasOnBoundary_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`
  \sa `::Kernel::HasOnBoundedSide_3`
  \sa `::Kernel::HasOnUnboundedSide_3`
  \sa `::Kernel::BoundedSide_3`

*/
class HasOnBoundary_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the boundary of `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the boundary of `t`.
  */
  bool operator()(const Kernel::Tetrahedron_3&t,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the boundary of `c`.
  */
  bool operator()(const Kernel::Iso_cuboid_3&c,
                  const Kernel::Point_3&p);

  /// @}

}; /* end Kernel::HasOnBoundary_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`
  \sa `::Kernel::HasOnUnboundedSide_2`
  \sa `::Kernel::HasOnBoundary_2`
  \sa `::Kernel::BoundedSide_2`

*/
class HasOnBoundedSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the bounded side of `c`.
  */
  bool operator()(const Kernel::Circle_2&c,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the bounded side of `i`.
  */
  bool operator()(const Kernel::Iso_rectangle_2&i,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the bounded side of `t`.
  */
  bool operator()(const Kernel::Triangle_2&t,
                  const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::HasOnBoundedSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`
  \sa `::Kernel::HasOnUnboundedSide_3`
  \sa `::Kernel::HasOnBoundary_3`
  \sa `::Kernel::BoundedSide_3`

*/
class HasOnBoundedSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the bounded side of `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the bounded side of `t`.
  */
  bool operator()(const Kernel::Tetrahedron_3&t,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the bounded side of `c`.
  */
  bool operator()(const Kernel::Iso_cuboid_3&c,
                  const Kernel::Point_3&p);

  /*!
    returns true iff the line segment `ab` is inside the union of the
    bounded sides of `s1` and `s2`.
  */
  bool operator()(const Kernel::Sphere_3& s1,
                  const Kernel::Sphere_3& s2,
                  const Kernel::Point_3& a,
                  const Kernel::Point_3& b);

  /// @}

}; /* end Kernel::HasOnBoundedSide_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`

*/
class HasOnNegativeSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff `p` lies on the negative side of `c`.
  */
  bool operator()(const Kernel::Circle_2&c,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the negative side of `l`
    (`l` is considered a half-space).
  */
  bool operator()(const Kernel::Line_2&l,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the negative side of `t`.
  */
  bool operator()(const Kernel::Triangle_2&t,
                  const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::HasOnNegativeSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class HasOnNegativeSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the negative side of `h`
    (`h` is considered a half-space).
  */
  bool operator()(const Kernel::Plane_3&h,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the negative side of `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the negative side of `t`.
  */
  bool operator()(const Kernel::Tetrahedron_3&t,
                  const Kernel::Point_3&p);

  /// @}

}; /* end Kernel::HasOnNegativeSide_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`

*/
class HasOnPositiveSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the positive side of `c`.
  */
  bool operator()(const Kernel::Circle_2&c,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the positive side of `l`
    (`l` is considered a half-space).
  */
  bool operator()(const Kernel::Line_2&l,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the positive side of `t`.
  */
  bool operator()(const Kernel::Triangle_2&t,
                  const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::HasOnPositiveSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class HasOnPositiveSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the positive side of `h`
    (`h` is considered a half-space).
  */
  bool operator()(const Kernel::Plane_3&h,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the positive side of `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the positive side of `t`.
  */
  bool operator()(const Kernel::Tetrahedron_3&t,
                  const Kernel::Point_3&p);

  /// @}

}; /* end Kernel::HasOnPositiveSide_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`
  \sa `::Kernel::HasOnBoundedSide_2`
  \sa `::Kernel::HasOnBoundary_2`
  \sa `::Kernel::BoundedSide_2`

*/
class HasOnUnboundedSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff `p` lies on the unbounded side of `c`.
  */
  bool operator()(const Kernel::Circle_2&c,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the unbounded side of `i`.
  */
  bool operator()(const Kernel::Iso_rectangle_2&i,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on the unbounded side of `t`.
  */
  bool operator()(const Kernel::Triangle_2&t,
                  const Kernel::Point_2&p);


  /// @}

}; /* end Kernel::HasOnUnboundedSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`
  \sa `::Kernel::HasOnBoundedSide_3`
  \sa `::Kernel::HasOnBoundary_3`
  \sa `::Kernel::BoundedSide_3`

*/
class HasOnUnboundedSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on the unbounded side of `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the unbounded side of `t`.
  */
  bool operator()(const Kernel::Tetrahedron_3&t,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on the unbounded side of `c`.
  */
  bool operator()(const Kernel::Iso_cuboid_3&c,
                  const Kernel::Point_3&p);


  /// @}

}; /* end Kernel::HasOnUnboundedSide_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class HasOn_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff `p` lies on `l`.
  */
  bool operator()(const Kernel::Line_2&l,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on `r`.
  */
  bool operator()(const Kernel::Ray_2&r,
                  const Kernel::Point_2&p);

  /*!
    returns true iff `p` lies on `s`.
  */
  bool operator()(const Kernel::Segment_2&s,
                  const Kernel::Point_2&p);

  /// @}

}; /* end Kernel::HasOn_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_3<Kernel>`
  \sa `CGAL::Line_3<Kernel>`
  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Point_3<Kernel>`
  \sa `CGAL::Ray_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Triangle_3<Kernel>`

*/
class HasOn_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `p` lies on `c`.
  */
  bool operator()(const Kernel::Circle_3&c,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on `l`.
  */
  bool operator()(const Kernel::Line_3&l,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on `r`.
  */
  bool operator()(const Kernel::Ray_3&r,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on `s`.
  */
  bool operator()(const Kernel::Segment_3&s,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `p` lies on `pl`.
  */
  bool operator()(const Kernel::Plane_3&pl,
                  const Kernel::Point_3&p);

  /*!
    returns true iff `l` lies on `pl`.
  */
  bool operator()(const Kernel::Plane_3&pl,
                  const Kernel::Line_3&l);

  /*!
    returns true iff `c` lies on `pl`.
  */
  bool operator()(const Kernel::Plane_3&pl,
                  const Kernel::Circle_3&c);

  /*!
    returns true iff `c` lies on `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Point_3&c);

  /*!
    returns true iff `c` lies on `s`.
  */
  bool operator()(const Kernel::Sphere_3&s,
                  const Kernel::Circle_3&c);

  /*!
    returns true iff `p` lies on `t`.
  */
  bool operator()(const Kernel::Triangle_3&t,
                  const Kernel::Point_3&p);

  /// @}

}; /* end Kernel::HasOn_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa \link intersection_grp `CGAL::intersection()` \endlink
*/
class Intersect_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the intersection region of two geometrical objects of type
    `Type1` and `Type2`, for all pairs `Type1` and `Type2`.
    For details see the reference manual page for \link intersection_grp `CGAL::intersection()` \endlink.
  */
  decltype(auto)
  operator()(Type1 obj1, Type2 obj2);

  /// @}

}; /* end Kernel::Intersect_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunctor, AdaptableTernaryFunctor}

  \sa intersection_linear_grp
*/
class Intersect_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    computes the intersection region of two geometrical
    objects of type `Type1` and `Type2`.
    For details see the reference manual page for \ref intersection_linear_grp.
  */
  decltype(auto)
  operator()(Type1 obj1, Type2 obj2);



  /// @}

}; /* end Kernel::Intersect_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Iso_rectangle_2<Kernel>`
  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`

*/
class IsDegenerate_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Circle_2&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Iso_rectangle_2&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Line_2&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Ray_2&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Segment_2&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Triangle_2&o);

  /// @}

}; /* end Kernel::IsDegenerate_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Circle_3<Kernel>`
  \sa `CGAL::Iso_cuboid_3<Kernel>`
  \sa `CGAL::Line_3<Kernel>`
  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Point_3<Kernel>`
  \sa `CGAL::Ray_3<Kernel>`
  \sa `CGAL::Segment_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`
  \sa `CGAL::Triangle_3<Kernel>`

*/
class IsDegenerate_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Circle_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Iso_cuboid_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Line_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Plane_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Ray_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Segment_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Sphere_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Tetrahedron_3&o);

  /*!
    returns true iff `o` is degenerate.
  */
  bool operator()(const Kernel::Triangle_3&o);

  /// @}

}; /* end Kernel::IsDegenerate_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class IsHorizontal_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff `o` is horizontal.
  */
  bool operator()(const Kernel::Line_2&o);

  /*!
    returns true iff `o` is horizontal.
  */
  bool operator()(const Kernel::Ray_2&o);

  /*!
    returns true iff `o` is horizontal.
  */
  bool operator()(const Kernel::Segment_2&o);

  /// @}

}; /* end Kernel::IsHorizontal_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}

  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Ray_2<Kernel>`
  \sa `CGAL::Segment_2<Kernel>`

*/
class IsVertical_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff `o` is vertical.
  */
  bool operator()(const Kernel::Line_2&o);

  /*!
    returns true iff `o` is vertical.
  */
  bool operator()(const Kernel::Ray_2&o);

  /*!
    returns true iff `o` is vertical.
  */
  bool operator()(const Kernel::Segment_2&o);


  /// @}

}; /* end Kernel::IsVertical_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `left_turn_grp`

*/
class LeftTurn_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns `true`, iff the three points `p`, `q`
    and `r` form a left turn.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);

  /// @}

}; /* end Kernel::LeftTurn_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `has_smaller_distance_to_point_grp`

*/
class LessDistanceToPoint_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the distance of `q` to `p` is
    smaller than the distance of `r` to `p`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);


  /// @}

}; /* end Kernel::LessDistanceToPoint_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `has_smaller_distance_to_point_grp`

*/
class LessDistanceToPoint_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the distance of `q` to `p` is
    smaller than the distance of `r` to `p`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q,
                  const Kernel::Point_3&r);

  /// @}

}; /* end Kernel::LessDistanceToPoint_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

*/
class LessRotateCCW_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff the three points `p`, `q`
    and `r` form a left turn or if they are collinear and
    the distance of `q` to `p` is
    larger than the distance of `r` to `p`, where `p` is the point
    passed to the object at construction.
    \pre `p` does not lie in the interior of the segment `rq`, i.e.\ `p` is an extreme point with respect to \f$ \{p,q,r\}\f$.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q,
                  const Kernel::Point_2&r);

  /// @}

}; /* end Kernel::LessRotateCCW_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `has_smaller_signed_distance_to_line_grp`

*/
class LessSignedDistanceToLine_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!

    returns `true` if the signed distance from `p` and the oriented line `l`
    is smaller than the signed distance of `q` and `l`.

  */
  bool operator()(const Kernel::Line_2& l,
                  const Kernel::Point_2& p,
                  const Kernel::Point_2& q);

  /*!

    returns `true` if the signed distance from `r` and the oriented line `l`
    defined by `p` and `q` is smaller than the signed distance of `s` and `l`.
    \pre `p != q`.

  */
  bool operator()(const Kernel::Point_2& p,
                  const Kernel::Point_2& q,
                  const Kernel::Point_2&r,
                  const Kernel::Point_2&s);

  /// @}

}; /* end Kernel::LessSignedDistanceToLine_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `has_smaller_signed_distance_to_plane_grp`

*/
class LessSignedDistanceToPlane_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true, iff the signed distance from point `q` to plane
    `p` is smaller than the signed distance from point `r` to `p`.
  */
  bool operator()(const Kernel::Plane_3& p,
                  const Kernel::Point_3& q,
                  const Kernel::Point_3& r);

  /*!
    returns true, iff the signed distance from point `q` to the plane
    `p` defined by `p1, p2, p3` is smaller than the signed distance
    from point `r` to `p`.
    \pre `p, q`, and `r` are not collinear.
  */
  bool operator()(const Kernel::Point_3& p1,
                  const Kernel::Point_3& p2,
                  const Kernel::Point_3& p3,
                  const Kernel::Point_3& q,
                  const Kernel::Point_3& r);

  /// @}

}; /* end Kernel::LessSignedDistanceToPlane_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `lexicographically_xyz_smaller_grp`

*/
class LessXYZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ x\f$-coordinate of `p` is smaller than the
    \f$ x\f$-coordinate of `q` or if they are the same and
    the \f$ y\f$-coordinate of `p` is smaller than the \f$ y\f$-coordinate of `q`, or,
    if both \f$ x\f$- and \f$ y\f$- coordinate are identical and
    the \f$ z\f$-coordinate of `p` is smaller than the \f$ z\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::LessXYZ_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `lexicographically_xy_smaller_grp`

*/
class LessXY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ x\f$-coordinate of `p` is smaller than the
    \f$ x\f$-coordinate of `q` or if they are the same and
    the \f$ y\f$-coordinate of `p` is smaller than the \f$ y\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);


  /// @}

}; /* end Kernel::LessXY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_xy_grp`

*/
class LessXY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ x\f$-coordinate of `p` is smaller than the
    \f$ x\f$-coordinate of `q` or if they are the same and
    the \f$ y\f$-coordinate of `p` is smaller than the \f$ y\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::LessXY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_x_grp`

*/
class LessX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns true iff the \f$ x\f$-coordinate of `p` is smaller than the
    \f$ x\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::LessX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_x_grp`

*/
class LessX_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ x\f$-coordinate of `p` is smaller than the
    \f$ x\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::LessX_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_yx_grp`

*/
class LessYX_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ y\f$-coordinate of `p` is smaller than the
    \f$ y\f$-coordinate of `q` or if they are the same and
    the \f$ x\f$-coordinate of `p` is smaller than the \f$ x\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);

  /// @}

}; /* end Kernel::LessYX_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_y_grp`

*/
class LessY_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ y\f$-coordinate of `p` is smaller than the
    \f$ y\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_2&p,
                  const Kernel::Point_2&q);

  ///@}

}; /* end Kernel::LessY_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `compare_y_grp`

*/
class LessY_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ y\f$-coordinate of `p` is smaller than the
    \f$ y\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::LessY_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \sa `compare_z_grp`

*/
class LessZ_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns true iff the \f$ z\f$-coordinate of `p` is smaller than the
    \f$ z\f$-coordinate of `q`.
  */
  bool operator()(const Kernel::Point_3&p,
                  const Kernel::Point_3&q);

  /// @}

}; /* end Kernel::LessZ_3 */




/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableUnaryFunction}
*/
class NonZeroCoordinateIndex_3
{
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns any of `0`, `1`, or `2` if the corresponding coordinate of the vector `v` is not
    equal to zero, and `-1` if `v` is the null vector.
  */
  int operator()(const Kernel::Vector_3& v);

  /// @}

};


/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableTernaryFunction}

  \sa `orientation_grp`

*/
class Orientation_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns \ref CGAL::LEFT_TURN, if `r` lies to the left of the oriented
    line `l` defined by `p` and `q`, returns \ref CGAL::RIGHT_TURN if `r`
    lies to the right of `l`, and returns \ref CGAL::COLLINEAR if `r` lies
    on `l`.
  */
  Orientation operator()(const Kernel::Point_2&p,
                         const Kernel::Point_2&q,
                         const Kernel::Point_2&r);

  /*!
    returns \ref CGAL::LEFT_TURN if `u` and `v` form a left turn,
    returns \ref CGAL::RIGHT_TURN if `u` and `v` form a right turn,
    and returns \ref CGAL::COLLINEAR if `u` and `v` are collinear.
    */
  Orientation operator()(const Kernel::Vector_2&u,
                         const Kernel::Vector_2&v);


  /// @}

}; /* end Kernel::Orientation_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `orientation_grp`

*/
class Orientation_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns \ref CGAL::POSITIVE, if `s` lies on the positive side of the oriented
    plane `h` defined by `p`, `q`, and `r`, returns \ref CGAL::NEGATIVE if `s`
    lies on the negative side of `h`, and returns \ref CGAL::COPLANAR if `s` lies
    on `h`.
  */
  Orientation operator()(const Kernel::Point_3&p,
                         const Kernel::Point_3&q,
                         const Kernel::Point_3&r,
                         const Kernel::Point_3&s);

  /*!
    returns \ref CGAL::POSITIVE if `u`, `v` and `w` are positively oriented,
    returns \ref CGAL::NEGATIVE if `u`, `v` and `w` are negatively oriented,
    and returns \ref CGAL::COPLANAR if `u`, `v` and `w` are coplanar.
  */
  Orientation operator()(const Kernel::Vector_3&u,
                         const Kernel::Vector_3&v,
                         const Kernel::Vector_3&w);

  /*!
    returns the orientation of the sphere `s`.
  */
  Orientation operator()(const Kernel::Sphere_3&s);


  /// @}

}; /* end Kernel::Orientation_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Circle_2<Kernel>`
  \sa `CGAL::Line_2<Kernel>`
  \sa `CGAL::Triangle_2<Kernel>`

*/
class OrientedSide_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or the constant \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented circle `c`.
  */
  Oriented_side operator()(const Kernel::Circle_2&c,
                           const Kernel::Point_2&p);

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or the constant \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented line `l`.
  */
  Oriented_side operator()(const Kernel::Line_2&l,
                           const Kernel::Point_2&p);

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or the constant \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented triangle `t`.
  */
  Oriented_side operator()(const Kernel::Triangle_2&t,
                           const Kernel::Point_2&p);

  /*!
  * returns \ref CGAL::ON_ORIENTED_BOUNDARY,
  * \ref CGAL::ON_NEGATIVE_SIDE, or the constant \ref CGAL::ON_POSITIVE_SIDE,
  * depending on the position of the circumcenter of `t` relative
  * to the oriented supporting line of `s`. The orientation of the
  * supporting line is the same as the orientation of `s`.
  */
  Oriented_side operator()(const Kernel::Segment_2& s,
                           const Kernel::Triangle_2& t);

  /// @}

}; /* end Kernel::OrientedSide_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableBinaryFunction}

  \sa `CGAL::Plane_3<Kernel>`
  \sa `CGAL::Sphere_3<Kernel>`
  \sa `CGAL::Tetrahedron_3<Kernel>`

*/
class OrientedSide_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented plane `h`.
  */
  Oriented_side operator()(const Kernel::Plane_3&h,
                           const Kernel::Point_3&p);

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented plane constructed
    from `plane_point` and `plane_vector`.
  */
  Oriented_side operator()(const Kernel::Point_3& plane_point,
                           const Kernel::Vector_3& plane_vector,
                           const Kernel::Point_3&p);

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented tetrahedron `t`.
  */
  Oriented_side operator()(const Kernel::Tetrahedron_3&t,
                           const Kernel::Point_3&p);

  /*!
    returns \ref CGAL::ON_ORIENTED_BOUNDARY,
    \ref CGAL::ON_NEGATIVE_SIDE, or \ref CGAL::ON_POSITIVE_SIDE,
    depending on the position of `p` relative to the oriented sphere `s`.
  */
  Oriented_side operator()(const Kernel::Sphere_3& s,
                           const Kernel::Point_3& p);


  /// @}

}; /* end Kernel::OrientedSide_3 */


/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableQuaternaryFunction}

\sa `CGAL::Weighted_point_2<Kernel>`
\sa `ComputePowerProduct_2` for the definition of orthogonality for power distances.
\sa `PowerSideOfOrientedPowerCircle_2`

*/
class PowerSideOfBoundedPowerCircle_2
{
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Let \f$ {z(p,q,r)}^{(w)}\f$ be the power circle of the weighted points
    \f$ (p,q,r)\f$. This method returns:

    - `ON_BOUNDARY` if `t` is orthogonal to \f$ {z(p,q,r)}^{(w)}\f$,

    - `ON_UNBOUNDED_SIDE` if `t` lies outside the bounded circle of
    center \f$ z(p,q,r)\f$ and radius \f$ \sqrt{ w_{z(p,q,r)}^2 + w_t^2 }\f$
    (which is equivalent to \f$ \Pi({t}^{(w)},{z(p,q,r)}^{(w)}) > 0\f$),

    - `ON_BOUNDED_SIDE` if `t` lies inside this bounded circle.

    The order of the points `p`, `q`, and `r` does not matter.

    \pre `p`, `q`, and `r` are not collinear.

    If all the points have a weight equal to 0, then
    `power_side_of_bounded_power_circle_2(p,q,r,t)` ==
      `side_of_bounded_circle(p,q,r,t)`.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_2 & p,
             const Kernel::Weighted_point_2 & q,
             const Kernel::Weighted_point_2 & r,
             const Kernel::Weighted_point_2 & t);

  /*!
    returns the sign of the power test of `t` with respect
    to the smallest circle orthogonal to `p` and `q`.

    \pre `p` and `q` have different bare points.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_2 & p,
             const Kernel::Weighted_point_2 & q,
             const Kernel::Weighted_point_2 & t);

  /*!
    returns the sign of the power test of `t` with respect
    to the smallest circle orthogonal to `p`.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_2 & p,
             const Kernel::Weighted_point_2 & t);
  /// @}
}; /* end Kernel::PowerSideOfBoundedPowerCircle_2 */

/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableQuinaryFunction}

\sa `CGAL::Weighted_point_3<Kernel>`
\sa `ComputePowerProduct_3` for the definition of orthogonality for power distances.
\sa `PowerSideOfOrientedPowerSphere_3`

*/
class PowerSideOfBoundedPowerSphere_3
{
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Let \f$ {z(p,q,r,s)}^{(w)}\f$ be the power sphere of the weighted points
    \f$ (p,q,r,s)\f$. This method returns:

    - `ON_BOUNDARY` if `t` is orthogonal to
    \f$ {z(p,q,r,s)}^{(w)}\f$,

    - `ON_UNBOUNDED_SIDE` if `t` lies outside the bounded sphere of
    center \f$ z(p,q,r,s)\f$ and radius \f$ \sqrt{ w_{z(p,q,r,s)}^2 + w_t^2 }\f$
    (which is equivalent to \f$ \Pi({t}^{(w)},{z(p,q,r,s)}^{(w)}) >0\f$),

    - `ON_BOUNDED_SIDE` if `t` lies inside this bounded sphere.

    The order of the points `p`, `q`, `r`, and `s` does not matter.

    \pre `p, q, r, s` are not coplanar.

    If all the points have a weight equal to 0, then
    `power_side_of_bounded_power_sphere_3(p,q,r,s,t)` ==
      `side_of_bounded_sphere(p,q,r,s,t)`.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_3 & p,
             const Kernel::Weighted_point_3 & q,
             const Kernel::Weighted_point_3 & r,
             const Kernel::Weighted_point_3 & s,
             const Kernel::Weighted_point_3 & t);

  /*!
    returns the sign of the power test of `t` with respect
    to the smallest sphere orthogonal to `p`, `q`, and `r`.

    \pre `p, q, r` are not collinear.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_3 & p,
             const Kernel::Weighted_point_3 & q,
             const Kernel::Weighted_point_3 & r,
             const Kernel::Weighted_point_3 & t);

  /*!
    returns the sign of the power test of `t` with respect
    to the smallest sphere orthogonal to `p` and `q`.

    \pre `p` and `q` have different bare points.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_3 & p,
             const Kernel::Weighted_point_3 & q,
             const Kernel::Weighted_point_3 & t);

  /*!
    returns the sign of the power test of `t` with respect
    to the smallest sphere orthogonal to `p`.
  */
  CGAL::Bounded_side
  operator()(const Kernel::Weighted_point_3 & p,
             const Kernel::Weighted_point_3 & t);
  /// @}

}; /* end Kernel::PowerSideOfBoundedPowerSphere_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `CGAL::Weighted_point_2<Kernel>`
  \sa `ComputePowerProduct_2` for the definition of power distance.
  \sa `PowerSideOfBoundedPowerCircle_2`
*/
class PowerSideOfOrientedPowerCircle_2
{
  public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the relative position of point `s` to the oriented power circle
    defined by `p`, `q`, and `r`.

    The order of the points `p`, `q` and `r` is important,
    since it determines the orientation of the implicitly
    constructed power circle.

    \pre the bare points corresponding to `p`, `q`, `r` are not collinear.
  */

  Oriented_side operator()(const Kernel::Weighted_point_2& p,
                           const Kernel::Weighted_point_2& q,
                           const Kernel::Weighted_point_2& r,
                           const Kernel::Weighted_point_2& s);
  /// @}

};

/*!
\ingroup PkgKernel23ConceptsFunctionObjects
\cgalConcept

\cgalRefines{AdaptableQuinaryFunction}

\sa `CGAL::Weighted_point_3<Kernel>`
\sa `ComputePowerProduct_3` for the definition of power distance.
\sa `PowerSideOfBoundedPowerSphere_3`

*/
class PowerSideOfOrientedPowerSphere_3
{
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    Let \f$ {z(p,q,r,s)}^{(w)}\f$ be the power sphere of the weighted points
    \f$ (p,q,r,s)\f$. Returns

    - `ON_ORIENTED_BOUNDARY` if `t` is orthogonal to
      \f$ {z(p,q,r,s)}^{(w)}\f$,

    - `ON_NEGATIVE_SIDE` if `t` lies outside the oriented sphere of
      center \f$ z(p,q,r,s)\f$ and radius \f$ \sqrt{ w_{z(p,q,r,s)}^2 + w_t^2 }\f$
      (which is equivalent to \f$ \Pi({t}^{(w)},{z(p,q,r,s)}^{(w)}) > 0 \f$),

    - `ON_POSITIVE_SIDE` if `t` lies inside this oriented sphere.

    The order of the points `p`, `q`, `r` and `s` is important,
    since it determines the orientation of the implicitly
    constructed power sphere.

    \pre `p, q, r, s` are not coplanar.

    If all the points have a weight equal to 0, then
    `power_side_of_oriented_power_sphere_3(p,q,r,s,t)` =
       `side_of_oriented_sphere(p,q,r,s,t)`.
  */
  Oriented_side operator()( const Kernel::Weighted_point_3& p,
                            const Kernel::Weighted_point_3& q,
                            const Kernel::Weighted_point_3& r,
                            const Kernel::Weighted_point_3& s,
                            const Kernel::Weighted_point_3& t) const;
  /// @}
};

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `side_of_bounded_circle_grp`

*/
class SideOfBoundedCircle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the relative position of point `t`
    to the circle defined by `p`, `q` and `r`. The order
    of the points `p`, `q` and `r` does not matter.
    \pre `p, q` and `r` are not collinear.
  */
  Bounded_side operator()(const Kernel::Point_2&p,
                          const Kernel::Point_2&q,
                          const Kernel::Point_2&r,
                          const Kernel::Point_2&t);

  /*!
    returns the position of the point `t` relative to the circle
    that has `pq` as its diameter.
  */
  Bounded_side operator()(const Kernel::Point_2&p,
                          const Kernel::Point_2&q,
                          const Kernel::Point_2&t);

  /// @}

}; /* end Kernel::SideOfBoundedCircle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuinaryFunction}

  \sa `side_of_bounded_sphere_grp`

*/
class SideOfBoundedSphere_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the relative position of point `t`
    to the sphere defined by `p`, `q`, `r`, and `s`. The order
    of the points `p`, `q`, `r`, and `s` does not matter.
    \pre `p, q, r` and `s` are not coplanar.
  */
  Bounded_side operator()(const Kernel::Point_3&p,
                          const Kernel::Point_3&q,
                          const Kernel::Point_3&r,
                          const Kernel::Point_3&s,
                          const Kernel::Point_3&t);

  /*!
    returns the position of the point `t` relative to the sphere
    passing through `p`, `q`, and `r` and whose center is in the plane defined
    by these three points.
  */
  Bounded_side operator()(const Kernel::Point_3&p,
                          const Kernel::Point_3&q,
                          const Kernel::Point_3&r,
                          const Kernel::Point_3&t);

  /*!
    returns the position of the point `t` relative to the sphere
    that has `pq` as its diameter.
  */
  Bounded_side operator()(const Kernel::Point_3&p,
                          const Kernel::Point_3&q,
                          const Kernel::Point_3&t);


  /// @}

}; /* end Kernel::SideOfBoundedSphere_3 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuaternaryFunction}

  \sa `side_of_oriented_circle_grp`

*/
class SideOfOrientedCircle_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{


  /*!
    returns the relative position of point `t`
    to the oriented circle defined by `p`, `q` and `r`.
    The order of the points `p`, `q` and `r` is important,
    since it determines the orientation of the implicitly
    constructed circle.

    If `p`, `q` and `r` are collinear, the circle degenerates in a line.
    \ref CGAL::ON_ORIENTED_BOUNDARY is returned if `t` is also collinear or if two
    points are identical,
    otherwise, `side_of_oriented_circle(r, q, t, p)` is returned.

  */
  Oriented_side operator()(const Kernel::Point_2&p,
                           const Kernel::Point_2&q,
                           const Kernel::Point_2&r,
                           const Kernel::Point_2&t);


  /// @}

}; /* end Kernel::SideOfOrientedCircle_2 */

/*!
  \ingroup PkgKernel23ConceptsFunctionObjects
  \cgalConcept

  \cgalRefines{AdaptableQuinaryFunction}

  \sa `side_of_oriented_sphere_grp`

*/
class SideOfOrientedSphere_3 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*!
    returns the relative position of point `t`
    to the oriented sphere defined by `p`, `q`, `r` and `s`.
    The order of the points `p`, `q`, `r`, and `s` is important,
    since it determines the orientation of the implicitly
    constructed sphere. If the points `p`, `q`, `r` and `s`
    are positive oriented, positive side is the bounded interior
    of the sphere.

    In case of degeneracies, \ref CGAL::ON_ORIENTED_BOUNDARY is returned
    if all points are coplanar. Otherwise, there is a cyclic permutation of the five points
    that puts four non coplanar points first, it is used to answer the predicate:
    e.g. `side_of_oriented_sphere(q, r, s, t, p)` is returned if `q`, `r`, `s`,
    and `t` are non coplanar.
  */
  Oriented_side operator()(const Kernel::Point_3&p,
                           const Kernel::Point_3&q,
                           const Kernel::Point_3&r,
                           const Kernel::Point_3&s,
                           const Kernel::Point_3&t);

  /// @}

}; /* end Kernel::SideOfOrientedSphere_3 */

} // end of Kernel namespace
