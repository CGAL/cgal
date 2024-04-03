namespace CGAL {

/*!
\ingroup kernel_classes2

The class `Aff_transformation_2` represents two-dimensional affine transformations.
The general form of an affine transformation is based on a homogeneous
representation of points. Thereby all transformations can be realized by
matrix multiplications.

Multiplying the transformation matrix by a scalar does not change the
represented transformation. Therefore, any transformation represented
by a matrix with rational entries can be represented by a
transformation matrix with integer entries as well. (Multiply the
matrix with the common denominator of the rational entries.) Hence, it
is sufficient to use the number type `Kernel::RT` to represent
the entries of the transformation matrix.

\cgal offers several specialized affine transformations. Different
constructors are provided to create them. They are parameterized with
a symbolic name to denote the transformation type, followed by
additional parameters. The symbolic name tags solve ambiguities in the
function overloading and they make the code more readable, i.e., what
type of transformation is created.

Since two-dimensional points have three homogeneous coordinates, we
have a \f$ 3\times 3\f$ matrix \f$ {(m_{ij})}_{i,\,j=0\ldots 2}\f$.

If the homogeneous representations are normalized (the homogenizing
coordinate is 1), then the upper left \f$ 2\times 2\f$ matrix realizes
linear transformations. In the matrix form of a translation, the
translation vector \f$ (v_0,\,v_1,\,1)\f$ appears in the last column of the
matrix. The entries \f$ m_{20}\f$ and \f$ m_{21}\f$ are always zero and
therefore do not appear in the constructors.

\cgalModels{Hashable if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`}

\sa `Identity_transformation`
\sa `Rotation`
\sa `Scaling`
\sa `Translation`
\sa `Reflection`
\sa `rational_rotation_approximation_grp`

\cgalHeading{Example}

\code
typedef Cartesian<double> K;
typedef Aff_transformation_2<K> Transformation;
typedef Point_2<K> Point;
typedef Vector_2<K> Vector;
typedef Direction_2<K> Direction;

Transformation rotate(ROTATION, sin(pi), cos(pi));
Transformation rational_rotate(ROTATION,Direction(1,1), 1, 100);
Transformation translate(TRANSLATION, Vector(-2, 0));
Transformation scale(SCALING, 3);

Point q(0, 1);
q = rational_rotate(q);

Point p(1, 1);

p = rotate(p);

p = translate(p);

p = scale(p);
\endcode

The same would have been achieved with

\code
Transformation transform = scale * (translate * rotate);
p = transform(Point(1.0, 1.0));
\endcode

\sa `CGAL::Aff_transformation_3<Kernel>`
\sa `CGAL::Identity_transformation`
\sa `CGAL::Reflection`
\sa `CGAL::Rotation`
\sa `CGAL::Scaling`
\sa `CGAL::Translation`

*/
template< typename Kernel >
class Aff_transformation_2 {
public:

/// \name Creation
/// @{

/*!
introduces an identity transformation.
*/
Aff_transformation_2(const Identity_transformation& );

/*!
introduces a translation by a vector `v`.
*/
Aff_transformation_2(const Translation,
const Vector_2<Kernel> &v);

/*!
approximates the rotation over the angle indicated by direction
`d`, such that the differences between the sines and cosines
of the rotation given by d and the approximating rotation
are at most \f$ num/den\f$ each.
\pre `num/den > 0` and  `d != 0`.
*/
Aff_transformation_2(const Rotation,
const Direction_2<Kernel> &d,
const Kernel::RT &num,
const Kernel::RT &den = RT(1));

/*!
introduces a rotation by the angle `rho`.
\pre <tt>sine\_rho<sup>2</sup> + cosine\_rho<sup>2</sup> == hw<sup>2</sup></tt>.
*/
Aff_transformation_2(const Rotation,
const Kernel::RT &sine_rho,
const Kernel::RT &cosine_rho,
const Kernel::RT &hw = RT(1));

/*!
introduces a scaling by a scale factor \f$ s/hw\f$.
*/
Aff_transformation_2(const Scaling,
const Kernel::RT &s,
const Kernel::RT &hw = RT(1));

/*!
introduces a reflection by a line `l`.
*/
Aff_transformation_2(const Reflection,
const Line_2<Kernel>& l);

/*!
introduces a general affine transformation in the
\f$3 \times 3\f$ matrix form
\f$
\small \mbox{\( \left(\begin{array}{ccc}
      m_{00} & m_{01} & m_{02}\\
      m_{10} & m_{11} & m_{12}\\
      0     &  0     & hw
    \end{array}\right) \)}
\f$.

The sub-matrix
\f$1\over hw\f$ \f$\small \mbox{\( \left(\begin{array}{cc}
                 m_{00} & m_{01}\\
                 m_{10} & m_{11}
                \end{array}\right) \) }
\f$ contains the scaling and rotation
information, the vector
\f$ \small
  \left(
    \begin{array}{c}
      m_{02}\\
      m_{12}
    \end{array}
  \right)
\f$
contains the translational part of the transformation.
*/
Aff_transformation_2(
const Kernel::RT &m00, const Kernel::RT &m01, const Kernel::RT &m02,
const Kernel::RT &m10, const Kernel::RT &m11, const Kernel::RT &m12,
const Kernel::RT &hw = RT(1));

/*!
introduces a general linear transformation
\f$\small \mbox{\(\left(\begin{array}{ccc}
                 m_{00} & m_{01} & 0\\
                 m_{10} & m_{11} & 0\\
                  0     &  0     & hw
              \end{array}\right)\)}\f$
i.e.\ there is no translational part.
*/
Aff_transformation_2(
const Kernel::RT &m00, const Kernel::RT &m01,
const Kernel::RT &m10, const Kernel::RT &m11,
const Kernel::RT &hw = RT(1));

/// @}

/*!
\name Operations

The main thing to do with transformations is to apply them on
geometric objects. Each class `Class_2<Kernel>` representing a
geometric object has a member function:
\code
Class_2<Kernel> transform(Aff_transformation_2<Kernel> t).
\endcode

The transformation classes provide a member function `transform()` for
points, vectors, directions, and lines. The same functionality is also
available through `operator()` overloads.
*/
/// @{

/*!

*/
Point_2<Kernel> transform(const Point_2<Kernel> &p) const;

/*!

*/
Vector_2<Kernel> transform(const Vector_2<Kernel> &p) const;

/*!

*/
Direction_2<Kernel> transform(const Direction_2<Kernel> &p) const;

/*!

*/
Line_2<Kernel> transform(const Line_2<Kernel> &p) const;

/*!

*/
Point_2<Kernel> operator()(const Point_2<Kernel> &p) const;

/*!

*/
Vector_2<Kernel> operator()(const Vector_2<Kernel> &p) const;

/*!

*/
Direction_2<Kernel> operator()(const Direction_2<Kernel> &p) const;

/*!

*/
Line_2<Kernel> operator()(const Line_2<Kernel> &p) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
composes two affine transformations.
*/
Aff_transformation_2<Kernel> operator*(const Aff_transformation_2<Kernel> &s) const;

/*!
gives the inverse transformation.
*/
Aff_transformation_2<Kernel> inverse() const;

/*!
compares two affine transformations.
*/
bool operator==(const Aff_transformation_2<Kernel> &s) const;

/*!
returns `true`, if the transformation is not reflecting,
i.e.\ the determinant of the involved linear transformation is
non-negative.
*/
bool is_even() const;

/*!
returns `true`, if the transformation is reflecting.
*/
bool is_odd() const;

/*!
returns `true`, if the object was constructed using the tag `CGAL::Scaling`, or is the result of the composition of only such scaling transformation objects.
*/
bool is_scaling() const;

/*!
returns `true`, if the object was constructed using the tag `CGAL::Translation`, or is the result of the composition of only such translation transformation objects.
*/
bool is_translation() const;

/*!
returns `true`, if the object was constructed using the tag `CGAL::Rotation`, or is the result of the composition of only such rotation transformation objects.
*/
bool is_rotation() const;

/*!
returns `true`, if the object was constructed using the tag `CGAL::Reflection`, or is the result of the composition of only such reflection transformation objects.
*/
bool is_reflection() const;

/// @}

/// \name Matrix Entry Access
/// @{

/*!

*/
Kernel::FT cartesian(int i, int j) const;

/*!
returns entry \f$ m_{ij}\f$ in a matrix representation in which \f$ m_{22}\f$ is 1.
*/
Kernel::FT m(int i, int j) const;

/*!

*/
Kernel::RT homogeneous(int i, int j) const;

/*!
returns entry \f$ m_{ij}\f$ in some fixed matrix representation.
*/
Kernel::RT hm(int i, int j) const;

/// @}

}; /* end Aff_transformation_2 */
} /* end namespace CGAL */
