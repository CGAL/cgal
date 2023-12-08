namespace CGAL {

/*!
\ingroup kernel_classes3

The class `Aff_transformation_3` represents three-dimensional affine transformations.
The general form of an affine transformation is based on a homogeneous
representation of points. Thereby all transformations can be realized by
matrix multiplication.

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

In three-dimensional space we have a \f$ 4\times 4\f$ matrix
\f$ {(m_{ij})}_{i,\,j=0\ldots 3}\f$. Entries \f$ m_{30}\f$, \f$ m_{31}\f$, and
\f$ m_{32}\f$ are always zero and therefore do not appear in the
constructors.

\cgalModels{Hashable if `Kernel` is a cartesian kernel and if `Kernel::FT` is `Hashable`}

\sa `CGAL::Aff_transformation_2<Kernel>`
\sa `CGAL::Identity_transformation`
\sa `CGAL::Reflection`
\sa `CGAL::Rotation`
\sa `CGAL::Scaling`
\sa `CGAL::Translation`

*/
template< typename Kernel >
class Aff_transformation_3 {
public:

/// \name Creation
/// @{

/*!
introduces an identity transformation.
*/
Aff_transformation_3(const Identity_transformation& );

/*!
introduces a translation by a vector `v`.
*/
Aff_transformation_3(const Translation,
const Vector_3<Kernel> &v);

/*!
introduces a scaling by a scale factor \f$ s/hw\f$.
*/
Aff_transformation_3(const Scaling,
const Kernel::RT &s,
const Kernel::RT &hw = RT(1));

/*!
introduces a general affine transformation of the matrix
form
\f$
\small \mbox{\(\left(\begin{array}{cccc}
                 m_{00} & m_{01} & m_{02} & m_{03}\\
                 m_{10} & m_{11} & m_{12} & m_{13}\\
                 m_{20} & m_{21} & m_{22} & m_{23}\\
                  0     &  0     &      0 & hw
              \end{array}\right)\)}
\f$.

The part \f$1\over hw\f$ \f$ \small \mbox{\(\left(\begin{array}{ccc}
                 m_{00} & m_{01} & m_{02}\\
                 m_{10} & m_{11} & m_{12}\\
                 m_{20} & m_{21} & m_{22}\\
               \end{array}\right)\)} \f$ defines the scaling
and rotational part of the transformation,
while the vector \f$1\over hw\f$ \f$\small \mbox{\(\left(\begin{array}{c}
                 m_{03}\\
                 m_{13}\\
                 m_{23}
              \end{array}\right)\)}\f$ contains the translational part.
*/
Aff_transformation_3(
const Kernel::RT &m00, const Kernel::RT &m01, const Kernel::RT &m02, const Kernel::RT &m03,
const Kernel::RT &m10, const Kernel::RT &m11, const Kernel::RT &m12, const Kernel::RT &m13,
const Kernel::RT &m20, const Kernel::RT &m21, const Kernel::RT &m22, const Kernel::RT &m23,
const Kernel::RT &hw = RT(1));

/*!
introduces a general linear transformation of the
matrix form \f$
\small \mbox{\(\left(\begin{array}{cccc}
                 m_{00} & m_{01} & m_{02} & 0\\
                 m_{10} & m_{11} & m_{12} & 0\\
                 m_{20} & m_{21} & m_{22} & 0\\
                  0     &  0     &      0 & hw
              \end{array}\right)\)}
\f$, i.e.\ an affine  transformation without translational part.
*/
Aff_transformation_3(
const Kernel::RT &m00, const Kernel::RT &m01, const Kernel::RT& m02,
const Kernel::RT &m10, const Kernel::RT &m11, const Kernel::RT& m12,
const Kernel::RT &m20, const Kernel::RT &m21, const Kernel::RT& m22,
const Kernel::RT &hw = RT(1));

/// @}

/*!
\name Operations

The main thing to do with transformations is to apply them on
geometric objects. Each class `Class_3<Kernel>` representing a
geometric object has a member function:

\code
Class_3<Kernel> transform(Aff_transformation_3<Kernel> t)
\endcode

The transformation classes provide a member function `transform()` for
points, vectors, directions, and planes. The same functionality is also
available through `operator()` overloads.
*/
/// @{

/*!

*/
Point_3<Kernel> transform(const Point_3<Kernel> &p) const;

/*!

*/
Vector_3<Kernel> transform(const Vector_3<Kernel> &p) const;

/*!

*/
Direction_3<Kernel> transform(const Direction_3<Kernel> &p) const;

/*!

*/
Plane_3<Kernel> transform(const Plane_3<Kernel> &p) const;

/*!

*/
Point_3<Kernel> operator()(const Point_3<Kernel> &p) const;

/*!

*/
Vector_3<Kernel> operator()(const Vector_3<Kernel> &p) const;

/*!

*/
Direction_3<Kernel> operator()(const Direction_3<Kernel> &p) const;

/*!

*/
Plane_3<Kernel> operator()(const Plane_3<Kernel> &p) const;

/*!
composes two affine transformations.
*/
Aff_transformation_3<Kernel>
operator*(const Aff_transformation_3<Kernel> &s) const;

/*!
gives the inverse transformation.
*/
Aff_transformation_3<Kernel> inverse() const;

/*!
compares two affine transformations.
*/
bool operator==(const Aff_transformation_3<Kernel> &s) const;


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


/// @}

/// \name Matrix Entry Access
/// @{


/*!

*/
Kernel::FT cartesian(int i, int j) const;

/*!
returns entry \f$ m_{ij}\f$ in a matrix representation in which \f$ m_{33}\f$ is 1.
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

}; /* end Aff_transformation_3 */
} /* end namespace CGAL */
