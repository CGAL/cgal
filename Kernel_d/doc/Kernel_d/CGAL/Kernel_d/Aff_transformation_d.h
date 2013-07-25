namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of the data type `Aff_transformation_d<Kernel>` is an 
affine transformation of \f$ d\f$-dimensional space. It is specified by a 
square matrix \f$ M\f$ of dimension \f$ d + 1\f$. All entries in the last row of 
`M` except the diagonal entry must be zero; the diagonal entry 
must be non-zero. A point \f$ p\f$ with homogeneous coordinates \f$ (p[0], 
\ldots, p[d])\f$ can be transformed into the point `p.transform(A)` 
= \f$ Mp\f$, where `A` is an affine transformation created from `M` 
by the constructors below. 

\cgalHeading{Implementation}

Affine Transformations are implemented by matrices of number type 
`RT` as a handle type. All operations like creation, 
initialization, input and output on a transformation \f$ t\f$ take time 
\f$ O(t.dimension()^2)\f$. `dimension()` takes constant time. 
The operations for inversion and composition have the cubic costs of 
the used matrix operations. The space requirement is 
\f$ O(t.dimension()^2)\f$. 

*/
template< typename Kernel >
class Aff_transformation_d {
public:

/// \name Types 
/// @{

/*!
the linear algebra layer. 
*/ 
typedef unspecified_type LA; 

/*!
the matrix type. 
*/ 
typedef unspecified_type Matrix; 

/// @} 

/// \name Creation 
/// @{

/*!
introduces some 
transformation. 
*/ 
Aff_transformation_d<Kernel>(); 

/*!
introduces the identity transformation in 
\f$ d\f$-dimensional space. 
*/ 
Aff_transformation_d<Kernel>(int d, 
Identity_transformation); 

/*!
introduces the 
transformation of \f$ d\f$-space specified by matrix \f$ M\f$.

\pre `M` is a square matrix of dimension \f$ d + 1\f$ where entries in the last row of `M` except the diagonal entry must be zero; the diagonal entry must be non-zero. 
*/ 
Aff_transformation_d<Kernel>(Matrix M); 

/*!
introduces the transformation of \f$ d\f$-space 
specified by a diagonal matrix with entries `set [start,end)` on 
the diagonal (a scaling of the space).

\pre `set [start,end)` is a vector of dimension \f$ d+1\f$. 
*/ 
template <typename Forward_iterator> 
Aff_transformation_d<Kernel>(Scaling, Forward_iterator start, 
Forward_iterator end); 

/*!
introduces the translation by vector \f$ v\f$. 
*/ 
Aff_transformation_d<Kernel>(Translation, Vector_d<Kernel> 
v); 

/*!
returns a scaling by a scale factor `num/den`.

\pre `den != 0 `. 
*/ 
Aff_transformation_d<Kernel>(int d, Scaling, RT num, RT 
den); 

/*!
returns a planar rotation 
with sine and cosine values `sin_num/den` and `cos_num/den` 
in the plane spanned by the base vectors \f$ b_{e1}\f$ and \f$ b_{e2}\f$ in 
\f$ d\f$-space. Thus the default use delivers a planar rotation in the 
\f$ x\f$-\f$ y\f$ plane.

\pre \f$ sin_num^2 + cos_num^2 = den^2\f$ and \f$ 0 \leq e_1 < e_2 < d\f$.
\pre `den != 0`. 

*/ 
Aff_transformation_d<Kernel>(int d, Rotation, RT sin_num, RT 
cos_num, RT den, int e1 = 0, int e2 = 1); 

/*!
returns a planar 
rotation within a two-dimensional linear subspace. The subspace is 
spanned by the base vectors \f$ b_{e1}\f$ and \f$ b_{e2}\f$ in \f$ d\f$-space. The 
rotation parameters are given by the \f$ 2\f$-dimensional direction 
`dir`, such that the difference between the sines and cosines of 
the rotation given by `dir` and the approximated rotation are at 
most `num/den` each. 
\pre `dir.dimension() == 2`, `!dir.is_degenerate()` and `num < den` is positive, `den != 0`, \f$ 0 \leq e_1 < e_2 < d\f$. 
*/ 
Aff_transformation_d<Kernel>(int d, Rotation, Direction_d<Kernel> 
dir, RT num, RT den, int e1 = 0, int e2 = 1); 

/// @} 

/// \name Operations 
/// @{

/*!
the dimension of the underlying space 
*/ 
int dimension() ; 

/*!
returns the transformation matrix 

*/ 
const Matrix& matrix() ; 

/*!
returns the inverse 
transformation.

\pre `t.matrix()` is invertible. 
*/ 
Aff_transformation_d<Kernel> inverse() ; 

/*!
composition of transformations. Note 
that transformations are not necessarily commutative. `t*s` is 
the transformation which transforms first by `t` and then by 
`s`. 
*/ 
Aff_transformation_d<Kernel> operator*(const 
Aff_transformation_d<Kernel>& s) ; 

/// @}

}; /* end Aff_transformation_d */
} /* end namespace CGAL */
