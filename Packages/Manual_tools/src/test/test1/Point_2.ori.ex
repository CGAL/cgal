
template < class R >
class CGAL_Point_2 {
public:


// SECTION: 2D Point
// ========================================================================
// 
// DEFINITION An object of the class CGAL_Point_2<R> is a point in the
// two-dimensional euclidean plane |E_2.
// 
// Remember that `R::RT' and `R::FT' denote a ring type and a field type.
// For the representation class `CGAL_Cartesian<T>' the two types are
// equivalent. For the representation class `CGAL_Homogeneous<T>' the ring
// type is `R::RT' == `T' and the field type is `R::FT' ==
// `CGAL_Quotient<T>'.
// 
// CREATION
// 
// New creation variable is: `p'
// 
// `#include <CGAL/Point_2.h>'

    CGAL_Point_2();
        // introduces an uninitialized variable `p'.

    CGAL_Point_2(const CGAL_Point_2<R> &q);
        // copy constructor.

    CGAL_Point_2(const R::RT &hx, const R::RT &hy, const R::RT &hw = R::RT(1));
        // introduces a point `p' initialized to (hx/hw,hy/hw). If the
        // third argument is not explicitely given it defaults to `R::RT(1
        // )'.

    CGAL_Point_2(const R::RT &hx, const R::RT &hy, const R::RT &hw);

    CGAL_Point_2(const R::RT &hx, const R::RT &hy);

// OPERATIONS

    CGAL_Point_2<R> & operator=(const CGAL_Point_2<R> &q);
        // Assignment.

    bool operator==(const CGAL_Point_2<R> &q) const;
        // Test for equality: Two points are equal, iff their x and y
        // coodinates are equal.

    bool operator!=(const CGAL_Point_2<R> &q) const;
        // Test for inequality.

// There are two sets of coordinate access functions, namely to the
// homogeneous and to the Cartesian coordinates. They can be used
// independently from the chosen representation type `R'.

    R::RT hx() const;
        // returns the homogeneous x coordinate.

    R::RT hy() const;
        // returns the homogeneous y coordinate.

    R::RT hw() const;
        // returns the homogenizing coordinate.

// Here come the Cartesian access functions. Note that you do not loose
// information with the homogeneous representation, because then the field
// type is a quotient.

    R::FT x() const;
        // returns the Cartesian x coordinate, that is hx/hw.

    R::FT y() const;
        // returns the Cartesian y coordinate, that is hy/hw.

// The following operations are for convenience and for making this point
// class compatible with code for higher dimensional points. Again they
// come in a Cartesian and homogeneous flavor.

    R::RT homogeneous(int i) const;
        // returns the i'th homogeneous coordinate of `p', starting with
        // 0. Precondition: 0<= i <= 2.

    R::FT cartesian(int i) const;
        // returns the i'th Cartesian coordinate of `p', starting with 0.
        // Precondition: 0<= i <= 1.

    R::FT operator[](int i) const;
        // returns `cartesian(i)'. Precondition: 0<= i <= 1.

    int dimension() const;
        // returns the dimension (the constant 2).

    CGAL_Bbox_2 bbox() const;
        // returns a bounding box containing `p'. Note that bounding boxes
        // are not parameterized with whatsoever.

    CGAL_Point_2<R> transform(const CGAL_Aff_transformation_2<R> &t) const;
        // returns the point obtained by applying t on `p'.

// The following operations can be applied on points:

    CGAL_Vector_2<R> operator-(const CGAL_Point_2<R> &p,  const CGAL_Point_2<R> &q);
        // returns the difference vector between `q' and `p'.

    CGAL_Point_2<R> operator+(const CGAL_Point_2<R> &p,  const CGAL_Vector_2<R> &v);
        // returns a point obtained by translating `p' by the vector `v'.

    CGAL_Point_2<R> operator-(const CGAL_Point_2<R> &p,  const CGAL_Vector_2<R> &v);
        // returns a point obtained by translating `p' by the vector `-v'.

// EXAMPLE
// 
// The following declaration creates two points with Cartesian double
// coordinates.
// 
// The variable p is uninitialized and should first be used on the left
// hand side of an assignment.
};

