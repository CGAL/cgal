
class Demo_Class1 {
public:


// SECTION: A Simple Class
// ========================================================================
// 
// Text\footnote{A footnote in a class environment.}.
// 
// New creation variable is: `p'
// 
// `#include< Demo_Class1.h>'

    Demo_Class1();
        // introduces a variable `p' initialized to the default. C++ code:
        // Demo_Class1();. Test ccStyle: `Underscore_within ccStyle'.

    Demo_Class1( const Demo_Class1 &);
        // copy constructor. C++ code: Demo_Class1(const Demo_Class1 &);

    FT x() const;
        // Cartesian x-coordinate\protect\footnotemark{}. C++ code: FT x()
        // const;

    const FT& y();
        // Cartesian y-coordinate. C++ code: const FT& y();

    Demo_Class1   transform( const CGAL_HAff_transformation<FT,RT> &t) const;  
        // Longish declarations forces the comment to start in the next
        // line.

    Demo_Class1   longish_function_name(   const CGAL_Aff_transformation<FT,RT> &t,  const Dummy_Type &q,  Long_Type_Name Variable_Also_Long) const;  
        // Even more longish declarations forces the parameters printed
        // one per line. This was the default formatting.

// \footnotetext{Another footnote in this class file, implemented with
// mark and text.}
};


Demo_Class& foo1( int& a, int* b);

Demo_Class* foo2( int& a, int* b);

Demo_Class &foo3( int &a, int *b);

Demo_Class *foo4( int &a, int *b);
class Demo_Class2 {
public:


// \section{Operator Test}
// 
// New creation variable is: `p'
// 
// Type casting through a conversion operator is the default behavior for
// the formatting routine if the return type before the operator keyword
// is empty.

    operator int () const;
        // Conversion operator.

    operator A<FT>() const;
        // Conversion operator.

// Sometimes, there is a choice between implementing an operator as a
// method or as a function. Both declarations will produce the same
// formatting, as demonstrated with the next two declarations.

    Demo_Class2  operator+(Demo_Class2 p, Demo_Class2 q);
        // Declaration via function.

    Demo_Class2  operator+(Demo_Class2 q);
        // Declaration via method.

// One can locally activate that the operator declaration is shown as it
// is written without operator formatting, const ...&, classname, or
// trailing const declarations for methods removal. This can be done with
// \ccTagFullDeclarations within a scope of braces {...}.

    Demo_Class2  operator+(const Demo_Class2& q) const;
        // Declaration via method.

// There is some laziness allowed in placing spaces around the operator
// characters. See the following examples:

    A  operator+(Demo_Class2 q);
        // C++ code: A operator+(Demo_Class2 q);

    A  operator +(Demo_Class2 q);
        // C++ code: A operator +(Demo_Class2 q);

    A  operator+ (Demo_Class2 q);
        // C++ code: A operator+ (Demo_Class2 q);

    A  operator + (Demo_Class2 q);
        // C++ code: A operator + (Demo_Class2 q);

// The keyword operator is reserved, but it can appear as a substring in
// another name. See the following examples that this style can handle
// such cases:

    A foo_operator(Demo_Class2 q);

    A noperator(Demo_Class2 q);

    A operatoro(Demo_Class2 q);

    A operator_(Demo_Class2 q);

    A operator0(Demo_Class2 q);

// A problem has occured in detecting the operator keyword if it was
// directly preceded by an & or * character. It is fixed as the following
// example demonstrates:

    Int &operator+=( Int a, Int b);
};

class ClassA {
public:


// \section{The List of All Operators}

    Ptr_Class  operator->(ClassA p);

    ClassA  operator[](ClassA p, int i);

    ClassA  operator()(ClassA p);

    ClassA  operator()(ClassA p, int i);

    ClassA  operator()(ClassA p, int i, int j);

    ClassA  operator()(ClassA p, int i, int j, int k);

    ClassA  operator()(ClassA p,  const A& a, B& b, C c, const D& d, ClassA  e);
        // all number and types of parameters are possible.

    ClassA  operator++(ClassA p);

    ClassA  operator++(ClassA p, int);
        // The postfix incr. operator has a hidden int parameter that the
        // formatting does not show.

    ClassA  operator--(ClassA p);

    ClassA  operator--(ClassA p, int);

    ClassA  operator~(ClassA p);

    ClassA  operator!(ClassA p);

    ClassA  operator-(ClassA p);

    ClassA  operator+(ClassA p);

    ClassA  operator&(ClassA p);

    ClassA  operator*(ClassA p);

    void*  operator new( size_t);
        // Hidden parameters are not shown. C++ code: \backslash method{
        // void* operator new( size_t);}.

    void  operator delete( void*, size_t);
        // Hidden parameters are not shown. C++ code: \backslash method{
        // void operator delete( void*, size_t);}

    void  operator delete[]( void*, size_t);
        // Hidden parameters are not shown again. C++ code: \backslash
        // method{ void operator delete[]( void*, size_t);}

    Member_Ptr  operator->*(ClassA p);

    ClassA  operator*(ClassA p, ClassA q);

    ClassA  operator/(ClassA p, ClassA q);

    ClassA  operator%(ClassA p, ClassA q);

    ClassA  operator+(ClassA p, ClassA q);

    ClassA  operator-(ClassA p, ClassA q);

    ClassA  operator<<(ClassA p, int i);

    ClassA  operator>>(ClassA p, int i);

    ClassA  operator<(ClassA p, ClassA q);

    ClassA  operator<=(ClassA p, ClassA q);

    ClassA  operator>(ClassA p, ClassA q);

    ClassA  operator>=(ClassA p, ClassA q);

    ClassA  operator==(ClassA p, ClassA q);

    ClassA  operator!=(ClassA p, ClassA q);

    ClassA  operator&(ClassA p, ClassA q);

    ClassA  operator^(ClassA p, ClassA q);

    ClassA  operator|(ClassA p, ClassA q);

    ClassA  operator&&(ClassA p, ClassA q);

    ClassA  operator||(ClassA p, ClassA q);

    ClassA  operator=(ClassA p, ClassA q);

    ClassA  operator*=(ClassA p, ClassA q);

    ClassA  operator/=(ClassA p, ClassA q);

    ClassA  operator%=(ClassA p, ClassA q);

    ClassA  operator+=(ClassA p, ClassA q);

    ClassA  operator-=(ClassA p, ClassA q);

    ClassA  operator<<=(ClassA p, ClassA q);

    ClassA  operator>>=(ClassA p, ClassA q);

    ClassA  operator&=(ClassA p, ClassA q);

    ClassA  operator|=(ClassA p, ClassA q);

    ClassA  operator^=(ClassA p, ClassA q);
};


Rep_Class::Nested_Class  foo(Rep_Class::Nested_Class p, Demo_Class q);
    // Declaration with scope.

Rep_Class :: Nested_Class  foo(Rep_Class :: Nested_Class p, Demo_Class q);
    // The same, surrounded by spaces.

Rep_Class::Nested_Class  TI::foo(Rep_Class::Nested_Class p, Demo_Class q);
    // Declaration with scope.

Rep_Class :: Nested_Class  CBP::foo(Rep_Class :: Nested_Class p, Demo_Class q);
    // The same, surrounded by spaces.
template < class FT<RT>  >
class ClassB {
public:


// SECTION: Demo Class Template
// ========================================================================
// 
// CREATION
// 
// New creation variable is: `p'

    ClassB();
        // default.

    ClassB( ClassB<FT<RT> > q);
        // copy.

    ClassB( A a, B *b);
        // arbitrary.

// OPERATIONS

    ClassB foo( ClassB q);
        // wrong, without template parameters.

    ClassB<FT<RT> > foo( ClassB<FT<RT> > q);
        // right, with template parameters.
};


enum Short { A, B, C};
    // Comment.

enum Funny_type_name1 { A_couple_of_entries1,  one_with_initialisation1 = 5, another1 = -3};
    // Comment.

enum Funny_type_name2 { A_couple_of_entries2,  one_with_initialisation2 = 5, another2 = -3};
    // Comment.

long int foo4;
    // Local variables are possible.

long int foo5 = 15;
    // Initialisation.

const long int foo6 = 15;
    // Make a constant.

typedef int integer;
    // Simple typedef.

typedef List< int> Integer_list;
    // Typedef including template parameters.

Intersection_type CGAL_intersection_type(  Polygon_2< R>, Polygon_2< R>);

enum global_enum { A_couple_of_entries3,  one_with_initialisation3 = 5, another3 = -3};

int CGAL_global_var;
template < class  R >
class CGAL_Point_2 {
public:


// \section{Default Parameters in Function Argument Lists}
// 
// The following example demonstrates the new ability to write default
// parameters with initializers in parantheses notion of C++.
// 
// New creation variable is: `p'

    CGAL_Point_2(const R::RT &x, const R::RT &y, const R::RT &w = R::RT(1.0));
        // blabla
};


int a_really_long_function_name( double paramter1, double  paramter2);
    // the default formatting. A bit more text is necessary to demonstrate
    // the right margin.

int a_really_long_function_name( double paramter1, double  paramter2);
    // the default formatting. A bit more text is necessary to demonstrate
    // the right margin.

int a_really_long_function_name( double paramter1, double  paramter2);
    // the alternative formatting. A bit more text is necessary to
    // demonstrate the right margin.

template<class A> int bar(A a);
    // A bit more text to demonstrate the right margin.

int a_really_long_function_name( double paramter1, double  paramter2);
    // the alternative formatting. A bit more text is necessary to
    // demonstrate the right margin.

template<class A> int bar(A a);
    // A bit more text to demonstrate the right margin.

int foo( int i, int j);
    // returns gnats(i,j).

int foo( int i, int j);
    // returns gnats(i,j).

int foo( int i, int j);
    // returns gnats(i,j).
class Gnu1 {
public:


// New creation variable is: `g'

    Gnu1( double d);
        // test.

    int foo( Gnats<T> gn);
        // blablabla.
};

class Gnu2 {
public:


// New creation variable is: `g'

    Gnu2( double d);
        // test.

    int foo( Gnats<T> gn);
        // blablabla.
};

template < class  C, class  C*  (C::*next)(), class  C* (C::*previous)() >
class CBP_Bidirectional_circulator {
public:


// New creation variable is: `circ'
// 
// NOTHING is interesting here.

    CBP_Bidirectional_circulator();
        // a const circulator `circ' with singular value.

    CBP_Bidirectional_circulator( const C* ptr);
        // a const circulator `circ' initialized to point to the element
        // `*ptr'.
};


int foo( double x);

int bar( double x);

int foo( double x);
    // Bla.

int bar( double x);
    // Blubb blubb.

int foo_baaaaaaaarrrr( double x);
    // Bla bal blabal blabal blabal blab.

int barfoooooooooooooooooo( double x);
    // Blubb blubblubb blubblubb blubblubb blubb.

int foo_baaaaaaaarrrr( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal blab.

int bar( double x);
    // Blubb blubb.

int foo_bar( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal blab.

int foo_bar( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal Blkj
    // fdkjbjlh flknlkj kjh kjh lkjh a bkuuhd kjwdhwkjdhh blab.

int bar( double x);
    // Blubb blubb.

int foo( double x);

int bar( double x);

int foo( double x);
    // Bla.

int bar( double x);
    // Blubb blubb.

int foo_baaaaaaaarrrr( double x);
    // Bla bal blabal blabal blabal blab.

int barfoooooooooooooooooo( double x);
    // Blubb blubblubb blubblubb blubblubb blubb.

int foo_baaaaaaaarrrr( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal blab.

int bar( double x);
    // Blubb blubb.

int foo_bar( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal blab.

int foo_bar( double x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh a bal blabal blabal blabal Blkj
    // fdkjbjlh flknlkj kjh kjh lkjh a bkuuhd kjwdhwkjdhh blab.

int bar( double x);
    // Blubb blubb.
class Base::Derived {
public:


// \section{Nested Classes}
// 
// First try, use the scope operator in the class name.
// 
// New creation variable is: `a'

    Base::Derived( double x);
        // a constructor.

    Base::Derived foo( Base::Derived b);
        // a function.
};


CGAL_Point<R> CGAL_f( CGAL_Vector<R> v);
    // the original CGAL prefix.

CGAL_Point<R> CGAL_f( CGAL_Vector<R> v);
    // assuming name spaces. ;-)

CGAL_Point<R> CGAL_f( CGAL_Vector<R> v);
    // make a Porsche from this declaration.
template < class T1 >
class CGAL_SegmentTree_2d {
public:


// The classname is CGAL_SegmentTree_2d<T1>::SegTree.
// 
// New creation variable is: `t'

    CGAL_SegmentTree_2d<T1>::SegTree (S2dList::const_iterator first);
        // constructor.

    CGAL_SegmentTree_2d<T1>::SegTree  foo (CGAL_SegmentTree_2d<T1>::SegTree same);
        // removal of own type in parameterlist of a member function.
};


CGAL_SegmentTree_2d<T1>::SegTree t(S2dList::const_iterator first);
    // fake constructor.
class ClassC {
public:


// New creation variable is: `pp'

    ClassC operator+(ClassC qq);
        // Declaration of an operator as a member function.
};


\tt k_th\ccFont _dim \tt k_th\ccFont _foo8;
    // a k_th-dimensional variable.
class \tt k_th\ccFont _dim_Class {
public:


// New creation variable is: `v'

    \tt k_th\ccFont _dim_Class();
        // The default constructor.

    \tt k_th\ccFont _dim_Class( \tt k_th\ccFont  _dim_param);
        // The custom constructor.

    \tt k_th\ccFont _dim_retvalue \tt k_th\ccFont  _dim_foo( \tt k_th\ccFont _dim_param);
        // a function.
};

template < class _T1, class  _C1, class  _T2, class   _C2, class \dots, class  _T\tt k\ccFont , class  _C\tt k\ccFont , class  _V >
class CGAL_Segment_tree_\tt k\ccFont  {
public:


// New creation variable is: `s'

    CGAL_Segment_tree_\tt k\ccFont <_T1, _C1, _T2,  _C2,\dots, _T\tt k\ccFont , _C\tt k\ccFont ,  _V>();
        // Creation of a k-dimensional Segment Tree, where the i-th
        // coordinate has type `_Ti' and a comparator of type `_Ci'.
        // Parameter type `Value' is set to `_V'.

    void CGAL_Segment_tree_\tt k \ccFont  window_query(GIS_interval_data_\tt k \ccFont   win, back_insert_iterator<GIS_interval_data_list_\tt k\ccFont> out);
        // All elements that intersect the associated d-dimensional half
        // open interval of `win' are placed in the list of `out'.
};

class ClassD {
public:


// \section{Customization Tags for the Style}
// 
// First, a declaration with all default substitution rule active.
// 
// New creation variable is: `pp'

    ClassD operator+(const ClassD& qq) const;
        // member function.

// Second, all rules switched off.

    ClassD operator+(const ClassD& qq) const;
        // member function.

// Third, back to the default.

    ClassD operator+(const ClassD& qq) const;
        // member function.
};

class Dummy {
public:


// \section{Test the removal or inlining of template declarations}
// 
// New creation variable is: `q'
// 
// Inlining of template declarations. The default is to format the
// template declaration in an extra line. Here the default with a struct
// and a global function.

    template <class T>   struct circulator_base1 {};
        // forward.

    template <class T>   Dummy();
        // constructor.

    template <class T>   T* value_type( const circulator_base<T>&);

// Here with the inlining on.

    template <class T>   struct circulator_base2 {};
        // forward.

    template <class T>   Dummy();
        // constructor.

    template <class T>   T* value_type( const circulator_base<T>&);

// Here with the removal on.

    template <class T>   struct circulator_base3 {};
        // forward.

    template <class T>   Dummy();
        // constructor.

    template <class T>   T* value_type( const circulator_base<T>&);
};

template < class T >
class Demo {
public:


// \section{A Function Pointer as an Operator Argument}
// 
// New creation variable is: `out'
// 
// \ccThreeToTwo

    Ascii_ostream& operator<<(ostream& (*f)(ostream&));
        // manipulators in streams.
};

template < class P, class  W, class  Key, class    Left, class  Right, class  Comp >
class CGAL_tree_point_interface1 {
public:


    Nested type required: CGAL_tree_point_interface1<P, W, Key,   Left, Right, Comp>::Tree_Point
        // complies to the container `Point'.
};

template < class Point, class  Window, class  Key, class    Point_func, class  Window_left_func, class    Window_right_func, class  Compare >
class CGAL_tree_point_interface2 {
public:


    Nested type required: CGAL_tree_point_interface2<Point, Window, Key,   Point_func, Window_left_func,   Window_right_func, Compare>::Tree_Point
        // complies to the container `Point'.
};

template < class Point, class  Window, class  Key, class    Point_func, class  Window_left_func, class    Window_right_func, class  Compare >
class CGAL_tree_point_interface3 {
public:


    Nested type required: CGAL_tree_point_interface3<Point, Window, Key,   Point_func, Window_left_func,   Window_right_func, Compare>::Tree_Point
        // complies to the container `Point'.
};

class Point_interface {
public:


// The following line must be done by hand if it should be formatted more
// nicely.
// 
// `typedef CGAL_tree_point_interface<Point, Window, Key, Point_func,
// Window_left_func, Window_right_func, Compare> Point_interface;'
// 
// Since here the specification is written in terms of the above `typedef'
// .

    Nested type required: Point_interface::Tree_Point
        // complies to the container `Point'.
};

template < class X >
class I {
public:


// New creation variable is: `i'

    I(I<X> j);
};


void math_foo();
    // Simple fraction {1 \over 2}.
class LinkA {
public:


// New creation variable is: `a'
// 
// A class LinkA with \ccHtmlNoClassIndex and \ccHtmlNoLinks to suppress
// cross linking:

     bool link1(const LinkA& a, const B& b);
        // a function is still indexed.

     bool link2(const LinkA& a, const B& b);
        // a method is never indexed.

    struct LocalX {};
        // a nested struct is not indexed.
};

class LinkB {
public:


// New creation variable is: `b'
// 
// A class LinkB without \ccHtmlNoClassIndex but the individual
// declarations are prefixed with \ccHtmlNoIndex and \ccHtmlNoLinks:

     bool link3(const LinkB& a, const B& b);
        // explicitly not indexed.

     bool link4(const LinkB& a, const B& b);
        // a method is never indexed.

    struct LocalY {};
        // a nested struct explicitly not indexed.
};


 bool link5(const LinkB& a, const B& b);
    // explicitly not indexed.

struct GlobalX {};
    // a global struct explicitly not indexed and not linked.

struct GlobalY {};
    // a global struct indexed and linked.
class CGAL_Geomview_stream {
public:


// New creation variable is: `G'

    template <class R>  CGAL_Geomview_stream&  operator<<(CGAL_Geomview_stream& G,  const CGAL_Point_2<R>& p);
        // Inserts the point `p' into the stream `G'. (global function)

    template <class R>  CGAL_Geomview_stream&  operator<<(const CGAL_Point_2<R>& p);
        // Inserts the point `p' into the stream `G'. (member function)
};

class FctClass {
public:


// New creation variable is: `f'
// 
// Use default settings.

    void method() const;
        // Normal Const Method.

    void operator+(FctClass b) const;
        // Normal Const Operator.

    void operator()();
        // Empty Fct. Operator.

    void operator()(int i);
        // Fct. Operator.

    void operator()() const;
        // Empty Fct. Operator.

    void operator()(int i, double d) const;
        // Fct. Operator.

    operator int();
        // Conversion Operator.

// Format full declarations.

    void method() const;
        // Normal Const Method.

    void operator+(FctClass b) const;
        // Normal Const Operator.

    void operator()();
        // Empty Fct. Operator.

    void operator()(int i);
        // Fct. Operator.

    void operator()() const;
        // Empty Fct. Operator.

    void operator()(int i, double d) const;
        // Fct. Operator.

    operator int();
        // Conversion Operator.
};

class Face {
public:


// SECTION: Removal of Own Class Name
// ========================================================================
// 
// New creation variable is: `f'

    Face(Face face);

    Face(Face_handle n);

    Face(Facehandle n);

    Face(handle_Face n);

    Face(Face5_handle n);

    Face(handle_5Face n);
};


CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak.

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal.

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal bur.

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal bur nic.

CGAL_Return<X,Y> foos_bar( int x);

CGAL_Retu_rn<X,Y> foos_bar( int x);

CGAL_Retu_1rn<X,Y> foos_bar( int x);

CGAL_Retu_12rn<X,Y> foos_bar( int x);

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int1 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int_1 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int_18 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

int foos_bar( double x);

CGAL_ReturnMMMMM<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_ReturnMMMMM<X,Y> foos_bar( int1 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_ReturnMMMMM<X,Y> foos_bar( int_1 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_ReturnMMMMM<X,Y> foos_bar( int_18 x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Retu_12rn<X,Y> foos_bar( double x);

CGAL_ReturnMMMMM<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( double<> x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( double<> x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak.

CGAL_Return<X,Y> foos_bar( double<> x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal.

CGAL_Return<X,Y> foos_bar( double<> x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal bur.

CGAL_Return<X,Y> foos_bar( double<> x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all bak bal bur nic.

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y,  double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y, double y,  double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Return<X,Y> foos_bar( int x, double y, double y, double y,  double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y, double y,   double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .

CGAL_Ret_12urn<X,Y> foos_bar( int x, double y, double y, double y,  double y, double y);
    // Blkj fdkjbjlh flknlkj kjh kjh lkjh all .
class IndentClass {
public:


// New creation variable is: `indent'

    IndentClass();
        // Short comment. Default constr.

// Now one with and one without comment.

    IndentClass();
        // Short comment. Default constr.

    IndentClass();

// Here with longer comment.

    IndentClass();
        // Short Short Short Short Short Short Short comment. Default
        // constr.

    IndentClass();
        // Short Short Short Short Short Short Short Short Short comment.
        // Default constr.

// Here with long parameter lists.

    IndentClass( double x, double y, double y, double y,   double y);
        // Short comment. Default constr.

    IndentClass( double x, double y, double y, double y,   double y);
        // Short comment. Default constr.

    IndentClass( double x, double y, double y, double y,   double y, double y);
        // Short comment. Default constr.

// Same, but with alternate parameter layout for long parameter names.
// Here with long parameter lists.

    IndentClass( double x, double y, double y, double y,   double y);
        // Short comment. Default constr.

    IndentClass( double x, double y, double y, double y,   double y);
        // Short comment. Default constr.

    IndentClass( double x, double y, double y, double y,   double y, double y);
        // Short comment. Default constr.

// Long template parameter list and short function.

    template < class Point, class Vector, class Traits,   class Magnificant> void floo();
        // thats it.
};


int foo( const Bla& ref) const;

int foo( const Bla& ref) const;
class NiceClass {
public:


    Nested type required: NiceClass::NType
        // comment.
};

template < class C >
class Circulator_from_container {
public:


// This is a ref class and should produce an index entry without the
// template parameter...
};

class Concept {
public:


// A reference concept is talked about here...
};

