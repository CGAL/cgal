namespace CGAL {

/*!
\ingroup nt_ralgebraic

An instance of this class represents an extension of the type `NT` by *one* square root of the type `ROOT`.

`NT` is required to be constructible from `ROOT`.

`NT` is required to be an `IntegralDomainWithoutDivision`.

`Sqrt_extension` is `RealEmbeddable` if NT is `RealEmbeddable`.

For example, let `Integer` be some type representing \f$ \Z\f$, then
`Sqrt_extension<Integer,Integer>` is able to represent \f$ \Z[\sqrt{\mathrm{root}}]\f$
for some arbitrary Integer \f$\mathrm{root}\f$. \cgalFootnote{\f$ R[a]\f$ denotes the extension of a ring \f$ R\f$ by an element \f$ a\f$. See also: <A HREF="http://mathworld.wolfram.com/ExtensionRing.html"><TT>http://mathworld.wolfram.com/ExtensionRing.html</TT></A>}
The value of \f$\mathrm{root}\f$ is set at
construction time, or set to zero if it is not specified.

Arithmetic operations among different extensions, say \f$ \Z[\sqrt{a}]\f$
and \f$ \Z[\sqrt{b}]\f$, are not supported.
The result would be in \f$ \Z[\sqrt{a},\sqrt{b}]\f$, which is not
representable by `Sqrt_extension<Integer,Integer>`.

\attention The user is responsible to check that arithmetic operations are carried out for elements from the same extensions only.

This is not tested by `Sqrt_extension` for efficiency reasons.
A violation of the precondition leads to undefined behavior.
Be aware that for efficiency reasons the given \f$\mathrm{root}\f$ is stored as it is given to
the constructor. In particular, an extension by a square root of a square is
considered as an extension.

Since elements of `Sqrt_extension` that lie in different extensions
are not interoperable with respect to any arithmetic operations, the full
value range of `Sqrt_extension` does not represent an algebraic structure.
However, each subset of the value range that represents the extension of
NT by a particular square root is a valid algebraic structure, since
this subset is closed under all provided arithmetic operations.
From there, `Sqrt_extension` can be used as if it were a model of an algebraic structure
concept, with the following correspondence:

<TABLE CELLSPACING=5 >
<TR>
<TD ALIGN=LEFT NOWRAP>
NT
<TD ALIGN=LEFT NOWRAP>
`Sqrt_extension`
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
<TR>
<TD ALIGN=LEFT NOWRAP>

<TR>
<TD ALIGN=LEFT NOWRAP>
`IntegralDomainWithoutDivision`
<TD ALIGN=LEFT NOWRAP>
`IntegralDomainWithoutDivision`
<TR>
<TD ALIGN=LEFT NOWRAP>
`IntegralDomain`
<TD ALIGN=LEFT NOWRAP>
`IntegralDomain`
<TR>
<TD ALIGN=LEFT NOWRAP>
`UniqueFactorizationDomain`
<TD ALIGN=LEFT NOWRAP>
`IntegralDomain`
<TR>
<TD ALIGN=LEFT NOWRAP>
`EuclideanRing`
<TD ALIGN=LEFT NOWRAP>
`IntegralDomain`
<TR>
<TD ALIGN=LEFT NOWRAP>
`Field`
<TD ALIGN=LEFT NOWRAP>
`Field`
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
</TABLE>

The extension of a `UniqueFactorizationDomain` or
`EuclideanRing` is just an `IntegralDomain`, since the extension in general destroys the unique factorization property. For instance consider \f$ \Z[\sqrt{10}]\f$, the extension of \f$ \Z\f$ by \f$ \sqrt{10}\f$: in \f$ \Z[\sqrt{10}]\f$ the element 10 has two different factorizations \f$ \sqrt{10} \cdot \sqrt{10}\f$ and \f$ 2 \cdot 5\f$. In particular, the factorization is not unique.

If `NT` is a model of `RealEmbeddable` the type `Sqrt_extension` is also considered as `RealEmbeddable`. However, by default it is not allowed to compare values from different extensions for efficiency reasons. In case such a comparison becomes necessary, use the member function compare with the according Boolean flag.
If such a comparison is a very frequent case, override the default of `DifferentExtensionComparable` by giving \cgalTagTrue as third template parameter. This effects the behavior of compare functions as well as the compare operators.

The fourth template argument, `FilterPredicates`, triggers an internal filter that may speed up comparisons and sign computations. In case `FilterPredicates` is set to \cgalTagTrue the type first computes a double interval containing the represented number and tries to perform the comparison or sign computation using this interval. Once computed, this interval is stored by the corresponding `Sqrt_extension` object for further usage. Note that this internal filter is switched off by default, since it may conflict with other filtering methods, such as `Lazy_exact_nt<Sqrt_extension>`.

In case `NT` is not `RealEmbeddable`, `DifferentExtensionComparable` as well as `FilterPredicates` have no effect.

\cgalModels `Assignable`
\cgalModels `CopyConstructible`
\cgalModels `DefaultConstructible`
\cgalModels `EqualityComparable`
\cgalModels `ImplicitInteroperable` with int
\cgalModels `ImplicitInteroperable` with NT
\cgalModels `Fraction` if NT is a `Fraction`
\cgalModels `RootOf_2`

\sa `IntegralDomainWithoutDivision`
\sa `IntegralDomain`
\sa `Field`
\sa `RealEmbeddable`
\sa \cgalTagTrue
\sa \cgalTagFalse

*/
template< typename NT, typename ROOT,
          typename DifferentExtensionComparable = Tag_false,
          typename FilterPredicates = Tag_false>
class Sqrt_extension {
public:

/// \name Creation
/// @{

/*!
Introduces a variable `ext` initialized with 0.
*/
Sqrt_extension ();

/*!
Copy constructor.
*/
Sqrt_extension (const Sqrt_extension& x);

/*!
Introduces a variable `ext` initialized with \f$ i\f$.
*/
Sqrt_extension (const int &i);

/*!
Introduces a variable `ext` initialized with \f$ x\f$.
*/
Sqrt_extension (const NT &x);

/*!
Introduces a variable `ext` initialized with \f$ x\f$.
\pre NT must constructible from NTX
*/
template <class NTX>
explicit Sqrt_extension(const NTX& x);

/*!
Constructor from int: `ext`\f$ = a0 +a1 \cdot sqrt(r)\f$. \pre \f$ r \neq0\f$
*/
Sqrt_extension (int a0, int a1, int r);

/*!
General constructor: `ext`\f$ = a0 + a1 \cdot sqrt(r)\f$. \pre \f$ r \neq0\f$
*/
Sqrt_extension (NT a0, NT a1, ROOT r);

/// @}



/// \name Operations
/// An object of type `Sqrt_extension` represent an expression of the
/// form: \f$ a0 + a1 \sqrt(\mathrm{root}) \f$.
/// @{

/*!
Const access operator for a0
*/
const NT & a0 () const ;

/*!
Const access operator for a1
*/
const NT & a1 () const ;

/*!
Const access operator for root
*/
const ROOT &         root () const;

/*!
Returns true in case root of `ext` is not zero.

Note that \f$ a1 == 0 \f$ does not imply \f$ \mathrm{root} == 0\f$.
*/
bool is_extended () const;

/*!
Simplifies the representation, in particular \f$\mathrm{root}\f$ is set to
zero if \f$ a1\f$ is zero, that is, `ext` becomes not extended.

Moreover, it propagates the simplify command to members
of `ext`. see also: `AlgebraicStructureTraits::Simplify`.

*/
void         simplify ();

/*!
returns true if `ext` represents the value zero.
*/
bool         is_zero () const;

/*!
Determines the sign of `ext` by (repeated) squaring.
\pre `Sqrt_extension` is `RealEmbeddable`.
*/
CGAL::Sign sign () const;

/*!
returns the absolute value of `ext`.
\pre `Sqrt_extension` is `RealEmbeddable`.
*/
Sqrt_extension abs () const;

/*!

Compares `ext` with y.

The optional bool `in_same_extension` indicates whether `ext`
and \f$ y\f$ are in the same extension of NT.

*/
CGAL::Comparison_result compare
(const Sqrt_extension& y, bool in_same_extension = !DifferentExtensionComparable::value) const;

/*!
\pre `(this->root()==0 or a.root()==0 or this->root() == a.root())`
*/
Sqrt_extension & operator+=(const Sqrt_extension& a);

/*!
\pre `(this->root()==0 or a.root()==0 or this->root() == a.root())`
*/
Sqrt_extension & operator-=(const Sqrt_extension& a);

/*!
\pre `(this->root()==0 or a.root()==0 or this->root() == a.root())`
*/
Sqrt_extension & operator*=(const Sqrt_extension& a);


/*!
\pre `(this->root()==0 or a.root()==0 or this->root() == a.root())`

In case `NT` is only an `IntegralDomain` operator/ implements integral
division. In case `NT` is a `Field` operator/ implements the field
division.
*/
Sqrt_extension & operator/=(const Sqrt_extension& a);

/// @}

}; /* end Sqrt_extension */

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
*/
Sqrt_extension operator+(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
*/
Sqrt_extension operator-(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
*/
Sqrt_extension operator*(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension

In case `NT` is only an `IntegralDomain` operator/ implements integral
division. In case `NT` is a `Field` operator/ implements the field
division.
*/
Sqrt_extension operator/(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
*/
bool operator==(const Sqrt_extension&a,const Sqrt_extension& b);



/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
*/
bool operator!=(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension

\attention Only exists when `Sqrt_extension` is `RealEmbeddable`.
*/
bool operator< (const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
\attention Only exists when `Sqrt_extension` is `RealEmbeddable`.
*/
bool operator<=(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
\attention Only exists when `Sqrt_extension` is `RealEmbeddable`.
*/
bool operator> (const Sqrt_extension&a,const Sqrt_extension& b);

/*!
\pre `(a.root()==0 or b.root()==0 or a.root() == b.root())`
\relates Sqrt_extension
\attention Only exists when `Sqrt_extension` is `RealEmbeddable`.
*/
bool operator>=(const Sqrt_extension&a,const Sqrt_extension& b);

/*!
writes `ext` to ostream `os`. The format depends on the `CGAL::IO::MODE` of `os`.

In case the mode is `CGAL::IO::ASCII` the format is `EXT[a0,a1,root]`.

In case the mode is `CGAL::IO::PRETTY` the format is human readable.

\attention `operator>>` must be defined for `ROOT` and `NT`.

\relates Sqrt_extension
*/
std::ostream& operator<<(std::ostream& os, const Sqrt_extension<NT,ROOT> &ext);

/*!
reads `ext` from istream `is` in format `EXT[a0,a1,root]`, the output format in mode `CGAL::IO::ASCII`

\attention `operator<<` must be defined exist for `ROOT` and `NT`.

\relates Sqrt_extension
*/
std::istream& operator>>(std::istream& is, const Sqrt_extension<NT,ROOT> &ext);

} /* end namespace CGAL */
