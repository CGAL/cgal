





/*!
  \ingroup PkgHandlesAndCirculatorsConcepts

A Circulator is similar to an Iterator, with the difference that it is designed for circular data structures.
<h1></h1>

\section circulatorConceptsIntro Circulators
Iterators in the \stl were tailored for linear sequences. The
specialization for circular data structures leads to slightly
different requirements which we will summarize in the
<I>circulators</I> concept. The main difference is that a circular
data structure has no natural past-the-end value. As a consequence, a
container supporting circulators will not have an
<TT>end()</TT>-member function. The semantic of a circulator range
differs from the semantic of an iterator range. For a circulator
`c` the range `[c, c)` denotes the sequence of all
elements in the data structure. For iterators, this range defines the
empty sequence. A separate test for an empty sequence has been added
to the circulator requirements: A comparison `c == nullptr`
for a circulator `c` is true for an empty sequence. As for \cpp,
we recommend the use of 0 instead of `nullptr`.

\subsection circulatorConceptsCirculators Circulators
Similar to \stl iterators, we distinguish between forward,
bidirectional, and random access circulators\cgalFootnote{Input
circulators are a contradiction, since any circulator is supposed to
return once to itself. Output circulators are not supported since they
would be indistinguishable from output iterators.}. Most requirements
for circulators are equal to those for iterators. We present the
changes, please refer to [\cgalCite{cgal:ms-strg-96}, chapter 18 or \cgalCite{cgal:ansi-is14882-98}]
for the iterator requirements.

<B>Past-the-end value:</B> There is no past-the-end value for circulators.

<B>Singular values:</B> There are no singular values for
circulators\cgalFootnote{Since circulators must be implemented as classes
anyway, there is no need to allow singular values for them. An
un-initalized circulator does not have a singular value, but is
supposed to refer to an empty sequence.}

<B>Empty sequence:</B> The comparison `c == nullptr` (or `c == 0`)
for a circulator `c` is `true` if `c`
denotes an empty sequence, and `false` otherwise.

<B>Dereferenceable values:</B> A circulator that does not denote an
empty sequence is dereferenceable.

<B>Reachability:</B> Each dereferenceable circulator can reach itself
with a finite and non-empty sequence of applications of `operator++`.

<B>Ranges:</B> For any circulator `c` the range `[c, c)` is a valid range.
If the circulator refers to an empty
sequence, the range `[c, c)` denotes the empty
range. Otherwise the circulator is dereferenceable and the range
`[c, c)` denotes the sequence of all elements in the data
structure. <I>Remark:</I> When a circulator is used in a place of an
iterator, as, for example, with an \stl algorithm, it will work as
expected with the only exception that, in \stl algorithms, the range
`[c, c)` denotes always the empty range.  This is not
a requirement, but a consequence of the requirements stated here and
the fact that the \stl requirements for iterator ranges are based on
the `operator++` and the `operator==`, which we use for
circulators as well. In principle, we face here the difference between
a `while` loop and a `do-while` loop.

<B>Types:</B> For a circulator of type `c` the following local
types are required:

<TABLE border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::value_type</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
value type the circulator refers to.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::reference</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
reference type used for the return type
of <TT>C::operator*()</TT>.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::pointer</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
pointer type used for the return type of
<TT>C::operator->()</TT>.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::size_type</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
unsigned integral type that can hold
the size of a sequence
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::difference_type</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
signed integral type that can hold
the distance between two circulators.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::iterator_category</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
circulator category.

</TABLE>

\subsection circulatorConceptsForward Forward Circulators

In the following, we assume that <TT>a</TT> and <TT>b</TT> are
circulators of type <TT>C</TT>, <TT>r</TT> is of type <TT>C&</TT> (is
assignable), and <TT>T</TT> denotes the value type of <TT>C</TT>. Let
<TT>D</TT> be the distance type of <TT>C</TT>. As for \cpp, we
recommend the use of 0 instead of <TT>nullptr</TT>.

<TABLE border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C()</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
a circulator equal to `nullptr` denoting an
empty sequence.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>a == nullptr</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
Returns <TT>true</TT> if <TT>a</TT> denotes an empty
sequence, <TT>false</TT> otherwise.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
For simplicity, <TT>nullptr == a</TT> is not required. The
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
behavior for comparisons with pointer-like
values different than <TT>nullptr</TT>
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
is undefined. A runtime assertion is
recommended.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>a != nullptr</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
Returns <TT>!(a == nullptr)</TT>.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>++r</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
Like for forward iterators, but a dereferenceable
circulator <TT>r</TT> will always
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
be dereferenceable after <TT>++r</TT> (no
past-the-end value). <I>Precondition:</I> <TT>r</TT>
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
does not denote an empty sequence.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>r++</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
Same as for <TT>++r</TT>.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::iterator_category</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
circulator category <TT>Forward_circulator_tag</TT>.

</TABLE>

\subsection circulatorConceptsBidirectional Bidirectional Circulators

The same requirements as for the forward circulators hold for
bidirectional iterators with the following change of the iterator
category:

<TABLE border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::iterator_category</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
circulator category <TT>Bidirectional_circulator_tag</TT>.

</TABLE>

\subsection circulatorConceptsRandomAccessCirculators Random Access Circulators

\anchor sectionMinCircleRequ

The same requirements as for the bidirectional circulators hold for
random access iterators with the following changes and extensions.

The idea of random access extends naturally to circulators using
equivalence classes modulus the length of the sequence. With this in
mind, the additional requirements for random access iterators hold
also for random access circulators. The only exception is that the
random access iterator is required to provide a total order on the
sequence, which a circulator cannot provide\cgalFootnote{One might define
an order by splitting the circle at a fixed point, e.g. the start
circulator provided from the data structure. This is what the adaptor
to iterators will do. Nonetheless, we do not require this for
circulators.}.

The difference of two circulators is not unique as for iterators. A
reasonable requirement demands that the result is in a certain range
`[1-size, size-1]`, where `size` is the size of the
sequence, and that whenever a circulator `a` is fixed that
the differences with all other circulators of the sequence form a
consistent ordering.

For the adaptor to iterators a minimal circulator
`d_min` is required for which the difference
`c - d_min` to all other circulators `c` is non negative.

<TABLE border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>b - a</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
limited range and consistent ordering
as explained above.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>a.min_circulator()</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
returns the minimal circulator from the
range `[a,a)`.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>C::iterator_category</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
circulator category <TT>Random_access_circulator_tag</TT>.

</TABLE>

\subsection circulatorConceptsConstCirculators Const Circulators

As with iterators, we distinguish between circulators and const
circulators. The expression <TT>*a = t</TT> with <TT>t</TT> of type
<TT>T</TT> is valid for mutable circulators. It is invalid for const
circulators.

\subsection subsec_CircContClass Circulators in Container Classes

For a container <TT>x</TT> of type <TT>X</TT> that supports
circulators <TT>c</TT> the following naming convention is recommended:

<TABLE  border="0"><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>X::Circulator</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
the type of the mutable circulator.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>X::Const_circulator</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
the type of the const circulator.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
<TT>c = x.begin()</TT>
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
the start circulator of the sequence.
It is of type <TT>X::Circulator</TT> for a
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
mutable container or <TT>X::Const_circulator</TT> for
a const container.
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

<TD ALIGN=LEFT VALIGN=TOP NOWRAP>


</TABLE>

If a container will support iterators and circulators, the member
function <TT>circulator_begin()</TT> is proposed. However, the support
of iterators and circulators simultaneously is not recommended, since
it would lead to fat interfaces. The natural choice should be
supported, the other concept will be available through adaptors.

\subsection subsec_Circ_ex Example

A generic <TT>contains</TT> function accepts a range of circulators
and a value. It returns `true` if the value is contained in the
sequence of items denoted by the range of circulators. As usual for
circular structures, a <TT>do</TT>-<TT>while</TT> loop is preferable,
such that for the specific input, <TT>c == d</TT>, all elements in the
sequence are reached. Note that the example simplifies if the sequence
is known to be non-empty, which is for example the common case in
polyhedral surfaces where vertices and facets have at least one
incident edge.

\code{.cpp}

template <class Circulator, class T>
bool contains( Circulator c, Circulator d, const T& value) {
  if (c != 0) {
    do {
      if (*c == value)
        return true;
    } while (++c != d);
  }
  return false;
}

\endcode
\cgalConcept
*/
class Circulator {
public:

}; /* end Circulator */

/*!
  \ingroup PkgHandlesAndCirculatorsConcepts
  See Subsection \ref circulatorConceptsForward in the page on Circulators.
\cgalConcept
*/
class ForwardCirculator {

}; /* end ForwardCirculator */

/*!
  \ingroup PkgHandlesAndCirculatorsConcepts
  See Subsection \ref circulatorConceptsBidirectional in the page on Circulators.
\cgalConcept
*/
class BidirectionalCirculator {

}; /* end ForwardCirculator */


/*!
  \ingroup PkgHandlesAndCirculatorsConcepts
  See Subsection \ref sectionMinCircleRequ in the page on Circulators.
\cgalConcept
*/
class RandomAccessCirculator {

}; /* end RandomAccessCirculator */
