@! ============================================================================
@! The CGAL Project
@! Implementation: 2D Smallest Enclosing Ellipse
@! ----------------------------------------------------------------------------
@! file  : web/Optimisation/Min_ellipse_2.aw
@! author: Bernd Gärtner, Sven Schönherr (sven@inf.fu-berlin.de)
@! ----------------------------------------------------------------------------
@! $Revision$
@! $Date$
@! ============================================================================
 
@p maximum_input_line_length = 180
@p maximum_output_line_length = 180

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amssymb}
@usepackage{cc_manual}
@article

\input{cprog.sty}
\setlength{\parskip}{1ex}

@! LaTeX macros
\newenvironment{pseudocode}[1]%
  {\vspace*{-0.5\baselineskip} \upshape \begin{tabbing}
     99 \= bla \= bla \= bla \= bla \= bla \= bla \= bla \= bla \kill
     #1 \+ \\}%
  {\end{tabbing}}
\newcommand{\keyword}[1]{\texttt{#1}}
\newcommand{\IF}{\keyword{IF} }
\newcommand{\THEN}{\keyword{THEN} \+ \\}
\newcommand{\ELSE}{\< \keyword{ELSE} \\}
\newcommand{\END}{\< \keyword{END} \- \\ }
\newcommand{\OR}{\keyword{OR} }
\newcommand{\FOR}{\keyword{FOR} }
\newcommand{\TO}{\keyword{TO} }
\newcommand{\DO}{\keyword{DO} \+ \\}
\newcommand{\RETURN}{\keyword{RETURN} }

\newcommand{\me}{\texttt{me}}
\newcommand{\SaveSpaceByHand}[1]{#1}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsrevision}{$Revision$}
\RCSdefDate{\rcsdate}{$Date$}

@t vskip 5 mm
@t title titlefont centre "CGAL -- 2D Smallest Enclosing Ellipse*"
@t vskip 1 mm
@t title smalltitlefont centre "Implementation Documentation"
@t vskip 8 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
\smallskip
\centerline{\rcsrevision\ , \rcsdate}
@t vskip 1 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Project
  No.~21957 (CGAL).}

@! ============================================================================
@! Introduction and Contents
@! ============================================================================

\section*{Introduction}

We provide an implementation of an optimisation algorithm for computing
the smallest (w.r.t.\ area) enclosing ellipse of a finite point set $P$
in the plane. The class template \ccc{CGAL_Min_ellipse_2} is implemented
as a semi-dynamic data structure, thus allowing to insert points while
maintaining the smallest enclosing circle. It is parameterized with a
traits class, that defines the abstract interface between the
optimisation algorithm and the primitives it uses. For ease of use, we
provide traits class adapters that interface the optimisation algorithm
with user supplied point classes.

This document is organized as follows. The algorithm is described in
Section~1. Section~2 contains the specifications as they appear in the
CGAL Reference Manual. Section~3 gives the implementations. In
Section~4 we provide a test program which performs some correctness
checks. Finally the product files are created in Section~5.

\tableofcontents

@! ============================================================================
@! The Algorithm
@! ============================================================================

\clearpage
\section{The Algorithm} \label{sec:algo}

The implementation is based on an algorithm by Welzl~\cite{w-sedbe-91a},
which we shortly describe now. The smallest (w.r.t.\ area) enclosing
ellipse of a finite point set $P$ in the plane, denoted by $me(P)$, is
built up incrementally, adding one point after another. Assume $me(P)$
has been constructed, and we would like to obtain $me(P \cup \{p\})$, $p$
some new point. There are two cases: if $p$ already lies inside $me(P)$,
then $me(P \cup \{p\}) = me(P)$. Otherwise $p$ must lie on the boundary
of $me(P \cup \{p\})$ (this is proved in~\cite{w-sedbe-91a} and not hard
to see), so we need to compute $me(P,\{p\})$, the smallest ellipse
enclosing $P$ with $p$ on the boundary. This is recursively done in the
same manner. In general, for point sets $P$,$B$, define $me(P,B)$ as the
smallest ellipse enclosing $P$ that has the points of $B$ on the boundary
(if defined). Although the algorithm finally delivers a ellipse
$me(P,\emptyset)$, it internally deals with ellipses that have a possibly
nonempty set $B$. Here is the pseudocode of Welzl's method.  To compute
$me(P)$, it is called with the pair $(P,\emptyset)$, assuming that
$P=\{p_1,\ldots,p_n\}$ is stored in a linked list.

\begin{pseudocode}{$\me(P,B)$:}
  $me := \me(\emptyset,B)$ \\
  \IF $|B| = 5$ \keyword{THEN} \RETURN $me$ \\
  \FOR $i := 1$ \TO $n$ \DO
      \IF $p_i \not\in me$ \THEN
          $me := \me(\{p_1,\ldots,p_{i-1}\}, B \cup \{p_i\})$ \\
          move $p_i$ to the front of $P$ \\
      \END
  \END
  \RETURN $me$ \\
\end{pseudocode}

Note the following: (a) $|B|$ is always bounded by 5, thus the computation of
$me(\emptyset,B)$ is easy. In our implementation, it is done by the private
member function \ccc{compute_ellipse}. (b) One can check that the method
maintains the invariant `$me(P,B)$ exists'. This justifies termination if
$|B| = 5$, because then $me(P,B)$ must be the unique ellipse with the points
of $B$ on the boundary, and $me(P,B)$ exists if and only if this ellipse
contains the points of $P$. Thus, no subsequent in-ellipse tests are
necessary anymore (for details see~\cite{w-sedbe-91a}).  (c) points which are
found to lie outside the current ellipse $me$ are considered `important' and
are moved to the front of the linked list that stores $P$. This is crucial
for the method's efficiency.

It can also be advisable to bring $P$ into random order before computation
starts. There are `bad' insertion orders which cause the method to be very
slow -- random shuffling gives these orders a very small probability.

@! ============================================================================
@! Specifications
@! ============================================================================

\clearpage
\section{Specifications}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

\renewcommand{\ccSection}{\ccSubsection}
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{}
\ccSetThreeColumns{typedef Traits::Ellipse}{}{%
  creates a variable \ccc{min_ellipse} of type
  \ccc{CGAL_Min_ellipse_2<Traits>}.}
\ccPropagateThreeToTwoColumns
\input{../../spec/Optimisation/Min_ellipse_2.tex}
\input{../../spec/Optimisation/Min_ellipse_2_adapterC2.tex}
\input{../../spec/Optimisation/Min_ellipse_2_adapterH2.tex}

@! ============================================================================
@! Implementations
@! ============================================================================

\clearpage
\section{Implementations}

@! ----------------------------------------------------------------------------
@! Class template CGAL_Min_ellipse_2<Traits>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont CGAL\_Min\_ellipse\_2<Traits>}

First, we declare the class template \ccc{CGAL_Min_ellipse_2}.

@macro<Min_ellipse_2 declaration> = @begin
    template < class _Traits >
    class CGAL_Min_ellipse_2;
@end

The actual work of the algorithm is done in the private member
functions \ccc{me} and \ccc{compute_ellipse}. The former directly
realizes the pseudocode of $\me(P,B)$, the latter solves the basic
case $\me(\emptyset,B)$, see Section~\ref{sec:algo}.

\emph{Workaround:} The GNU compiler (g++ 2.7.2[.?]) does not accept types
with scope operator as argument type or return type in class template
member functions. Therefore, all member functions are implemented in
the class interface.

The class interface looks as follows.

@macro <Min_ellipse_2 interface> = @begin
    template < class _Traits >
    class CGAL_Min_ellipse_2 {
      public:
	@<Min_ellipse_2 public interface>

      private:
	// private data members
	@<Min_ellipse_2 private data members>

	// copying and assignment not allowed!
	CGAL_Min_ellipse_2( CGAL_Min_ellipse_2<_Traits> const&);
	CGAL_Min_ellipse_2<_Traits>&
            operator = ( CGAL_Min_ellipse_2<_Traits> const&);

    @<dividing line>

    // Class implementation
    // ====================

      public:
	// Access functions and predicates
	// -------------------------------
	@<Min_ellipse_2 access functions `number_of_...'>

	@<Min_ellipse_2 predicates `is_...'>

	@<Min_ellipse_2 access functions>

	@<Min_ellipse_2 predicates>

      private:
	// Privat member functions
	// -----------------------
	@<Min_ellipse_2 private member function `compute_ellipse'>

	@<Min_ellipse_2 private member function `me'>

      public:
	// Constructors
	// ------------
	@<Min_ellipse_2 constructors>

	// Destructor
	// ----------
	@<Min_ellipse_2 destructor>

	// Modifiers
	// ---------
	@<Min_ellipse_2 modifiers>

	// Validity check
	// --------------
	// Not yet implemented!
	@!<Min_ellipse_2 validity check>

	// Miscellaneous
	// -------------
	@<Min_ellipse_2 miscellaneous>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Public Interface}

The functionality is described and documented in the specification
section, so we do not comment on it here.

@macro <Min_ellipse_2 public interface> = @begin
    // types
    typedef           _Traits                      Traits;
    typedef typename  _Traits::Point               Point;
    typedef typename  _Traits::Ellipse             Ellipse;
    typedef typename  list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                Support_point_iterator;

    /**************************************************************************
    WORKAROUND: The GNU compiler (g++ 2.7.2[.*]) does not accept types
    with scope operator as argument type or return type in class template
    member functions. Therefore, all member functions are implemented in
    the class interface.

    // creation
    CGAL_Min_ellipse_2( const Point*  first,
			const Point*  last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits());
    CGAL_Min_ellipse_2( list<Point>::const_iterator first,
			list<Point>::const_iterator last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits());
    CGAL_Min_ellipse_2( istream_iterator<Point,ptrdiff_t> first,
			istream_iterator<Point,ptrdiff_t> last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits())
    CGAL_Min_ellipse_2( Traits const& traits = Traits());
    CGAL_Min_ellipse_2( Point  const& p,
			Traits const& traits = Traits());
    CGAL_Min_ellipse_2( Point  const& p,
			Point  const& q,
			Traits const& traits = Traits());
    CGAL_Min_ellipse_2( Point  const& p1,
			Point  const& p2,
			Point  const& p3,
			Traits const& traits = Traits());
    ~CGAL_Min_ellipse_2( );

    // access functions
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;

    Point_iterator  points_begin( ) const;
    Point_iterator  points_end  ( ) const;

    Support_point_iterator  support_points_begin( ) const;
    Support_point_iterator  support_points_end  ( ) const;

    Point const&  support_point( int i) const;

    Ellipse const&  ellipse( ) const;

    // predicates
    CGAL_Bounded_side  bounded_side( Point const& p) const;
    bool  has_on_bounded_side      ( Point const& p) const;
    bool  has_on_boundary          ( Point const& p) const;
    bool  has_on_unbounded_side    ( Point const& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;

    // modifiers
    void  insert( Point const& p);
    void  insert( const Point* first,
                  const Point* last );
    void  insert( list<Point>::const_iterator first,
                  list<Point>::const_iterator last );
    void  insert( istream_iterator<Point,ptrdiff_t> first,
                  istream_iterator<Point,ptrdiff_t> last );
    void  clear( );

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;

    // miscellaneous
    Traits const&  traits( ) const;
    **************************************************************************/
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Data Members}

First, the traits class object is stored.

@macro <Min_ellipse_2 private data members> += @begin
    Traits       tco;				// traits class object
@end

The points of $P$ are internally stored as a linked list that allows to
bring points to the front of the list in constant time. We use the
sequence container \ccc{list} from STL~\cite{sl-stl-95}.

@macro <Min_ellipse_2 private data members> += @begin
    list<Point>  points;			// doubly linked list of points
@end

The support set $S$ of at most five support points is stored in an
array \ccc{support_points}, the actual number of support points is
given by \ccc{n_support_points}. During the computations, the set of
support points coincides with the set $B$ appearing in the pseudocode
for $\me(P,B)$, see Section~\ref{sec:algo}.

\emph{Workaround:} The array of support points is allocated dynamically,
because the SGI compiler (mipspro CC 7.1) does not accept a static
array here.

@macro <Min_ellipse_2 private data members> += @begin
    int          n_support_points;		// number of support points
    Point*       support_points;		// array of support points
@end

Finally, the actual ellipse is stored in a variable \ccc{ellipse}
provided by the traits class object, by the end of computation equal to
$me(P)$. During computation, \ccc{tco.ellipse} equals the ellipse $me$
appearing in the pseudocode for $\me(P,B)$, see Section~\ref{sec:algo}.

@! ----------------------------------------------------------------------------
\subsubsection{Constructors and Destructor}

We provide several different constructors, which can be put into two
groups. The constructors in the first group, i.e. the more important
ones, build the smallest enclosing ellipse $me(P)$ from a point set $P$,
given by a begin iterator and a past-the-end iterator. Usually this
would be done by a single template member function, but since most
compilers do not support this yet, we provide specialized constructors
for C~arrays (using pointers as iterators), for STL sequence containers
\ccc{vector<Point>} and \ccc{list<Point>} and for the STL input stream
iterator \ccc{istream_iterator<Point>}.  Actually, the constructors for
a C~array and a \ccc{vector<point>} are the same, since the random
access iterator of \ccc{vector<Point>} is implemented as \ccc{Point*}.

All three constructors of the first group copy the points into the
internal list \ccc{points}. If randomization is demanded, the points
are copied to a vector and shuffled at random, before being copied to
\ccc{points}. Finally the private member function $me$ is called to
compute $me(P)=me(P,\emptyset)$.

@macro <Min_ellipse_2 constructors> += @begin
    // STL-like constructor for C array and vector<Point>
    CGAL_Min_ellipse_2( const Point*  first,
			const Point*  last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// range not empty?
	if ( ( last-first) > 0) {

	    // store points
	    if ( randomize) {

	        // shuffle points at random
	        vector<Point> v( first, last);
		random_shuffle( v.begin(), v.end(), random);
		copy( v.begin(), v.end(), back_inserter( points)); }
	    else
		copy( first, last, back_inserter( points)); }

	// compute me
	me( points.end(), 0);
    }

    // STL-like constructor for list<Point>
    CGAL_Min_ellipse_2( list<Point>::const_iterator first,
			list<Point>::const_iterator last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// compute number of points
	list<Point>::size_type n = 0;
	CGAL__distance( first, last, n);
	if ( n > 0) {

	    // store points
	    if ( randomize) {

	        // shuffle points at random
		vector<Point> v;
		v.reserve( n);
		copy( first, last, back_inserter( v));
		random_shuffle( v.begin(), v.end(), random);
		copy( v.begin(), v.end(), back_inserter( points)); }
	    else
		copy( first, last, back_inserter( points)); }

	// compute me
	me( points.end(), 0);
    }

    // STL-like constructor for input stream iterator istream_iterator<Point>
    CGAL_Min_ellipse_2( istream_iterator<Point,ptrdiff_t>  first,
			istream_iterator<Point,ptrdiff_t>  last,
			bool          randomize = false,
			CGAL_Random&  random    = CGAL_random,
			Traits const& traits    = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// range not empty?
	if ( first != last) {

	    // store points
	    if ( randomize) {

	        // shuffle points at random
		vector<Point> v;
		copy( first, last, back_inserter( v));
		random_shuffle( v.begin(), v.end(), random);
		copy( v.begin(), v.end(), back_inserter( points)); }
	    else
		copy( first, last, back_inserter( points)); }

	// compute me
	me( points.end(), 0);
    }

@end

The remaining constructors are actually specializations of the previous
ones, building the smallest enclosing ellipse for up to five points.
The idea is the following: recall that for any point set $P$ there
exists $S \subseteq P$, $|S| \leq 5$ with $me(S) = me(P)$ (in fact,
such a set $S$ is determined by the algorithm). Once $S$ has been
computed (or given otherwise), $me(P)$ can easily be reconstructed from
$S$ in constant time. To make this reconstruction more convenient, a
constructor is available for each size of $|S|$, ranging from 0 to 5.
For $|S|=0$, we get the default constructor, building $me(\emptyset)$.

@macro <Min_ellipse_2 constructors> += @begin

    // default constructor
    inline
    CGAL_Min_ellipse_2( Traits const& traits = Traits())
	: tco( traits), n_support_points( 0)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// initialize ellipse
	tco.ellipse.set();

	CGAL_optimisation_postcondition( is_empty());
    }

    // constructor for one point
    inline
    CGAL_Min_ellipse_2( Point const& p, Traits const& traits = Traits())
	: tco( traits), points( 1, p), n_support_points( 1)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// initialize ellipse
	support_points[ 0] = p;
	tco.ellipse.set( p);

	CGAL_optimisation_postcondition( is_degenerate());
    }

    // constructor for two points
    inline
    CGAL_Min_ellipse_2( Point const& p,
			Point const& q,
			Traits const& traits = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// store points
	points.push_back( p);
	points.push_back( q);

	// compute me
	me( points.end(), 0);
    }

    // constructor for three points
    inline
    CGAL_Min_ellipse_2( Point const& p1,
			Point const& p2,
			Point const& p3,
			Traits const& traits = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// store points
	points.push_back( p1);
	points.push_back( p2);
	points.push_back( p3);

	// compute me
	me( points.end(), 0);
    }

    // constructor for four points
    inline
    CGAL_Min_ellipse_2( Point const& p1,
			Point const& p2,
			Point const& p3,
			Point const& p4,
			Traits const& traits = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// store points
	points.push_back( p1);
	points.push_back( p2);
	points.push_back( p3);
	points.push_back( p4);

	// compute me
	me( points.end(), 0);
    }

    // constructor for five points
    inline
    CGAL_Min_ellipse_2( Point const& p1,
			Point const& p2,
			Point const& p3,
			Point const& p4,
			Point const& p5,
			Traits const& traits = Traits())
	: tco( traits)
    {
	// allocate support points' array
	support_points = new Point[ 5];

	// store points
	points.push_back( p1);
	points.push_back( p2);
	points.push_back( p3);
	points.push_back( p4);
	points.push_back( p3);

	// compute me
	me( points.end(), 0);
    }
@end

The destructor only frees the memory of the support points' array.

@macro <Min_ellipse_2 destructor> = @begin
    inline
    ~CGAL_Min_ellipse_2( )
    {
	// free support points' array
	delete[] support_points;
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Access Functions}

These functions are used to retrieve information about the current
status of the \ccc{CGAL_Min_ellipse_2<Traits>} object. They are all very
simple (and therefore \ccc{inline}) and mostly rely on corresponding
access functions of the data members of \ccc{CGAL_Min_ellipse_2<Traits>}.

First, we define the \ccc{number_of_...} methods.

@macro <Min_ellipse_2 access functions `number_of_...'> = @begin
    // #points and #support points
    inline
    int
    number_of_points( ) const
    {
	return( points.size());
    }

    inline
    int
    number_of_support_points( ) const
    {
	return( n_support_points);
    }
@end

Then, we have the access functions for points and support points.

@macro <Min_ellipse_2 access functions> += @begin
    // access to points and support points
    inline
    Point_iterator
    points_begin( ) const
    {
	return( points.begin());
    }

    inline
    Point_iterator
    points_end( ) const
    {
	return( points.end());
    }

    inline
    Support_point_iterator
    support_points_begin( ) const
    {
	return( support_points);
    }

    inline
    Support_point_iterator
    support_points_end( ) const
    {
	return( support_points+n_support_points);
    }

    // random access for support points    
    inline
    Point const&
    support_point( int i) const
    {
	CGAL_optimisation_precondition( (i >= 0) &&
					(i <  number_of_support_points()));
	return( support_points[ i]);
    }
@end

Finally, the access function \ccc{ellipse}.

@macro <Min_ellipse_2 access functions> += @begin
    // ellipse
    inline
    Ellipse const&
    ellipse( ) const
    {
	return( tco.ellipse);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The predicates \ccc{is_empty} and \ccc{is_degenerate} are used in
preconditions and postconditions of some member functions. Therefore we
define them \ccc{inline} and put them in a separate macro.

@macro <Min_ellipse_2 predicates `is_...'> = @begin
    // is_... predicates
    inline
    bool
    is_empty( ) const
    {
	return( number_of_support_points() == 0);
    }

    inline
    bool
    is_degenerate( ) const
    {
	return( number_of_support_points() <  2);
    }
@end

The remaining predicates perform in-ellipse tests, based on the
corresponding predicates of class \ccc{Ellipse}.

@macro <Min_ellipse_2 predicates> = @begin
    // in-ellipse test predicates
    inline
    CGAL_Bounded_side
    bounded_side( Point const& p) const
    {
	return( tco.ellipse.bounded_side( p));
    }

    inline
    bool
    has_on_bounded_side( Point const& p) const
    {
	return( tco.ellipse.has_on_bounded_side( p));
    }

    inline
    bool
    has_on_boundary( Point const& p) const
    {
	return( tco.ellipse.has_on_boundary( p));
    }

    inline
    bool
    has_on_unbounded_side( Point const& p) const
    {
	return( tco.ellipse.has_on_unbounded_side( p));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers}

There is another way to build up $me(P)$, other than by supplying
the point set $P$ at once. Namely, $me(P)$ can be built up
incrementally, adding one point after another. If you look at the
pseudocode in the introduction, this comes quite naturally. The
modifying method \ccc{insert}, applied with point $p$ to a
\ccc{CGAL_Min_ellipse_2<Traits>} object representing $me(P)$,
computes $me(P \cup \{p\})$, where work has to be done only if $p$
lies outside $me(P)$. In this case, $me(P \cup \{p\}) = me(P,\{p\})$
holds, so the private member function \ccc{me} is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudocode in
Section~\ref{sec:algo}.

@macro <Min_ellipse_2 modifiers> += @begin
    void
    insert( Point const& p)
    {
	// p not in current ellipse?
	if ( has_on_unbounded_side( p)) {

	    // p new support point
	    support_points[ 0] = p;

	    // recompute me
	    me( points.end(), 1);

	    // store p as the first point in list
	    points.push_front( p); }
	else

	    // append p to the end of the list
	    points.push_back( p);
    }
@end

Inserting a range of points would usually be done by a single template
member function, but since most compilers do not support this yet, we
provide specialized \ccc{insert} functions for C~arrays (using pointers
as iterators), for STL sequence containers \ccc{vector<Point>} and
\ccc{list<Point>} and for the STL input stream iterator
\ccc{istream_iterator<Point>}. Actually, the \ccc{insert} function for
a C~array and a \ccc{vector<point>} are the same, since the random
access iterator of \ccc{vector<Point>} is implemented as \ccc{Point*}.

The following \ccc{insert} functions perform for each point \ccc{p} in
the range $[\mbox{\ccc{first}},\mbox{\ccc{last}})$ a call \ccc{insert(p)}.

@macro <Min_ellipse_2 modifiers> += @begin
    inline
    void
    insert( const Point* first, const Point* last)
    {
	for ( ; first != last; ++first)
	    insert( *first);
    }

    inline
    void
    insert( list<Point>::const_iterator first,
	    list<Point>::const_iterator last )
    {
	for ( ; first != last; ++first)
	    insert( *first);
    }

    inline
    void
    insert( istream_iterator<Point,ptrdiff_t>  first,
	    istream_iterator<Point,ptrdiff_t>  last )
    {
	for ( ; first != last; ++first)
	    insert( *first);
    }
@end

The member function \ccc{clear} deletes all points from a
\ccc{CGAL_Min_ellipse_2<Traits>} object and resets it to the
empty ellipse.

@macro <Min_ellipse_2 modifiers> += @begin
    void  clear( )
    {
	points.erase( points.begin(), points.end());
	n_support_points = 0;
	tco.ellipse.set();
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Miscellaneous}

The member function \ccc{traits} returns a const reference to the
traits class object.

@macro <Min_ellipse_2 miscellaneous> = @begin
    inline
    Traits const&
    traits( ) const
    {
	return( tco);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Min_ellipse_2 I/O operators declaration> = @begin
    template < class _Traits >
    ostream& operator << ( ostream& os, CGAL_Min_ellipse_2<_Traits> const& me);

    template < class _Traits >
    istream& operator >> ( istream& is, CGAL_Min_ellipse_2<_Traits>      & me);
@end

@macro <Min_ellipse_2 I/O operators> = @begin
    template < class _Traits >
    ostream&
    operator << ( ostream& os, CGAL_Min_ellipse_2<_Traits> const& min_ellipse)
    {
	typedef typename  CGAL_Min_ellipse_2<_Traits>::Point  Point;

	switch ( CGAL_get_mode( os)) {

	  case CGAL_IO::PRETTY:
	    os << endl;
	    os << "CGAL_Min_ellipse_2( |P| = "<< min_ellipse.number_of_points()
	       << ", |S| = " << min_ellipse.number_of_support_points() << endl;
	    os << "  P = {" << endl;
	    os << "    ";
	    copy( min_ellipse.points_begin(), min_ellipse.points_end(),
		  ostream_iterator<Point>( os, ",\n    "));
	    os << "}" << endl;
	    os << "  S = {" << endl;
	    os << "    ";
	    copy( min_ellipse.support_points_begin(),
		  min_ellipse.support_points_end(),
		  ostream_iterator<Point>( os, ",\n    "));
	    os << "}" << endl;
	    os << "  ellipse = " << min_ellipse.ellipse() << endl;
	    os << ")" << endl;
	    break;

	  case CGAL_IO::ASCII:
	    copy( min_ellipse.points_begin(), min_ellipse.points_end(),
		  ostream_iterator<Point>( os, "\n"));
	    break;

	  case CGAL_IO::BINARY:
	    copy( min_ellipse.points_begin(), min_ellipse.points_end(),
		  ostream_iterator<Point>( os));
	    break;

	  default:
	    CGAL_optimisation_assertion_msg( false,
					     "CGAL_get_mode( os) invalid!");
	    break; }

	return( os);
    }

    template < class Traits >
    istream&
    operator >> ( istream& is, CGAL_Min_ellipse_2<Traits>& min_ellipse)
    {
	switch ( CGAL_get_mode( is)) {

	  case CGAL_IO::PRETTY:
	    cerr << endl;
	    cerr << "Stream must be in ascii or binary mode" << endl;
	    break;

	  case CGAL_IO::ASCII:
	  case CGAL_IO::BINARY:
	    typedef typename  CGAL_Min_ellipse_2<Traits>::Point  Point;
	    typedef           istream_iterator<Point,ptrdiff_t>  Is_it;
	    min_ellipse.clear();
	    min_ellipse.insert( Is_it( is), Is_it());
	    break;

	  default:
	    CGAL_optimisation_assertion_msg( false, "CGAL_IO::mode invalid!");
	    break; }

	return( is);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Member Function {\ccFont compute\_ellipse}}

This is the method for computing the basic case $\me(\emptyset,B)$,
the set $B$ given by the first \ccc{n_support_points} in the array
\ccc{support_points}. It is realized by a simple case analysis,
noting that $|B| \leq 5$.

@macro <Min_ellipse_2 private member function `compute_ellipse'> = @begin
    // compute_ellipse
    inline
    void
    compute_ellipse( )
    {
	switch ( n_support_points) {
	  case 5:
	    tco.ellipse.set( support_points[ 0],
			     support_points[ 1],
			     support_points[ 2],
			     support_points[ 3],
			     support_points[ 4]);
            break;
	  case 4:
	    tco.ellipse.set( support_points[ 0],
			     support_points[ 1],
			     support_points[ 2],
			     support_points[ 3]);
            break;
	  case 3:
	    tco.ellipse.set( support_points[ 0],
			     support_points[ 1],
			     support_points[ 2]);
            break;
	  case 2:
	    tco.ellipse.set( support_points[ 0], support_points[ 1]);
	    break;
	  case 1:
	    tco.ellipse.set( support_points[ 0]);
	    break;
	  case 0:
	    tco.ellipse.set( );
	    break;
	  default:
	    CGAL_optimisation_assertion( ( n_support_points >= 0) &&
					 ( n_support_points <= 5) ); }
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Member Function {\ccFont me}}

This function computes the general ellipse $me(P,B)$, where $P$ contains
the points in the range $[$\ccc{points.begin()}$,$\ccc{last}$)$ and $B$
is given by the first \ccc{n_sp} support points in the array
\ccc{support_points}. The function is directly modelled after the
pseudocode above.

@macro <Min_ellipse_2 private member function `me'> = @begin
    void
    me( Point_iterator const& last, int n_sp)
    {
	// compute ellipse through support points
	n_support_points = n_sp;
	compute_ellipse();
	if ( n_sp == 5) return;

	// test first n points
	list<Point>::iterator  point_iter( points.begin());
	for ( ; last != point_iter; ) {
	    Point const& p( *point_iter);

	    // p not in current ellipse?
	    if ( has_on_unbounded_side( p)) {

		// recursive call with p as additional support point
		support_points[ n_sp] = p;
		me( point_iter, n_sp+1);

		// move current point to front
		if ( point_iter != points.begin()) {		// p not first?
		    points.push_front( p);
		    points.erase( point_iter++); }
		else
		    ++point_iter; }
	    else
		++point_iter; }
    }
@end

@! ============================================================================
@! Tests
@! ============================================================================

\clearpage
\section{Test}

We test \ccc{CGAL_Min_ellipse_2} with the traits class implementation
for optimisation algorithms, using exact arithmetic, i.e.\ Cartesian
representation with number type \ccc{CGAL_Quotient<CGAL_Gmpz>} and
homogeneous representation with number type \ccc{CGAL_Gmpz}.

@macro <Min_ellipse_2 test (includes and typedefs)> = @begin
    #include <CGAL/Gmpz.h>
    #include <CGAL/Cartesian.h>
    #include <CGAL/Homogeneous.h>
    #include <CGAL/Optimisation_traits_2.h>
    #include <CGAL/Min_ellipse_2.h>
    #include <CGAL/Random.h>
    #include <CGAL/IO/Verbose_ostream.h>
    #include <algo.h>
    #include <assert.h>
    #include <fstream.h>

    typedef  CGAL_Gmpz				 Rt;
    typedef  CGAL_Quotient< CGAL_Gmpz >		 Ft;

    typedef  CGAL_Cartesian< Ft >		 RepC;
    typedef  CGAL_Homogeneous< Rt >		 RepH;
    typedef  CGAL_Optimisation_traits_2< RepC >	 TraitsC;
    typedef  CGAL_Optimisation_traits_2< RepH >	 TraitsH;
@end

The command line option \ccc{-verbose} enables verbose output.

@macro <Min_ellipse_2 test (verbose option)> = @begin
    bool  verbose = false;
    if ( ( argc > 1) && ( strcmp( argv[ 1], "-verbose") == 0)) {
        verbose = true;
	--argc;
	++argv; }
@end

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

We call each function of class \ccc{CGAL_Min_ellipse_2<Traits>} at least
once to ensure code coverage.

@macro <Min_ellipse_2 test (code coverage)> = @begin
    cover_Min_ellipse_2( verbose, TraitsC(), Rt());
    cover_Min_ellipse_2( verbose, TraitsH(), Rt());
@end

@macro <Min_ellipse_2 test (code coverage test function)> = @begin
    template < class Traits, class RT >
    void
    cover_Min_ellipse_2( bool verbose, Traits const&, RT const&)
    {
	typedef  CGAL_Min_ellipse_2< Traits >	Min_ellipse;
	typedef  Min_ellipse::Point		Point;
	typedef  Min_ellipse::Ellipse		Ellipse;

	CGAL_Verbose_ostream verr( verbose);

	// generate `n' points at random
	const int    n = 20;
	CGAL_Random  random_x, random_y;
	Point	     random_points[ n];
	int	     i;
	verr << n << " random points from [0,128)^2:" << endl;
	for ( i = 0; i < n; ++i)
	    random_points[ i] = Point( RT( random_x( 128)),
				       RT( random_y( 128)));
	if ( verbose)
	    for ( i = 0; i < n; ++i)
		cerr << i << ": " << random_points[ i] << endl;

	// cover code
	verr << endl << "default constructor...";
	{
	    Min_ellipse  me;
	    bool  is_valid = me.is_valid( verbose);
	    bool  is_empty = me.is_empty();
	    assert( is_valid); 
	    assert( is_empty);
	}

	verr << endl << "one point constructor...";
	{
	    Min_ellipse  me( random_points[ 0]);
	    bool  is_valid      = me.is_valid( verbose);
	    bool  is_degenerate = me.is_degenerate();
	    assert( is_valid);
	    assert( is_degenerate);
	}

	verr << endl << "two points constructor...";
	{
	    Min_ellipse  me( random_points[ 1],
			     random_points[ 2]);
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == 2);
	}

	verr << endl << "three points constructor...";
	{    
	    Min_ellipse  me( random_points[ 3],
			     random_points[ 4],
			     random_points[ 5]);
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == 3);
	}

	verr << endl << "four points constructor...";
	{    
	    Min_ellipse  me( random_points[ 6],
			     random_points[ 7],
			     random_points[ 8],
			     random_points[ 9]);
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == 4);
	}

	verr << endl << "five points constructor...";
	{    
	    Min_ellipse  me( random_points[ 10],
			     random_points[ 11],
			     random_points[ 12],
			     random_points[ 13],
			     random_points[ 14]);
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == 5);
	}

	verr << endl << "Point* constructor...";
	Min_ellipse  me( random_points, random_points+9);
	{
	    Min_ellipse  me2( random_points, random_points+9, true);
	    bool  is_valid  = me .is_valid( verbose);
	    bool  is_valid2 = me2.is_valid( verbose);
	    assert( is_valid);
	    assert( is_valid2);
	    assert( me .number_of_points() == 9);
	    assert( me2.number_of_points() == 9);
	    assert( me.ellipse() == me2.ellipse());
	}

	verr << endl << "list<Point>::const_iterator constructor...";
	{
	    Min_ellipse  me1( me.points_begin(), me.points_end());
	    Min_ellipse  me2( me.points_begin(), me.points_end(), true);
	    bool  is_valid1 = me1.is_valid( verbose);
	    bool  is_valid2 = me2.is_valid( verbose);
	    assert( is_valid1);
	    assert( is_valid2);
	    assert( me1.number_of_points() == 9);
	    assert( me2.number_of_points() == 9);
	    assert( me.ellipse() == me1.ellipse());
	    assert( me.ellipse() == me2.ellipse());
	}

	verr << endl << "#points already called above.";

	verr << endl << "points access already called above.";

	verr << endl << "support points access...";
	{
	    Point  support_point;
	    Min_ellipse::Support_point_iterator
		iter( me.support_points_begin());
            for ( i = 0; i < me.number_of_support_points(); ++i, ++iter) {
		support_point = me.support_point( i);
		assert( support_point == *iter); }
	    Min_ellipse::Support_point_iterator
		end_iter( me.support_points_end());
	    assert( iter == end_iter);
	}

	verr << endl << "ellipse access already called above...";

	verr << endl << "in-ellipse predicates...";
	{
	    Point              p;
	    CGAL_Bounded_side  bounded_side;
	    bool               has_on_bounded_side;
	    bool	       has_on_boundary;
	    bool	       has_on_unbounded_side;
	    for ( i = 0; i < 9; ++i) {
		p = random_points[ i];
		bounded_side          = me.bounded_side( p);
		has_on_bounded_side   = me.has_on_bounded_side( p);
		has_on_boundary       = me.has_on_boundary( p);
		has_on_unbounded_side = me.has_on_unbounded_side( p);
	    assert( bounded_side != CGAL_ON_UNBOUNDED_SIDE);
	    assert( has_on_bounded_side || has_on_boundary);
	    assert( ! has_on_unbounded_side); }
	}

	verr << endl << "is_... predicates already called above.";

	verr << endl << "single point insert...";
	me.insert( random_points[ 9]);
	{
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == 10);
	}

	verr << endl << "Point* insert...";
	me.insert( random_points+10, random_points+n);
	{
	    bool  is_valid = me.is_valid( verbose);
	    assert( is_valid);
	    assert( me.number_of_points() == n);
	}

	verr << endl << "list<Point>::const_iterator insert...";
	{
	    Min_ellipse  me2;
	    me2.insert( me.points_begin(), me.points_end());
	    bool  is_valid = me2.is_valid( verbose);
	    assert( is_valid);
	    assert( me2.number_of_points() == n);
	    
	    verr << endl << "clear...";
	    me2.clear();	
	          is_valid = me2.is_valid( verbose);
	    bool  is_empty = me2.is_empty();
	    assert( is_valid); 
	    assert( is_empty);
    	}

	verr << endl << "validity check already called several times.";

	verr << endl << "traits class access...";
	{
	    Traits  traits( me.traits());
	}

	verr << endl << "I/O...";
	{
	    verr << endl << "  writing `test_Min_ellipse_2.ascii'...";
	    ofstream os( "test_Min_ellipse_2.ascii");
	    CGAL_set_ascii_mode( os);
	    os << me;
	}
	{
	    verr << endl << "  writing `test_Min_ellipse_2.pretty'...";
	    ofstream os( "test_Min_ellipse_2.pretty");
	    CGAL_set_pretty_mode( os);
	    os << me;
	}
	{
	    verr << endl << "  writing `test_Min_ellipse_2.binary'...";
	    ofstream os( "test_Min_ellipse_2.binary");
	    CGAL_set_binary_mode( os);
	    os << me;
	}
	{
	    verr << endl << "  reading `test_Min_ellipse_2.ascii'...";
	    Min_ellipse me_in;
	    ifstream is( "test_Min_ellipse_2.ascii");
	    CGAL_set_ascii_mode( is);
	    is >> me_in;
	    bool    is_valid = me_in.is_valid( verbose);
	    assert( is_valid);
	    assert( me_in.number_of_points() == n);
	    assert( me_in.ellipse() == me.ellipse());
	}
	verr << endl;
    }
@end

@! ----------------------------------------------------------------------------
@! External Test Sets
@! ----------------------------------------------------------------------------

\subsection{External Test Sets}

In addition, some data files can be given as command line arguments.
A data file contains pairs of \ccc{int}s, namely the x- and
y-coordinates of a set of points. The first number in the file is the
number of points. A short description of the test set is given at the
end of each file.

@macro <Min_ellipse_2 test (external test sets)> = @begin
    while ( argc > 1) {

	typedef  CGAL_Min_ellipse_2< TraitsH >	Min_ellipse;
	typedef  Min_ellipse::Point		Point;
	typedef  Min_ellipse::Ellipse		Ellipse;

	CGAL_Verbose_ostream verr( verbose);

	// read points from file
	verr << endl << "input file: `" << argv[ 1] << "'" << flush;

	list<Point>  points;
	int          n, x, y;
	ifstream     in( argv[ 1]);
	in >> n;
	assert( in);
	for ( int i = 0; i < n; ++i) {
	    in >> x >> y;
	    assert( in);
	    points.push_back( Point( x, y)); }

	// compute and check min_ellipse
	Min_ellipse  me2( points.begin(), points.end());
	bool  is_valid = me2.is_valid( verbose);
	assert( is_valid);

	// next file
	--argc;
	++argv; }
@end

@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@! ----------------------------------------------------------------------------
@! Min_ellipse_2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2.h}

@file <../../include/CGAL/Optimisation/Min_ellipse_2.h> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Optimisation/Min_ellipse_2.h")

    #ifndef CGAL_MIN_ELLIPSE_2_H
    #define CGAL_MIN_ELLIPSE_2_H

    // Class declaration
    // =================
    @<Min_ellipse_2 declaration>

    // Class interface
    // ===============
    // includes
    #ifndef CGAL_RANDOM_H
    #  include <CGAL/Random.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/optimisation_assertions.h>
    #endif
    #ifndef CGAL_OPTIMISATION_MISC_H
    #  include <CGAL/optimisation_misc.h>
    #endif
    #include <list.h>
    #include <vector.h>
    #include <algo.h>
    #include <iostream.h>

    @<Min_ellipse_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators declaration>

    #ifdef CGAL_INCLUDE_TEMPLATE_CODE
    #  include <CGAL/Min_ellipse_2.C>
    #endif

    #endif // CGAL_MIN_ELLIPSE_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2.C}

@file <../../include/CGAL/Optimisation/Min_ellipse_2.C> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Optimisation/Min_ellipse_2.C")

    // Class implementation (continued)
    // ================================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators>

    @<end of file line>
@end

@i ../file_header.awi
 
@macro <Min_ellipse_2 header>(1) many = @begin
    @<file header>("2D Smallest Enclosing Ellipse",@1,
		   "Optimisation/Min_ellipse_2",
		   "Bernd Gärtner, Sven Schönherr",
		   "$Revision$","$Date$")
@end

@! ----------------------------------------------------------------------------
@! test_Min_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{test\_Min\_ellipse\_2.C}

@file <../../test/Optimisation/test_Min_ellipse_2.C> = @begin
    @<Min_ellipse_2 header>("test/Optimisation/test_Min_ellipse_2.C")

    @<Min_ellipse_2 test (includes and typedefs)>

    // code coverage test function
    // ---------------------------
    @<Min_ellipse_2 test (code coverage test function)>

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        // command line options
	// --------------------
	// option `-verbose'
	@<Min_ellipse_2 test (verbose option)>

	// code coverage
	// -------------
	@<Min_ellipse_2 test (code coverage)>

	// external test sets
	// -------------------
	@<Min_ellipse_2 test (external test sets)>
    }

    @<end of file line>
@end

@! ============================================================================
@! Bibliography
@! ============================================================================

\clearpage
\bibliographystyle{plain}
\bibliography{geom,cgal}

@! ===== EOF ==================================================================

\clearpage

@section{OLD}
@subsection{old}

Most important, a \ccc{CGAL_Min_ellipse_2} object at any time
keeps the actual ellipse $\me(\emptyset,B)$. Unlike in the class
\ccc{CGAL_Min_ellipse_2}, this is not a \ccc{CGAL_Ellipse_2}, on
the one hand, because this type does not exist yet, on the other
hand because the representation of $\me(\emptyset,B)$ is different
for any number of support points. For $|B|=3,5$, we store an explicit
ellipse (although in different formats, see below), for $|B|=4$, the
representation is implicit, because an explicit one might contain
irrational coordinates. For $|B|\leq 2$, we have no representation
at all. In this case, the ellipse is degenerate, and the in-ellipse test
(which is the only routine we need the ellipse for) can be decided directly
from the set $B$. 

Although at any point in time, only one representation is
valid and will be accessed, all three representations (for $|B|=3,4,5$)
are simultaneously stored in a \ccc{CGAL_Min_ellipse_2} object. Since
none of them is really space-consuming, this is not a problem. 
(Unions, which would be just the right concept in this situation, are not
allowed to be used with classes as members, even if the respective
constructors do absolutely nothing). Each representation has its own
class. 

@macro <private classes> zero += @begin
    @<Ellipse_3 class>
    @<Ellipse_4 class>
    @<Ellipse_5 class>
@end

Here are the actual representations. 

@macro <private data members> zero += @begin
    Ellipse_3 ellipse_3;
    Ellipse_4 ellipse_4;
    Ellipse_5 ellipse_5;
@end   
    

@! ----------------------------------------------------------------------------
@subsubsection{ Ellipse representation classes}
@! ----------------------------------------------------------------------------
\label{ellipse_rep}
We have three classes,
\ccc{Ellipse_3}, \ccc{Ellipse_4} and \ccc{Ellipse_5}, being in charge
of sets $B$ with 3,4, or 5 points. All classes have a method \ccc{set} to
compute some representation of $\me(\emptyset,B)$, suitable to do in-ellipse
tests `$p\in \me(\emptyset,B)$'. To perform these tests, each class has a
\ccc{bounded_side} method.

@macro <Ellipse_3 class> zero = @begin
    class Ellipse_3 {
	public:
	    @<Ellipse_3 data members>

	    // default constructor
	    // -------------------
	    Ellipse_3 () {}

	    // set method
	    // ----------
	    void set (const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2,
	              const CGAL_Point_2<R>& p3)
	    {
	    	@<Ellipse_3::set body>
	    }

	    // bounded_side method
	    // -------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p) const
	    {
	    	@<Ellipse_3::bounded_side body>
	    }
    };
@end 

@macro <Ellipse_4 class> zero = @begin
    class Ellipse_4 {
	public:
	    @<Ellipse_4 data members>

	    // default constructor
	    // -------------------
	    Ellipse_4 () {}

	    // set method
	    // ----------
	    void set (const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2,
	              const CGAL_Point_2<R>& p3, const CGAL_Point_2<R>& p4)
	    {
	    	@<Ellipse_4::set body>
	    }

	    // bounded_side method, non-const
	    // ------------------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p)
	    {
	    	@<Ellipse_4::bounded_side body>
	    }
    };
@end 

The class \ccc{Ellipse_5} does not have a method to compute the ellipse 
from its five support points. The reason is that by the time an 
\ccc{Ellipse_5} object is set up, the ellipse through the five points 
has already been computed and stored in the \ccc{Ellipse_4} representation, 
see implementation section below. Thus it suffices to `steal' this ellipse.
To this end, a reference to the \ccc{Ellipse_4} object is passed to the
\ccc{set} method. Because the ellipse in question is computed by the
\ccc{bounded_side} method of the \ccc{Ellipse_4} class, the latter method
is not a \ccc{const} method, unlike all other \ccc{bounded_side} methods.

@macro <Ellipse_5 class> zero = @begin
    class Ellipse_5 {
	public:
	    @<Ellipse_5 data members>

 	    // default constructor
	    // -------------------
	    Ellipse_5 () {}

	    // set method
	    // ----------
 	    void set (Ellipse_4& ellipse_4)
	    {
	    	@<Ellipse_5::set body>	
	    }    
	
	    // bounded_side method
	    // -------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p) const
	    {
	 	@<Ellipse_5::bounded_side body>
	    }	   
    };
@end 	    

@! ----------------------------------------------------------------------------
@section{Class \texttt{CGAL\_Min\_ellipse\_2}: Implementation}
@! ----------------------------------------------------------------------------

The implementation consists of several parts, each of which is
described in the sequel. The actual work is hidden in the member
functions \ccc{set} and \ccc{bounded_side} of the ellipse
representation classes. The remaining functions -- in particular
the private member functions \ccc{me} and \ccc{compute_ellipse} -- are
implemented completely similar to the corresponding functions of the
class \ccc{CGAL_Min_ellipse_2}. An exception are the predicates for
the in-ellipse tests. In \ccc{CGAL_Min_ellipse_2}, the in-ellipse predicates
were directly mapped to the corresponding predicates over the 
\ccc{CGAL_ellipse_2} object stored in the class. In our case,
the \ccc{bounded_side} method is mapped to the corresponding one
of the respective ellipse representation class (and evaluated directly
if the current ellipse $me$ is defined by less than three points); the predicates (\ccc{has_on_bounded_side}, \ldots) are then implemented by 
directly referring to the \ccc{bounded_side} method.  

@! ----------------------------------------------------------------------------
@subsection{ Update operations}
@! ----------------------------------------------------------------------------

There is another way to build up $\me(P)$, other than by supplying the
point set $P$ at once. Namely, $\me(P)$ can be built up incrementally,
adding one point after another. If you look at the pseudocode in the
introduction, this comes quite naturally. The method \ccc{insert},
applied with point $p$ to a \ccc{CGAL_Min_ellipse_2} object
representing $\me(P)$, computes $\me(P \cup \{p\})$, where work has to
be done only if $p$ lies outside $\me(P)$.  In this case, $\me(P \cup
\{p\}) = \me(P,\{p\})$ holds, so the method $\me$ is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudocode and
the `main' constructor above.

@macro <update operations> zero += @begin
    // update operations
    // -----------------
    template < class R >
    void
    CGAL_Min_ellipse_2<R>::
    insert( const CGAL_Point_2<R>& p)
    {
	// store point
	int old_n = number_of_points();
	points.push_back( _Point( p));  // NOTE: p is not linked with list yet!

	// p not in current ellipse?
	if ( has_on_unbounded_side( p)) {

	    // p new support point
	    support_points[ 0] = p;

	    // recompute me
	    me( old_n, 1); }

	// make p the first point in list
	if ( old_n > 0)                                  // old list not empty?
	    points[ i_first_point].pred = old_n;
	points[ old_n].succ = i_first_point;
	i_first_point = old_n;
    }
@end

The operation \ccc{reserve} does nothing but tell the
\ccc{CGAL_Min_ellipse_2} that some number $n$ of points might
eventually be inserted, allowing the object to allocate storage for
them at once. Inserting the points without doing this might lead to
overhead caused by moving the array \ccc{points} around in memory
while it grows.

@macro <update operations> += @begin

    template < class R > inline
    void
    CGAL_Min_ellipse_2<R>::
    reserve( int n)
    {
	points.reserve( n);
    }
@end

@! ----------------------------------------------------------------------------
@subsection{ Access Operations and Predicates}
@! ----------------------------------------------------------------------------

Next we have predicates testing the relative position of a point
w.r.t. the ellipse. They rely on the \ccc{bounded_side} method
which for less than 3 support points is evaluated directly, while
for at least three support points, the corresponding methods
of the ellipse representation classes are used. A degenerate ellipse 
defined by less than two support points has no bounded side. In
this case, the ellipse itself is the boundary, everything else
belongs to the unbounded side. Note that the \ccc{bounded_side} method
is a \ccc{const} method although the underlying method of the class
\ccc{Ellipse_4} is not. The difference is that, `from the outside', the
\ccc{CGAL_Min_ellipse_2} object appears unchanged after calling
\ccc{bounded_side} (logical constness), while internally, 
we explicitly make use of the changes to an \ccc{Ellipse_4} object
caused by a call to \ccc{bounded_side}, see the discussion in 
subsection \ref{ellipse_rep}.

@macro <access operations and predicates> zero += @begin
    // in-ellipse predicates
    // -------------------
    template < class R > inline
    CGAL_Bounded_side
    CGAL_Min_ellipse_2<R>::
    bounded_side( const CGAL_Point_2<R>& p) const
    {
	switch (n_support_points) {
	  case 5: 
	    return ellipse_5.bounded_side (p);
	    break;
	  case 4:
	    // cast away ellipse_4's constness 
	    return ( const_cast( Ellipse_4&, ellipse_4)).bounded_side (p);
	    break;
	  case 3:
	    return ellipse_3.bounded_side (p);
	    break;
	  case 2: // Min_ellipse is line segment spanned by support points
	    if (CGAL_between (support_points [0], p, support_points [1]))
		return CGAL_ON_BOUNDARY;
	    else
		return CGAL_ON_UNBOUNDED_SIDE;
	    break;
	  case 1: // Min_ellipse is equal to its support point
	    if (p==support_points [0])
		return CGAL_ON_BOUNDARY;
	    else
		return CGAL_ON_UNBOUNDED_SIDE;
	    break;
	  case 0: // Min_ellipse is empty
	    return CGAL_ON_UNBOUNDED_SIDE;
	  default:
 	    // we should never get here ...
	    CGAL_Min_ellipse_2_assertion 
		( (number_of_support_points() >= 0) &&
		  (number_of_support_points() <= 5) );
	    // ... and we definitely never get here
	    return CGAL_ON_BOUNDARY; 	
	    break; 
	}
    }
@end
 
@! ----------------------------------------------------------------------------
@subsection{ Private Member Function \texttt{compute\_ellipse}}
@! ----------------------------------------------------------------------------

This is the method for computing $me(\emptyset,B)$ the set $B$ given
by the first \ccc{n_support_points} indices in the array
\ccc{i_support_points}. It is realized by a case analysis,
noting that $|B| \leq 5$. If $|B|\leq 2$, nothing is done,
in the other cases the \ccc{set} methods for \ccc{Ellipse_3}, 
\ccc{Ellipse_4}, or \ccc{Ellipse_5} are called. 

@macro <private member functions> zero += @begin
    template < class R >
    void
    CGAL_Min_ellipse_2<R>::
    compute_ellipse( )
    {
	switch ( n_support_points) {
	  case 5:
	    // `steal' ellipse from Ellipse_4 represemtation
	    ellipse_5.set (ellipse_4);	    
            break;
	  case 4:
	    ellipse_4.set (support_points[0], support_points[1],
			   support_points[2], support_points[3]);
	    break;
	  case 3:
	    ellipse_3.set (support_points[0], support_points[1],
			   support_points[2]); 
	    break;
	  default:
	    // do not compute any representation of a degenerate ellipse
	    break; 
	}
    }
@end

@! ----------------------------------------------------------------------------
@subsection{ Class \texttt{Conic<R>}}
@! ----------------------------------------------------------------------------

In their implementations, the ellipse representation classes rely on a
concept more general than ellipses, namely on {\em conics}. Ellipses are
special conics, in addition there are {\em hyperbolas} and {\em parabolas}. 
Conics play a particularly important role in the \ccc{Ellipse_4} class.

A conic (in linear form) is the set of points $q=(x,y)^T$ satisfying
\begin{equation}
{\cal C}(q) := rx^2 + sy^2 + 2txy + 2ux + 2vy + w = 0. 
\end{equation}
$r,s,t,u,v,w$ being real parameters. 

Alternatively, ${\cal C}$ can be written as 
$${\cal C}(q) = (x,y)\left(\begin{array}{cc}r & t \\ t & s \end{array}
\right)(x,y)^T + 2(u,v)(x,y)^T + w,$$
and the determinant of the matrix 
$$M=\left(\begin{array}{cc}r & t \\ t & s \end{array}
\right)$$ tells the type of ${\cal C}$. We also denote $\det(M)$ 
by $\det({\cal C})$. 

If $\det({\cal C})>0$, ${\cal C}$ is an ellipse, for $\det({\cal C})<0$ 
we get a hyperbola, and if $\det({\cal C})=0$, ${\cal C}$ is a parabola. 

The {\em linear combination} $\lambda{\cal C}_1 +\mu{\cal C}_2$ of two
conics ${\cal C}_1$ and ${\cal C}_2$ is defined as the conic ${\cal C}$
with $${\cal C}(q) := \lambda{\cal C}_1(q) +\mu{\cal C}_2(q)=0.$$

We introduce a conic class that supports 
\begin{itemize}
\item definition of a conic directly from values $r,s,t,u,v,w$,
\item computation of the linear combination of two conics, 
\item conic evaluation (i.e. computation of the value 
${\cal C}(p)$ for given $p\in\R^2$), 
\item computation of $\det({\cal C})$,
\item normalization, i.e. scaling of the conic by $-1$ -- if 
necessary  -- to guarantee $r\geq 0$. 
If $\det(M) > 0$, normalization guarantees that $M$ 
is positive definite afterwards (which is important
for doing in-ellipse tests, see subsequent subsections).
\end{itemize}

@macro <conic class> zero = @begin
    template <class R>
    class Conic {
	public:
	    R::FT  r;
	    R::FT  s;
	    R::FT  t;
	    R::FT  u;
	    R::FT  v;
	    R::FT  w;

	    // default constructor
	    // -------------------
	    Conic () {} 
	    
	    // member functions
	    // ----------------
	    // definition of conic by values r, s, t, u, v, w
	    void set (R::FT _r, R::FT _s, R::FT _t, R::FT _u,
		      R::FT _v, R::FT _w) 
	    {
		r=_r; s=_s; t=_t; u=_u; v=_v; w=_w;
	    }

	    // definition of conic by linear combination of conics
	    void set (R::FT lambda, Conic& conic1, 
		      R::FT mu, Conic& conic2)
	    {
		r = lambda*conic1.r + mu*conic2.r;
		s = lambda*conic1.s + mu*conic2.s;
		t = lambda*conic1.t + mu*conic2.t;
		u = lambda*conic1.u + mu*conic2.u;
		v = lambda*conic1.v + mu*conic2.v;
		w = lambda*conic1.w + mu*conic2.w;
	    }

	    // evaluation of conic at point p
	    R::FT eval (const CGAL_Point_2<R>& p) const
	    {
		R::FT x = p.x(), y = p.y();
		return r*x*x + s*y*y + R::FT(2)*(t*x*y + u*x + v*y) + w;
	    }

	    // determinant of matrix M = [[r,t],[t,s]]
	    R::FT det () const
	    {
		return r*s-t*t;
	    }

	    // scaling of conic such that r becomes nonnegative
	    void normalize ()
	    {
	   	if (r<R::FT(0)) {
		   r=-r; s=-s; t=-t; u=-u; v=-v; w=-w;
	  	}
	    } 
    };
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_3 Class}
@! ----------------------------------------------------------------------------

This class stores $\me(\emptyset,B)$ directly as an ellipse 
in {\em center form}, i.e. in the form $$(q-c)^TM(q-c)-z = 0,$$
$c\in \R^2$ the center. 
For support points $p_1,p_2,p_3$, the ellipse is given by
$$c=\frac{1}{3}\sum_{i=1}^3p_i,\quad M^{-1} = 
\frac{2}{3}\sum_{i=1}^3(p_i-c)(p_i-c)^T, \quad z = 1,$$
see \cite{gs-seefe-97}.

Let $$M^{-1} = \frac{2}{3}
\left(\begin{array}{cc} m_{11} & m_{12} \\ m_{12} & m_{22}
\end{array}\right).$$ Then 
$$M = \frac{3}{2}\left(\begin{array}{cc} m_{11} & m_{12} \\ m_{12} & m_{22}
\end{array}\right)^{-1} = \frac{3}{2(m_{11}m_{22}-m_{12}^2)}
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right).$$

This implies 
\begin{eqnarray*}
(q-c)^TM(q-c)-z 
&=& \frac{3}{2(m_{11}m_{22}-m_{12}^2)}(q-c)^T
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right)(q-c)-1,
\end{eqnarray*}
which is equal in sign to 
\begin{equation}
\label{ellipse_3_test}
(q-c)^T
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right)(q-c)-\frac{2}{3}(m_{11}m_{22}-m_{12}^2).
\end{equation}
(For this, we need
the fact that $M$ is positive definite, hence $\det(M),\det(M')$ and
$\det(M'^{-1})$ are positive.)
If $z'$ denotes $2(m_{11}m_{22}-m_{12}^2)/3$, an
\ccc{Ellipse_3} object stores the values 
$c, z', m_{11},2m_{12},m_{22}$.
This enables the subsequent
\ccc{bounded_side} routine to perform the in-ellipse test by 
just evaluating the sign of (\ref{ellipse_3_test}).

@macro <Ellipse_3 data members> zero = @begin
    R::FT  cx;
    R::FT  cy;
    R::FT  z;
    R::FT  m11;
    R::FT  m12;
    R::FT  m22;
@end
  
@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{set} method}
@! ----------------------------------------------------------------------------

Using the formulas above, the private data members are
straightforward to set.

@macro <Ellipse_3::set body> zero = @begin
    R::FT  x1 = p1.x();
    R::FT  y1 = p1.y();
    R::FT  x2 = p2.x(); 
    R::FT  y2 = p2.y();
    R::FT  x3 = p3.x(); 
    R::FT  y3 = p3.y();

    // center c
    cx = (x1+x2+x3)/R::FT(3);
    cy = (y1+y2+y3)/R::FT(3);

    R::FT  x1_cx = x1-cx;
    R::FT  x2_cx = x2-cx;
    R::FT  x3_cx = x3-cx;
    R::FT  y1_cy = y1-cy;
    R::FT  y2_cy = y2-cy;
    R::FT  y3_cy = y3-cy;

    // matrix M
    m11 = (x1_cx)*(x1_cx) + (x2_cx)*(x2_cx) + (x3_cx)*(x3_cx);
    m12 = (x1_cx)*(y1_cy) + (x2_cx)*(y2_cy) + (x3_cx)*(y3_cy);
    m22 = (y1_cy)*(y1_cy) + (y2_cy)*(y2_cy) + (y3_cy)*(y3_cy);

    // z'
    z = R::FT(2)*(m11*m22-m12*m12)/R::FT(3);

    // assert positive definiteness
    CGAL_Min_ellipse_2_assertion( (z > R::FT(0)) && (m11 > R::FT(0)) );
@end

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{bounded\_side} method}
@! ----------------------------------------------------------------------------

This is quite obvious now by (\ref{ellipse_3_test}). 
Since $M$ is positive definite, we have
$p=(x,y)^T$ in/on/outside the ellipse iff
$$(x-c_x,y-c_y)^T \left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right) (x-c_x,y-c_y) - z' \left\{\begin{array}{c}<\\=\\>
\end{array}\right. 0.$$

@macro <Ellipse_3::bounded_side body> zero = @begin
    R::FT  x = p.x();
    R::FT  y = p.y();
    R::FT  x_cx = x-cx;
    R::FT  y_cy = y-cy;
    R::FT  discr = m22*(x_cx)*(x_cx) + m11*(y_cy)*(y_cy)
		   - R::FT(2)*m12*(x_cx)*(y_cy) - z;	    
    return static_cast( CGAL_Bounded_side, CGAL_sign( discr));
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_4 Class}
@! ----------------------------------------------------------------------------

This is by far the most complicated one among the three representation 
classes. The \ccc{set} method computes some implicit representation of 
the ellipse, derived from the four support points. The \ccc{bounded_side}
predicate works over this implicit representation and is not
straightforward anymore.  

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{set} method}
@! ----------------------------------------------------------------------------

Subsequent computations assume that $p_1,\ldots p_4$ are in convex position 
and enumerated in clockwise or counterclockwise order (we would like 
to work without the latter assumption but don't know how). While convex
position is guaranteed by the ambient algorithm, a (counter)clockwise
odering can be achieved in a preprocessing step, as follows.

We first check whether $p_1$ and $p_3$ lie on different sides of the
line through $p_2,p_4$. If so, the segment $\overline{p_2p_4}$ is a
diagonal of the convex quadrilateral defined by the four points, so
the given order is ok. Otherwise, we check whether $\overline{p_3p_2}$
defines a diagonal. From the previous test we already know the relative
position of $p_4$ w.r.t. this segment, so it remains to check the position
of $p_1$. If $p_1,p_4$ lie on different sides, $\overline{p_3p_2}$ is a
diagonal. To achieve correct enumeration of the points, we just need to
swap points $p_1$ and $p_2$. If $p_1,p_4$ lie on the same side of
$\overline{p_3p_2}$, then $\overline{p_1p_2}$ must be a diagonal. For
the correct order, we swap $p_2$ and $p_3$.

@macro <Ellipse_4::set body> zero = @begin
    // copy p_i's to allow swapping
    CGAL_Point_2<R>  q1 = p1;
    CGAL_Point_2<R>  q2 = p2;
    CGAL_Point_2<R>  q3 = p3;
    CGAL_Point_2<R>  q4 = p4;

    // preprocessing: bring q_i's into (counter)clockwise order
    CGAL_Orientation side1_24 = CGAL_orientation (q2, q4, q1),
	     	     side3_24 = CGAL_orientation (q2, q4, q3);
    if (side1_24 == side3_24) { 
	CGAL_Orientation side1_32 = CGAL_orientation (q3, q2, q1);
	if (side1_32 != side3_24)  // side3_24 == side4_32
	    swap (q1, q2); 	// swap from STL algobase
	else 
	    swap (q2, q3);
    }

    // assert correct ordering
    CGAL_Min_ellipse_2_assertion(
        ( CGAL_orientation( q1, q2, q3) == CGAL_orientation( q2, q3, q4)) &&
	( CGAL_orientation( q2, q3, q4) == CGAL_orientation( q3, q4, q1)) &&
	( CGAL_orientation( q3, q4, q1) == CGAL_orientation( q4, q1, q2)) &&
 	( CGAL_orientation( q1, q2, q3) != CGAL_COLLINEAR               ) );

    // do the actual computations
    @<Ellipse_4 setup>
@end 

In all terminology that follows, we refer to \cite{gs-epsee-97,gs-seefe-97}. 
The implicit representation of
$\me(\emptyset,\{p_1,p_2,p_3,p_4\})$ first of all consists of two
special conics ${\cal C}_1,{\cal C}_2$ through the for points. 
Assuming that $q_i=(x_i,y_i)$, these conics are given by parameters 
$r_1,s_1,t_1,u_1,v_1,w_1$ and $r_2,s_2,t_2,u_2,v_2,w_2$, where
\begin{eqnarray*}
  r_1 & = & (y_1-y_2)(y_3-y_4), \\
  s_1 & = & (x_1-x_2)(x_3-x_4), \\
  2t_1 & = & -(x_1-x_2)(y_3-y_4)-(y_1-y_2)(x_3-x_4), \\
  2u_1 & = & (x_1y_2-x_2y_1)(y_3-y_4) + (x_3y_4-x_4y_3)(y_1-y_2), \\
  2v_1 & = & (x_1y_2-x_2y_1)(x_4-x_3) + (x_3y_4-x_4y_3)(x_2-x_1), \\
  w_1 & = & (x_1y_2-x_2y_1)(x_3y_4-x_4y_3)
\end{eqnarray*}
and
\begin{eqnarray*}
  r_2 & = & (y_2-y_3)(y_4-y_1), \\
  s_2 & = & (x_2-x_3)(x_4-x_1), \\
  2t_2 & = & -(x_2-x_3)(y_4-y_1)-(y_2-y_3)(x_4-x_1), \\
  2u_2 & = & (x_2y_3 - x_3y_2)(y_4-y_1) + (x_4y_1-x_1y_4)(y_2-y_3), \\
  2v_2 & = & (x_2y_3 - x_3y_2)(x_1-x_4) + (x_4y_1-x_1y_4)(x_3-x_2), \\
  w_2 & = & (x_2y_3 - x_3y_2)(x_4y_1-x_1y_4).
\end{eqnarray*}

${\cal C}_1$ and ${\cal C}_2$ are stored as data members in an
\ccc{Ellipse_4} object.

@macro <Ellipse_4 data members> zero += @begin
    Conic<R>  conic1;
    Conic<R>  conic2;
@end

@macro <Ellipse_4 setup> zero += @begin
    // point coordinates
    R::FT  x1 = q1.x();
    R::FT  y1 = q1.y();
    R::FT  x2 = q2.x();
    R::FT  y2 = q2.y();
    R::FT  x3 = q3.x();
    R::FT  y3 = q3.y();
    R::FT  x4 = q4.x();
    R::FT  y4 = q4.y();

    // auxiliary values
    R::FT  y1_y2 = y1-y2;
    R::FT  y3_y4 = y3-y4;
    R::FT  x1_x2 = x1-x2;
    R::FT  x3_x4 = x3-x4;
    R::FT  y2_y3 = y2-y3;
    R::FT  y4_y1 = y4-y1;
    R::FT  x2_x3 = x2-x3;
    R::FT  x4_x1 = x4-x1;
    R::FT  x1y2_x2y1 = x1*y2 - x2*y1;
    R::FT  x3y4_x4y3 = x3*y4 - x4*y3;
    R::FT  x2y3_x3y2 = x2*y3 - x3*y2;
    R::FT  x4y1_x1y4 = x4*y1 - x1*y4;

    // define conics
    conic1.set (  y1_y2 * y3_y4, 
		  x1_x2 * x3_x4, 
	        -(x1_x2 * y3_y4 + y1_y2 * x3_x4)/R::FT(2),
		 (x1y2_x2y1 * y3_y4 + x3y4_x4y3 * y1_y2)/R::FT(2),
		-(x1y2_x2y1 * x3_x4 + x3y4_x4y3 * x1_x2)/R::FT(2),
    		  x1y2_x2y1 * x3y4_x4y3 ); 
    conic2.set (  y2_y3 * y4_y1,
		  x2_x3 * x4_x1,
		-(x2_x3 * y4_y1 + y2_y3 * x4_x1)/R::FT(2),
		 (x2y3_x3y2 * y4_y1 + x4y1_x1y4 * y2_y3)/R::FT(2),
		-(x2y3_x3y2 * x4_x1 + x4y1_x1y4 * x2_x3)/R::FT(2),
		  x2y3_x3y2 * x4y1_x1y4 );
@end

In addition, we store values that are needed to perform the in-ellipse
test w.r.t. $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, given a fifth point $p$.
Some of these values are independent of $p$, they are computed here. 

In \cite{gs-epsee-97,gs-seefe-97} it is shown that the unique conic
${\cal C}_0$ through $\{p_1,p_2,p_3,p_4,p\}$ is given by
$$\lambda_0 {\cal C}_1 + \mu_0 {\cal C}_2, $$ where
$\lambda_0,\mu_0$ are given by
$\lambda_0 = {\cal C}_2(p), \mu_0 = -{\cal C}_1(p).$

Furthermore, to do the in-ellipse test with $p$ over $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, the sign of a certain 
derivative needs to be computed. For this, we define
\begin{equation}
\left(\begin{array}{c}\lambda(\tau)\\ \mu(\tau)\end{array}\right) = 
\left(\begin{array}{c}\lambda_0\\ \mu_0\end{array}\right) + {\tau}
\left(\begin{array}{r} -r2 \\ r1 \end{array}\right)
\end{equation}
and then consider the conic
$$
{\cal C}(\tau) := \lambda(\tau) {\cal C}_1 + \mu(\tau) {\cal C}_2,
$$
given by parameters $r(\tau),s(\tau), t(\tau), u(\tau),v(\tau),w(\tau)$.
To evaluate the derivative in question, we need the values
$$r(0)=r_0, s(0)=s_0, t(0)=t_0, u(0)=u_0, v(0)=v_0, w(0)=w_0$$
and 
$$r'(0), s'(0), t'(0), u'(0), v'(0), w'(0),$$
primed functions being derivatives w.r.t. $\tau$
(see subsection \ref{ellipse_4::bounded_side}).
While the former values depend on the fifth point $p$, the primed values 
do not, so they can be computed and stored already here. We get
\begin{eqnarray*}
r'(0) &=& r_1 r_2 - r_2 r_1 = 0, \\
s'(0) &=& r_1 s_2 - r_2 s_1, \\
t'(0) &=& r_1 t_2 - r_2 t_1, \\
u'(0) &=& r_1 u_2 - r_2 u_1, \\
v'(0) &=& r_1 v_2 - r_2 v_1, \\
w'(0) &=& r_1 w_2 - r_2 w_1.
\end{eqnarray*}

@macro <Ellipse_4 data members> += @begin
    R::FT  ds;	  
    R::FT  dt;
    R::FT  du;
    R::FT  dv;	
    R::FT  dw;
@end

@macro <Ellipse_4 setup> += @begin
    ds = conic1.r*conic2.s - conic2.r*conic1.s;
    dt = conic1.r*conic2.t - conic2.r*conic1.t;
    du = conic1.r*conic2.u - conic2.r*conic1.u;
    dv = conic1.r*conic2.v - conic2.r*conic1.v;
    dw = conic1.r*conic2.w - conic2.r*conic1.w;
@end

Finally, we need some ellipse through the four points. Such an ellipse
can be obtained as $$E = \lambda {\cal C}_1 + \mu {\cal C}_2,$$
with \begin{eqnarray*}
\lambda &=& 2\gamma-\beta, \\
\mu &=& 2\alpha-\beta,
\end{eqnarray*}
$\alpha = \det({\cal C}_1), \gamma = \det({\cal C}_2), 
\beta = r_1s_2+r_2s_1-2t_1t_2$.
This ellipse is stored as another conic in \ccc{Ellipse_4}.

@macro <Ellipse_4 data members> += @begin
    Conic<R> ellipse;
@end

@macro <Ellipse_4 setup> += @begin
    R::FT  beta   = conic1.r*conic2.s + conic2.r*conic1.s 
		    - R::FT(2)*conic1.t*conic2.t;
    R::FT  lambda = R::FT(2)*conic2.det()-beta;
    R::FT  mu     = R::FT(2)*conic1.det()-beta;
    ellipse.set (lambda, conic1, mu, conic2);
    ellipse.normalize();
    CGAL_Min_ellipse_2_assertion( ellipse.det() > R::FT(0));
@end

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{bounded\_side} method}
@! ----------------------------------------------------------------------------
\label{ellipse_4::bounded_side}
Given a fifth point $p$, the first action is to compute the unique conic
${\cal C}_0$ through $p_1,p_2,p_3,p_4,p$, as already described in the 
previous subsubsection. Since ${\cal C}_0$ is possibly being
recycled for usage in the \ccc{Ellipse_5} class later, we need to 
store it. Thus, ${\cal C}_0$ is the only data member of \ccc{Ellipse_4}
that is not independent of the query point $p$. (For this reason, the
\ccc{bounded_side} method is not a \ccc{const} method; moreover, the change 
to the \ccc{Ellipse_4} object caused by calling the method does not fall 
under logical constness, because the change becomes quite visible when 
${\cal C}_0$ is later passed to the \ccc{set} method of \ccc{Ellipse_5}.)

@macro <Ellipse_4 data members> += @begin
    Conic<R> conic0;
@end

Depending on the type of ${\cal C}_0$ (given by its determinant), 
different actions are taken.

@macro <Ellipse_4::bounded_side body> zero = @begin
    R::FT  lambda_0 =  conic2.eval(p);
    R::FT  mu_0     = -conic1.eval(p);
    conic0.set (lambda_0, conic1, mu_0, conic2);
    R::FT d = conic0.det();
    if (d > R::FT(0)) {
	@<handle ellipse case>
    }
    else {
	@<handle hyperbola/parabola case>
    }
@end

In the hyperbola/parabola case, we test $p$ against the ellipse
$E$ precomputed by the \ccc{Ellipse_4} \ccc{set} method, and just return
the result (which is the same for all ellipses, in particular for
$\me(\emptyset,B)$).

@macro <handle hyperbola/parabola case> zero = @begin
    R::FT discr = ellipse.eval(p); 
    return static_cast (CGAL_Bounded_side, CGAL_sign (discr));
@end

To handle the ellipse case, we need to evaluate two derivatives and 
combine the results. Let us review these two tasks.

(i) compute $$\rho := \frac{\partial}{\partial\tau}E^{\tau}(p),$$
where 
\begin{eqnarray*}
E^{\tau}(p) &=&  
(\lambda_0-\tau r_2){\cal C}_1(p) + (\mu_0+\tau r_1){\cal C}_2(p) \\
&=& \tau (r_1{\cal C}_2(p) - r_2{\cal C}_1(p)).
\end{eqnarray*}
This means, $$\rho = \lambda_0r_1 + \mu_0r_2 = r_0.$$ 

@macro <handle ellipse case> zero += @begin
    R::FT rho = conic0.r;
@end

(ii) in the terminology of the previous subsubsection, 
the second task is to decide the sign of 
$$\frac{\partial}{\partial\tau}\det
\left(\frac{M(\tau)}{z(\tau)}\right)\mid_{\tau=0},$$
where 
$$M(\tau) = \left(\begin{array}{cc} r(\tau) & t(\tau) \\ t(\tau) & s(\tau)
\end{array}\right)$$ and 

$$z(\tau)=\left(u(\tau), v(\tau)\right)M(\tau)^{-1}
\left(u(\tau), v(\tau)\right)^T - w(\tau).$$

For the sake of readability, we omit the parameter $\tau$ 
in the sequel. Noting that 
$$M^{-1} = \frac{1}{\det(M)}\left(\begin{array}{cc}s & - t \\ 
- t & r\end{array}\right),$$
we get
$$z = \frac{1}{\det(M)} (u^2s-2uvt+v^2r) - w.$$

Let us introduce the following abbreviations.
\begin{eqnarray*}
d &:=& \det(M), \\
Z &:=& u^2s-2uvt+v^2r.
\end{eqnarray*}
With primes ($d',Z'$ etc.) we denote derivatives w.r.t. $\tau$.
Now we can write $$\frac{\partial}{\partial \tau}\det(M/z) = (d/z^2)' =
\frac{d'z^2 - 2dzz'}{z^4} = \frac{d'z -2dz'}{z^3}.$$

Since $d(0),z(0)>0$ (recall that ${\cal C}_0$ is an ellipse), this is
equal in sign to $$f := d(d'z-2dz'),$$
at least for $\tau=0$ which is the value we are interested in. 
Now we have
\begin{eqnarray*}
d'z &=& d'(\frac{1}{d}Z-w) = \frac{d'}{d}Z-d'w, \\
dz' &=& d(\frac{Z'd - Zd'}{d^2} - w') = \frac{Z'd - Zd'}{d} - dw'.
\end{eqnarray*}
Hence 
\begin{eqnarray*} 
f    &=& d'Z-dd'w - 2(Z'd -Zd' - d^2w') \\
     &=& 3d'Z + d(2dw'-d'w - 2Z').
\end{eqnarray*}
Now it's time to look at the explicit expressions for $d,d',Z,Z'$. 
Rewriting $Z$ as
\begin{eqnarray*} 
Z &=& u (us - vt) + v (vr-ut) \\
  &=& uZ_1+vZ_2,
\end{eqnarray*} we get 
\begin{eqnarray*}
d   &=& rs-t^2, \\
d'  &=& r's + rs' - 2tt', \\
Z_1' &=& u's + us' - v't -vt', \\
Z_2' &=& v'r + vr'  -u't -ut', \\
Z'   &=& u'Z_1 + uZ_1' + v'Z_2 + v Z_2'.
\end{eqnarray*}
Evaluated for $\tau=0$, all these values can be computed from 
$r(0),s(0),t(0),u(0),v(0),w(0)$ (the defining values of the conic
${\cal C}_0$) and their corresponding primed values which have already
been computed by the \ccc{set} method. $d(0)=\det({\cal C}_0)$ is already
defined, the other values follow (recall that $r'(0)=0$). For all this,
${\cal C}_0$ is assumed to be normalized.

@macro <handle ellipse case> += @begin
    conic0.normalize();
    R::FT  dd  = conic0.r*ds - R::FT(2)*conic0.t*dt;
    R::FT  Z1  = conic0.u*conic0.s - conic0.v*conic0.t;
    R::FT  Z2  = conic0.v*conic0.r - conic0.u*conic0.t;
    R::FT  dZ1 = du*conic0.s + conic0.u*ds - dv*conic0.t - conic0.v*dt;
    R::FT  dZ2 = dv*conic0.r - du*conic0.t - conic0.u*dt;
    R::FT  Z   = conic0.u*Z1 + conic0.v*Z2;
    R::FT  dZ  = du*Z1 + conic0.u*dZ1 +dv*Z2 + conic0.v*dZ2;
    R::FT  f   = R::FT(3)*dd*Z + d*(R::FT(2)*d*dw-dd*conic0.w-R::FT(2)*dZ);
@end

The sign of $f(0)$ is the one we are interested in. We deduce 
from \cite{gs-epsee-97,gs-seefe-97} that
$$p \left\{\begin{array}{c} {\rm inside} \\ {\rm on} \\ {\rm outside}
\end{array}\right\} E^*
\Leftrightarrow \rho~f(0) \left\{\begin{array}{c}< 0 \\ = 0 \\ > 0
\end{array}\right.,$$
$E^{*}$ the smallest ellipse through $p_1,p_2,p_3,p_4$. 

@macro <handle ellipse case> += @begin
    return static_cast (CGAL_Bounded_side, CGAL_sign (rho)*CGAL_sign (f));
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_5 Class}
@! ----------------------------------------------------------------------------

This class is the simplest one among the three representation classes, since
it does not need to compute its ellipse. Recall that the current support set
$B$ attains cardinality five only if previously, a point
$p$ was found to lie outside $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, 
$\{p_1,p_2,p_3,p_4\}$ the support set at that time. However, during this test
with $p$, the unique conic (which is then an ellipse)  through 
$\{p_1,p_2,p_3,p_4,p\}$ has already been computed by the \ccc{Ellipse_4}
representation. In fact, it is the conic ${\cal C}_0$ addressed in the
previous subsection, and it suffices to store a pointer to it in the
\ccc{Ellipse_5} object, which is then initialized by the \ccc{set} method.

@macro <Ellipse_5 data members> zero = @begin
    Conic<R>* ellipse;
@end

@macro <Ellipse_5::set body> zero = @begin
    ellipse = &(ellipse_4.conic0);
@end

The \ccc{bounded_side} method is then straightforward, noting that the
ellipse has already been normalized by \ccc{ellipse_4}.

@macro <Ellipse_5::bounded_side body> zero = @begin
    R::FT discr = (*ellipse).eval(p);
    return static_cast (CGAL_Bounded_side, CGAL_sign (discr));
@end 
