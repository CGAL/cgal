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

\newcommand{\linebreakByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

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
maintaining the smallest enclosing ellipse. It is parameterized with a
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
\input{../../doc_tex/basic/Optimisation/Min_ellipse_2.tex}
\input{../../doc_tex/basic/Optimisation/Optimisation_ellipse_2.tex}
\input{../../doc_tex/basic/Optimisation/Min_ellipse_2_adapterC2.tex}
\input{../../doc_tex/basic/Optimisation/Min_ellipse_2_adapterH2.tex}

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
        @<Min_ellipse_2 validity check>

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
    Traits       tco;                           // traits class object
@end

The points of $P$ are internally stored as a linked list that allows to
bring points to the front of the list in constant time. We use the
sequence container \ccc{list} from STL~\cite{sl-stl-95}.

@macro <Min_ellipse_2 private data members> += @begin
    list<Point>  points;                        // doubly linked list of points
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
    int          n_support_points;              // number of support points
    Point*       support_points;                // array of support points
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
        points.push_back( p5);

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
\subsubsection{Validity Check}

A \ccc{CGAL_Min_ellipse_2<Traits>} object can be checked for validity.
This means, it is checked whether (a) the ellipse contains all points
of its defining set $P$, (b) the ellipse is the smallest ellipse
spanned by its support set, and (c) the support set is minimal, i.e.\
no support point is redundant. (\emph{Note:} (b) and (c) are not yet
implemented. Instead we check if the support set lies on the boundary
of the ellipse.) The function \ccc{is_valid} is mainly intended for
debugging user supplied traits classes but also for convincing the
anxious user that the traits class implementation is correct. If
\ccc{verbose} is \ccc{true}, some messages concerning the performed
checks are written to standard error stream. The second parameter
\ccc{level} is not used, we provide it only for consistency with
interfaces of other classes.

@macro <Min_ellipse_2 validity check> = @begin
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        CGAL_Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL_Min_ellipse_2<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;

        // containment check (a)
        @<Min_ellipse_2 containment check>

        // support set checks (b)+(c) (not yet implemented)
        @!<Min_ellipse_2 support set checks>

        // alternative support set check
        @<Min_ellipse_2 support set check>

        verr << "  object is valid!" << endl;
        return( true);
    }
@end

The containment check (a) is easy to perform, just a loop over all
points in \ccc{points}.

@macro <Min_ellipse_2 containment check> = @begin
    verr << "  a) containment check..." << flush;
    Point_iterator point_iter;
    for ( point_iter  = points_begin();
          point_iter != points_end();
          ++point_iter)
        if ( has_on_unbounded_side( *point_iter)) 
            return( CGAL__optimisation_is_valid_fail( verr,
                        "ellipse does not contain all points"));
    verr << "passed." << endl;
@end

The alternative support set check is easy to perform, just a loop over all
support points in \ccc{support_points}.

@macro <Min_ellipse_2 support set check> = @begin
    verr << "  +) support set check..." << flush;
    Support_point_iterator support_point_iter;
    for ( support_point_iter  = support_points_begin();
          support_point_iter != support_points_end();
          ++support_point_iter)
        if ( ! has_on_boundary( *support_point_iter)) 
            return( CGAL__optimisation_is_valid_fail( verr,
                        "ellipse does not have all \
                         support points on the boundary"));
    verr << "passed." << endl;
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
                if ( point_iter != points.begin()) {            // p not first?
                    points.push_front( p);
                    points.erase( point_iter++); }
                else
                    ++point_iter; }
            else
                ++point_iter; }
    }
@end

@! ----------------------------------------------------------------------------
@! Class template CGAL_Optimisation_ellipse_2<R>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont CGAL\_Optimisation\_ellipse\_2<R>}

First, we declare the class template \ccc{CGAL_Optimisation_ellipse_2},

@macro<Optimisation_ellipse_2 declaration> = @begin
    template < class _R >
    class CGAL_Optimisation_ellipse_2;

    class ostream;
    class istream;
    class CGAL_Window_stream;
@end

\emph{Workaround:} The GNU compiler (g++ 2.7.2[.?]) does not accept types
with scope operator as argument type or return type in class template
member functions. Therefore, all member functions are implemented in
the class interface.

The class interface looks as follows.

@macro <Optimisation_ellipse_2 interface> = @begin
    template < class _R >
    class CGAL_Optimisation_ellipse_2 {
        friend  ostream&  operator << CGAL_NULL_TMPL_ARGS (
            ostream&, CGAL_Optimisation_ellipse_2<_R> const&);
        friend  istream&  operator << CGAL_NULL_TMPL_ARGS (
            istream&, CGAL_Optimisation_ellipse_2<_R> &);
        friend  CGAL_Window_stream& operator << CGAL_NULL_TMPL_ARGS (
            CGAL_Window_stream&, CGAL_Optimisation_ellipse_2<_R> const&);
      public:
        @<Optimisation_ellipse_2 public interface>

      private:
        // private data members
        @<Optimisation_ellipse_2 private data members>

    @<dividing line>

    // Class implementation
    // ====================

      public:
        // Set functions
        // -------------
        @<Optimisation_ellipse_2 set functions>

        // Access functions
        // ----------------
        @<Optimisation_ellipse_2 access functions>

        // Equality tests
        // --------------
        @<Optimisation_ellipse_2 equality tests>

        // Predicates
        // ----------
        @<Optimisation_ellipse_2 predicates>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Public Interface}

The functionality is described and documented in the specification
section, so we do not comment on it here.

@macro <Optimisation_ellipse_2 public interface> = @begin
    // types
    typedef           _R               R;
    typedef           typename _R::RT  RT;
    typedef           typename _R::FT  FT;
    typedef           CGAL_Point_2<R>  Point;
    typedef           CGAL_Conic_2<R>  Conic;

    /**************************************************************************
    WORKAROUND: The GNU compiler (g++ 2.7.2[.*]) does not accept types
    with scope operator as argument type or return type in class template
    member functions. Therefore, all member functions are implemented in
    the class interface.

    // creation
    void  set( );
    void  set( Point const& p);
    void  set( Point const& p, Point const& q);
    void  set( Point const& p1, Point const& p2, Point const& p3);
    void  set( Point const& p1, Point const& p2,
               Point const& p3, Point const& p4);
    void  set( Point const& p1, Point const& p2,
               Point const& p3, Point const& p4, Point const& p5);

    // access functions    
    int  number_of_boundary_points()

    // equality tests
    bool  operator == ( CGAL_Optimisation_ellipse_2<R> const& e) const;
    bool  operator != ( CGAL_Optimisation_ellipse_2<R> const& e) const;

    // predicates
    CGAL_Bounded_side  bounded_side( Point const& p) const;
    bool  has_on_bounded_side      ( Point const& p) const;
    bool  has_on_boundary          ( Point const& p) const;
    bool  has_on_unbounded_side    ( Point const& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
    **************************************************************************/
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Data Members}

The representation of the ellipse depends on the number of given
boundary points, stored in \ccc{n_boundary_points}.

@macro <Optimisation_ellipse_2 private data members> += @begin
    int    n_boundary_points;                   // number of boundary points
@end

In the degenerate cases with zero to two boundary points, the given
points are stored directly in \ccc{boundary_point1} and
\ccc{boundary_point2}, resp.

@macro <Optimisation_ellipse_2 private data members> += @begin
    Point  boundary_point1, boundary_point2;    // two boundary points
@end

Given three or five points, the ellipse is represented as a conic,
using the class \ccc{CGAL_Conic_2<R>}. The case with four boundary
points is the most complicated one, since in general a direct
representation with one conic has irrational
coordinates~\cite{gs-seefe-97a}. Therefore the ellipse is represented
implicitly as a linear combination of two conics.

@macro <Optimisation_ellipse_2 private data members> += @begin
    Conic  conic1, conic2;                      // two conics
@end

Finally, in the case of four boundary points, we need the gradient
vector of the linear combination for the volume derivative in the
in-ellipse test.

@macro <Optimisation_ellipse_2 private data members> += @begin
    RT     dr, ds, dt, du, dv, dw;              // the gradient vector
@end

@! ----------------------------------------------------------------------------
\subsubsection{Set Functions}

We provide set functions taking zero, one, two, three, four or five
boundary points. They all set the variable to the smallest ellipse
through the given points. \emph{Note:} The set function taking five
boundary points only uses the fifth point from its input together with
the two internally represented conics to compute the ellipse. The
algorithm in Section~\ref{sec:algo} garantees that this set function
is only called an ellipse that already have the first four points as
its boundary points.

@macro <Optimisation_ellipse_2 set functions> = @begin
    inline
    void
    set( )
    {
        n_boundary_points = 0;
    }
    
    inline
    void
    set( Point const& p)
    {
        n_boundary_points = 1;
        boundary_point1   = p;
    }
    
    inline
    void
    set( Point const& p, Point const& q)
    {
        n_boundary_points = 2;
        boundary_point1   = p;
        boundary_point2   = q;
    }
    
    inline
    void
    set( Point const& p1, Point const& p2, Point const& p3)
    {
        n_boundary_points = 3;
        conic1.set_ellipse( p1, p2, p3);
    }
    
    inline
    void
    set( Point const& p1, Point const& p2, Point const& p3, Point const& p4)
    {
        n_boundary_points = 4;
        Conic::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);
        dr = RT( 0);
        ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
        dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
        du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
        dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
        dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
    }
    
    inline
    void
    set( Point const&, Point const&,
         Point const&, Point const&, Point const& p5)
    {
        n_boundary_points = 5;
        conic1.set( conic1, conic2, p5);
        conic1.analyse();
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Access Functions}

@macro <Optimisation_ellipse_2 access functions> = @begin
    inline
    int
    number_of_boundary_points( ) const
    {
        return( n_boundary_points);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Equality Tests}

@macro <Optimisation_ellipse_2 equality tests> = @begin
    bool
    operator == ( CGAL_Optimisation_ellipse_2<R> const& e) const
    {
        if ( n_boundary_points != e.n_boundary_points)
            return( false);

        switch ( n_boundary_points) {
          case 0:
            return( true);
          case 1:
            return( boundary_point1 == e.boundary_point1);
          case 2:
            return(    (    ( boundary_point1 == e.boundary_point1)
                         && ( boundary_point2 == e.boundary_point2))
                    || (    ( boundary_point1 == e.boundary_point2)
                         && ( boundary_point2 == e.boundary_point1)));
          case 3:
          case 5:
            return( conic1 == e.conic1);
          case 4:
            return(    (    ( conic1 == e.conic1)
                         && ( conic2 == e.conic2))
                    || (    ( conic1 == e.conic2)
                         && ( conic2 == e.conic1)));
          default:
            CGAL_optimisation_assertion(    ( n_boundary_points >= 0)
                                         && ( n_boundary_points <= 5)); }
        // keeps g++ happy
        return( false);
    }
    
    inline
    bool
    operator != ( CGAL_Optimisation_ellipse_2<R> const& e) const
    {
        return( ! operator == ( e));
    }
@end
    
@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The following predicates perform in-ellipse tests and check for
emptyness and degeneracy, resp. The way to evaluate the in-ellipse
test depends on the number of boundary points and is realised by a
case analysis. Again, the case with four points is the most difficult
one.

@macro <Optimisation_ellipse_2 predicates> = @begin
    inline
    CGAL_Bounded_side
    bounded_side( Point const& p) const
    {
        switch ( n_boundary_points) {
          case 0:
            return( CGAL_ON_UNBOUNDED_SIDE);
          case 1:
            return( ( p == boundary_point1) ?
                           CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
          case 2:
            return(    ( p == boundary_point1)
                    || ( p == boundary_point2)
                    || ( CGAL_are_ordered_along_line(
                             boundary_point1, p, boundary_point2)) ?
                         CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
          case 3:
          case 5:
            return( conic1.convex_side( p));
          case 4: {
            Conic c;
            c.set( conic1, conic2, p);
            c.analyse();
            if ( ! c.is_ellipse()) {
                c.set_ellipse( conic1, conic2);
                c.analyse();
                return( c.convex_side( p)); }
            else {
                int tau_star = -c.vol_derivative( dr, ds, dt, du, dv, dw);
                return( CGAL_static_cast( CGAL_Bounded_side,
                                          CGAL_sign( tau_star))); } }
          default:
            CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                         ( n_boundary_points <= 5) ); }
        // keeps g++ happy
        return( CGAL_Bounded_side( 0));
    }

    inline
    bool
    has_on_bounded_side( Point const& p) const
    {
        return( bounded_side( p) == CGAL_ON_BOUNDED_SIDE);
    }

    inline
    bool
    has_on_boundary( Point const& p) const
    {
        return( bounded_side( p) == CGAL_ON_BOUNDARY);
    }

    inline
    bool
    has_on_unbounded_side( Point const& p) const
    {
        return( bounded_side( p) == CGAL_ON_UNBOUNDED_SIDE);
    }

    inline
    bool
    is_empty( ) const
    {
        return( n_boundary_points == 0);
    }

    inline
    bool
    is_degenerate( ) const
    {
        return( n_boundary_points < 3);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Optimisation_ellipse_2 I/O operators declaration> = @begin
    template < class _R >
    ostream&
    operator << ( ostream& os, CGAL_Optimisation_ellipse_2<_R> const& e);

    template < class _R >
    istream&
    operator >> ( istream& is, CGAL_Optimisation_ellipse_2<_R>      & e);
@end

@macro <Optimisation_ellipse_2 I/O operators> = @begin
    template < class _R >
    ostream&
    operator << ( ostream& os, CGAL_Optimisation_ellipse_2<_R> const& e)
    {
        const char* const  empty       = "";
        const char* const  pretty_head = "CGAL_Optimisation_ellipse_2( ";
        const char* const  pretty_sep  = ", ";
        const char* const  pretty_tail = ")";
        const char* const  ascii_sep   = " ";

        const char*  head = empty;
        const char*  sep  = empty;
        const char*  tail = empty;

        switch ( CGAL_get_mode( os)) {
          case CGAL_IO::PRETTY:
            head = pretty_head;
            sep  = pretty_sep;
            tail = pretty_tail;
            break;
          case CGAL_IO::ASCII:
            sep  = ascii_sep;
            break;
          case CGAL_IO::BINARY:
            break;
          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL_get_mode( os) invalid!");
            break; }

        os << head << e.n_boundary_points;
        switch ( e.n_boundary_points) {
          case 0:
            break;
          case 1:
            os << sep << e.boundary_point1;
            break;
          case 2:
            os << sep << e.boundary_point1
               << sep << e.boundary_point2;
            break;
          case 3:
          case 5:
            os << sep << e.conic1;
            break;
          case 4:
            os << sep << e.conic1
               << sep << e.conic2;
            break; }
        os << tail;

        return( os);
    }

    template < class _R >
    istream&
    operator >> ( istream& is, CGAL_Optimisation_ellipse_2<_R>& e)
    {
        switch ( CGAL_get_mode( is)) {

          case CGAL_IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL_IO::ASCII:
          case CGAL_IO::BINARY:
            CGAL_read( is, e.n_boundary_points);
            switch ( e.n_boundary_points) {
              case 0:
                break;
              case 1:
                is >> e.boundary_point1;
                break;
              case 2:
                is >> e.boundary_point1
                   >> e.boundary_point2;
                break;
              case 3:
              case 5:
                is >> e.conic1;
                break;
              case 4:
                is >> e.conic1
                   >> e.conic2;
                break; }
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL_get_mode( is) invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
@! Class template CGAL_Min_ellipse_2_adapterC2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont CGAL\_Min\_ellipse\_2\_adapterC2<PT,DA>}

First, we declare the class templates \ccc{CGAL_Min_ellipse_2},
\ccc{CGAL_Min_ellipse_2_adapterC2} and
\ccc{CGAL__Min_ellipse_2_adapterC2__Ellipse}.

@macro<Min_ellipse_2_adapterC2 declarations> = @begin
    template < class _Traits >
    class CGAL_Min_ellipse_2;

    template < class _PT, class _DA >
    class CGAL_Min_ellipse_2_adapterC2;

    template < class _PT, class _DA >
    class CGAL__Min_ellipse_2_adapterC2__Ellipse;
@end

The actual work of the adapter is done in the nested class
\ccc{Ellipse}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{ellipse} containing the current ellipse is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current ellipse,
it is declared \ccc{friend}.

@macro <Min_ellipse_2_adapterC2 interface and implementation> = @begin
    template < class _PT, class _DA >
    class CGAL_Min_ellipse_2_adapterC2 {
      public:
        // types
        typedef  _PT  PT;
        typedef  _DA  DA;

        // nested types
        typedef  PT                                             Point;
        typedef  CGAL__Min_ellipse_2_adapterC2__Ellipse<PT,DA>  Ellipse;

      private:
        DA      dao;                                    // data accessor object
        Ellipse ellipse;                                // current ellipse
        friend class CGAL_Min_ellipse_2< CGAL_Min_ellipse_2_adapterC2<PT,DA> >;

      public:
        // creation
        @<Min_ellipse_2_adapterC2 constructors>

        // operations
        @<Min_ellipse_2_adapterC2 operations>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Constructors}

@macro <Min_ellipse_2_adapterC2 constructors> = @begin
    CGAL_Min_ellipse_2_adapterC2( DA const& da = DA())
        : dao( da), ellipse( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_ellipse_2_adapterC2 operations> = @begin
    CGAL_Orientation
    orientation( Point const& p, Point const& q, Point const& r) const
    {
        typedef  typename _DA::FT  FT;

        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;
        
        dao.get( p, px, py);
        dao.get( q, qx, qy);
        dao.get( r, rx, ry);

        return( CGAL_static_cast( CGAL_Orientation,
                    CGAL_sign( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Ellipse}

@macro <Min_ellipse_2_adapterC2 nested type `Ellipse'> = @begin
    template < class _PT, class _DA >
    class CGAL__Min_ellipse_2_adapterC2__Ellipse {
      public:
        // typedefs
        typedef  _PT  PT;
        typedef  _DA  DA;

        typedef           CGAL_ConicCPA2< PT, DA>  CT;
        typedef  typename _DA::FT                  FT;

      private:
        // data members
        int  n_boundary_points;                 // number of boundary points
        PT   boundary_point1, boundary_point2;  // two boundary points
        CT   conic1, conic2;                    // two conics
        FT   dr, ds, dt, du, dv, dw;            // the gradient vector

      public:
        // types
        typedef  PT  Point;

        // creation
        CGAL__Min_ellipse_2_adapterC2__Ellipse( DA const& da)
            : conic1( da), conic2( da)
        { }

        void
        set( )
        {
            n_boundary_points = 0;
        }
    
        void
        set( Point const& p)
        {
            n_boundary_points = 1;
            boundary_point1   = p;
        }
    
        void
        set( Point const& p, Point const& q)
        {
            n_boundary_points = 2;
            boundary_point1   = p;
            boundary_point2   = q;
        }
    
        void
        set( Point const& p1, Point const& p2, Point const& p3)
        {
            n_boundary_points = 3;
            conic1.set_ellipse( p1, p2, p3);
        }
    
        void
        set( Point const& p1, Point const& p2,
             Point const& p3, Point const& p4)
        {
            n_boundary_points = 4;
            CT::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);
            dr = FT( 0);
            ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
            dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
            du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
            dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
            dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
        }
    
        void
        set( Point const&, Point const&,
             Point const&, Point const&, Point const& p5)
        {
            n_boundary_points = 5;
            conic1.set( conic1, conic2, p5);
            conic1.analyse();
        }

        // predicates
        CGAL_Bounded_side
        bounded_side( Point const& p) const
        {
            switch ( n_boundary_points) {
              case 0:
                return( CGAL_ON_UNBOUNDED_SIDE);
              case 1:
                return( ( p == boundary_point1) ?
                               CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
              case 2:
                return(    ( p == boundary_point1)
                        || ( p == boundary_point2)
                        || ( CGAL_are_ordered_along_lineC2( boundary_point1, p,
                                               boundary_point2, conic1.da())) ?
                                    CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
              case 3:
              case 5:
                return( conic1.convex_side( p));
              case 4: {
                CT c( conic1.da());
                c.set( conic1, conic2, p);
                c.analyse();
                if ( ! c.is_ellipse()) {
                    c.set_ellipse( conic1, conic2);
                    c.analyse();
                    return( c.convex_side( p)); }
                else {
                    int tau_star = -c.vol_derivative( dr, ds, dt, du, dv, dw);
                    return( CGAL_static_cast( CGAL_Bounded_side,
                                              CGAL_sign( tau_star))); } }
              default:
                CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                             ( n_boundary_points <= 5) ); }
            // keeps g++ happy
            return( CGAL_Bounded_side( 0));
        }

        bool
        has_on_bounded_side( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_BOUNDED_SIDE);
        }

        bool
        has_on_boundary( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_BOUNDARY);
        }

        bool
        has_on_unbounded_side( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_UNBOUNDED_SIDE);
        }

        bool
        is_empty( ) const
        {
            return( n_boundary_points == 0);
        }

        bool
        is_degenerate( ) const
        {
            return( n_boundary_points < 3);
        }

        // additional operations for checking
        bool
        operator == (
            CGAL__Min_ellipse_2_adapterC2__Ellipse<_PT,_DA> const& e) const
        {
            if ( n_boundary_points != e.n_boundary_points)
                return( false);

            switch ( n_boundary_points) {
              case 0:
                return( true);
              case 1:
                return( boundary_point1 == e.boundary_point1);
              case 2:
                return(    (    ( boundary_point1 == e.boundary_point1)
                             && ( boundary_point2 == e.boundary_point2))
                        || (    ( boundary_point1 == e.boundary_point2)
                             && ( boundary_point2 == e.boundary_point1)));
              case 3:
              case 5:
                return( conic1 == e.conic1);
              case 4:
                return(    (    ( conic1 == e.conic1)
                             && ( conic2 == e.conic2))
                        || (    ( conic1 == e.conic2)
                             && ( conic2 == e.conic1)));
              default:
                CGAL_optimisation_assertion(    ( n_boundary_points >= 0)
                                             && ( n_boundary_points <= 5)); }
            // keeps g++ happy
            return( false);
        }

        // I/O
        friend
        ostream&
        operator << ( ostream& os,
                      CGAL__Min_ellipse_2_adapterC2__Ellipse<_PT,_DA> const& e)
        {
            const char* const  empty       = "";
            const char* const  pretty_head =
                                     "CGAL_Min_ellipse_2_adapterC2::Ellipse( ";
            const char* const  pretty_sep  = ", ";
            const char* const  pretty_tail = ")";
            const char* const  ascii_sep   = " ";

            const char*  head = empty;
            const char*  sep  = empty;
            const char*  tail = empty;

            switch ( CGAL_get_mode( os)) {
              case CGAL_IO::PRETTY:
                head = pretty_head;
                sep  = pretty_sep;
                tail = pretty_tail;
                break;
              case CGAL_IO::ASCII:
                sep  = ascii_sep;
                break;
              case CGAL_IO::BINARY:
                break;
              default:
                CGAL_optimisation_assertion_msg( false,
                                                "CGAL_get_mode( os) invalid!");
                break; }

            os << head << e.n_boundary_points;
            switch ( e.n_boundary_points) {
              case 0:
                break;
              case 1:
                os << sep << e.boundary_point1;
                break;
              case 2:
                os << sep << e.boundary_point1
                   << sep << e.boundary_point2;
                break;
              case 3:
              case 5:
                os << sep << e.conic1;
                break;
              case 4:
                os << sep << e.conic1
                   << sep << e.conic2;
                break; }
            os << tail;

            return( os);
        }

        friend
        istream&
        operator >> ( istream& is,
                      CGAL__Min_ellipse_2_adapterC2__Ellipse<_PT,_DA>& e)
        {
            switch ( CGAL_get_mode( is)) {

              case CGAL_IO::PRETTY:
                cerr << endl;
                cerr << "Stream must be in ascii or binary mode" << endl;
                break;

              case CGAL_IO::ASCII:
              case CGAL_IO::BINARY:
                CGAL_read( is, e.n_boundary_points);
                switch ( e.n_boundary_points) {
                  case 0:
                    break;
                  case 1:
                    is >> e.boundary_point1;
                    break;
                  case 2:
                    is >> e.boundary_point1
                       >> e.boundary_point2;
                    break;
                  case 3:
                  case 5:
                    is >> e.conic1;
                    break;
                  case 4:
                    is >> e.conic1
                       >> e.conic2;
                    break; }
                break;

              default:
                CGAL_optimisation_assertion_msg( false,
                                                 "CGAL_IO::mode invalid!");
                break; }

            return( is);
        }
    };
@end

@! ----------------------------------------------------------------------------
@! Class template CGAL_Min_ellipse_2_adapterH2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont CGAL\_Min\_ellipse\_2\_adapterH2<PT,DA>}

First, we declare the class templates \ccc{Min_ellipse_2},
\ccc{CGAL_Min_ellipse_2_adapterH2} and
\ccc{CGAL__Min_ellipse_2_adapterH2__Ellipse}.

@macro<Min_ellipse_2_adapterH2 declarations> = @begin
    template < class _Traits >
    class CGAL_Min_ellipse_2;

    template < class _PT, class _DA >
    class CGAL_Min_ellipse_2_adapterH2;

    template < class _PT, class _DA >
    class CGAL__Min_ellipse_2_adapterH2__Ellipse;
@end

The actual work of the adapter is done in the nested class
\ccc{Ellipse}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{ellipse} containing the current ellipse is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current ellipse,
it is declared \ccc{friend}.

@macro <Min_ellipse_2_adapterH2 interface and implementation> = @begin
    template < class _PT, class _DA >
    class CGAL_Min_ellipse_2_adapterH2 {
      public:
        // types
        typedef  _PT  PT;
        typedef  _DA  DA;

        // nested types
        typedef  PT                                             Point;
        typedef  CGAL__Min_ellipse_2_adapterH2__Ellipse<PT,DA>  Ellipse;

      private:
        DA      dao;                                    // data accessor object
        Ellipse ellipse;                                // current ellipse
        friend class CGAL_Min_ellipse_2< CGAL_Min_ellipse_2_adapterH2<PT,DA> >;

      public:
        // creation
        @<Min_ellipse_2_adapterH2 constructors>

        // operations
        @<Min_ellipse_2_adapterH2 operations>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Constructors}

@macro <Min_ellipse_2_adapterH2 constructors> = @begin
    CGAL_Min_ellipse_2_adapterH2( DA const& da = DA())
        : dao( da), ellipse( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_ellipse_2_adapterH2 operations> = @begin
    CGAL_Orientation
    orientation( Point const& p, Point const& q, Point const& r) const
    {
        typedef  typename _DA::RT  RT;

        RT  phx;
        RT  phy;
        RT  phw;
        RT  qhx;
        RT  qhy;
        RT  qhw;
        RT  rhx;
        RT  rhy;
        RT  rhw;
        
        dao.get( p, phx, phy, phw);
        dao.get( q, qhx, qhy, qhw);
        dao.get( r, rhx, rhy, rhw);

        return( CGAL_static_cast( CGAL_Orientation,
                    CGAL_sign( ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                             - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Ellipse}

@macro <Min_ellipse_2_adapterH2 nested type `Ellipse'> = @begin
    template < class _PT, class _DA >
    class CGAL__Min_ellipse_2_adapterH2__Ellipse {
      public:
        // typedefs
        typedef  _PT  PT;
        typedef  _DA  DA;

        typedef           CGAL_ConicHPA2< PT, DA>  CT;
        typedef  typename _DA::RT                  RT;

      private:
        // data members
        int  n_boundary_points;                 // number of boundary points
        PT   boundary_point1, boundary_point2;  // two boundary points
        CT   conic1, conic2;                    // two conics
        RT   dr, ds, dt, du, dv, dw;            // the gradient vector

      public:
        // types
        typedef  PT  Point;

        // creation
        CGAL__Min_ellipse_2_adapterH2__Ellipse( DA const& da)
            : conic1( da), conic2( da)
        { }

        void
        set( )
        {
            n_boundary_points = 0;
        }
    
        void
        set( Point const& p)
        {
            n_boundary_points = 1;
            boundary_point1   = p;
        }
    
        void
        set( Point const& p, Point const& q)
        {
            n_boundary_points = 2;
            boundary_point1   = p;
            boundary_point2   = q;
        }
    
        void
        set( Point const& p1, Point const& p2, Point const& p3)
        {
            n_boundary_points = 3;
            conic1.set_ellipse( p1, p2, p3);
        }
    
        void
        set( Point const& p1, Point const& p2,
             Point const& p3, Point const& p4)
        {
            n_boundary_points = 4;
            CT::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);
            dr = RT( 0);
            ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
            dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
            du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
            dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
            dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
        }
    
        void
        set( Point const&, Point const&,
             Point const&, Point const&, Point const& p5)
        {
            n_boundary_points = 5;
            conic1.set( conic1, conic2, p5);
            conic1.analyse();
        }

        // predicates
        CGAL_Bounded_side
        bounded_side( Point const& p) const
        {
            switch ( n_boundary_points) {
              case 0:
                return( CGAL_ON_UNBOUNDED_SIDE);
              case 1:
                return( ( p == boundary_point1) ?
                               CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
              case 2:
                return(    ( p == boundary_point1)
                        || ( p == boundary_point2)
                        || ( CGAL_are_ordered_along_lineH2( boundary_point1, p,
                                               boundary_point2, conic1.da())) ?
                                    CGAL_ON_BOUNDARY : CGAL_ON_UNBOUNDED_SIDE);
              case 3:
              case 5:
                return( conic1.convex_side( p));
              case 4: {
                CT c( conic1.da());
                c.set( conic1, conic2, p);
                c.analyse();
                if ( ! c.is_ellipse()) {
                    c.set_ellipse( conic1, conic2);
                    c.analyse();
                    return( c.convex_side( p)); }
                else {
                    int tau_star = -c.vol_derivative( dr, ds, dt, du, dv, dw);
                    return( CGAL_static_cast( CGAL_Bounded_side,
                                              CGAL_sign( tau_star))); } }
              default:
                CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                             ( n_boundary_points <= 5) ); }
            // keeps g++ happy
            return( CGAL_Bounded_side( 0));
        }

        bool
        has_on_bounded_side( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_BOUNDED_SIDE);
        }

        bool
        has_on_boundary( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_BOUNDARY);
        }

        bool
        has_on_unbounded_side( Point const& p) const
        {
            return( bounded_side( p) == CGAL_ON_UNBOUNDED_SIDE);
        }

        bool
        is_empty( ) const
        {
            return( n_boundary_points == 0);
        }

        bool
        is_degenerate( ) const
        {
            return( n_boundary_points < 3);
        }

        // additional operations for checking
        bool
        operator == (
            CGAL__Min_ellipse_2_adapterH2__Ellipse<_PT,_DA> const& e) const
        {
            if ( n_boundary_points != e.n_boundary_points)
                return( false);

            switch ( n_boundary_points) {
              case 0:
                return( true);
              case 1:
                return( boundary_point1 == e.boundary_point1);
              case 2:
                return(    (    ( boundary_point1 == e.boundary_point1)
                             && ( boundary_point2 == e.boundary_point2))
                        || (    ( boundary_point1 == e.boundary_point2)
                             && ( boundary_point2 == e.boundary_point1)));
              case 3:
              case 5:
                return( conic1 == e.conic1);
              case 4:
                return(    (    ( conic1 == e.conic1)
                             && ( conic2 == e.conic2))
                        || (    ( conic1 == e.conic2)
                             && ( conic2 == e.conic1)));
              default:
                CGAL_optimisation_assertion(    ( n_boundary_points >= 0)
                                             && ( n_boundary_points <= 5)); }
            // keeps g++ happy
            return( false);
        }

        // I/O
        friend
        ostream&
        operator << ( ostream& os,
                      CGAL__Min_ellipse_2_adapterH2__Ellipse<_PT,_DA> const& e)
        {
            const char* const  empty       = "";
            const char* const  pretty_head =
                                     "CGAL_Min_ellipse_2_adapterH2::Ellipse( ";
            const char* const  pretty_sep  = ", ";
            const char* const  pretty_tail = ")";
            const char* const  ascii_sep   = " ";

            const char*  head = empty;
            const char*  sep  = empty;
            const char*  tail = empty;

            switch ( CGAL_get_mode( os)) {
              case CGAL_IO::PRETTY:
                head = pretty_head;
                sep  = pretty_sep;
                tail = pretty_tail;
                break;
              case CGAL_IO::ASCII:
                sep  = ascii_sep;
                break;
              case CGAL_IO::BINARY:
                break;
              default:
                CGAL_optimisation_assertion_msg( false,
                                                "CGAL_get_mode( os) invalid!");
                break; }

            os << head << e.n_boundary_points;
            switch ( e.n_boundary_points) {
              case 0:
                break;
              case 1:
                os << sep << e.boundary_point1;
                break;
              case 2:
                os << sep << e.boundary_point1
                   << sep << e.boundary_point2;
                break;
              case 3:
              case 5:
                os << sep << e.conic1;
                break;
              case 4:
                os << sep << e.conic1
                   << sep << e.conic2;
                break; }
            os << tail;

            return( os);
        }

        friend
        istream&
        operator >> ( istream& is,
                      CGAL__Min_ellipse_2_adapterH2__Ellipse<_PT,_DA>& e)
        {
            switch ( CGAL_get_mode( is)) {

              case CGAL_IO::PRETTY:
                cerr << endl;
                cerr << "Stream must be in ascii or binary mode" << endl;
                break;

              case CGAL_IO::ASCII:
              case CGAL_IO::BINARY:
                CGAL_read( is, e.n_boundary_points);
                switch ( e.n_boundary_points) {
                  case 0:
                    break;
                  case 1:
                    is >> e.boundary_point1;
                    break;
                  case 2:
                    is >> e.boundary_point1
                       >> e.boundary_point2;
                    break;
                  case 3:
                  case 5:
                    is >> e.conic1;
                    break;
                  case 4:
                    is >> e.conic1
                       >> e.conic2;
                    break; }
                break;

              default:
                CGAL_optimisation_assertion_msg( false,
                                                 "CGAL_IO::mode invalid!");
                break; }

            return( is);
        }
    };
@end

@! ============================================================================
@! Tests
@! ============================================================================

\clearpage
\section{Test}

We test \ccc{CGAL_Min_ellipse_2} with the traits class implementation
for optimisation algorithms, using exact arithmetic, i.e.\ Cartesian
representation with number type \ccc{CGAL_Quotient<CGAL_Gmpz>} or
\ccc{CGAL_Quotient<integer>} and homogeneous representation with
number type \ccc{CGAL_Gmpz} or \ccc{integer}.

@macro <Min_ellipse_2 test (includes and typedefs)> = @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Homogeneous.h>
    #include <CGAL/Optimisation_traits_2.h>
    #include <CGAL/Min_ellipse_2.h>
    #include <CGAL/Min_ellipse_2_adapterC2.h>
    #include <CGAL/Min_ellipse_2_adapterH2.h>
    #include <CGAL/IO/Verbose_ostream.h>
    #include <assert.h>
    #include <string.h>
    #include <fstream.h>

    #ifdef CGAL_USE_LEDA
    #  include <CGAL/leda_integer.h>
       typedef  leda_integer                     Rt;
       typedef  CGAL_Quotient< leda_integer >    Ft;
    #else
    #  include <CGAL/Gmpz.h>
       typedef  CGAL_Gmpz                        Rt;
       typedef  CGAL_Quotient< CGAL_Gmpz >       Ft;
    #endif

    typedef  CGAL_Cartesian< Ft >                RepC;
    typedef  CGAL_Homogeneous< Rt >              RepH;
    typedef  CGAL_Optimisation_traits_2< RepC >  TraitsC;
    typedef  CGAL_Optimisation_traits_2< RepH >  TraitsH;
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
        typedef  CGAL_Min_ellipse_2< Traits >   Min_ellipse;
        typedef  Min_ellipse::Point             Point;
        typedef  Min_ellipse::Ellipse           Ellipse;

        CGAL_Verbose_ostream verr( verbose);

        // generate `n' points at random
        const int    n = 20;
        CGAL_Random  random_x, random_y;
        Point        random_points[ n];
        int          i;
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
            bool               has_on_boundary;
            bool               has_on_unbounded_side;
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
@! Adapters
@! ----------------------------------------------------------------------------

\subsection{Traits Class Adapters}

We define two point classes (one with Cartesian, one with homogeneous
representation) and corresponding data accessors.

@macro <Min_ellipse_2 test (point classes)> = @begin
    // 2D Cartesian point class
    class MyPointC2 {
      public:
        typedef  ::Ft  FT;
      private:
        FT _x;
        FT _y;
      public:
        MyPointC2( ) { }
        MyPointC2( FT const& x, FT const& y) : _x( x), _y( y) { }

        FT const&  x( ) const { return( _x); }
        FT const&  y( ) const { return( _y); }

        bool
        operator == ( MyPointC2 const& p) const
        {
            return( ( _x == p._x) && ( _y == p._y));
        }

        friend
        ostream&
        operator << ( ostream& os, MyPointC2 const& p)
        {
            return( os << p._x << ' ' << p._y);
        }

        friend
        istream&
        operator >> ( istream& is, MyPointC2& p)
        {
            return( is >> p._x >> p._y);
        }
    };

    // 2D Cartesian point class data accessor
    class MyPointC2DA {
      public:
        typedef  ::Ft  FT;

        FT const&  get_x( MyPointC2 const& p) const { return( p.x()); }
        FT const&  get_y( MyPointC2 const& p) const { return( p.y()); }

        void
        get( MyPointC2 const& p, FT& x, FT& y) const
        {
            x = get_x( p);
            y = get_y( p);
        }

        void
        set( MyPointC2& p, FT const& x, FT const& y) const
        {
            p = MyPointC2( x, y);
        }
    };


    // 2D homogeneous point class
    class MyPointH2 {
      public:
        typedef  ::Rt  RT;
      private:
        RT _hx;
        RT _hy;
        RT _hw;
      public:
        MyPointH2( ) { }
        MyPointH2( RT const& hx, RT const& hy, RT const& hw = RT( 1))
            : _hx( hx), _hy( hy), _hw( hw) { }

        RT const&  hx( ) const { return( _hx); }
        RT const&  hy( ) const { return( _hy); }
        RT const&  hw( ) const { return( _hw); }

        bool
        operator == ( MyPointH2 const& p) const
        {
            return( ( _hx*p._hw == p._hx*_hw) && ( _hy*p._hw == p._hy*_hw));
        }

        friend
        ostream&
        operator << ( ostream& os, MyPointH2 const& p)
        {
            return( os << p._hx << ' ' << p._hy << ' ' << p._hw);
        }

        friend
        istream&
        operator >> ( istream& is, MyPointH2& p)
        {
            return( is >> p._hx >> p._hy >> p._hw);
        }
    };

    // 2D homogeneous point class data accessor
    class MyPointH2DA {
      public:
        typedef  ::Rt  RT;

        RT const&  get_hx( MyPointH2 const& p) const { return( p.hx()); }
        RT const&  get_hy( MyPointH2 const& p) const { return( p.hy()); }
        RT const&  get_hw( MyPointH2 const& p) const { return( p.hw()); }

        void
        get( MyPointH2 const& p, RT& hx, RT& hy, RT& hw) const
        {
            hx = get_hx( p);
            hy = get_hy( p);
            hw = get_hw( p);
        }

        void
        set( MyPointH2& p, RT const& hx, RT const& hy, RT const& hw) const
        {
            p = MyPointH2( hx, hy, hw);
        }
    };
@end

To test the traits class adapters we use the code coverage test function.

@macro <Min_ellipse_2 test (adapters test)> = @begin
    typedef  CGAL_Min_ellipse_2_adapterC2< MyPointC2, MyPointC2DA >  AdapterC2;
    typedef  CGAL_Min_ellipse_2_adapterH2< MyPointH2, MyPointH2DA >  AdapterH2;
    cover_Min_ellipse_2( verbose, AdapterC2(), Rt());
    cover_Min_ellipse_2( verbose, AdapterH2(), Rt());
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

        typedef  CGAL_Min_ellipse_2< TraitsH >  Min_ellipse;
        typedef  Min_ellipse::Point             Point;
        typedef  Min_ellipse::Ellipse           Ellipse;

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

@file <include/CGAL/Min_ellipse_2.h> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Min_ellipse_2.h")

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
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/optimisation_basic.h>
    #endif
    #ifndef CGAL_PROTECT_LIST_H
    #  include <list.h>
    #endif
    #ifndef CGAL_PROTECT_VECTOR_H
    #include <vector.h>
    #endif
    #ifndef CGAL_PROTECT_ALGO_H
    #include <algo.h>
    #endif
    #ifndef CGAL_PROTECT_IOSTREAM_H
    #include <iostream.h>
    #endif

    @<Min_ellipse_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators declaration>

    #ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #  include <CGAL/Min_ellipse_2.C>
    #endif

    #endif // CGAL_MIN_ELLIPSE_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2.C}

@file <include/CGAL/Min_ellipse_2.C> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Min_ellipse_2.C")

    // Class implementation (continued)
    // ================================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_ellipse_2.h
@! ----------------------------------------------------------------------------

\subsection{Optimisation\_ellipse\_2.h}

@file <include/CGAL/Optimisation_ellipse_2.h> = @begin
    @<Optimisation_ellipse_2 header>("include/CGAL/Optimisation_ellipse_2.h")

    #ifndef CGAL_OPTIMISATION_ELLIPSE_2_H
    #define CGAL_OPTIMISATION_ELLIPSE_2_H

    // Class declaration
    // =================
    @<Optimisation_ellipse_2 declaration>

    // Class interface
    // ===============
    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_CONIC_2_H
    #  include <CGAL/Conic_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/optimisation_assertions.h>
    #endif

    @<Optimisation_ellipse_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Optimisation_ellipse_2 I/O operators declaration>

    #ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #  include <CGAL/Optimisation_ellipse_2.C>
    #endif

    #endif // CGAL_OPTIMISATION_ELLIPSE_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{Optimisation\_ellipse\_2.C}

@file <include/CGAL/Optimisation_ellipse_2.C> = @begin
    @<Optimisation_ellipse_2 header>("include/CGAL/Optimisation_ellipse_2.C")

    // Class implementation (continued)
    // ================================
    // includes
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/optimisation_assertions.h>
    #endif

    // I/O
    // ---
    @<Optimisation_ellipse_2 I/O operators>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_adapterC2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_adapterC2.h}

@file <include/CGAL/Min_ellipse_2_adapterC2.h> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Min_ellipse_2_adapterC2.h")

    #ifndef CGAL_MIN_ELLIPSE_2_ADAPTERC2_H
    #define CGAL_MIN_ELLIPSE_2_ADAPTERC2_H

    // Class declarations
    // ==================
    @<Min_ellipse_2_adapterC2 declarations>

    // Class interface and implementation
    // ==================================
    // includes
    #ifndef CGAL_CONICCPA2_H
    #  include <CGAL/ConicCPA2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/optimisation_assertions.h>
    #endif

    template < class PT, class DA >
    bool
    CGAL_are_ordered_along_lineC2( PT const& p, PT const& q, PT const& r,
                                   DA const& da)
    {
        typedef  typename DA::FT  FT;

        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;
        
        da.get( p, px, py);
        da.get( q, qx, qy);
        da.get( r, rx, ry);

        // p,q,r collinear?
        if ( ! CGAL_is_zero( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx)))
            return( false);

        // p,q,r vertical?
        if ( px != rx)
            return(    ( ( px < qx) && ( qx < rx))
                    || ( ( rx < qx) && ( qx < px)));
        else
            return(    ( ( py < qy) && ( qy < ry))
                    || ( ( ry < qy) && ( qy < py)));
    }

    @<Min_ellipse_2_adapterC2 interface and implementation>

    // Nested type `Ellipse'
    @<Min_ellipse_2_adapterC2 nested type `Ellipse'>

    #endif // CGAL_MIN_ELLIPSE_2_ADAPTERC2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_adapterH2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_adapterH2.h}

@file <include/CGAL/Min_ellipse_2_adapterH2.h> = @begin
    @<Min_ellipse_2 header>("include/CGAL/Min_ellipse_2_adapterH2.h")

    #ifndef CGAL_MIN_ELLIPSE_2_ADAPTERH2_H
    #define CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

    // Class declarations
    // ==================
    @<Min_ellipse_2_adapterH2 declarations>

    // Class interface and implementation
    // ==================================
    // includes
    #ifndef CGAL_CONICHPA2_H
    #  include <CGAL/ConicHPA2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/optimisation_assertions.h>
    #endif

    template < class PT, class DA >
    bool
    CGAL_are_ordered_along_lineH2( PT const& p, PT const& q, PT const& r,
                                   DA const& da)
    {
        typedef  typename DA::RT  RT;

        RT  phx;
        RT  phy;
        RT  phw;
        RT  qhx;
        RT  qhy;
        RT  qhw;
        RT  rhx;
        RT  rhy;
        RT  rhw;
        
        da.get( p, phx, phy, phw);
        da.get( q, qhx, qhy, qhw);
        da.get( r, rhx, rhy, rhw);

        // p,q,r collinear?
        if ( ! CGAL_is_zero(   ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                             - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw)))
            return( false);

        // p,q,r vertical?
        if ( phx*rhw != rhx*phw)
            return(    ( ( phx*qhw < qhx*phw) && ( qhx*rhw < rhx*qhw))
                    || ( ( rhx*qhw < qhx*rhw) && ( qhx*phw < phx*qhw)));
        else
            return(    ( ( phy*qhw < qhy*phw) && ( qhy*rhw < rhy*qhw))
                    || ( ( rhy*qhw < qhy*rhw) && ( qhy*phw < phy*qhw)));
    }

    @<Min_ellipse_2_adapterH2 interface and implementation>

    // Nested type `Ellipse'
    @<Min_ellipse_2_adapterH2 nested type `Ellipse'>

    #endif // CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{test\_Min\_ellipse\_2.C}

@file <test/Optimisation/test_Min_ellipse_2.C> = @begin
    @<Min_ellipse_2 header>("test/Optimisation/test_Min_ellipse_2.C")

    @<Min_ellipse_2 test (includes and typedefs)>

    // code coverage test function
    // ---------------------------
    @<Min_ellipse_2 test (code coverage test function)>

    // point classes for adapters test
    // -------------------------------
    @<Min_ellipse_2 test (point classes)>

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

        // adapters test
        // -------------
        @<Min_ellipse_2 test (adapters test)>

        // external test sets
        // -------------------
        @<Min_ellipse_2 test (external test sets)>
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! File Header
@! ----------------------------------------------------------------------------

\subsection*{File Header}

@i ../file_header.awi
 
@macro <Min_ellipse_2 header>(1) many = @begin
    @<file header>("2D Smallest Enclosing Ellipse",@1,
                   "Optimisation/Min_ellipse_2",
                   "Sven Schönherr <sven@@inf.fu-berlin.de>",
                   "Bernd Gärtner",
                   "ETH Zurich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
                   "$Revision$","$Date$")
@end

@macro <Optimisation_ellipse_2 header>(1) many = @begin
    @<file header>("2D Optimisation Ellipse",@1,
                   "Optimisation/Min_ellipse_2",
                   "Sven Schönherr <sven@@inf.fu-berlin.de>",
                   "Bernd Gärtner",
                   "ETH Zurich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
                   "$Revision$","$Date$")
@end
@! ============================================================================
@! Bibliography
@! ============================================================================

\clearpage
\bibliographystyle{plain}
\bibliography{geom,cgal}

@! ===== EOF ==================================================================
