@! ============================================================================
@! The CGAL Library
@! Implementation: 2D Smallest Enclosing Ellipse
@! ----------------------------------------------------------------------------
@! file  : web/Min_ellipse_2.aw
@! author: Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Min_ellipse_2 WIP $
@! $Id$
@! $Date$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amssymb}
@usepackage{path}
@usepackage{cc_manual,cc_manual_index}
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

\newcommand{\linebreakByHand}{\ccTexHtml{\linebreak[4]}{}}
\newcommand{  \newlineByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsRevision}{$Id$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 5 mm
@t title titlefont centre "2D Smallest Enclosing Ellipse*"
@t vskip 1 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
\smallskip
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Chapter: Geometric Optimisation $+ \\
    \verb+$CGAL_Package: Min_ellipse_2 WIP+\cgalWIP\verb+$+ \\
  \end{tabular}
\end{center}
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
in the plane. The class template \ccc{Min_ellipse_2} is implemented
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
nonempty set $B$. Here is the pseudo-code of Welzl's method.  To compute
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

@! \clearpage
@! \section{Specifications}
@! 
@! \emph{Note:} Below some references are undefined, they refer to sections
@! in the \cgal\ Reference Manual.
@! 
@! \renewcommand{\ccFont}{\tt}
@! \renewcommand{\ccEndFont}{}
@! \newcommand{\cgalColumnLayout}{\ccTexHtml{%
@!   \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
@!   \ccPropagateThreeToTwoColumns}{}}
@! \newcommand{\cgalSetMinEllipseLayout}{%
@!   \ccSetThreeColumns{Support_point_iterator}{}{returns
@!     \ccc{ON_BOUNDED_SIDE}, \ccc{ON_BOUNDARY}, or \ccc{ON_UNBOUNDED_SIDE}}
@! %  \ccSetThreeColumns{Support_point_iterator}{}{creates a variable
@! %    \ccc{min_ellipse} of type \ccc{CGAL_Min_ellipse_2<Traits>}.}
@!   \ccPropagateThreeToTwoColumns}
@! \newcommand{\cgalSetOptTraitsAdaptLayout}{\ccTexHtml{%
@!     \ccSetThreeColumns{CGAL_Oriented_side}{}{returns constants
@!       \ccc{CGAL_LEFT_TURN}, \ccc{CGAL_COLLINEAR}}
@!     \ccPropagateThreeToTwoColumns}{}}
@! \input{../../doc_tex/basic/Optimisation/Min_ellipse_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Optimisation_ellipse_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_ellipse_2_traits_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_ellipse_2_adapterC2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_ellipse_2_adapterH2.tex}

@! ============================================================================
@! Implementations
@! ============================================================================

\clearpage
\section{Implementations}

@! ----------------------------------------------------------------------------
@! Class template Min_ellipse_2<Traits>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_ellipse\_2<Traits>}

First, we declare the class template \ccc{Min_ellipse_2}.

@macro<Min_ellipse_2 declaration> = @begin
    template < class Traits_ >
    class Min_ellipse_2;
@end

The actual work of the algorithm is done in the private member
functions \ccc{me} and \ccc{compute_ellipse}. The former directly
realizes the pseudo-code of $\me(P,B)$, the latter solves the basic
case $\me(\emptyset,B)$, see Section~\ref{sec:algo}.

The class interface looks as follows.

@macro <Min_ellipse_2 interface> = @begin
    template < class Traits_ >
    class Min_ellipse_2 {
      public:
        @<Min_ellipse_2 public interface>

      private:
        // private data members
        @<Min_ellipse_2 private data members>

        // copying and assignment not allowed!
        Min_ellipse_2( const Min_ellipse_2<Traits_>&);
        Min_ellipse_2<Traits_>& operator = ( const Min_ellipse_2<Traits_>&);

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
        // Private member functions
        // ------------------------
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
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Ellipse                  Ellipse;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;

    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.

    // creation
    template < class InputIterator >
    Min_ellipse_2( InputIterator first,
                   InputIterator last,
                   bool          randomize = false,
                   Random&       random    = default_random,
                   const Traits& traits    = Traits());
    
    Min_ellipse_2( const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p,
                   const Point&  q,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Point&  p4,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Point&  p4,
                   const Point&  p5,
                   const Traits& traits = Traits());
    ~Min_ellipse_2( );

    // access functions
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;

    Point_iterator  points_begin( ) const;
    Point_iterator  points_end  ( ) const;

    Support_point_iterator  support_points_begin( ) const;
    Support_point_iterator  support_points_end  ( ) const;

    const Point&  support_point( int i) const;

    const Ellipse&  ellipse( ) const;

    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;

    // modifiers
    void  insert( const Point& p);
    void  insert( const Point* first,
                  const Point* last );
    void  insert( std::list<Point>::const_iterator first,
                  std::list<Point>::const_iterator last );
    void  insert( std::istream_iterator<Point,std::ptrdiff_t> first,
                  std::istream_iterator<Point,std::ptrdiff_t> last );
    void  clear( );

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;

    // miscellaneous
    const Traits&  traits( ) const;
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
    std::list<Point>  points;                   // doubly linked list of points
@end

The support set $S$ of at most five support points is stored in an
array \ccc{support_points}, the actual number of support points is
given by \ccc{n_support_points}. During the computations, the set of
support points coincides with the set $B$ appearing in the pseudo-code
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
appearing in the pseudo-code for $\me(P,B)$, see Section~\ref{sec:algo}.

@! ----------------------------------------------------------------------------
\subsubsection{Constructors and Destructor}

We provide several different constructors, which can be put into two
groups. The constructors in the first group, i.e. the more important
ones, build the smallest enclosing ellipse $me(P)$ from a point set
$P$, given by a begin iterator and a past-the-end iterator, realized
as a member template.

All constructors of the first group copy the points into the internal
list \ccc{points}. If randomization is demanded, the points are copied
to a vector and shuffled at random, before being copied to \ccc{points}.
Finally the private member function $me$ is called to compute
$me(P)=me(P,\emptyset)$.

@macro <Min_ellipse_2 constructors> += @begin
    // STL-like constructor (member template)
    template < class InputIterator >
    Min_ellipse_2( InputIterator first,
                   InputIterator last,
                   bool          randomize
#if !defined(_MSC_VER) || _MSC_VER > 1300
                                           = false
#endif
                                                  ,
                   Random&       random    = default_random,
                   const Traits& traits    = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];

        // range of points not empty?
        if ( first != last) {

            // store points
            if ( randomize) {

                // shuffle points at random
                std::vector<Point> v( first, last);
                std::random_shuffle( v.begin(), v.end(), random);
                std::copy( v.begin(), v.end(),
                           std::back_inserter( points)); }
            else
                std::copy( first, last, std::back_inserter( points)); }

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
    Min_ellipse_2( const Traits& traits = Traits())
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
    Min_ellipse_2( const Point& p, const Traits& traits = Traits())
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
    Min_ellipse_2( const Point& p1, const Point& p2,
                   const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];

        // store points
        points.push_back( p1);
        points.push_back( p2);

        // compute me
        me( points.end(), 0);

        CGAL_optimisation_postcondition( is_degenerate());
    }

    // constructor for three points
    inline
    Min_ellipse_2( const Point& p1, const Point& p2, const Point& p3,
                   const Traits& traits = Traits())
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
    Min_ellipse_2( const Point& p1, const Point& p2,
                   const Point& p3, const Point& p4,
                   const Traits& traits = Traits())
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
    Min_ellipse_2( const Point& p1, const Point& p2, const Point& p3,
                   const Point& p4, const Point& p5,
                   const Traits& traits = Traits())
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
    ~Min_ellipse_2( )
    {
        // free support points' array
        delete[] support_points;
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Access Functions}

These functions are used to retrieve information about the current
status of the \ccc{Min_ellipse_2<Traits>} object. They are all very
simple (and therefore \ccc{inline}) and mostly rely on corresponding
access functions of the data members of \ccc{Min_ellipse_2<Traits>}.

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
    const Point&
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
    const Ellipse&
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
        return( number_of_support_points() <  3);
    }
@end

The remaining predicates perform in-ellipse tests, based on the
corresponding predicates of class \ccc{Ellipse}.

@macro <Min_ellipse_2 predicates> = @begin
    // in-ellipse test predicates
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( tco.ellipse.bounded_side( p));
    }

    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_bounded_side( p));
    }

    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( tco.ellipse.has_on_boundary( p));
    }

    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_unbounded_side( p));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers}

There is another way to build up $me(P)$, other than by supplying
the point set $P$ at once. Namely, $me(P)$ can be built up
incrementally, adding one point after another. If you look at the
pseudo-code in the introduction, this comes quite naturally. The
modifying method \ccc{insert}, applied with point $p$ to a
\ccc{Min_ellipse_2<Traits>} object representing $me(P)$,
computes $me(P \cup \{p\})$, where work has to be done only if $p$
lies outside $me(P)$. In this case, $me(P \cup \{p\}) = me(P,\{p\})$
holds, so the private member function \ccc{me} is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudo-code in
Section~\ref{sec:algo}.

@macro <Min_ellipse_2 modifiers> += @begin
    void
    insert( const Point& p)
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

Inserting a range of points is done by a single member template. In case
a compiler does not support this yet, we provide specialized \ccc{insert}
functions for C~arrays (using pointers as iterators), for STL sequence
containers \ccc{vector<Point>} and \ccc{list<Point>} and for the STL
input stream iterator \ccc{istream_iterator<Point>}.  Actually, the
\ccc{insert} function for a C~array and a \ccc{vector<point>} are the
same, since the random access iterator of \ccc{vector<Point>} is
implemented as \ccc{Point*}.

The following \ccc{insert} functions perform a call \ccc{insert(p)} for
each point \ccc{p} in the range $[\mbox{\ccc{first}},\mbox{\ccc{last}})$.

@macro <Min_ellipse_2 modifiers> += @begin
    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
    {
        for ( ; first != last; ++first)
            insert( *first);
    }
@end

The member function \ccc{clear} deletes all points from a
\ccc{Min_ellipse_2<Traits>} object and resets it to the
empty ellipse.

@macro <Min_ellipse_2 modifiers> += @begin
    void
    clear( )
    {
        points.erase( points.begin(), points.end());
        n_support_points = 0;
        tco.ellipse.set();
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Validity Check}

A \ccc{Min_ellipse_2<Traits>} object can be checked for validity.  This
means, it is checked whether (a) the ellipse contains all points of its
defining set $P$, (b) the ellipse is the smallest ellipse spanned by its
support set, and (c) the support set is minimal, i.e.\ no support point is
redundant. (\emph{Note:} (b) and (c) are not yet implemented. Instead we
check if the support set lies on the boundary of the ellipse.) If
\ccc{verbose} is \ccc{true}, some messages concerning the performed checks
are written to standard error stream.  The second parameter \ccc{level} is
not used, we provide it only for consistency with interfaces of other
classes.

@macro <Min_ellipse_2 validity check> = @begin
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        using namespace std;
        
        CGAL::Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL::Min_ellipse_2<Traits>::" << endl;
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
            return( CGAL::_optimisation_is_valid_fail( verr,
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
            return( CGAL::_optimisation_is_valid_fail( verr,
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
    const Traits&
    traits( ) const
    {
        return( tco);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Min_ellipse_2 I/O operators declaration> = @begin
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os, const Min_ellipse_2<Traits_>& me);

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is,       Min_ellipse_2<Traits_>& me);
@end

@macro <Min_ellipse_2 I/O operators> = @begin
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Min_ellipse_2<Traits_>& min_ellipse)
    {
        using namespace std;

        typedef  Min_ellipse_2<Traits_>::Point  Point;
        typedef  ostream_iterator<Point>        Os_it;

        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << endl;
            os << "CGAL::Min_ellipse_2( |P| = "<<min_ellipse.number_of_points()
               << ", |S| = " << min_ellipse.number_of_support_points() << endl;
            os << "  P = {" << endl;
            os << "    ";
            copy( min_ellipse.points_begin(), min_ellipse.points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S = {" << endl;
            os << "    ";
            copy( min_ellipse.support_points_begin(),
                  min_ellipse.support_points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  ellipse = " << min_ellipse.ellipse() << endl;
            os << ")" << endl;
            break;

          case CGAL::IO::ASCII:
            copy( min_ellipse.points_begin(), min_ellipse.points_end(),
                  Os_it( os, "\n"));
            break;

          case CGAL::IO::BINARY:
            copy( min_ellipse.points_begin(), min_ellipse.points_end(),
                  Os_it( os));
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is, CGAL::Min_ellipse_2<Traits_>& min_ellipse)
    {
        using namespace std;
        
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            typedef  Min_ellipse_2<Traits_>::Point  Point;
            typedef  istream_iterator<Point>       Is_it;
            min_ellipse.clear();
            min_ellipse.insert( Is_it( is), Is_it());
            break;

          default:
            CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Graphical Output}

@macro<Min_ellipse_2 graphical output operator> = @begin
    #ifdef CGAL_MIN_ELLIPSE_2_H
    #ifndef CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2
    #define CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2

    template< class Traits_ >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream &ws,
                  const CGAL::Min_ellipse_2<Traits_>& min_ellipse)
    {
        typedef CGAL::Min_ellipse_2<Traits_>::Point_iterator  Point_iterator;

        Point_iterator  first( min_ellipse.points_begin());
        Point_iterator  last ( min_ellipse.points_end());
        for ( ; first != last; ++first)
            ws << *first;
        return( ws << min_ellipse.ellipse());
    }
    #endif // CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2
    #endif // CGAL_MIN_ELLIPSE_2_H
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
\ccc{support_points}. The function is directly modeled after the
pseudo-code above.

@macro <Min_ellipse_2 private member function `me'> = @begin
    void
    me( const Point_iterator& last, int n_sp)
    {
        // compute ellipse through support points
        n_support_points = n_sp;
        compute_ellipse();
        if ( n_sp == 5) return;

        // test first n points
        typename std::list<Point>::iterator  point_iter = points.begin();
        for ( ; last != point_iter; ) {
            const Point& p = *point_iter;

            // p not in current ellipse?
            if ( has_on_unbounded_side( p)) {

                // recursive call with p as additional support point
                support_points[ n_sp] = p;
                me( point_iter, n_sp+1);

                // move current point to front
                points.splice( points.begin(), points, point_iter++); }

            else
                ++point_iter; }
    }
@end

@! ----------------------------------------------------------------------------
@! Class template Optimisation_ellipse_2<K>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Optimisation\_ellipse\_2<K>}

First, we declare the class template \ccc{Optimisation_ellipse_2},

@macro<Optimisation_ellipse_2 declaration> = @begin
    template < class K_ >
    class Optimisation_ellipse_2;
@end

The class interface looks as follows.

@macro <Optimisation_ellipse_2 interface> = @begin
    template < class K_ >
    class Optimisation_ellipse_2 {
        /*
        friend  std::ostream&  operator << <> (
            std::ostream&, const Optimisation_ellipse_2<K_>&);
        friend  std::istream&  operator >> <> (
            std::istream&, Optimisation_ellipse_2<K_> &);
        friend  CGAL::Window_stream& operator << <> (
            CGAL::Window_stream&, const Optimisation_ellipse_2<K_>&);
        */
      public:
        @<Optimisation_ellipse_2 public interface>

      /* private: */
        // private data members
        @<Optimisation_ellipse_2 private data members>

    @<dividing line>

    // Class implementation
    // ====================

      public:
        // Constructor
        // -----------
        @<Optimisation_ellipse_2 constructor>
        
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
    typedef           K_                K;
    typedef  typename K_::RT            RT;
    typedef  typename K_::FT            FT;
    typedef           CGAL::Point_2<K>  Point;
    typedef           CGAL::Conic_2<K>  Conic;

    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.

    // creation
    Optimisation_ellipse_2( );

    void  set( );
    void  set( const Point& p);
    void  set( const Point& p,  const Point& q);
    void  set( const Point& p1, const Point& p2, const Point& p3);
    void  set( const Point& p1, const Point& p2,
               const Point& p3, const Point& p4);
    void  set( const Point& p1, const Point& p2,
               const Point& p3, const Point& p4, const Point& p5);

    // access functions    
    int  number_of_boundary_points()

    // equality tests
    bool  operator == ( const Optimisation_ellipse_2<K>& e) const;
    bool  operator != ( const Optimisation_ellipse_2<K>& e) const;

    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;

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
points, if any, are stored directly in \ccc{boundary_point1} and
\ccc{boundary_point2}, resp.

@macro <Optimisation_ellipse_2 private data members> += @begin
    Point  boundary_point1, boundary_point2;    // two boundary points
@end

Given three or five points, the ellipse is represented as a conic,
using the class \ccc{Conic_2<K>}. The case with four boundary
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
\subsubsection{Constructor}

Only a default constructor is needed.

@macro <Optimisation_ellipse_2 constructor> = @begin
    inline
    Optimisation_ellipse_2( )
    { }
@end
    
@! ----------------------------------------------------------------------------
\subsubsection{Set Functions}

We provide set functions taking zero, one, two, three, four or five
boundary points. They all set the variable to the smallest ellipse
through the given points. \emph{Note:} The set function taking five
boundary points only uses the fifth point from its input together with
the two internally represented conics to compute the ellipse. The
algorithm in Section~\ref{sec:algo} guarantees that this set function
is only called if the current ellipse already has the first four
points as its boundary points.

@macro <Optimisation_ellipse_2 set functions> = @begin
    inline
    void
    set( )
    {
        n_boundary_points = 0;
    }
    
    inline
    void
    set( const Point& p)
    {
        n_boundary_points = 1;
        boundary_point1   = p;
    }
    
    inline
    void
    set( const Point& p, const Point& q)
    {
        n_boundary_points = 2;
        boundary_point1   = p;
        boundary_point2   = q;
    }
    
    inline
    void
    set( const Point& p1, const Point& p2, const Point& p3)
    {
        n_boundary_points = 3;
        conic1.set_ellipse( p1, p2, p3);
    }
    
    inline
    void
    set( const Point& p1, const Point& p2, const Point& p3, const Point& p4)
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
    set( const Point&, const Point&,
         const Point&, const Point&, const Point& p5)
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

    Conic_2< Cartesian< double > >
    to_double( ) const
    {
        CGAL_optimisation_precondition( ! is_degenerate());

        double t = 0.0;
        
        if ( n_boundary_points == 4)
            t = conic1.vol_minimum( dr, ds, dt, du, dv, dw);
        
        Conic_2<K> c( conic1);
        Conic_2< Cartesian<double> > e;
        e.set( CGAL::to_double( c.r()) + t*CGAL::to_double( dr),
               CGAL::to_double( c.s()) + t*CGAL::to_double( ds),
               CGAL::to_double( c.t()) + t*CGAL::to_double( dt),
               CGAL::to_double( c.u()) + t*CGAL::to_double( du),
               CGAL::to_double( c.v()) + t*CGAL::to_double( dv),
               CGAL::to_double( c.w()) + t*CGAL::to_double( dw));

        return( e);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Equality Tests}

@macro <Optimisation_ellipse_2 equality tests> = @begin
    bool
    operator == ( const Optimisation_ellipse_2<K>& e) const
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
    operator != ( const Optimisation_ellipse_2<K>& e) const
    {
        return( ! operator == ( e));
    }
@end
    
@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The following predicates perform in-ellipse tests and check for
emptiness and degeneracy, resp. The way to evaluate the in-ellipse
test depends on the number of boundary points and is realised by a
case analysis. Again, the case with four points is the most difficult
one.

@macro <Optimisation_ellipse_2 predicates> = @begin
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        switch ( n_boundary_points) {
          case 0:
            return( CGAL::ON_UNBOUNDED_SIDE);
          case 1:
            return( ( p == boundary_point1) ?
                           CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
          case 2:
            return(    ( p == boundary_point1)
                    || ( p == boundary_point2)
                    || ( CGAL::are_ordered_along_line(
                             boundary_point1, p, boundary_point2)) ?
                         CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
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
                int tau_star = c.vol_derivative( dr, ds, dt, du, dv, dw);
                return( CGAL::Bounded_side( CGAL_NTS sign( tau_star))); } }
          default:
            CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                         ( n_boundary_points <= 5) ); }
        // keeps g++ happy
        return( CGAL::Bounded_side( 0));
    }

    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDED_SIDE);
    }

    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDARY);
    }

    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_UNBOUNDED_SIDE);
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
    template < class K_ >
    std::ostream&
    operator << ( std::ostream&, const CGAL::Optimisation_ellipse_2<K_>&);

    template < class K_ >
    std::istream&
    operator >> ( std::istream&, CGAL::Optimisation_ellipse_2<K_>&);
@end

@macro <Optimisation_ellipse_2 I/O operators> = @begin
    template < class K_ >
    std::ostream&
    operator << ( std::ostream& os, const CGAL::Optimisation_ellipse_2<K_>& e)
    {
        const char* const  empty       = "";
        const char* const  pretty_head = "CGAL::Optimisation_ellipse_2( ";
        const char* const  pretty_sep  = ", ";
        const char* const  pretty_tail = ")";
        const char* const  ascii_sep   = " ";

        const char*  head = empty;
        const char*  sep  = empty;
        const char*  tail = empty;

        switch ( CGAL::get_mode( os)) {
          case CGAL::IO::PRETTY:
            head = pretty_head;
            sep  = pretty_sep;
            tail = pretty_tail;
            break;
          case CGAL::IO::ASCII:
            sep  = ascii_sep;
            break;
          case CGAL::IO::BINARY:
            break;
          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
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

    template < class K_ >
    std::istream&
    operator >> ( std::istream& is, CGAL::Optimisation_ellipse_2<K_>& e)
    {
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            CGAL::read( is, e.n_boundary_points);
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
                                             "CGAL::get_mode( is) invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Graphical Output}

@macro<Optimisation_ellipse_2 graphical output operator> = @begin
    #ifdef CGAL_OPTIMISATION_ELLIPSE_2_H
    #ifndef CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2
    #define CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2

    template< class Traits_ >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream &ws,
                  const CGAL::Optimisation_ellipse_2<Traits_>& oe)
    {
        switch ( oe.n_boundary_points) {
          case 0:
            break;
          case 1:
            ws << oe.boundary_point1;
            break;
          case 2: {
            double  px1( CGAL::to_double( oe.boundary_point1.x()));
            double  py1( CGAL::to_double( oe.boundary_point1.y()));
            double  px2( CGAL::to_double( oe.boundary_point2.x()));
            double  py2( CGAL::to_double( oe.boundary_point2.y()));
            ws.draw_segment( px1, py1, px2, py2); }
            break;
          case 3:
          case 4:
          case 5:
            ws << oe.to_double();
            break;
          default:
            CGAL_optimisation_assertion( ( oe.n_boundary_points >= 0) &&
                                         ( oe.n_boundary_points <= 5) ); }
        return( ws);
    }

    #endif // CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2
    #endif // CGAL_OPTIMISATION_ELLIPSE_2_H
@end


@! ----------------------------------------------------------------------------
@! Class template Min_ellipse_2_traits_2<K>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_ellipse\_2\_traits\_2<K>}

First, we declare the class templates \ccc{Min_ellipse_2} and
\ccc{Min_ellipse_2_traits_2}.

@macro<Min_ellipse_2_traits_2 declarations> = @begin
    template < class Traits_ >
    class Min_ellipse_2;

    template < class K_ >
    class Min_ellipse_2_traits_2;
@end

Since the actual work of the traits class is done in the nested type
\ccc{Ellipse}, we implement the whole class template in its interface.

The variable \ccc{ellipse} containing the current ellipse is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current ellipse,
it is declared \ccc{friend}.

@macro <Min_ellipse_2_traits_2 interface and implementation> = @begin
    template < class K_ >
    class Min_ellipse_2_traits_2 {
      public:
        // types
        typedef  K_                               K;
        typedef  CGAL::Point_2<K>                 Point;
        typedef  CGAL::Optimisation_ellipse_2<K>  Ellipse;

    private:
        // data members
        Ellipse  ellipse;                                // current ellipse

        // friends
        friend  class CGAL::Min_ellipse_2< CGAL::Min_ellipse_2_traits_2<K> >;

      public:
        // creation (use default implementations)
        // Min_ellipse_2_traits_2( );
        // Min_ellipse_2_traits_2( Min_ellipse_2_traits_2<K> const&);
    };
@end   

@! ----------------------------------------------------------------------------
@! Class template Min_ellipse_2_adapterC2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_ellipse\_2\_adapterC2<PT,DA>}

First, we declare the class templates \ccc{Min_ellipse_2},
\ccc{Min_ellipse_2_adapterC2} and
\ccc{_Min_ellipse_2_adapterC2__Ellipse}.

@macro<Min_ellipse_2_adapterC2 declarations> = @begin
    template < class Traits_ >
    class Min_ellipse_2;

    template < class PT_, class DA_ >
    class Min_ellipse_2_adapterC2;

    template < class PT_, class DA_ >
    class _Min_ellipse_2_adapterC2__Ellipse;
@end

The actual work of the adapter is done in the nested class
\ccc{Ellipse}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{ellipse} containing the current ellipse is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current ellipse,
it is declared \ccc{friend}.

@macro <Min_ellipse_2_adapterC2 interface and implementation> = @begin
    template < class PT_, class DA_ >
    class Min_ellipse_2_adapterC2 {
      public:
        // types
        typedef  PT_  PT;
        typedef  DA_  DA;

        // nested types
        typedef  PT                                             Point;
        typedef  _Min_ellipse_2_adapterC2__Ellipse<PT,DA>  Ellipse;

      private:
        DA      dao;                                    // data accessor object
        Ellipse ellipse;                                // current ellipse
        friend class Min_ellipse_2< Min_ellipse_2_adapterC2<PT,DA> >;

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
    Min_ellipse_2_adapterC2( const DA& da = DA())
        : dao( da), ellipse( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_ellipse_2_adapterC2 operations> = @begin
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::FT  FT;

        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;
        
        dao.get( p, px, py);
        dao.get( q, qx, qy);
        dao.get( r, rx, ry);

        return( static_cast< CGAL::Orientation>(
                   CGAL_NTS sign( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Ellipse}

@macro <Min_ellipse_2_adapterC2 nested type `Ellipse'> = @begin
    template < class PT_, class DA_ >
    std::ostream&  operator << ( std::ostream& os,
        const _Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& c);

    template < class PT_, class DA_ >
    std::istream&  operator >> ( std::istream& is,
        _Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& c);

    template < class PT_, class DA_ >
    class _Min_ellipse_2_adapterC2__Ellipse {
      public:
        // typedefs
        typedef  PT_  PT;
        typedef  DA_  DA;

        typedef           ConicCPA2< PT, DA>  CT;
        typedef  typename DA_::FT                  FT;

      private:
        // data members
        int  n_boundary_points;                 // number of boundary points
        PT   boundary_point1, boundary_point2;  // two boundary points
        CT   conic1, conic2;                    // two conics
        FT   dr, ds, dt, du, dv, dw;            // the gradient vector

        friend
        std::ostream&  operator << <> ( std::ostream& os,
            const _Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& c);

        friend
        std::istream&  operator >> <> ( std::istream& is,
            _Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& c);

      public:
        // types
        typedef  PT  Point;

        // creation
        _Min_ellipse_2_adapterC2__Ellipse( const DA& da)
            : conic1( da), conic2( da)
        { }

        void
        set( )
        {
            n_boundary_points = 0;
        }
    
        void
        set( const Point& p)
        {
            n_boundary_points = 1;
            boundary_point1   = p;
        }
    
        void
        set( const Point& p, const Point& q)
        {
            n_boundary_points = 2;
            boundary_point1   = p;
            boundary_point2   = q;
        }
    
        void
        set( const Point& p1, const Point& p2, const Point& p3)
        {
            n_boundary_points = 3;
            conic1.set_ellipse( p1, p2, p3);
        }
    
        void
        set( const Point& p1, const Point& p2,
             const Point& p3, const Point& p4)
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
        set( const Point&, const Point&,
             const Point&, const Point&, const Point& p5)
        {
            n_boundary_points = 5;
            conic1.set( conic1, conic2, p5);
            conic1.analyse();
        }

        // predicates
        CGAL::Bounded_side
        bounded_side( const Point& p) const
        {
            switch ( n_boundary_points) {
              case 0:
                return( CGAL::ON_UNBOUNDED_SIDE);
              case 1:
                return( ( p == boundary_point1) ?
                               CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
              case 2:
                return(    ( p == boundary_point1)
                        || ( p == boundary_point2)
                        || ( CGAL::are_ordered_along_lineC2( boundary_point1, p,
                                               boundary_point2, conic1.da())) ?
                                    CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
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
                    int tau_star = c.vol_derivative( dr, ds, dt, du, dv, dw);
                    return( CGAL::Bounded_side( CGAL_NTS sign( tau_star))); } }
              default:
                CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                             ( n_boundary_points <= 5) ); }
            // keeps g++ happy
            return( CGAL::Bounded_side( 0));
        }

        bool
        has_on_bounded_side( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_BOUNDED_SIDE);
        }

        bool
        has_on_boundary( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_BOUNDARY);
        }

        bool
        has_on_unbounded_side( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_UNBOUNDED_SIDE);
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
            const CGAL::_Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& e) const
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

        bool
        operator != (
            const CGAL::_Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& e) const
        {
            return( ! ( *this == e));
        }
    };

    // I/O
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const CGAL::_Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& e)
    {
        const char* const  empty       = "";
        const char* const  pretty_head =
                                 "CGAL::Min_ellipse_2_adapterC2::Ellipse( ";
        const char* const  pretty_sep  = ", ";
        const char* const  pretty_tail = ")";
        const char* const  ascii_sep   = " ";

        const char*  head = empty;
        const char*  sep  = empty;
        const char*  tail = empty;

        switch ( CGAL::get_mode( os)) {
          case CGAL::IO::PRETTY:
            head = pretty_head;
            sep  = pretty_sep;
            tail = pretty_tail;
            break;
          case CGAL::IO::ASCII:
            sep  = ascii_sep;
            break;
          case CGAL::IO::BINARY:
            break;
          default:
            CGAL_optimisation_assertion_msg( false,
                                            "CGAL::get_mode( os) invalid!");
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

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream& is,
                  CGAL::_Min_ellipse_2_adapterC2__Ellipse<PT_,DA_>& e)
    {
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            CGAL::read( is, e.n_boundary_points);
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
                                             "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
@! Class template Min_ellipse_2_adapterH2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_ellipse\_2\_adapterH2<PT,DA>}

First, we declare the class templates \ccc{Min_ellipse_2},
\ccc{Min_ellipse_2_adapterH2} and
\ccc{_Min_ellipse_2_adapterH2__Ellipse}.

@macro<Min_ellipse_2_adapterH2 declarations> = @begin
    template < class Traits_ >
    class Min_ellipse_2;

    template < class PT_, class DA_ >
    class Min_ellipse_2_adapterH2;

    template < class PT_, class DA_ >
    class _Min_ellipse_2_adapterH2__Ellipse;
@end

The actual work of the adapter is done in the nested class
\ccc{Ellipse}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{ellipse} containing the current ellipse is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current ellipse,
it is declared \ccc{friend}.

@macro <Min_ellipse_2_adapterH2 interface and implementation> = @begin
    template < class PT_, class DA_ >
    class Min_ellipse_2_adapterH2 {
      public:
        // types
        typedef  PT_  PT;
        typedef  DA_  DA;

        // nested types
        typedef  PT                                             Point;
        typedef  _Min_ellipse_2_adapterH2__Ellipse<PT,DA>  Ellipse;

      private:
        DA      dao;                                    // data accessor object
        Ellipse ellipse;                                // current ellipse
        friend class Min_ellipse_2< Min_ellipse_2_adapterH2<PT,DA> >;

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
    Min_ellipse_2_adapterH2( const DA& da = DA())
        : dao( da), ellipse( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_ellipse_2_adapterH2 operations> = @begin
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::RT  RT;

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

        return( static_cast< CGAL::Orientation>(
                 CGAL_NTS sign( ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                              - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Ellipse}

@macro <Min_ellipse_2_adapterH2 nested type `Ellipse'> = @begin
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream&,
                  const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream&,
                  CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

    template < class PT_, class DA_ >
    class _Min_ellipse_2_adapterH2__Ellipse {
      public:
        // typedefs
        typedef  PT_  PT;
        typedef  DA_  DA;

        typedef           CGAL::ConicHPA2< PT, DA>  CT;
        typedef  typename DA_::RT                   RT;

      private:
        // data members
        int  n_boundary_points;                 // number of boundary points
        PT   boundary_point1, boundary_point2;  // two boundary points
        CT   conic1, conic2;                    // two conics
        RT   dr, ds, dt, du, dv, dw;            // the gradient vector

        friend  std::ostream&  operator << <> ( std::ostream&,
            const _Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

        friend  std::istream&  operator >> <> ( std::istream&,
            _Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

      public:
        // types
        typedef  PT  Point;

        // creation
        _Min_ellipse_2_adapterH2__Ellipse( const DA& da)
            : conic1( da), conic2( da)
        { }

        void
        set( )
        {
            n_boundary_points = 0;
        }
    
        void
        set( const Point& p)
        {
            n_boundary_points = 1;
            boundary_point1   = p;
        }
    
        void
        set( const Point& p, const Point& q)
        {
            n_boundary_points = 2;
            boundary_point1   = p;
            boundary_point2   = q;
        }
    
        void
        set( const Point& p1, const Point& p2, const Point& p3)
        {
            n_boundary_points = 3;
            conic1.set_ellipse( p1, p2, p3);
        }
    
        void
        set( const Point& p1, const Point& p2,
             const Point& p3, const Point& p4)
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
        set( const Point&, const Point&,
             const Point&, const Point&, const Point& p5)
        {
            n_boundary_points = 5;
            conic1.set( conic1, conic2, p5);
            conic1.analyse();
        }

        // predicates
        CGAL::Bounded_side
        bounded_side( const Point& p) const
        {
            switch ( n_boundary_points) {
              case 0:
                return( CGAL::ON_UNBOUNDED_SIDE);
              case 1:
                return( ( p == boundary_point1) ?
                               CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
              case 2:
                return(    ( p == boundary_point1)
                        || ( p == boundary_point2)
                        || ( CGAL::are_ordered_along_lineH2( boundary_point1, p,
                                               boundary_point2, conic1.da())) ?
                                    CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
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
                    int tau_star = c.vol_derivative( dr, ds, dt, du, dv, dw);
                    return( CGAL::Bounded_side( CGAL_NTS sign( tau_star))); } }
              default:
                CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                             ( n_boundary_points <= 5) ); }
            // keeps g++ happy
            return( CGAL::Bounded_side( 0));
        }

        bool
        has_on_bounded_side( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_BOUNDED_SIDE);
        }

        bool
        has_on_boundary( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_BOUNDARY);
        }

        bool
        has_on_unbounded_side( const Point& p) const
        {
            return( bounded_side( p) == CGAL::ON_UNBOUNDED_SIDE);
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
            const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e) const
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

        bool
        operator != (
            const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e) const
        {
            return( ! ( *this == e));
        }
    };

    // I/O
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e)
    {
        const char* const  empty       = "";
        const char* const  pretty_head =
                                 "CGAL::Min_ellipse_2_adapterH2::Ellipse( ";
        const char* const  pretty_sep  = ", ";
        const char* const  pretty_tail = ")";
        const char* const  ascii_sep   = " ";

        const char*  head = empty;
        const char*  sep  = empty;
        const char*  tail = empty;

        switch ( CGAL::get_mode( os)) {
          case CGAL::IO::PRETTY:
            head = pretty_head;
            sep  = pretty_sep;
            tail = pretty_tail;
            break;
          case CGAL::IO::ASCII:
            sep  = ascii_sep;
            break;
          case CGAL::IO::BINARY:
            break;
          default:
            CGAL_optimisation_assertion_msg( false,
                                            "CGAL::get_mode( os) invalid!");
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

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream& is,
                  CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e)
    {
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            CGAL::read( is, e.n_boundary_points);
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
                                             "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ============================================================================
@! Tests
@! ============================================================================

\clearpage
\section{Test}

We test \ccc{Min_ellipse_2} with the traits class implementation
for optimisation algorithms, using exact arithmetic, i.e.\ Cartesian
representation with number type \ccc{Quotient<Gmpz>} or
\ccc{Quotient<integer>} and homogeneous representation with
number type \ccc{Gmpz} or \ccc{integer}.

@macro <Min_ellipse_2 test (includes and typedefs)> = @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Homogeneous.h>
    #include <CGAL/Min_ellipse_2.h>
    #include <CGAL/Min_ellipse_2_traits_2.h>
    #include <CGAL/Min_ellipse_2_adapterC2.h>
    #include <CGAL/Min_ellipse_2_adapterH2.h>
    #include <CGAL/IO/Verbose_ostream.h>
    #include <cassert>
    #include <cstring>
    #include <fstream>

    #ifdef CGAL_USE_LEDA_FOR_OPTIMISATION_TEST
    #  include <CGAL/leda_integer.h>
       typedef  leda_integer                       Rt;
       typedef  CGAL::Quotient< leda_integer >     Ft;
    #else
    #  include <CGAL/Gmpz.h>
       typedef  CGAL::Gmpz                         Rt;
       typedef  CGAL::Quotient< CGAL::Gmpz >       Ft;
    #endif

    typedef  CGAL::Cartesian< Ft >                 KerC;
    typedef  CGAL::Homogeneous< Rt >               KerH;
    typedef  CGAL::Min_ellipse_2_traits_2< KerC >  TraitsC;
    typedef  CGAL::Min_ellipse_2_traits_2< KerH >  TraitsH;
@end

The command line option \ccc{-verbose} enables verbose output.

@macro <Min_ellipse_2 test (verbose option)> = @begin
    bool  verbose = false;
    if ( ( argc > 1) && ( std::strcmp( argv[ 1], "-verbose") == 0)) {
        verbose = true;
        --argc;
        ++argv; }
@end

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

We call each function of class \ccc{Min_ellipse_2<Traits>} at least
once to ensure code coverage.

@macro <Min_ellipse_2 test (code coverage)> = @begin
    cover_Min_ellipse_2( verbose, TraitsC(), Rt());
    cover_Min_ellipse_2( verbose, TraitsH(), Rt());
@end

@macro <Min_ellipse_2 test (code coverage test function)> = @begin
    template < class Traits, class RT >
    void
    cover_Min_ellipse_2( bool verbose, const Traits&, const RT&)
    {
        using namespace std;

        typedef  CGAL::Min_ellipse_2< Traits >  Min_ellipse;
        typedef  typename Min_ellipse::Point    Point;
        typedef  typename Min_ellipse::Ellipse  Ellipse;

        CGAL::Verbose_ostream verr( verbose);

        // generate `n' points at random
        const int     n = 20;
        CGAL::Random  random_x, random_y;
        Point         random_points[ n];
        int           i;
        verr << n << " random points from [0,128)^2:" << endl;
        for ( i = 0; i < n; ++i) {
            random_points[ i] = Point( RT( random_x( 128)),
                                       RT( random_y( 128)));
            verr << i << ": " << random_points[ i] << endl;
        }

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
            bool  is_valid      = me.is_valid( verbose);
            bool  is_degenerate = me.is_degenerate();
            assert( is_valid);
            assert( is_degenerate);
        }

        verr << endl << "three points constructor...";
        {    
            Min_ellipse  me( random_points[ 3],
                             random_points[ 4],
                             random_points[ 5]);
            bool  is_valid = me.is_valid( verbose);
            int   num_pts  = me.number_of_points();
            assert( is_valid);
            assert( num_pts == 3);
        }

        verr << endl << "four points constructor...";
        {    
            Min_ellipse  me( random_points[ 6],
                             random_points[ 7],
                             random_points[ 8],
                             random_points[ 9]);
            bool  is_valid = me.is_valid( verbose);
            int   num_pts  = me.number_of_points();
            assert( is_valid);
            assert( num_pts == 4);
        }

        verr << endl << "five points constructor...";
        {    
            Min_ellipse  me( random_points[ 10],
                             random_points[ 11],
                             random_points[ 12],
                             random_points[ 13],
                             random_points[ 14]);
            bool  is_valid = me.is_valid( verbose);
            int   num_pts  = me.number_of_points();
            assert( is_valid);
            assert( num_pts == 5);
        }

        verr << endl << "Point* constructor...";
        Min_ellipse  me( random_points, random_points+9, false);
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
            Min_ellipse  me1( me.points_begin(), me.points_end(), false);
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
            typedef  typename Min_ellipse::Support_point_iterator
                                              Support_point_iterator;
            Point                   support_point;
            Support_point_iterator  iter( me.support_points_begin());
            for ( i = 0; i < me.number_of_support_points(); ++i, ++iter) {
                support_point = me.support_point( i);
                assert( support_point == *iter); }
            Support_point_iterator  end_iter( me.support_points_end());
            assert( iter == end_iter);
        }

        verr << endl << "ellipse access already called above...";

        verr << endl << "in-ellipse predicates...";
        {
            Point               p;
            CGAL::Bounded_side  bounded_side;
            bool                has_on_bounded_side;
            bool                has_on_boundary;
            bool                has_on_unbounded_side;
            for ( i = 0; i < 9; ++i) {
                p = random_points[ i];
                bounded_side          = me.bounded_side( p);
                has_on_bounded_side   = me.has_on_bounded_side( p);
                has_on_boundary       = me.has_on_boundary( p);
                has_on_unbounded_side = me.has_on_unbounded_side( p);
            assert( bounded_side != CGAL::ON_UNBOUNDED_SIDE);
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
            CGAL::set_ascii_mode( os);
            os << me;
        }
        {
            verr << endl << "  writing `test_Min_ellipse_2.pretty'...";
            ofstream os( "test_Min_ellipse_2.pretty");
            CGAL::set_pretty_mode( os);
            os << me;
        }
        {
            verr << endl << "  writing `test_Min_ellipse_2.binary'...";
            ofstream os( "test_Min_ellipse_2.binary");
            CGAL::set_binary_mode( os);
            os << me;
        }
        {
            verr << endl << "  reading `test_Min_ellipse_2.ascii'...";
            Min_ellipse me_in;
            ifstream is( "test_Min_ellipse_2.ascii");
            CGAL::set_ascii_mode( is);
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
    class MyPointC2;

    std::ostream&  operator << ( std::ostream&, const MyPointC2&);
    std::istream&  operator >> ( std::istream&,       MyPointC2&);
    
    class MyPointC2 {
      public:
        typedef  ::Ft  FT;
      private:
        FT x_;
        FT y_;
      public:
        MyPointC2( ) { }
        MyPointC2( const FT& x, const FT& y) : x_( x), y_( y) { }

        const FT&  x( ) const { return( x_); }
        const FT&  y( ) const { return( y_); }

        bool
        operator == ( const MyPointC2& p) const
        {
            return( ( x_ == p.x_) && ( y_ == p.y_));
        }

        bool
        operator != ( const MyPointC2& p) const
        {
            return( ( x_ != p.x_) || ( y_ != p.y_));
        }

        friend
        std::ostream&  operator << ( std::ostream& os, const MyPointC2& p);

        friend
        std::istream&  operator >> ( std::istream& is,       MyPointC2& p);
    };

    std::ostream&
    operator << ( std::ostream& os, const MyPointC2& p)
    {
        return( os << p.x_ << ' ' << p.y_);
    }

    std::istream&
    operator >> ( std::istream& is, MyPointC2& p)
    {
        return( is >> p.x_ >> p.y_);
    }

    // 2D Cartesian point class data accessor
    class MyPointC2DA {
      public:
        typedef  ::Ft  FT;

        MyPointC2DA( ) { }

        const FT&  get_x( const MyPointC2& p) const { return( p.x()); }
        const FT&  get_y( const MyPointC2& p) const { return( p.y()); }

        void
        get( const MyPointC2& p, FT& x, FT& y) const
        {
            x = get_x( p);
            y = get_y( p);
        }

        void
        set( MyPointC2& p, const FT& x, const FT& y) const
        {
            p = MyPointC2( x, y);
        }
    };


    // 2D homogeneous point class
    class MyPointH2;

    std::ostream&  operator << ( std::ostream&, const MyPointH2&);
    std::istream&  operator >> ( std::istream&,       MyPointH2&);
    
    class MyPointH2 {
      public:
        typedef  ::Rt  RT;
      private:
        RT hx_;
        RT hy_;
        RT hw_;
      public:
        MyPointH2( ) { }
        MyPointH2( const RT& hx, const RT& hy, const RT& hw = RT( 1))
            : hx_( hx), hy_( hy), hw_( hw) { }

        const RT&  hx( ) const { return( hx_); }
        const RT&  hy( ) const { return( hy_); }
        const RT&  hw( ) const { return( hw_); }

        bool
        operator == ( const MyPointH2& p) const
        {
            return( ( hx_*p.hw_ == p.hx_*hw_) && ( hy_*p.hw_ == p.hy_*hw_));
        }

        bool
        operator != ( const MyPointH2& p) const
        {
            return( ( hx_*p.hw_ != p.hx_*hw_) || ( hy_*p.hw_ != p.hy_*hw_));
        }

        friend
        std::ostream&  operator << ( std::ostream& os, const MyPointH2& p);

        friend
        std::istream&  operator >> ( std::istream& is,       MyPointH2& p);
    };

    std::ostream&
    operator << ( std::ostream& os, const MyPointH2& p)
    {
        return( os << p.hx_ << ' ' << p.hy_ << ' ' << p.hw_);
    }

    std::istream&
    operator >> ( std::istream& is, MyPointH2& p)
    {
        return( is >> p.hx_ >> p.hy_ >> p.hw_);
    }
    
    // 2D homogeneous point class data accessor
    class MyPointH2DA {
      public:
        typedef  ::Rt  RT;

        MyPointH2DA( ) { }

        const RT&  get_hx( const MyPointH2& p) const { return( p.hx()); }
        const RT&  get_hy( const MyPointH2& p) const { return( p.hy()); }
        const RT&  get_hw( const MyPointH2& p) const { return( p.hw()); }

        void
        get( const MyPointH2& p, RT& hx, RT& hy, RT& hw) const
        {
            hx = get_hx( p);
            hy = get_hy( p);
            hw = get_hw( p);
        }

        void
        set( MyPointH2& p, const RT& hx, const RT& hy, const RT& hw) const
        {
            p = MyPointH2( hx, hy, hw);
        }
    };
@end

To test the traits class adapters we use the code coverage test function.

@macro <Min_ellipse_2 test (adapters test)> = @begin
    typedef CGAL::Min_ellipse_2_adapterC2< MyPointC2, MyPointC2DA >  AdapterC2;
    typedef CGAL::Min_ellipse_2_adapterH2< MyPointH2, MyPointH2DA >  AdapterH2;
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

        typedef  CGAL::Min_ellipse_2< TraitsH >  Min_ellipse;
        typedef  Min_ellipse::Point              Point;
        typedef  Min_ellipse::Ellipse            Ellipse;

        CGAL::Verbose_ostream verr( verbose);

        // read points from file
        verr << std::endl << "input file: `" << argv[ 1] << "'" << std::flush;

        std::list<Point>  points;
        int               n, x, y;
        std::ifstream     in( argv[ 1]);
        in >> n;
        assert( in);
        for ( int i = 0; i < n; ++i) {
            in >> x >> y;
            assert( in);
            points.push_back( Point( x, y)); }

        // compute and check min_ellipse
        Min_ellipse  me2( points.begin(), points.end(), false);
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

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Min_ellipse_2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2.h}

@file <include/CGAL/Min_ellipse_2.h> = @begin
    @<file header>(
        "include/CGAL/Min_ellipse_2.h",
        "2D Smallest Enclosing Ellipse")

    #ifndef CGAL_MIN_ELLIPSE_2_H
    #define CGAL_MIN_ELLIPSE_2_H

    // includes
    #include <CGAL/Optimisation/basic.h>
    #include <CGAL/Random.h>
    #include <list>
    #include <vector>
    #include <algorithm>
    #include <iostream>

    @<namespace begin>("CGAL")
    
    // Class declaration
    // =================
    @<Min_ellipse_2 declaration>

    // Class interface
    // ===============
    @<Min_ellipse_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators declaration>

    @<namespace end>("CGAL")
    
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
    @<file header>(
        "include/CGAL/Min_ellipse_2.C",
        "2D Smallest Enclosing Ellipse")

    @<namespace begin>("CGAL")
    
    // Class implementation (continued)
    // ================================
    // I/O
    // ---
    @<Min_ellipse_2 I/O operators>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_ellipse_2.h
@! ----------------------------------------------------------------------------

\subsection{Optimisation\_ellipse\_2.h}

@file <include/CGAL/Optimisation_ellipse_2.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation_ellipse_2.h",
        "2D Optimisation Ellipse")

    #ifndef CGAL_OPTIMISATION_ELLIPSE_2_H
    #define CGAL_OPTIMISATION_ELLIPSE_2_H

    // the following include is needed by `to_double()'
    #ifndef CGAL_CARTESIAN_H
    #  include <CGAL/Cartesian.h>
    #endif

    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_CONIC_2_H
    #  include <CGAL/Conic_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif
    #ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
    #  include <CGAL/IO/forward_decl_window_stream.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declaration
    // =================
    @<Optimisation_ellipse_2 declaration>

    // Class interface
    // ===============
    @<Optimisation_ellipse_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Optimisation_ellipse_2 I/O operators declaration>

    @<namespace end>("CGAL")
    
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
    @<file header>(
        "include/CGAL/Optimisation_ellipse_2.C",
        "2D Optimisation Ellipse")

    @<namespace begin>("CGAL")
    
    // Class implementation (continued)
    // ================================

    // I/O
    // ---
    @<Optimisation_ellipse_2 I/O operators>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_traits_2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_traits\_2.h}

@file <include/CGAL/Min_ellipse_2_traits_2.h> = @begin
    @<file header>(
        "include/CGAL/Min_ellipse_2_traits_2.h",
        "default traits class for 2D Smallest Enclosing Ellipse")

    #ifndef CGAL_MIN_ELLIPSE_2_TRAITS_2_H
    #define CGAL_MIN_ELLIPSE_2_TRAITS_2_H

    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ELLIPSE_2_H
    #  include <CGAL/Optimisation_ellipse_2.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_ellipse_2_traits_2 declarations>

    // Class interface and implementation
    // ==================================
    @<Min_ellipse_2_traits_2 interface and implementation>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_ELLIPSE_2_TRAITS_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_adapterC2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_adapterC2.h}

@file <include/CGAL/Min_ellipse_2_adapterC2.h> = @begin
    @<file header>(
        "include/CGAL/Min_ellipse_2_adapterC2.h",
        "traits class adapter for 2D Smallest Enclosing Ellipse")

    #ifndef CGAL_MIN_ELLIPSE_2_ADAPTERC2_H
    #define CGAL_MIN_ELLIPSE_2_ADAPTERC2_H

    // includes
    #ifndef CGAL_CONICCPA2_H
    #  include <CGAL/ConicCPA2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_ellipse_2_adapterC2 declarations>

    // Class interface and implementation
    // ==================================
    template < class PT, class DA >
    bool
    are_ordered_along_lineC2( const PT& p, const PT& q, const PT& r,
                              const DA& da)
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
        if ( ! CGAL_NTS is_zero( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx)))
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

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_ELLIPSE_2_ADAPTERC2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_adapterH2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_adapterH2.h}

@file <include/CGAL/Min_ellipse_2_adapterH2.h> = @begin
    @<file header>(
        "include/CGAL/Min_ellipse_2_adapterH2.h",
        "traits class adapter for 2D Smallest Enclosing Ellipse")

    #ifndef CGAL_MIN_ELLIPSE_2_ADAPTERH2_H
    #define CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

    // includes
    #ifndef CGAL_CONICHPA2_H
    #  include <CGAL/ConicHPA2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_ellipse_2_adapterH2 declarations>

    // Class interface and implementation
    // ==================================
    template < class PT, class DA >
    bool
    are_ordered_along_lineH2( const PT& p, const PT& q, const PT& r,
                              const DA& da)
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
        if ( ! CGAL_NTS is_zero(  ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
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

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_ellipse_2_Window_stream.h
@! ----------------------------------------------------------------------------

\subsection{Min\_ellipse\_2\_Window\_stream.h}

@file <include/CGAL/IO/Min_ellipse_2_Window_stream.h> = @begin
    @<file header>(
        "include/CGAL/IO/Min_ellipse_2_Window_stream.h",
        "graphical output to `leda_window' for Min_ellipse_2 algo.")

    // Each of the following operators is individually 
    // protected against multiple inclusion.

    // Window_stream I/O operators
    // ===========================
    // includes
    #ifndef CGAL_CONIC_2_WINDOW_STREAM_H
    #  include <CGAL/IO/Conic_2_Window_stream.h>
    #endif

    // Optimisation_ellipse_2
    // ----------------------
    @<Optimisation_ellipse_2 graphical output operator>

    // Min_ellipse_2
    // -------------
    @<Min_ellipse_2 graphical output operator>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_ellipse_2.C
@! ----------------------------------------------------------------------------

\subsection{test\_Min\_ellipse\_2.C}

@file <test/Min_ellipse_2/test_Min_ellipse_2.C> = @begin
    @<file header>(
        "test/Min_ellipse_2/test_Min_ellipse_2.C",
        "test program for 2D Smallest Enclosing Ellipse")

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

        return( 0);
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! File Header
@! ----------------------------------------------------------------------------

\subsection*{File Header}

@i share/file_header.awi
 
And here comes the specific file header for the product files of this
web file.

@macro <file header>(2) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Min_ellipse_2",
        "Geometric Optimisation",
        "Min_ellipse_2",
        "$Id$","$Date$",
        "Sven Schönherr <sven@@inf.ethz.ch>, Bernd Gärtner",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ============================================================================
@! Bibliography
@! ============================================================================

\clearpage
\bibliographystyle{plain}
\bibliography{geom,../doc_tex/basic/Optimisation/cgal}

@! ===== EOF ==================================================================
