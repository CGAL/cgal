@! ============================================================================
@! The CGAL Library
@! Implementation: Smallest Enclosing Sphere in Arbitrary Dimension
@! ----------------------------------------------------------------------------
@! file  : web/Min_sphere_d.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Min_sphere_d_new WIP $
@! $Revision$
@! $Date$
@! ============================================================================
 
@documentclass[twoside,fleqn]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amsmath}
@usepackage{amssymb}
@usepackage{path}
@usepackage{cc_manual,cc_manual_index}
@article

\input{cprog.sty}
\setlength{\skip\footins}{3ex}

\pagestyle{headings}

@! LaTeX macros
\newcommand{\remark}[2]{[\textbf{#1:} \emph{#2}]}

\newcommand{\linebreakByHand}{\ccTexHtml{\linebreak[4]}{}}
\newcommand{  \newlineByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

\renewcommand{\sectionmark}[1]{\markboth{\uppercase{#1}}{}}

\newcommand{\subsectionRef}[2]{
  \addtocounter{subsection}{1}
  \addcontentsline{toc}{subsection}{\protect\numberline{\thesubsection}#1: #2}
  \markright{\thesubsection~~#1: #2}}

@! settings for `cc_manual.sty'
\ccDefGlobalScope{CGAL::}
\renewcommand{\ccRefPageEnd}{\clearpage}

\newcommand{\cgalColumnLayout}{%
  \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
  \ccPropagateThreeToTwoColumns}
\newcommand{\cgalMinSphereLayout}{%
  \ccSetThreeColumns{}{min_sphere.center()\,}{returns
    \ccGlobalScope\ccc{ON_BOUNDED_SIDE},
    \ccGlobalScope\ccc{ON_BOUNDARY},}
  \ccPropagateThreeToTwoColumns}

\newcommand{\ccRequirements}{\ccHeading{Requirements}}
\newcommand{\ccRequire}{\ccCommentHeading{Requirements}}

@! ============================================================================
@! Title
@! ============================================================================

\thispagestyle{empty}

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 20 mm
@t title titlefont centre "Smallest Enclosing Sphere"
@t vskip 0 mm
@t title titlefont centre "in Arbitrary Dimension*"
@t vskip 10 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
\begin{center}
  \textbf{ETH Z{\"u}rich}
\end{center}
@t vskip 10 mm
{\small
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Package: Min_sphere_d WIP+\cgalWIP\verb+$+ \\
    \verb+$CGAL_Chapter: Geometric Optimisation $+ \\
  \end{tabular}
\end{center}
}
@t vskip 30 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Project
  No.~28155 (GALIA), and by a grant from the Swiss Federal Office for
  Education and Sciences for this project.}
\renewcommand{\thefootnote}{\arabic{footnote}}

@! --------
@! Abstract
@! --------

\begin{abstract}
  We provide an implementation for computing the smallest (w.r.t.~volume)
  enclosing sphere of a finite point set in arbitrary dimension. The
  problem is formulated as a quadratic program and a dedicated
  solver~\cite{gs-eegqp-00} is used to obtain the solution.
\end{abstract}

@! --------
@! Contents
@! --------

\clearpage

\newlength{\defaultparskip}
\setlength{\defaultparskip}{\parskip}
\setlength{\parskip}{1ex}

\tableofcontents

\setlength{\parskip}{\defaultparskip}


@! ============================================================================
@! Introduction
@! ============================================================================

\clearpage
\markright{\uppercase{Introduction}}
\section{Introduction}

We consider the problem of finding the unique sphere of smallest volume
enclosing a finite set of points in $d$-dimensional Euclidean space $\E_d$.
It can be formulated as an optimization problem with linear constraints and
a convex quadratic objective function~\cite{gs-eegqp-00}.

@! ----------------------------------------------------------------------------
@! Smallest Enclosing Sphere as a Quadratic Programming Problem
@! ----------------------------------------------------------------------------

\subsection{Smallest Enclosing Sphere as a Quadratic Programming Problem}

If the point set is given as $P = \{p_1,\dots,p_n\}$, we want to find a
point~$p^*$ such that $max_{i=1}^n \|p_i-p^*\|$ is minimized. The
point~$p^*$ then is the center of the smallest enclosing sphere. Define the
$d\!\times\!n$-matrix $C := (p_1,\dots,p_n)$ and consider the quadratic
programming problem
%
\begin{equation} \label{eq:MS_as_QP}
  \begin{array}{lll}
    \text{(MS)} & \text{minimize}   & x^T C^T C\, x
                                        - \sum_{i=1}^n p_i^Tp_i\, x_i \\[0.8ex]
                & \text{subject to} & \sum_{i=1}^n x_i = 1, \\[0.5ex]
                &                   & x \geq 0.
  \end{array}
\end{equation}
%
Let $x^* = (x^*_1,\dots,x^*_n)$ be its optimal solution, then the point
\[
   p^* = \sum_{i=1}^n p_i x^*_i
\]
is the center of the smallest enclosing sphere of $P$. The squared radius
is the negative value of the objective function at $x^*$.


@! ============================================================================
@! Reference Pages
@! ============================================================================

\clearpage
\section{Reference Pages} \label{sec:reference_pages}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

@p maximum_input_line_length = 102

@! ----------------------------------------------------------------------------
@! Class: Min_sphere_d
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{CGAL::Min\_sphere\_d\texttt{<}Traits\texttt{>}}
%\input{../doc_tex/basic/Optimisation/Optimisation_ref/Min_sphere_d.tex}

@! ----------------------------------------------------------------------------
@! Concept: OptimisationDTraits
@! ----------------------------------------------------------------------------

\subsectionRef{Concept}{Optimisation\_d\_traits}
\input{../../Optimisation_basic/doc_tex/basic/Optimisation/Optimisation_ref/OptimisationDTraits.tex}

@p maximum_input_line_length = 80


@! ============================================================================
@! Implementation
@! ============================================================================

\clearpage
\section{Implementation} \label{sec:implementation}

@! ----------------------------------------------------------------------------
@! The Class Template CGAL::Min_sphere_d<Traits>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Min\_sphere\_d\texttt{<}Traits\texttt{>}}

The class template \ccc{Min_sphere_d} expects a model of the concept
\ccc{OptimisationDTraits} (see Section~\ref{ccRef_OptimisationDTraits}.2)
as its template argument.

@macro <Min_sphere_d declarations> += @begin
    template < class Traits_ >
    class Min_sphere_d;
@end

The interface consists of the public types and member functions described
in Section~\ref{ccRef_CGAL::Min_sphere_d<Traits>}.1 and of some private
types, private member functions, and data members.

@macro <Min_sphere_d interface> = @begin
    template < class Traits_ >
    class Min_sphere_d {
      public:
        // self
        typedef  Traits_                    Traits;
        typedef  Min_sphere_d<Traits>       Self;

        // types from the traits class
        typedef  typename Traits::Point_d   Point;

        typedef  typename Traits::Rep_tag   Rep_tag;
        
        typedef  typename Traits::RT        RT;
        typedef  typename Traits::FT        FT;

        typedef  typename Traits::Access_dimension_d
                                            Access_dimension_d;
        typedef  typename Traits::Access_coordinates_begin_d
                                            Access_coordinates_begin_d;

        typedef  typename Traits::Construct_point_d
                                            Construct_point_d;

        typedef  typename Traits::ET        ET;
        typedef  typename Traits::NT        NT;
        
      private:
        @<Min_sphere_d Solver type>
        @<Min_sphere_d QP_solver types>
        @<Min_sphere_d private types>

      public:
        @<Min_sphere_d types>
        
        @<Min_sphere_d member functions>
        
      private:
        @<Min_sphere_d data members>

        @<Min_sphere_d private member functions>
    };
@end

@! ----------------------------------------------------------------------------
\subsubsection{Data Members}

Mainly, we have to store the given input points, the center and the squared
radius of the smallest enclosing sphere, and an instance of the quadratic
programming solver. Additional variables, that are used in the member
functions described below, are introduced when they appear for the first
time.

We start with the traits class object.

@macro <Min_sphere_d data members> += @begin

    Traits                   tco;       // traits class object
@end

The inputs points are kept in a vector to have random access to them.
Their dimension is stored separately.

@macro <Min_sphere_d standard includes> += @begin
    #ifndef CGAL_PROTECT_VECTOR
    #  include <vector>
    #  define CGAL_PROTECT_VECTOR
    #endif
@end

@macro <Min_sphere_d private types> += @begin
    // private types
    typedef  std::vector<Point>         Point_vector;
@end

@macro <Min_sphere_d data members> += @begin

    Point_vector             points;    // input points
    int                      d;         // dimension of input points
@end

The center and the squared radius of the smallest enclosing sphere are
stored with rational representation, i.e.~numerators and denominators are
kept separately. The vector \ccc{center_coords} contains $d+1$ entries,
the numerators of the $d$ coordinates and the common denominator.

@macro <Min_sphere_d private types> += @begin
    typedef  std::vector<ET>            ET_vector;
@end

@macro <Min_sphere_d data members> += @begin
    
    ET_vector                center_coords;     // center of small.encl.sphere

    ET                       sqr_rad_numer;     // squared radius of
    ET                       sqr_rad_denom;     // smallest enclosing sphere
@end

We store an instance of the quadratic programming solver described
in~\cite{s-qpego1-00}. The details are given in
Section~\ref{sec:using_qp_solver} below, here it suffice to know that there
is a variable \ccc{solver} of type \ccc{Solver}.

@macro <Min_sphere_d private types: quadratic programming solver> zero = @begin
    typedef  ...                        Solver;
@end

@macro <Min_sphere_d data members> += @begin

    Solver                   solver;    // quadratic programming solver
@end

@! ----------------------------------------------------------------------------
\subsubsection{Creation}

Two constructors are provided. If the user wants to get some verbose output
(of the underlying QP solver), he can override the default arguments of
\ccc{verbose} and \ccc{stream}.

@macro <Min_sphere_d standard includes> += @begin
    #ifndef CGAL_PROTECT_IOSTREAM
    #  include <iostream>
    #  define CGAL_PROTECT_IOSTREAM
    #endif
@end

@macro <Min_sphere_d member functions> += @begin
    // creation
    Min_sphere_d( const Traits&  traits  = Traits(),
                  int            verbose = 0,
                  std::ostream&  stream  = std::cout)
      : tco( traits), d( -1), solver( verbose, stream)
        {
            @<Min_sphere_d QP-solver set-up>
        }
@end

The second constructor expects a set of points given via an iterator range.
It calls the \ccc{set} member function described in
Subsection~\ref{sec:modifiers} to store the points and to compute the
smallest enclosing sphere of the given point set.

@macro <Min_sphere_d member functions> += @begin

    template < class InputIterator >
    Min_sphere_d( InputIterator  first,
                  InputIterator  last,
                  const Traits&  traits = Traits(),
                  int            verbose = 0,
                  std::ostream&  stream  = std::cout)
      : tco( traits), solver( verbose, stream)
        {
            @<Min_sphere_d QP-solver set-up>
            set( first, last);
        }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Access}

The following types and member functions give access to the set of points
contained in the smallest enclosing sphere.

@macro <Min_sphere_d types> += @begin
    // public types
    typedef  typename Point_vector::const_iterator
                                        Point_iterator;
@end

@macro <Min_sphere_d member functions> += @begin

    // access to point set
    int  ambient_dimension( ) const { return d; }

    int  number_of_points( ) const { return points.size(); }

    Point_iterator  points_begin( ) const { return points.begin(); }
    Point_iterator  points_end  ( ) const { return points.end  (); }
@end

To access the support points, we exploit the following fact. A point~$p_i$
is a support point, iff its corresponding variable $x_i$ (of the QP solver)
is basic. Thus the number of support points is equal to the number of basic
variables, if the smallest enclosing sphere is not empty.

@macro <Min_sphere_d member functions> += @begin

    // access to support points
    int
    number_of_support_points( ) const
        { return is_empty() ? 0 : solver.number_of_basic_variables(); }
@end

If $i$ is the index of the $k$-th basic variable, then $p_i$ is the $k$-th
support point. To access a point given its index, we use the following
function class.

@macro <Min_sphere_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_ACCESS_BY_INDEX_H
    #  include <CGAL/_QP_solver/Access_by_index.h>
    #endif
@end

@macro <Min_sphere_d private types> += @begin

    typedef  CGAL::Access_by_index<typename std::vector<Point>::const_iterator>
                                        Point_by_index;
@end

The indices of the basic variables can be accessed with the following
iterator.

@macro <Min_sphere_d QP_solver types> += @begin
    // types from the QP solver
    typedef  typename Solver::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;
@end

Combining the function class with the index iterator gives the support
point iterator.

@macro <Min_sphere_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
    #  include <CGAL/_QP_solver/Join_random_access_iterator.h>
    #endif
@end

@macro <Min_sphere_d types> += @begin

    typedef  CGAL::Join_random_access_iterator_1<
                 Basic_variable_index_iterator, Point_by_index >
                                        Support_point_iterator;
@end

@macro <Min_sphere_d member functions> += @begin

    Support_point_iterator
    support_points_begin() const
        { return Support_point_iterator(
                     solver.basic_variables_index_begin(),
                     Point_by_index( points.begin())); }

    Support_point_iterator
    support_points_end() const
        { return Support_point_iterator(
                     is_empty() ? solver.basic_variables_index_begin()
                                : solver.basic_variables_index_end(),
                     Point_by_index( points.begin())); }
@end

The following types and member functions give access to the center and the
squared radius of the smallest enclosing sphere.

@macro <Min_sphere_d types> += @begin

    typedef  typename ET_vector::const_iterator
                                        Coordinate_iterator;
@end

@macro <Min_sphere_d member functions> += @begin

    // access to center (rational representation)
    Coordinate_iterator
    center_coordinates_begin( ) const { return center_coords.begin(); }

    Coordinate_iterator
    center_coordinates_end  ( ) const { return center_coords.end  (); }

    // access to squared radius (rational representation)
    ET  squared_radius_numerator  ( ) const { return sqr_rad_numer; }
    ET  squared_radius_denominator( ) const { return sqr_rad_denom; }
@end

For convinience, we also provide member functions for accessing the center
as a single point of type \ccc{Point} and the squared radius as a single
number of type \ccc{FT}. Both functions only work, if an implicit
conversion from number type \ccc{ET} to number type \ccc{RT} is available,
e.g.~if both types are the same.

@macro <Min_sphere_d member functions> += @begin

    // access to center and squared radius
    // NOTE: an implicit conversion from ET to RT must be available!
    Point
    center( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return tco.construct_point_d_object()( ambient_dimension(),
                                                 center_coordinates_begin(),
                                                 center_coordinates_end()); }

    FT
    squared_radius( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_radius_numerator  ()) /
                 FT( squared_radius_denominator()); }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

We use the private member function \ccc{sqr_dist} to compute the
squared distance of a given point to the center of the smallest
enclosing sphere.

@macro <Min_sphere_d CGAL includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_H
    #  include <CGAL/function_objects.h>
    #endif
    #ifndef CGAL_IDENTITY_H
    #  include <CGAL/_QP_solver/identity.h>
    #endif
@end

@macro <Min_sphere_d private member functions> += @begin

    // squared distance to center
    ET
    sqr_dist( const Point& p) const
        { return std::inner_product(
              center_coords.begin(), center_coords.end()-1,
              tco.access_coordinates_begin_d_object()( p), ET( 0),
              std::plus<ET>(),
              CGAL::compose1_2( 
                  CGAL::compose2_1( std::multiplies<ET>(),
                      CGAL::identity<ET>(), CGAL::identity<ET>()),
                  CGAL::compose2_2( std::minus<ET>(),
                      CGAL::identity<ET>(),
                      std::bind2nd( std::multiplies<ET>(),
                                    center_coords.back())))); }
@end

Now the implementation of the sidedness predicates is straight forward.

@macro <Min_sphere_d member functions> += @begin

    // predicates
    CGAL::Bounded_side
    bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          return CGAL::Bounded_side( CGAL::NTS::sign(
              sqr_rad_numer - sqr_dist( p))); }

    bool
    has_on_bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          return ( sqr_dist( p) < sqr_rad_numer); }

    bool
    has_on_boundary( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          return ( sqr_dist( p) == sqr_rad_numer); }

    bool
    has_on_unbounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          return( sqr_dist( p) > sqr_rad_numer); }
@end

The smallest enclosing sphere is \emph{empty}, if it contains no points,
and it is \emph{degenerate}, if it has less than two support points.

@macro <Min_sphere_d member functions> += @begin

    bool  is_empty     ( ) const { return number_of_points() == 0; }
    bool  is_degenerate( ) const { return number_of_support_points() < 2; }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers} \label{sec:modifiers}

These private member functions are used by the following \ccc{set} and
\ccc{insert} member functions to set and check the dimension of the input
points, respectively.

@macro <Min_sphere_d private member functions> += @begin

    // set dimension of input points
    void
    set_dimension( )
        { d = ( points.size() == 0 ? -1 :
                    tco.access_dimension_d_object()( points[ 0])); }

    // check dimension of input points
    bool
    check_dimension( unsigned int  offset = 0)
        { return ( std::find_if( points.begin()+offset, points.end(),
                                 CGAL::compose1_1( std::bind2nd(
                                     std::not_equal_to<int>(), d),
                                     tco.access_dimension_d_object()))
                   == points.end()); }
@end

The \ccc{set} member function copies the input points into the internal
variable \ccc{points} and calls the private member function
\ccc{compute_min_sphere} (described in Section~\ref{sec:using_qp_solver})
to compute the smallest enclosing sphere.

@macro <Min_sphere_d member functions> += @begin

    // modifiers
    template < class InputIterator >
    void
    set( InputIterator first, InputIterator last)
        { if ( points.size() > 0) points.erase( points.begin(), points.end());
          std::copy( first, last, std::back_inserter( points));
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension(),
              "Not all points have the same dimension.");
          compute_min_sphere(); }
@end

The \ccc{insert} member functions append the given point(s) to the point
set and recompute the smallest enclosing sphere.

@macro <Min_sphere_d member functions> += @begin

    void
    insert( const Point& p)
        { CGAL_optimisation_precondition( is_empty() ||
              ( tco.access_dimension_d_object()( p) == d));
          points.push_back( p);
          compute_min_sphere(); }

    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
        { CGAL_optimisation_precondition_code( int old_n = points.size());
          points.insert( points.end(), first, last);
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension( old_n),
              "Not all points have the same dimension.");
          compute_min_sphere(); }
@end

The \ccc{clear} member function deletes all points and resets the smallest
enclosing sphere to the empty sphere.

@macro <Min_sphere_d member functions> += @begin

    void
    clear( )
        { points.erase( points.begin(), points.end());
          compute_min_sphere(); }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Validity Check}

A \ccc{Min_sphere_d<Traits>} object can be checked for validity. This
means, it is checked whether (a) the sphere contains all points of its
defining set $P$, and (b) the sphere is the smallest sphere spanned by its
support set $S$ and the support set is minimal, i.e.~no support point is
redundant. The function \ccc{is_valid} is mainly intended for debugging
user supplied traits classes but also for convincing the anxious user that
the traits class implementation is correct. If \ccc{verbose} is \ccc{true},
some messages concerning the performed checks are written to standard error
stream. The second parameter \ccc{level} is not used, we provide it only
for consistency with interfaces of other classes.

@macro <Min_sphere_d member functions> += @begin

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
@end

@macro <Min_sphere_d validity check> = @begin
    // validity check
    template < class Traits_ >
    bool
    Min_sphere_d<Traits_>::
    is_valid( bool verbose, int level) const
    {
        CGAL_USING_NAMESPACE_STD
        
        CGAL::Verbose_ostream verr( verbose);
        verr << "CGAL::Min_sphere_d<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;
        
        // containment check (a)
        // ---------------------
        @<Min_sphere_d validity check: containment>

        // support set check (b)
        // ---------------------
        @<Min_sphere_d validity check: support set>

        verr << "  object is valid!" << endl;
        return( true);
    }
@end

The containment check (a) is easy to perform, just a loop over all
points in $|P|$.

@macro <Min_sphere_d validity check: containment> = @begin
    verr << "  (a) containment check..." << flush;

    Point_iterator  point_it = points_begin();
    for ( ; point_it != points_end(); ++point_it) {
        if ( has_on_unbounded_side( *point_it)) 
            return CGAL::_optimisation_is_valid_fail( verr,
                       "sphere does not contain all points");
    }

    verr << "passed." << endl;
@end

To validate the support set, we check whether all support points lie on the
boundary of the smallest enclosing sphere, and if the center is strictly
contained in the convex hull of the support set. Since the center is a
linear combination of the support points, it suffice to check if the
coefficients are positive and at most $1$.

@macro <Min_sphere_d validity check: support set> = @begin
    verr << "  (b) support set check..." << flush;

    // all support points on boundary?
    Support_point_iterator  support_point_it = support_points_begin();
    for ( ; support_point_it != support_points_end(); ++support_point_it) {
        if ( ! has_on_boundary( *support_point_it)) 
            return CGAL::_optimisation_is_valid_fail( verr,
                "sphere does not have all support points on its boundary");
    }

    // center strictly in convex hull of support points?
    typename Solver::Basic_variable_numerator_iterator
        num_it = solver.basic_variables_numerator_begin();
    for ( ; num_it != solver.basic_variables_numerator_end(); ++num_it) {
        if ( ! (    CGAL::NTS::is_positive( *num_it)
                 && *num_it <= solver.variables_common_denominator()))
            return CGAL::_optimisation_is_valid_fail( verr,
              "center does not lie strictly in convex hull of support points");
    }

    verr << "passed." << endl;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Miscellaneous}

The member function \ccc{traits} returns a const reference to the
traits class object.

@macro <Min_sphere_d member functions> += @begin

    // traits class access
    const Traits&  traits( ) const { return tco; }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Min_sphere_d I/O operators declaration> = @begin
    // I/O operators
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os, const Min_sphere_d<Traits_>& min_sphere);

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is,       Min_sphere_d<Traits_>& min_sphere);
@end

@macro <Min_sphere_d I/O operators> = @begin
    // output operator
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Min_sphere_d<Traits_>& min_sphere)
    {
        CGAL_USING_NAMESPACE_STD

        typedef  Min_sphere_d<Traits_>::Point  Point;
        typedef  ostream_iterator<Point>       Os_it;
        typedef  typename Traits_::ET          ET;
        typedef  ostream_iterator<ET>          Et_it;

        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Min_sphere_d( |P| = " << min_sphere.number_of_points()
               << ", |S| = " << min_sphere.number_of_support_points() << endl;
            os << "  P = {" << endl;
            os << "    ";
            copy( min_sphere.points_begin(), min_sphere.points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S = {" << endl;
            os << "    ";
            copy( min_sphere.support_points_begin(),
                  min_sphere.support_points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  center = ( ";
            copy( min_sphere.center_coordinates_begin(),
                  min_sphere.center_coordinates_end(),
                  Et_it( os, " "));
            os << ")" << endl;
            os << "  squared radius = "
               << min_sphere.squared_radius_numerator() << " / "
               << min_sphere.squared_radius_denominator() << endl;
            os << ")" << endl;
            break;

          case CGAL::IO::ASCII:
            copy( min_sphere.points_begin(), min_sphere.points_end(),
                  Os_it( os, "\n"));
            break;

          case CGAL::IO::BINARY:
            copy( min_sphere.points_begin(), min_sphere.points_end(),
                  Os_it( os));
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    // input operator
    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is, CGAL::Min_sphere_d<Traits_>& min_sphere)
    {
        CGAL_USING_NAMESPACE_STD
        
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            typedef  CGAL::Min_sphere_d<Traits_>::Point  Point;
            typedef  istream_iterator<Point>             Is_it;
            min_sphere.set( Is_it( is), Is_it());
            break;

          default:
            CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end


@! ----------------------------------------------------------------------------
@! Using the Quadratic Programming Solver
@! ----------------------------------------------------------------------------

\subsection{Using the Quadratic Programming Solver}
\label{sec:using_qp_solver}

We use the solver described in~\cite{s-qpego1-00} to determine the solution
of the quadratic programming problem~(\ref{eq:MS_as_QP}).

@macro <Min_sphere_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_QP_SOLVER_H
    #  include <CGAL/_QP_solver/QP_solver.h>
    #endif
@end


@! ----------------------------------------------------------------------------
\subsubsection{Representing the Quadratic Program}

We need a model of the concept \ccc{QP_representation}, which defines the
number types and iterators used by the QP solver.

@macro <Min_sphere_d declarations> += @begin
    
    template < class ET_, class NT_, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_min_sphere_d;
@end

@macro <Min_sphere_d QP representation> = @begin
    template < class ET_, class NT_, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_min_sphere_d {
        typedef  ET_                    ET;
        typedef  NT_                    NT;

        @<Min_sphere_d QP representation: iterator types>

        typedef  CGAL::Tag_false        Is_lp;
    };
@end

The matrix $A$ has only one row filled with $1$s, the vector $b$ has
exactly one $1$-entry. We use the class template
\ccc{Const_value_iterator<T>} to represent $A$ and $b$.

@macro <Min_sphere_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_CONST_VALUE_ITERATOR_H
    #  include <CGAL/_QP_solver/Const_value_iterator.h>
    #endif
@end

@macro <Min_sphere_d QP representation: iterator types> += @begin
    typedef  CGAL::Const_value_iterator< CGAL::Const_value_iterator<NT> >
                                        A_iterator;
    typedef  CGAL::Const_value_iterator<NT>
                                        B_iterator;
@end

The vector $c$ is stored in the data member \ccc{c_vector}.

@macro <Min_sphere_d data members> += @begin

    std::vector<NT>          c_vector;  // vector `c' of QP
@end

@macro <Min_sphere_d QP representation: iterator types> += @begin
    typedef  typename std::vector<NT>::const_iterator
                                        C_iterator;
@end

Because of its size ($n\!\times\!n$), the matrix $D$ is represented
implicitly. By~(\ref{eq:MS_as_QP}) we have $D = C^T C$, i.e.~$D_{i,j} =
p_i^T p_j$. Row~$i$ of $D$ is determined by $p_i$ and an iterator to the
point set. The entry in column~$j$ then is the inner product of $p_i$ and
$p_j$.

@macro <Min_sphere_d inner-product function class> = @begin
    template < class NT, class Point,
               class Access_coord, class Access_dim >
    class QP_rep_inner_product
      : public std::unary_function< Point, NT > {
        Point         p_i;
        Access_coord  da_coord;
        Access_dim    da_dim;
      public:
        QP_rep_inner_product( ) { }
        QP_rep_inner_product( const Point& p,
                              const Access_coord& ac,
                              const Access_dim&   ad)
            : p_i( p), da_coord( ac), da_dim( ad) { }

        NT  operator( ) ( const Point& p_j) const
            { return std::inner_product( da_coord( p_i),
                                         da_coord( p_i)+da_dim( p_i),
                                         da_coord( p_j), NT( 0),
                                         std::plus<NT>(),
                                         std::multiplies<NT>()); }
    };
@end

@macro <Min_sphere_d declarations> += @begin
    
    template < class NT, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_row_of_d;
@end

@macro <Min_sphere_d row-of-D function class> = @begin
    template < class NT, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    class QP_rep_row_of_d {
        Point_iterator  pts_it;
        Access_coord    da_coord;
        Access_dim      da_dim;
    public:
        typedef  CGAL::QP_rep_inner_product<
                     NT, Point, Access_coord, Access_dim >
                                        Inner_product;
        typedef  CGAL::Join_random_access_iterator_1<
                     Point_iterator, Inner_product >
                                        Row_of_d;

        typedef  Point                  argument_type;
        typedef  Row_of_d               result_type;
                                
        QP_rep_row_of_d( ) { }
        QP_rep_row_of_d( const Point_iterator& it,
                         const Access_coord&   ac,
                         const Access_dim&     ad)
            : pts_it( it), da_coord( ac), da_dim( ad) { }

        Row_of_d  operator( ) ( const Point& p_i) const
            { return Row_of_d( pts_it, Inner_product( p_i, da_coord, da_dim));}
    };
@end

@macro <Min_sphere_d QP representation: iterator types> += @begin
    typedef  CGAL::Join_random_access_iterator_1<
                 Point_iterator, QP_rep_row_of_d<
                     NT, Point, Point_iterator,
                     Access_coord, Access_dim > >
                                        D_iterator;
@end

Now we are able to define the fully specialized type of the QP solver.

@macro <Min_sphere_d Solver type> = @begin
    // QP solver
    typedef  CGAL::QP_rep_min_sphere_d<
                 ET, NT, Point, typename std::vector<Point>::const_iterator,
                 Access_coordinates_begin_d, Access_dimension_d >
                                        QP_rep;
    typedef  CGAL::QP_solver< QP_rep >  Solver;
    typedef  typename Solver::Pricing_strategy
                                        Pricing_strategy;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Computing the Smallest Enclosing Sphere}

We set up the quadratic program, solve it, and compute center and squared
radius of the smallest enclosing sphere.

@macro <Min_sphere_d private member functions> += @begin
    
    // compute smallest enclosing sphere
    void
    compute_min_sphere( )
    {
        if ( is_empty()) {
            center_coords.resize( 1);
            sqr_rad_numer = -ET( 1);
            return;
        }
        
        // set up and solve QP
        @<Min_sphere_d compute_min_sphere: set up and solve QP>

        // compute center and squared radius
        @<Min_sphere_d compute_min_sphere: compute center and ...>
    }
@end

@macro <Min_sphere_d compute_min_sphere: set up and solve QP> = @begin
    c_vector.resize( points.size());
    int i;
    for ( i = 0; i < number_of_points(); ++i) {
        c_vector[ i] = -std::inner_product(
            tco.access_coordinates_begin_d_object()( points[ i]),
            tco.access_coordinates_begin_d_object()( points[ i])+d,
            tco.access_coordinates_begin_d_object()( points[ i]), NT( 0),
            std::plus<NT>(), std::multiplies<NT>());
    }
    typedef  typename QP_rep::A_iterator A_it;
    typedef  typename QP_rep::B_iterator B_it;
    typedef  typename QP_rep::D_iterator D_it;
    B_it  const_one( 1);
    solver.set( points.size(), 1, d+1,
                A_it( const_one),
                const_one, c_vector.begin(),
                D_it( points.begin(),
                                             CGAL::QP_rep_row_of_d< NT,
                                             Point,
                                             Point_iterator,
                                             Access_coordinates_begin_d,
                                             Access_dimension_d >(
                        points.begin(),
                        tco.access_coordinates_begin_d_object(),
                        tco.access_dimension_d_object())));
    solver.init();
    solver.solve();
@end

@macro <Min_sphere_d compute_min_sphere: compute center and ...> = @begin
    center_coords.resize( ambient_dimension()+1);
    std::fill( center_coords.begin(), center_coords.end(), ET( 0));
    for ( i = 0; i < solver.number_of_basic_variables(); ++i) {
        ET  value = solver.basic_variables_numerator_begin()[ i];
        int index = solver.basic_variables_index_begin()[ i];
        for ( int j = 0; j < d; ++j)
            center_coords[ j] += value
            * tco.access_coordinates_begin_d_object()( points[ index])[ j];
    }
    center_coords[ d] = solver.variables_common_denominator();
    sqr_rad_numer     = -solver.solution_numerator();
    sqr_rad_denom     = center_coords[ d] * center_coords[ d];
@end

@! ----------------------------------------------------------------------------
\subsubsection{Choosing the Pricing Strategy}

@macro <Min_sphere_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_PARTIAL_EXACT_PRICING_H
    #  include <CGAL/_QP_solver/Partial_exact_pricing.h>
    #endif
    #ifndef CGAL_PARTIAL_FILTERED_PRICING_H
    #  include <CGAL/_QP_solver/Partial_filtered_pricing.h>
    #endif
@end

@macro <Min_sphere_d data members> += @begin
    
    typename Solver::Pricing_strategy*  // pricing strategy
                             strategyP; // of the QP solver
@end

@macro <Min_sphere_d QP-solver set-up> many = @begin
    set_pricing_strategy( NT());
@end

@macro <Min_sphere_d private member functions> += @begin
    
    #ifdef _MSC_VER
    
    template < class NT >
    void  set_pricing_strategy( NT)
    { }
    
    #else
    
    template < class NT >
    void  set_pricing_strategy( NT)
        { strategyP = new CGAL::Partial_filtered_pricing<QP_rep>;
          solver.set_pricing_strategy( *strategyP); }

    void  set_pricing_strategy( ET)
        { strategyP = new CGAL::Partial_exact_pricing<QP_rep>;
          solver.set_pricing_strategy( *strategyP); }

    #endif
@end


@! ============================================================================
@! Test Programs
@! ============================================================================

\clearpage
\section{Test Programs} \label{sec:test_programs}

We provide three programs for testing in two-, three- and
$d$-dimensional space, namely \ccc{test_Min_sphere_d_2},
\ccc{test_Min_sphere_d_3}, and \ccc{test_Min_sphere_d_d}.  Verbose
output can be enabled by giving a number between 0 and 3 at the
command line.
  
@macro <Min_sphere_d test: command line argument> many = @begin
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr ( verbose >= 0); verr  << "";
    CGAL::Verbose_ostream  verr0( verbose == 0); verr0 << "";
    CGAL::Verbose_ostream  verrX( verbose >  0); verrX << "";
@end

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

The function \ccc{test_Min_sphere_d}, invoked with a set of points and a
traits class model, calls each function of \ccc{Min_sphere_d} at least once
to ensure code coverage. If \ccc{verbose} is set to $-1$, the function is
``silent'', otherwise some diagnosing output is written to the standard
error stream.

@macro <Min_sphere_d test function> = @begin
    #define COVER(text,code) \
                verr0.out().width( 26); verr0 << text << "..." << flush; \
                verrX.out().width(  0); verrX << "==> " << text << endl \
                  << "----------------------------------------" << endl; \
                { code } verr0 << "ok."; verr << endl;
    
    template < class ForwardIterator, class Traits >
    void
    test_Min_sphere_d( ForwardIterator first, ForwardIterator last,
                       const Traits& traits, int verbose)
    {
        CGAL_USING_NAMESPACE_STD
        
        typedef  CGAL::Min_sphere_d< Traits >  Min_sphere;
        typedef  typename Traits::Point_d      Point;

        CGAL::Verbose_ostream  verr ( verbose >= 0);
        CGAL::Verbose_ostream  verr0( verbose == 0);
        CGAL::Verbose_ostream  verrX( verbose >  0);
        CGAL::set_pretty_mode( verr.out());

        bool  is_valid_verbose = ( verbose > 0);

        // constructors
        COVER( "default constructor",
            Min_sphere  ms( traits, verbose, verr.out());
            assert( ms.is_valid( is_valid_verbose));
            assert( ms.is_empty());
        )

        COVER( "point set constructor",
            Min_sphere  ms( first, last, traits, verbose, verr.out());
            verrX << endl;
            assert( ms.is_valid( is_valid_verbose));
        )

        Min_sphere  min_sphere( first, last);
        COVER( "ambient dimension",
            Min_sphere  ms;
            assert( ms.ambient_dimension() == -1);
            verrX << min_sphere.ambient_dimension() << endl;
        )

        COVER( "(number of) points",
            verrX << min_sphere.number_of_points() << endl;
            typename Min_sphere::Point_iterator
                point_it = min_sphere.points_begin();
            for ( ; point_it != min_sphere.points_end(); ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_sphere.points_end() - min_sphere.points_begin())
                    == min_sphere.number_of_points());
        )
        
        COVER( "(number of) support points",
            verrX << min_sphere.number_of_support_points() << endl;
            typename Min_sphere::Support_point_iterator
                point_it = min_sphere.support_points_begin();
            for ( ; point_it != min_sphere.support_points_end(); ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_sphere.support_points_end()
                      - min_sphere.support_points_begin())
                    == min_sphere.number_of_support_points());
        )

        COVER( "center and squared radius",
            verrX << "center (as range):";
            typename Min_sphere::Coordinate_iterator  coord_it;
            for ( coord_it  = min_sphere.center_coordinates_begin();
                  coord_it != min_sphere.center_coordinates_end();
                  ++coord_it) {
                verrX << ' ' << *coord_it;
            }
            verrX << endl;
            verrX << "squared radius numerator  : "
                  << min_sphere.squared_radius_numerator()   << endl;
            verrX << "squared radius denominator: "
                  << min_sphere.squared_radius_denominator() << endl;
        )

        COVER( "predicates",
            CGAL::Bounded_side  bounded_side;
            bool                has_on_bounded_side;
            bool                has_on_boundary;
            bool                has_on_unbounded_side;
            Point               p;
            typename Min_sphere::Point_iterator
                point_it = min_sphere.points_begin();
            for ( ; point_it != min_sphere.points_end(); ++point_it) {
                p = *point_it;
                bounded_side          = min_sphere.bounded_side( p);
                has_on_bounded_side   = min_sphere.has_on_bounded_side( p);
                has_on_boundary       = min_sphere.has_on_boundary( p);
                has_on_unbounded_side = min_sphere.has_on_unbounded_side( p);
                verrX.out().width( 2);
                verrX << bounded_side          << "  "
                      << has_on_bounded_side   << ' '
                      << has_on_boundary       << ' '
                      << has_on_unbounded_side << endl;
                assert( bounded_side != CGAL::ON_UNBOUNDED_SIDE);
                assert( has_on_bounded_side || has_on_boundary);
                assert( ! has_on_unbounded_side);
            }
        )        

        COVER( "clear",
            min_sphere.clear();
            verrX << "min_sphere is" << ( min_sphere.is_empty() ? "" : " not")
                  << " empty." << endl;
            assert( min_sphere.is_empty());
        )

        COVER( "insert (single point)",
            min_sphere.insert( *first);
            assert( min_sphere.is_valid( is_valid_verbose));
            assert( min_sphere.is_degenerate());
        )

        COVER( "insert (point set)",
            min_sphere.insert( first, last);
            assert( min_sphere.is_valid( is_valid_verbose));
        )

        COVER( "traits class access",
            min_sphere.traits();
        )

        COVER( "I/O",
            verrX << min_sphere;
        )
    }
@end


@! ----------------------------------------------------------------------------
@! Traits Class Models
@! ----------------------------------------------------------------------------

\subsection{Traits Class Models}

We perform the tests with the traits class models
\ccc{Optimisation_d_traits_2}, \ccc{Optimisation_d_traits_3}, and
\ccc{Optimisation_d_traits_d} based on the two-, three-, and
$d$-dimensional \cgal~kernel, respectively. All three traits class models
are used twice, firstly with one exact number type (the ``default'' use)
and secondly with three different number types (the ``advanced'' use).
Since the current implementation of the underlying linear programming
solver can only handle input points with Cartesian representation, we use
\cgal's Cartesian kernel for testing.  (The homogeneous kernel is used in
the additional test with other implementations described in
Section~\ref{sec:other_implementations}.)

Some of the following macros are parameterized with the dimension,
e.g.~with $2$, $3$, or $d$.

@macro <Min_sphere_d test: includes>(1) many += @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Homogeneous.h>
    #include <CGAL/Min_sphere_d_new.h>
    #include <CGAL/Optimisation_d_traits_@1.h>
@end

We use the number type \ccc{leda_integer} from \leda{} for the first
variant.

@macro <Min_sphere_d test: typedefs>(1) many += @begin
    // test variant 1 (needs LEDA)
    #ifdef CGAL_USE_LEDA
    # include <CGAL/leda_integer.h>
      typedef  CGAL::Cartesian<leda_integer>       K_1;
      typedef  CGAL::Optimisation_d_traits_@1<K_1>  Traits_1;
    # define TEST_VARIANT_1 \
        "Optimisation_d_traits_@1< Cartesian<leda_integer> >"
      CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_integer)
    #endif
@end

The second variant uses points with \ccc{int} coordinates. The exact number
type used by the underlying quadratic programming solver is
\ccc{GMP::Double}, i.e.~an arbitrary precise floating-point type based on
\textsc{Gmp}'s integers. To speed up the pricing, we use \ccc{double}
arithmetic.

@macro <Min_sphere_d test: typedefs> += @begin

    // test variant 2 (needs GMP)
    #ifdef CGAL_USE_GMP
    # include <CGAL/_QP_solver/Double.h>
      typedef  CGAL::Cartesian< int >                                 K_2;
      typedef  CGAL::Optimisation_d_traits_@1<K_2,GMP::Double,double>  Traits_2;
    # define TEST_VARIANT_2 \
        "Optimisation_d_traits_@1< Cartesian<int>, GMP::Double, double >"
    #endif
@end

The test sets consist of $100$ points with $20$-bit random integer
coordinates. In $2$- and $3$-space we use \cgal's point generators to build
the test sets with points lying almost (due to rounding errors) on a circle
or sphere, respectively.

@macro <Min_sphere_d test: includes> += @begin
    
    #include <CGAL/Random.h>
    #include <vector>
@end

@macro <Min_sphere_d test: includes (2/3D)>(1) many = @begin
    
    #include <CGAL/point_generators_@1.h>
    #include <CGAL/copy_n.h>
    #include <iterator>
@end

@macro <Min_sphere_d test: generate point set>(3) = @begin
    std::vector<K_@1::Point_@2>  points_@1;
    points_@1.reserve( 100);
    CGAL::copy_n( CGAL::Random_points_on_@3_@2<K_@1::Point_@2>( 0x100000),
                  100, std::back_inserter( points_@1));
@end

The traits class model with $d$-dimensional points is tested with $d = 5$
(variant 1) and $d = 10$ (variant 2). The points are distributed uniformly
in a $d$-cube.

@macro <Min_sphere_d test: generate point set (dD)>(1) = @begin
    std::vector<K_@1::Point_d>  points_@1;
    points_@1.reserve( 100);
    {
        int d = 5*@1;
        std::vector<int>  coords( d);
        int  i, j;
        for ( i = 0; i < 100; ++i) {
            for ( j = 0; j < d; ++j)
                coords[ j] = CGAL::default_random( 0x100000);
            points_@1.push_back( K_@1::Point_d( d, coords.begin(),
                                                   coords.end()));
        }
    }
@end

Finally we call the test function (described in the last section).

@macro <Min_sphere_d test: includes> += @begin

    #include "test_Min_sphere_d.h"
@end

@macro <Min_sphere_d test: call test function>(1) many = @begin
    CGAL::test_Min_sphere_d( points_@1.begin(), points_@1.end(),
                             Traits_@1(), verbose);
@end

Each of the two test variants is compiled and executed only if the
respective number type is available.

@macro <Min_sphere_d test: test variant output>(1) many = @begin
    verr << endl
         << "==================================="
         << "===================================" << endl
         << "Testing `Min_sphere_d' with traits class model" << endl
         << "==> " << TEST_VARIANT_@1 << endl
         << "==================================="
         << "===================================" << endl
         << endl;
@end

@macro <Min_sphere_d test: test variant>(3) many = @begin
    #ifdef TEST_VARIANT_@1
        @<Min_sphere_d test: test variant output>(@1)

        // generate point set
        @<Min_sphere_d test: generate point set>(@1,@2,@3)

        // call test function
        @<Min_sphere_d test: call test function>(@1)
    #endif
@end

@macro <Min_sphere_d test: test variant (dD)>(1) many = @begin
    #ifdef TEST_VARIANT_@1

        @<Min_sphere_d test: test variant output>(@1)

        // generate point set
        @<Min_sphere_d test: generate point set (dD)>(@1)

        // call test function
        @<Min_sphere_d test: call test function>(@1)

    #endif
@end

@! ----------------------------------------------------------------------------
@! Other Implementations
@! ----------------------------------------------------------------------------

\subsection{Other Implementations}
\label{sec:other_implementations}

As an additional correctness check, we compare the results of our
implementation to those of other implementations available in \cgal. In 2D,
we use the class template \ccc{Min_circle_2<Traits>} parameterized with the
traits class model \ccc{Min_circle_2_traits_2<R>}.

@macro <Min_sphere_d test: includes (2D)> = @begin
    #include <CGAL/Min_circle_2.h>
    #include <CGAL/Min_circle_2_traits_2.h>
@end

In 3D and dD, we use the ``other'' class template
\ccc{Min_sphere_d<Traits>} parameterized with the traits class models
\ccc{Optimisation_d_traits_3<R>} and \ccc{Optimisation_d_traits_d<R>},
respectively. To avoid name conflicts, we ``define'' the prefix
\ccc{OTHER_}.

@macro <Min_sphere_d test: includes (3/dD)> many = @begin
    #define  Min_sphere_d               OTHER_Min_sphere_d
    #undef  CGAL_MIN_SPHERE_D_H
    #undef  CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #define CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #include <CGAL/Min_sphere_d.h>
    #undef  Min_sphere_d
@end

The traits class models of the other algorithms use the representation class
model \ccc{Homogeneous<leda_integer>}.

@macro <Min_sphere_d test: typedefs (2D)> = @begin
    
    // comparing (needs LEDA)
    #ifdef CGAL_USE_LEDA
      typedef  CGAL::Homogeneous<leda_integer>   K_3;
      typedef  CGAL::Min_circle_2_traits_2<K_3>  Traits_3;
      typedef  CGAL::Min_sphere_d<Traits_1>      Min_sphere_d;
      typedef  CGAL::Min_circle_2<Traits_3>      O_Min_sphere_d;
    #endif
@end

@macro <Min_sphere_d test: typedefs (3/dD)>(1) many = @begin
    
    // comparing (needs LEDA)
    #ifdef CGAL_USE_LEDA
      typedef  CGAL::Homogeneous<leda_integer>     K_3;
      typedef  CGAL::Optimisation_d_traits_@1<K_3>  Traits_3;
      typedef  CGAL::Min_sphere_d<Traits_1>        Min_sphere_d;
      typedef  CGAL::Min_sphere_d<Traits_3>        O_Min_sphere_d;
    #endif
@end

In order to reuse the points from the first test variant (see above), we
have to convert them to points with homogeneous representation.

@macro <Min_sphere_d test: convert point set (2D)> = @begin
    std::vector<K_3::Point_2>  points_3;
    points_3.reserve( points_1.size());
    {
        unsigned int i;
        for ( i = 0; i < points_1.size(); ++i) {
            points_3.push_back( K_3::Point_2( points_1[ i][ 0],
                                              points_1[ i][ 1]));
        }
    }
@end

@macro <Min_sphere_d test: convert point set (3D)> = @begin
    std::vector<K_3::Point_3>  points_3;
    points_3.reserve( points_1.size());
    {
        unsigned int i;
        for ( i = 0; i < points_1.size(); ++i) {
            points_3.push_back( K_3::Point_3( points_1[ i][ 0],
                                              points_1[ i][ 1],
                                              points_1[ i][ 2]));
        }
    }
@end

@macro <Min_sphere_d test: convert point set (dD)> = @begin
    std::vector<K_3::Point_d>  points_3;
    points_3.reserve( points_1.size());
    {
        int          d = points_1[ 0].dimension();
        unsigned int i;
        for ( i = 0; i < points_1.size(); ++i) {
            points_3.push_back( K_3::Point_d( d, points_1[ i].begin(),
                                                 points_1[ i].end()));
        }
    }
@end

Both algorithms to compare compute the smallest enclosing sphere and are
checked for validity.

@macro <Min_sphere_d test: compute smallest enclosing spheres>(1) many = @begin
    Min_sphere_d  ms( points_1.begin(), points_1.end(),
                      Traits_1(), verbose);
    verrX << endl << ms << endl;
    assert( ms.is_valid( verbose > 0));
        
    O_Min_sphere_d  o_ms( points_3.begin(), points_3.end(), @1);
    verrX << endl << o_ms << endl;
    assert( o_ms.is_valid( verbose > 0));
    verrX << endl;
@end

Finally we check whether center and squared radius are the same.

@macro <Min_sphere_d test: check center and squared radius>(1) many = @begin
    COVER( "center",
        O_Min_sphere_d::Point  o_ms_center = o_ms@1.center();
    
        verrX << "center (as point): " << ms.center()
              << "  [NOTE: coordinates are truncated!]" << endl;
        
        int           d     = points_1[ 0].dimension();
        leda_integer  den   = ms.center_coordinates_begin()[ d];
        leda_integer  o_den = o_ms_center.homogeneous( d);
        for ( int j = 0; j < d; ++j) {
            assert( ms.center_coordinates_begin()[ j]*o_den
                    == o_ms_center.homogeneous( j)*den);
        }
        verrX << "centers are equal." << endl;
    );
        
    COVER( "squared radius",
        verrX << "squared radius: " << ms.squared_radius()
              << "  [NOTE: value is truncated!]" << endl;
        
        assert( CGAL::Quotient<leda_integer>(
                    ms.squared_radius_numerator(),
                    ms.squared_radius_denominator())
                == o_ms@1.squared_radius());
        verrX << "squared radii are equal." << endl;
    );
@end

\ldots

@macro <Min_sphere_d test: additional test output>(1) many = @begin
    verr << endl
         << "==================================="
         << "===================================" << endl
         << "Comparing `Min_sphere_d' with `@1'" << endl
         << "==================================="
         << "===================================" << endl
         << endl;
@end

@macro <Min_sphere_d test: additional test (2D)> = @begin
    #ifdef CGAL_USE_LEDA

        @<Min_sphere_d test: additional test output>(Min_circle_2)

        // convert point set
        @<Min_sphere_d test: convert point set (2D)>

        // compute smallest enclosing spheres
        @<Min_sphere_d test: compute smallest enclosing spheres>("false")

        // check center and squared radius
        @<Min_sphere_d test: check center and squared radius>(".circle()")

    #endif
@end
    
@macro <Min_sphere_d test: additional test (3D)> = @begin
    #ifdef CGAL_USE_LEDA

        @<Min_sphere_d test: additional test output>(OTHER_Min_sphere_d)

        // convert point set
        @<Min_sphere_d test: convert point set (3D)>

        // compute smallest enclosing spheres
        @<Min_sphere_d test: compute smallest enclosing spheres>("Traits_3()")

        // check center and squared radius
        @<Min_sphere_d test: check center and squared radius>("")

    #endif
@end
    
@macro <Min_sphere_d test: additional test (dD)> = @begin
    #ifdef CGAL_USE_LEDA

        @<Min_sphere_d test: additional test output>(OTHER_Min_sphere_d)

        // convert point set
        @<Min_sphere_d test: convert point set (dD)>

        // compute smallest enclosing spheres
        @<Min_sphere_d test: compute smallest enclosing spheres>("Traits_3()")

        // check center and squared radius
        @<Min_sphere_d test: check center and squared radius>("")

    #endif
@end
    
@! ============================================================================
@! Files
@! ============================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Min_sphere_d_new.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Min\_sphere\_d.h}

@file <include/CGAL/Min_sphere_d_new.h> = @begin
    @<file header>(
        "include/CGAL/Min_sphere_d_new.h",
        "Smallest enclosing sphere in arbitrary dimension")

    #ifndef CGAL_MIN_SPHERE_D_H
    #define CGAL_MIN_SPHERE_D_H

    // includes
    // --------
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<Min_sphere_d CGAL includes>
    @<Min_sphere_d CGAL/QP_solver includes>
    @<Min_sphere_d standard includes>

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_sphere_d declarations>
    
    // Class interfaces
    // ================
    @<Min_sphere_d interface>

    @<Min_sphere_d inner-product function class>

    @<Min_sphere_d row-of-D function class>

    @<Min_sphere_d QP representation>

    // Function declarations
    // =====================
    @<Min_sphere_d I/O operators declaration>

    @<dividing line>
    
    // Class implementation
    // ====================

    @<Min_sphere_d validity check>
    
    @<Min_sphere_d I/O operators>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_SPHERE_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_sphere_d.h
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d.h}

@file <test/Min_sphere_d_new/test_Min_sphere_d.h> = @begin
    @<file header>(
        "test/Min_sphere_d_new/test_Min_sphere_d.h",
        "test function for smallest enclosing sphere")

    #ifndef CGAL_TEST_MIN_SPHERE_D_H
    #define CGAL_TEST_MIN_SPHERE_D_H

    // includes
    #ifndef CGAL_IO_VERBOSE_OSTREAM_H
    #  include <CGAL/IO/Verbose_ostream.h>
    #endif
    #include <cassert>
    
    @<namespace begin>("CGAL")

    @<Min_sphere_d test function>
    
    @<namespace end>("CGAL")

    #endif // CGAL_TEST_MIN_SPHERE_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_sphere_d_2.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_2.C}

@file <test/Min_sphere_d_new/test_Min_sphere_d_2.C> = @begin
    @<file header>(
        "test/Min_sphere_d_new/test_Min_sphere_d_2.C",
        "test program for smallest enclosing sphere (2D traits class)")

    // includes
    // --------
    @<Min_sphere_d test: includes>(2)
    @<Min_sphere_d test: includes (2D)>
    @<Min_sphere_d test: includes (2/3D)>(2)
    
    // typedefs
    // --------
    @<Min_sphere_d test: typedefs>(2)
    @<Min_sphere_d test: typedefs (2D)>
    
    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        CGAL_USING_NAMESPACE_STD

        // command line arguments
        @<Min_sphere_d test: command line argument>

        // code coverage
        // -------------
        @<Min_sphere_d test: test variant>(1,2,circle)

        @<Min_sphere_d test: test variant>(2,2,circle)

        // additional tests
        // ----------------
        @<Min_sphere_d test: additional test (2D)>

        return 0;
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_sphere_d_3.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_3.C}

@file <test/Min_sphere_d_new/test_Min_sphere_d_3.C> = @begin
    @<file header>(
        "test/Min_sphere_d_new/test_Min_sphere_d_3.C",
        "test program for smallest enclosing sphere (3D traits class)")

    // includes
    // --------
    @<Min_sphere_d test: includes>(3)
    @<Min_sphere_d test: includes (3/dD)>
    @<Min_sphere_d test: includes (2/3D)>(3)

    // typedefs
    // --------
    @<Min_sphere_d test: typedefs>(3)
    @<Min_sphere_d test: typedefs (3/dD)>(3)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        CGAL_USING_NAMESPACE_STD
    
        // command line arguments
        @<Min_sphere_d test: command line argument>
    
        // code coverage
        // -------------
        @<Min_sphere_d test: test variant>(1,3,sphere)

        @<Min_sphere_d test: test variant>(2,3,sphere)

        // additional tests
        // ----------------
        @<Min_sphere_d test: additional test (3D)>

        return 0;
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_sphere_d_d.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_d.C}

@file <test/Min_sphere_d_new/test_Min_sphere_d_d.C> = @begin
    @<file header>(
        "test/Min_sphere_d_new/test_Min_sphere_d_d.C",
        "test program for smallest enclosing sphere (dD traits class)")

    // includes
    // --------
    @<Min_sphere_d test: includes>(d)
    @<Min_sphere_d test: includes (3/dD)>

    // typedefs
    // --------
    @<Min_sphere_d test: typedefs>(d)
    @<Min_sphere_d test: typedefs (3/dD)>(d)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        CGAL_USING_NAMESPACE_STD

        // command line arguments
        @<Min_sphere_d test: command line argument>

        // code coverage
        // -------------
        @<Min_sphere_d test: test variant (dD)>(1)

        @<Min_sphere_d test: test variant (dD)>(2)

        // additional tests
        // ----------------
        @<Min_sphere_d test: additional test (dD)>


        return 0;
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
        "Min_sphere_d_new",
        "Geometric Optimisation",
        "Min_sphere_d",
        "$Revision$","$Date$",
        "Bernd Gärtner, Sven Schönherr <sven@@inf.ethz.ch>",
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
