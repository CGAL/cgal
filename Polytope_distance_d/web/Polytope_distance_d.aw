@! ============================================================================
@! The CGAL Library
@! Implementation: Distance of Polytopes in Arbitrary Dimension
@! ----------------------------------------------------------------------------
@! file  : web/Polytope_distance_d.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Polytope_distance_d WIP $
@! $Id$
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
\newcommand{\cgalPolytopeDistanceLayout}{\cgalColumnLayout}

\newcommand{\ccRequirements}{\ccHeading{Requirements}}
\newcommand{\ccRequire}{\ccCommentHeading{Requirements}}


@! ============================================================================
@! Title
@! ============================================================================

\thispagestyle{empty}

\RCSdef{\rcsRevision}{$Id$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 20 mm
@t title titlefont centre "Distance of Convex Polytopes"
@t vskip 0 mm
@t title titlefont centre "in Arbitrary Dimension*"
@t vskip 10 mm
@t title smalltitlefont centre "Sven Schönherr"
\begin{center}
  \textbf{ETH Z{\"u}rich}
\end{center}
@t vskip 10 mm
{\small
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Package: Polytope_distance_d WIP+\cgalWIP\verb+$+ \\
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
  We provide an implementation for computing the (squared) distance
  between to convex polytopes in arbitrary dimension. The problem is
  formulated as a quadratic program and a dedicated
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

We consider the problem of finding the distance between to convex
polytopes, given as the convex hulls of two finite point sets in
$d$-dimensional Euclidean space $\E_d$.  It can be formulated as an
optimization problem with linear constraints and a convex quadratic
objective function~\cite{gs-eegqp-00}.

@! ----------------------------------------------------------------------------
@! Polytope Distance as a Quadratic Programming Problem
@! ----------------------------------------------------------------------------

\subsection{Polytope Distance as a Quadratic Programming Problem}

If the point sets are given as $P = \{p_1,\dots,p_r\}$ and $Q =
\{q_1,\dots,q_s\}$, $n=r+s$, we want to find points $p^*$ and $q^*$ in the
convex hull of $P$ and $Q$, respectively, such that $|p^*-q^*|$ is
minimized.  Define the $d\!\times\!n$-matrix $C :=
(p_1,\dots,p_r,-q_1,\dots,-q_s)$ and consider the quadratic programming
problem
%
\begin{equation} \label{eq:PD_as_QP}
  \begin{array}{lll}
    \text{(PD)} & \text{minimize}   & x^T C^T C\, x \\[0.8ex]
                & \text{subject to} & \sum_{i=1}^r x_i     = 1, \\[0.5ex]
                &                   & \sum_{i=1}^s x_{i+r} = 1, \\[0.5ex]
                &                   & x \geq 0.
  \end{array}
\end{equation}
%
Let $x^* = (x^*_1,\dots,x^*_n)$ be its optimal solution, then the points
%
\begin{align*}
  p^* &= \sum_{i=1}^r x_i p_i, \\
  q^* &= \sum_{i=1}^s x_{i+r} q_i
\end{align*}
%
realize the distance between the two polytopes. The squared distance is the
value of the objective function at $x^*$.


@! ============================================================================
@! Reference Pages
@! ============================================================================

\clearpage
\section{Reference Pages} \label{sec:reference_pages}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

@p maximum_input_line_length = 102

@! ----------------------------------------------------------------------------
@! Class: Polytope_distance_d
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{CGAL::Polytope\_distance\_d\texttt{<}Traits\texttt{>}}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/Polytope_distance_d.tex}

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
@! The Class Template CGAL::Polytope_distance_d<Traits>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Polytope\_distance\_d\texttt{<}Traits\texttt{>}}

The class template \ccc{Polytope_distance_d} expects a model of the concept
\ccc{OptimisationDTraits} (see
Section~\ref{ccRef_OptimisationDTraits}.2) as its template argument.

@macro <Poly_dist_d declarations> += @begin
    template < class Traits_ >
    class Polytope_distance_d;
@end

The interface consists of the public types and member functions described
in Section~\ref{ccRef_CGAL::Polytope_distance_d<Traits>}.1 and of some
private types, private member functions, and data members.

@macro <Poly_dist_d interface> = @begin
    template < class Traits_ >
    class Polytope_distance_d {
      public:
        // self
        typedef  Traits_                    Traits;
        typedef  Polytope_distance_d<Traits>
                                            Self;

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
        @<Poly_dist_d Solver type>
        @<Poly_dist_d private types>

      public:
        @<Poly_dist_d types>
        
        @<Poly_dist_d member functions>
        
      private:
        @<Poly_dist_d data members>

        @<Poly_dist_d private member functions>
    };
@end

@! ----------------------------------------------------------------------------
\subsubsection{Data Members}

Mainly, we have to store the given input points, the two points realizing
the distance, and an instance of the quadratic programming solver.
Additional variables, that are used in the member functions described
below, are introduced when they appear for the first time.

We start with the traits class object.

@macro <Poly_dist_d data members> += @begin

    Traits                   tco;       // traits class object
@end

The inputs points are kept in a vector to have random access to them.
Their dimension is stored separately.

@macro <Poly_dist_d standard includes> += @begin
    #include <vector>
@end

@macro <Poly_dist_d private types> += @begin
    // private types
    typedef  std::vector<Point>         Point_vector;
@end

@macro <Poly_dist_d data members> += @begin

    Point_vector             p_points;  // points of P
    Point_vector             q_points;  // points of Q
    int                      d;         // dimension of input points
@end

The two points realizing the distance between the polytopes are stored with
rational representation, i.e.~numerators and denominators are kept
separately. Vectors \ccc{p_coords} and \ccc{q_coords} each contain $d+1$
entries, the numerators of the $d$ coordinates and the common denominator.

@macro <Poly_dist_d private types> += @begin
    typedef  std::vector<ET>            ET_vector;
@end

@macro <Poly_dist_d data members> += @begin
    
    ET_vector                p_coords;          // realizing point of P
    ET_vector                q_coords;          // realizing point of Q
@end

We store an instance of the quadratic programming solver described
in~\cite{s-qpego1-00}. The details are given in
Section~\ref{sec:using_qp_solver} below, here it suffice to know that there
is a variable \ccc{solver} of type \ccc{Solver}.

@macro <Poly_dist_d private types: quadratic programming solver> zero = @begin
    typedef  ...                        Solver;
@end

@macro <Poly_dist_d data members> += @begin

    Solver                   solver;    // quadratic programming solver
@end

@! ----------------------------------------------------------------------------
\subsubsection{Creation}

Two constructors are provided. If the user wants to get some verbose output
(of the underlying QP solver), he can override the default arguments of
\ccc{verbose} and \ccc{stream}.

@macro <Poly_dist_d standard includes> += @begin
    #include <iostream>
@end

@macro <Poly_dist_d member functions> += @begin
    // creation
    Polytope_distance_d( const Traits&  traits  = Traits(),
                         int            verbose = 0,
                         std::ostream&  stream  = std::cout)
      : tco( traits), d( -1), solver( verbose, stream)
        {
            @<Poly_dist_d QP-solver set-up>
        }
@end

The second constructor expects two sets of points given via iterator ranges.
It calls the \ccc{set} member function described in
Subsection~\ref{sec:modifiers} to store the points and to compute the
(squared) distance of the given polytopes.

@macro <Poly_dist_d member functions> += @begin

    template < class InputIterator1, class InputIterator2 >
    Polytope_distance_d( InputIterator1 p_first,
                         InputIterator1 p_last,
                         InputIterator2 q_first,
                         InputIterator2 q_last,
                         const Traits&  traits = Traits(),
                         int            verbose = 0,
                         std::ostream&  stream  = std::cout)
      : tco( traits), solver( verbose, stream)
        {
            @<Poly_dist_d QP-solver set-up>
            set( p_first, p_last, q_first, q_last);
        }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Access}

The following types and member functions give access to the points of the
given polytopes.

@macro <Poly_dist_d types> += @begin
    // public types
    typedef  typename Point_vector::const_iterator
                                        Point_iterator;
@end

@macro <Poly_dist_d member functions> += @begin

    // access to point sets
    int  ambient_dimension( ) const { return d; }

    int  number_of_points( ) const { return p_points.size()+q_points.size();}
    
    int  number_of_points_p( ) const { return p_points.size(); }
    int  number_of_points_q( ) const { return q_points.size(); }

    Point_iterator  points_p_begin( ) const { return p_points.begin(); }
    Point_iterator  points_p_end  ( ) const { return p_points.end  (); }

    Point_iterator  points_q_begin( ) const { return q_points.begin(); }
    Point_iterator  points_q_end  ( ) const { return q_points.end  (); }
@end

To access the support points, we exploit the following fact. A point~$p_i$
(or $q_i$) is a support point, iff its corresponding variable $x_i$ (of the
QP solver) is basic. Thus the number of support points is equal to the
number of basic variables, if the distance is finite.

@macro <Poly_dist_d member functions> += @begin

    // access to support points
    int
    number_of_support_points( ) const
        { return is_finite() ? solver.number_of_basic_variables() : 0; }
@end

Before we can access the support points of $P$ and $Q$, we have to divide
the set of basic variables into two sets corresponding to the support
points of $P$ and $Q$, respectively. The indices of the support points are
stored in \ccc{p_support_indices} and \ccc{q_support_indices}, while the
actual split-up is done in the private member function
\ccc{compute_distance} described below in
Section~\ref{sec:using_qp_solver}.

@macro <Poly_dist_d private types> += @begin

    typedef  std::vector<int>           Index_vector;
@end

@macro <Poly_dist_d data members> += @begin

    Index_vector             p_support_indices;
    Index_vector             q_support_indices;
@end

@macro <Poly_dist_d member functions> += @begin

    int  number_of_support_points_p() const { return p_support_indices.size();}
    int  number_of_support_points_q() const { return q_support_indices.size();}
@end

To access a point given its index, we use the following function class.

@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_ACCESS_BY_INDEX_H
    #  include <CGAL/_QP_solver/Access_by_index.h>
    #endif
@end

@macro <Poly_dist_d private types> += @begin

    typedef  CGAL::Access_by_index<typename std::vector<Point>::const_iterator>
                                        Point_by_index;
@end

Combining the function class with the index iterator gives the support
point iterator.

@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
    #  include <CGAL/_QP_solver/Join_random_access_iterator.h>
    #endif
@end

@macro <Poly_dist_d types> += @begin

    typedef typename Index_vector::const_iterator IVCI;

    typedef  CGAL::Join_random_access_iterator_1<
                 IVCI, Point_by_index >
                                        Support_point_iterator;
@end

@macro <Poly_dist_d member functions> += @begin

    Support_point_iterator
    support_points_p_begin() const
        { return Support_point_iterator(
                     p_support_indices.begin(),
                     Point_by_index( p_points.begin())); }

    Support_point_iterator
    support_points_p_end() const
        { return Support_point_iterator(
                     is_finite() ? p_support_indices.end()
                                 : p_support_indices.begin(),
                     Point_by_index( p_points.begin())); }

    Support_point_iterator
    support_points_q_begin() const
        { return Support_point_iterator(
                     q_support_indices.begin(),
                     Point_by_index( q_points.begin())); }

    Support_point_iterator
    support_points_q_end() const
        { return Support_point_iterator(
                     is_finite() ? q_support_indices.end()
                                 : q_support_indices.begin(),
                     Point_by_index( q_points.begin())); }
@end

The following types and member functions give access to the squared
distance and the realizing points.

@macro <Poly_dist_d types> += @begin

    typedef  typename ET_vector::const_iterator
                                        Coordinate_iterator;
@end

@macro <Poly_dist_d member functions> += @begin

    // access to realizing points (rational representation)
    Coordinate_iterator
    realizing_point_p_coordinates_begin( ) const { return p_coords.begin(); }

    Coordinate_iterator
    realizing_point_p_coordinates_end  ( ) const { return p_coords.end  (); }

    Coordinate_iterator
    realizing_point_q_coordinates_begin( ) const { return q_coords.begin(); }

    Coordinate_iterator
    realizing_point_q_coordinates_end  ( ) const { return q_coords.end  (); }

    // access to squared distance (rational representation)
    ET  squared_distance_numerator  ( ) const
        { return solver.solution_numerator(); }
    
    ET  squared_distance_denominator( ) const
        { return solver.solution_denominator(); }
@end

For convinience, we also provide member functions for accessing the squared
distance as a single number of type \ccc{FT} and the realizing points as
single points of type \ccc{Point}. These functions only work, if an
implicit conversion from number type \ccc{ET} to number type \ccc{RT} is
available, e.g.~if both types are the same.

@macro <Poly_dist_d member functions> += @begin

    // access to realizing points and squared distance
    // NOTE: an implicit conversion from ET to RT must be available!
    Point
    realizing_point_p( ) const
        { CGAL_optimisation_precondition( is_finite());
          return tco.construct_point_d_object()( ambient_dimension(),
                     realizing_point_p_coordinates_begin(),
                     realizing_point_p_coordinates_end  ()); }

    Point
    realizing_point_q( ) const
        { CGAL_optimisation_precondition( is_finite());
          return tco.construct_point_d_object()( ambient_dimension(),
                     realizing_point_q_coordinates_begin(),
                     realizing_point_q_coordinates_end  ()); }

    FT
    squared_distance( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_distance_numerator  ()) /
                 FT( squared_distance_denominator()); }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The distance of the two polytopes is \emph{finite} if none of the point
sets is empty, it is \emph{zero} if the polytopes intersect, and it is
\emph{degenerate}, if it is not finite or zero.

@macro <Poly_dist_d member functions> += @begin

    bool  is_finite( ) const
        { return ( number_of_points_p() > 0) && ( number_of_points_q() > 0); }
    
    bool  is_zero( ) const
        { return CGAL_NTS is_zero( squared_distance_numerator()); }
    
    bool  is_degenerate( ) const { return ( ! is_finite()); }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers} \label{sec:modifiers}

These private member functions are used by the following \ccc{set} and
\ccc{insert} member functions to set and check the dimension of the input
points, respectively.

@macro <Poly_dist_d CGAL includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_H
    #  include <CGAL/function_objects.h>
    #endif
@end

@macro <Poly_dist_d private member functions> += @begin

    // set dimension of input points
    void
    set_dimension( )
        { d = ( p_points.size() > 0 ?
                    tco.access_dimension_d_object()( p_points[ 0]) :
                q_points.size() > 0 ?
                    tco.access_dimension_d_object()( q_points[ 0]) :
                -1); }
        
    // check dimension of input points
    template < class InputIterator >
    bool
    check_dimension( InputIterator first, InputIterator last)
        { return ( std::find_if( first, last,
                                 CGAL::compose1_1( std::bind2nd(
                                     std::not_equal_to<int>(), d),
                                     tco.access_dimension_d_object()))
                   == last); }
@end

The \ccc{set} member function copies the input points into the internal
variables \ccc{p_points} and/or \ccc{q_points} and calls the private member
function \ccc{compute_distance} (described in
Section~\ref{sec:using_qp_solver}) to compute the (squared) distance of the
polytopes.

@macro <Poly_dist_d member functions> += @begin

    // modifiers
    template < class InputIterator1, class InputIterator2 >
    void
    set( InputIterator1 p_first, InputIterator1 p_last,
         InputIterator2 q_first, InputIterator2 q_last)
        { if ( p_points.size() > 0)
              p_points.erase( p_points.begin(), p_points.end());
          if ( q_points.size() > 0)
              q_points.erase( q_points.begin(), q_points.end());
          std::copy( p_first, p_last, std::back_inserter( p_points));
          std::copy( q_first, q_last, std::back_inserter( q_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
                 check_dimension( p_points.begin(), p_points.end())
              && check_dimension( q_points.begin(), q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
          
    template < class InputIterator >
    void
    set_p( InputIterator p_first, InputIterator p_last)
        { if ( p_points.size() > 0)
              p_points.erase( p_points.begin(), p_points.end());
          std::copy( p_first, p_last, std::back_inserter( p_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( p_points.begin(), p_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }

    template < class InputIterator >
    void
    set_q( InputIterator q_first, InputIterator q_last)
        { if ( q_points.size() > 0)
              q_points.erase( q_points.begin(), q_points.end());
          std::copy( q_first, q_last, std::back_inserter( q_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( q_points.begin(), q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
@end

The \ccc{insert} member functions append the given point(s) to the point
set(s) and recompute the (squared) distance.

@macro <Poly_dist_d member functions> += @begin

    void
    insert_p( const Point& p)
        { CGAL_optimisation_precondition( ( ! is_finite()) || 
              ( tco.access_dimension_d_object()( p) == d));
          p_points.push_back( p);
          compute_distance(); }

    void
    insert_q( const Point& q)
        { CGAL_optimisation_precondition( ( ! is_finite()) || 
              ( tco.access_dimension_d_object()( q) == d));
          q_points.push_back( q);
          compute_distance(); }

    template < class InputIterator1, class InputIterator2 >
    void
    insert( InputIterator1 p_first, InputIterator1 p_last,
            InputIterator2 q_first, InputIterator2 q_last)
        { CGAL_optimisation_precondition_code( int old_r = p_points.size());
          CGAL_optimisation_precondition_code( int old_s = q_points.size());
          p_points.insert( p_points.end(), p_first, p_last);
          q_points.insert( q_points.end(), q_first, q_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
                 check_dimension( p_points.begin()+old_r, p_points.end())
              && check_dimension( q_points.begin()+old_s, q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }

    template < class InputIterator >
    void
    insert_p( InputIterator p_first, InputIterator p_last)
        { CGAL_optimisation_precondition_code( int old_r = p_points.size());
          p_points.insert( p_points.end(), p_first, p_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( p_points.begin()+old_r, p_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }

    template < class InputIterator >
    void
    insert_q( InputIterator q_first, InputIterator q_last)
        { CGAL_optimisation_precondition_code( int old_s = q_points.size());
          q_points.insert( q_points.end(), q_first, q_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( q_points.begin()+old_s, q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
@end

The \ccc{clear} member function deletes all points and resets the (squared)
distance to infinity.

@macro <Poly_dist_d member functions> += @begin

    void
    clear( )
        { p_points.erase( p_points.begin(), p_points.end());
          q_points.erase( q_points.begin(), q_points.end());
          compute_distance(); }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Validity Check}

A \ccc{Polytope_distance_d<Traits>} object can be checked for validity.
This means, it is checked whether the polytopes are separated by two
hyperplanes orthogonal to vector $p-q$ and passing through $p$ and $q$,
respectively.  The function \ccc{is_valid} is mainly intended for debugging
user supplied traits classes but also for convincing the anxious user that
the traits class implementation is correct. If \ccc{verbose} is \ccc{true},
some messages concerning the performed checks are written to standard error
stream. The second parameter \ccc{level} is not used, we provide it only
for consistency with interfaces of other classes.

@macro <Poly_dist_d member functions> += @begin

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
@end

@macro <Poly_dist_d validity check> = @begin
    // validity check
    template < class Traits_ >
    bool
    Polytope_distance_d<Traits_>::
    is_valid( bool verbose, int level) const
    {
        using namespace std;
        
        CGAL::Verbose_ostream verr( verbose);
        verr << "CGAL::Polytope_distance_d<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P+Q| = " << number_of_points_p()
             <<          '+' << number_of_points_q()
             <<   ", |S| = " << number_of_support_points_p()
             <<          '+' << number_of_support_points_q() << endl;

        if ( is_finite()) {
            
            // compute normal vector
            ET_vector  normal( d), diff( d);
            ET  et_0 = 0, den = solver.variables_common_denominator();
            int i, j;
            for ( j = 0; j < d; ++j) normal[ j] = p_coords[ j] - q_coords[ j];

            // check P
            // -------
            @<Poly_dist_d validity check: check polytope>(P,p,>)

            // check Q
            // -------
            @<Poly_dist_d validity check: check polytope>(Q,q,<)
        }
        
        verr << "  object is valid!" << endl;
        return( true);
    }
@end

@macro <Poly_dist_d validity check: check polytope>(3) many = @begin
    verr << "  checking @1..." << flush;

    // check point set
    for ( i = 0; i < number_of_points_@2(); ++i) {
        for ( j = 0; j < d; ++j) {
            diff[ j] = @2_coords[ j] - den
              * tco.access_coordinates_begin_d_object()( @2_points[ i])[ j];
        }
        if ( std::inner_product( diff.begin(), diff.end(),
                                 normal.begin(), et_0) @3 et_0)
            return CGAL::_optimisation_is_valid_fail( verr,
                       "polytope @1 is not separated by its hyperplane");
    }

    verr << "passed." << endl;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Miscellaneous}

The member function \ccc{traits} returns a const reference to the
traits class object.

@macro <Poly_dist_d member functions> += @begin

    // traits class access
    const Traits&  traits( ) const { return tco; }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Poly_dist_d I/O operators declaration> = @begin
    // I/O operators
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Polytope_distance_d<Traits_>& poly_dist);

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is,
                  Polytope_distance_d<Traits_>& poly_dist);
@end

@macro <Poly_dist_d I/O operators> = @begin
    // output operator
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Polytope_distance_d<Traits_>& poly_dist)
    {
        using namespace std;

        typedef  Polytope_distance_d<Traits_>::Point  Point;
        typedef  ostream_iterator<Point>       Os_it;
        typedef  typename Traits_::ET          ET;
        typedef  ostream_iterator<ET>          Et_it;

        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Polytope_distance_d( |P+Q| = "
               << poly_dist.number_of_points_p() << '+'
               << poly_dist.number_of_points_q() << ", |S| = "
               << poly_dist.number_of_support_points_p() << '+'
               << poly_dist.number_of_support_points_q() << endl;
            os << "  P = {" << endl;
            os << "    ";
            copy( poly_dist.points_p_begin(), poly_dist.points_p_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  Q = {" << endl;
            os << "    ";
            copy( poly_dist.points_q_begin(), poly_dist.points_q_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S_P = {" << endl;
            os << "    ";
            copy( poly_dist.support_points_p_begin(),
                  poly_dist.support_points_p_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S_Q = {" << endl;
            os << "    ";
            copy( poly_dist.support_points_q_begin(),
                  poly_dist.support_points_q_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  p = ( ";
            copy( poly_dist.realizing_point_p_coordinates_begin(),
                  poly_dist.realizing_point_p_coordinates_end(),
                  Et_it( os, " "));
            os << ")" << endl;
            os << "  q = ( ";
            copy( poly_dist.realizing_point_q_coordinates_begin(),
                  poly_dist.realizing_point_q_coordinates_end(),
                  Et_it( os, " "));
            os << ")" << endl;
            os << "  squared distance = "
               << poly_dist.squared_distance_numerator() << " / "
               << poly_dist.squared_distance_denominator() << endl;
            break;

          case CGAL::IO::ASCII:
            os << poly_dist.number_of_points_p() << endl;
            copy( poly_dist.points_p_begin(),
                  poly_dist.points_p_end(),
                  Os_it( os, "\n"));
            os << poly_dist.number_of_points_q() << endl;
            copy( poly_dist.points_q_begin(),
                  poly_dist.points_q_end(),
                  Os_it( os, "\n"));
            break;

          case CGAL::IO::BINARY:
            os << poly_dist.number_of_points_p() << endl;
            copy( poly_dist.points_p_begin(),
                  poly_dist.points_p_end(),
                  Os_it( os));
            os << poly_dist.number_of_points_q() << endl;
            copy( poly_dist.points_q_begin(),
                  poly_dist.points_q_end(),
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
    operator >> ( std::istream& is,
                  CGAL::Polytope_distance_d<Traits_>& poly_dist)
    {
        using namespace std;
        /*        
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            typedef  CGAL::Polytope_distance_d<Traits_>::Point  Point;
            typedef  istream_iterator<Point>             Is_it;
            poly_dist.set( Is_it( is), Is_it());
            break;

          default:
            CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
            break; }
        */
        return( is);
    }
@end


@! ----------------------------------------------------------------------------
@! Using the Quadratic Programming Solver
@! ----------------------------------------------------------------------------

\subsection{Using the Quadratic Programming Solver}
\label{sec:using_qp_solver}

We use the solver described in~\cite{s-qpego1-00} to determine the solution
of the quadratic programming problem~(\ref{eq:PD_as_QP}).

@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_QP_SOLVER_H
    #  include <CGAL/_QP_solver/QP_solver.h>
    #endif
@end


@! ----------------------------------------------------------------------------
\subsubsection{Representing the Quadratic Program}

We need a model of the concept \ccc{QP_representation}, which defines the
number types and iterators used by the QP solver.

@macro <Poly_dist_d declarations> += @begin
    
    template < class ET_, class NT_, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_poly_dist_d;
@end

@macro <Poly_dist_d QP representation> = @begin
    template < class ET_, class NT_, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_poly_dist_d {
        typedef  ET_                    ET;
        typedef  NT_                    NT;

        @<Poly_dist_d QP representation: iterator types>

        typedef  CGAL::Tag_false        Is_lp;
    };
@end

The matrix $A$ has only two rows where each column contains a~$0$ and
a~$1$. We store matrix $A$ as a vector of vectors in \ccc{a_matrix}.

@macro <Poly_dist_d private types> += @begin

    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;
@end

@macro <Poly_dist_d data members> += @begin

    NT_matrix                a_matrix;  // matrix `A' of QP
@end

@macro <Poly_dist_d declarations> += @begin

    template < class NT >
    struct QP_rep_row_of_a {
        typedef  std::vector<NT>                         argument_type;
        typedef  typename argument_type::const_iterator  result_type;
      
        result_type
        operator ( ) ( const argument_type& v) const { return v.begin(); }
    };
@end
    
@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
    #  include <CGAL/_QP_solver/Join_random_access_iterator.h>
    #endif
@end


@macro <Poly_dist_d QP representation: iterator types> += @begin
    typedef  std::vector< std::vector<NT> >
                                        NT_matrix;

    typedef  CGAL::Join_random_access_iterator_1<
                 typename NT_matrix::const_iterator,
                 QP_rep_row_of_a<NT> >  A_iterator;
@end

The vector $b$ has exactly two $1$-entries, while vector $c$ has only
$0$-entries. We use the class template \ccc{Const_value_iterator<T>} to
represent $b$ and $c$.

@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_CONST_VALUE_ITERATOR_H
    #  include <CGAL/_QP_solver/Const_value_iterator.h>
    #endif
@end

@macro <Poly_dist_d QP representation: iterator types> += @begin
    typedef  CGAL::Const_value_iterator<NT>
                                        B_iterator;
    typedef  CGAL::Const_value_iterator<NT>
                                        C_iterator;
@end

Because of its size ($n\!\times\!n$), the matrix $D$ is represented
implicitly. By~(\ref{eq:PD_as_QP}) we have $D = C^T C$, i.e.~$D_{i,j}$ is
the inner product of the $i$-th and $j$-th column of $C$. Row~$i$ of $D$ is
determined by $p_i$ or $-q_{i-|P|}$, respectively, and iterators to both
point sets. Then the entry in column~$j$ is the inner product of the stored
point and $p_j$ or $-q_{j-|P|}$, respectively.

@macro <Poly_dist_d declarations> += @begin
    
    template < class Point, class Point_iterator >
    struct QP_rep_signed_point_iterator;

    template < class NT, class Point,
               class Access_coord, class Access_dim >
    class QP_rep_signed_inner_product;

    template < class NT, class Point, class Point_iterator,
               class Access_coord, class Access_dim >
    struct QP_rep_row_of_d;
@end

@macro <Poly_dist_d signed-point iterator> = @begin
    template < class Point, class PointIterator >
    struct QP_rep_signed_point_iterator {
      public:
        typedef  std::pair<Point,CGAL::Sign>      value_type;
        typedef  ptrdiff_t                        difference_type;
        typedef  value_type*                      pointer;
        typedef  value_type&                      reference;
        typedef  std::random_access_iterator_tag  iterator_category;

        typedef  QP_rep_signed_point_iterator<Point,PointIterator>  Self;
        typedef  value_type                                         Val;
        typedef  difference_type                                    Dist;
        typedef  pointer                                            Ptr;

        // forward operations
        QP_rep_signed_point_iterator(
            const PointIterator& it_p = PointIterator(), Dist n_p = 0,
            const PointIterator& it_q = PointIterator())
            : p_it( it_p), q_it( it_q), n( n_p), curr( 0) { }

        bool   operator == ( const Self& it) const { return (curr == it.curr);}
        bool   operator != ( const Self& it) const { return (curr != it.curr);}

        Val    operator *  ( ) const
            { return ( curr < n) ? std::make_pair( *p_it, CGAL::POSITIVE)
                                 : std__make_pair( *q_it, CGAL::NEGATIVE); }

        Self&  operator ++ (    )
                   { if ( ++curr <= n) ++p_it; else ++q_it; return *this; }
        Self   operator ++ ( int)
                   { Self tmp = *this; operator++(); return tmp; }

        // bidirectional operations
        Self&  operator -- (    )
                   { if ( --curr <  n) --p_it; else --q_it; return *this; }
        Self   operator -- ( int)
                   { Self tmp = *this; operator--(); return tmp; }

        // random access operations
        Self&  operator += ( Dist i)
                   {
                       if ( curr+i <= n) {
                           curr += i;
                           p_it += i;
                       } else {
                           if ( curr < n) p_it += n-curr;
                           curr += i;
                           q_it += curr-n;
                       }
                       return *this;
                   }

        Self&  operator -= ( Dist i)
                   {
                       if ( curr-i < n) {
                           if ( curr > n) q_it -= curr-n;
                           curr -= i;
                           p_it  -= n-curr;
                       } else {
                           curr -= i;
                           q_it  -= i;
                       }
                       return *this;
                   }

        Self   operator +  ( Dist i) const { Self tmp = *this; return tmp+=i; }
        Self   operator -  ( Dist i) const { Self tmp = *this; return tmp-=i; }

        Dist   operator -  ( const Self& it) const { return curr - it.curr; }

        Val    operator [] ( int i) const
            { return ( curr+i < n)
                  ? std::make_pair( p_it[ i  ], CGAL::POSITIVE)
                  : std::make_pair( q_it[ i-n], CGAL::NEGATIVE); }

        bool   operator <  ( const Self&) const { return ( curr <  it.curr); }
        bool   operator >  ( const Self&) const { return ( curr >  it.curr); }
        bool   operator <= ( const Self&) const { return ( curr <= it.curr); }
        bool   operator >= ( const Self&) const { return ( curr >= it.curr); }

      private:
        PointIterator  p_it;
        PointIterator  q_it;
        Dist           n;
        Dist           curr;
    };
@end

@macro <Poly_dist_d signed-inner-product function class> = @begin
    template < class NT, class Point,
               class Access_coord, class Access_dim >
    class QP_rep_signed_inner_product {
        Point         p_i;
        CGAL::Sign    s_i;
        Access_coord  da_coord;
        Access_dim    da_dim;
      public:
        typedef  std::pair<Point,CGAL::Sign>  argument_type;
        typedef  NT                           result_type;
        
        QP_rep_signed_inner_product( ) { }
        QP_rep_signed_inner_product( const argument_type& p_signed,
                                     const Access_coord&  ac,
                                     const Access_dim&    ad)
            : p_i( p_signed.first), s_i( p_signed.second),
              da_coord( ac), da_dim( ad) { }

        NT  operator( ) ( const argument_type& p_signed) const
            { NT ip = std::inner_product( da_coord( p_i),
                                          da_coord( p_i)+da_dim( p_i),
                                          da_coord( p_signed.first), NT( 0),
                                          std::plus<NT>(),
                                          std::multiplies<NT>());
              return ( s_i*p_signed.second == CGAL::POSITIVE) ? ip : -ip; }
    };
@end

@macro <Poly_dist_d row-of-D function class> = @begin
    template < class NT, class Point, class Signed_point_iterator,
               class Access_coord, class Access_dim >
    class QP_rep_row_of_d {
        Signed_point_iterator  signed_pts_it;
        Access_coord           da_coord;
        Access_dim             da_dim;
    public:
        typedef  CGAL::QP_rep_signed_inner_product<
                     NT, Point, Access_coord, Access_dim >
                                        Signed_inner_product;
        typedef  CGAL::Join_random_access_iterator_1<
                     Signed_point_iterator, Signed_inner_product >
                                        Row_of_d;

        typedef  std::pair<Point,CGAL::Sign>
                                        argument_type;
        typedef  Row_of_d               result_type;
                                
        QP_rep_row_of_d( ) { }
        QP_rep_row_of_d( const Signed_point_iterator& it,
                         const Access_coord&          ac,
                         const Access_dim&            ad)
            : signed_pts_it( it), da_coord( ac), da_dim( ad) { }

        Row_of_d  operator( ) ( const argument_type& p_signed) const
        { return Row_of_d( signed_pts_it,
                           Signed_inner_product( p_signed, da_coord, da_dim));}
    };
@end

@macro <Poly_dist_d QP representation: iterator types> += @begin
    
    typedef  CGAL::QP_rep_signed_point_iterator< Point, Point_iterator>
                                        Signed_point_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 Signed_point_iterator,
                 QP_rep_row_of_d< NT, Point, Signed_point_iterator,
                                  Access_coord, Access_dim > >
                                        D_iterator;
@end

Now we are able to define the fully specialized type of the QP solver.

@macro <Poly_dist_d Solver type> = @begin
    // QP solver
    typedef  CGAL::QP_rep_poly_dist_d<
                 ET, NT, Point, typename std::vector<Point>::const_iterator,
                 Access_coordinates_begin_d, Access_dimension_d >
                                        QP_rep;
    typedef  CGAL::QP_solver< QP_rep >  Solver;
    typedef  typename Solver::Pricing_strategy
                                        Pricing_strategy;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Computing the Distance of the Polytopes}

We set up the quadratic program, solve it, and compute the two points
realizing the distance.

@macro <Poly_dist_d private member functions> += @begin
    
    // compute (squared) distance
    void
    compute_distance( )
    {
        // clear support points
        p_support_indices.erase( p_support_indices.begin(),
                                 p_support_indices.end());
        q_support_indices.erase( q_support_indices.begin(),
                                 q_support_indices.end());
        if ( ( p_points.size() == 0) || ( q_points.size() == 0)) return;
        
        // set up and solve QP
        @<Poly_dist_d compute_distance: set up and solve QP>

        // compute support and realizing points
        @<Poly_dist_d compute_distance: compute support and realizing points>
    }
@end

@macro <Poly_dist_d compute_distance: set up and solve QP> = @begin
    int i, j;
    NT  nt_0 = 0, nt_1 = 1;

    // matrix A
    a_matrix.erase( a_matrix.begin(), a_matrix.end());
    a_matrix.insert( a_matrix.end(),
                     number_of_points(), NT_vector( 2, nt_0));
    for ( j = 0; j < number_of_points_p(); ++j) a_matrix[ j][ 0] = nt_1;
    for (      ; j < number_of_points  (); ++j) a_matrix[ j][ 1] = nt_1;

    // set-up
    typedef  QP_rep_signed_point_iterator< Point, Point_iterator >
                                        Signed_point_iterator;
    Signed_point_iterator
        signed_pts_it(p_points.begin(), p_points.size(), q_points.begin());
    QP_rep_row_of_d< NT, Point, Signed_point_iterator,
                     Access_coordinates_begin_d, Access_dimension_d >
        row_of_d( signed_pts_it,
                  tco.access_coordinates_begin_d_object(),
                  tco.access_dimension_d_object());
        
    typedef  typename QP_rep::A_iterator A_it;
    typedef  typename QP_rep::B_iterator B_it;
    typedef  typename QP_rep::C_iterator C_it;
    typedef  typename QP_rep::D_iterator D_it;
    
    solver.set( number_of_points(), 2, d+2,
                A_it( a_matrix.begin()), B_it( 1), C_it( 0),
                D_it( signed_pts_it, row_of_d));

    // solve
    solver.init();
    solver.solve();
@end

@macro <Poly_dist_d compute_distance: compute support and realizing points> = @begin
    ET  et_0 = 0;
    int r    = number_of_points_p();
    p_coords.resize( ambient_dimension()+1);
    q_coords.resize( ambient_dimension()+1);
    std::fill( p_coords.begin(), p_coords.end(), et_0);
    std::fill( q_coords.begin(), q_coords.end(), et_0);
    for ( i = 0; i < solver.number_of_basic_variables(); ++i) {
        ET  value = solver.basic_variables_numerator_begin()[ i];
        int index = solver.basic_variables_index_begin()[ i];
        if ( index < r) {
            for ( int j = 0; j < d; ++j) {
                p_coords[ j]
                    += value * tco.access_coordinates_begin_d_object()(
                                                   p_points[ index  ])[ j];
            }
            p_support_indices.push_back( index);
        } else {
            for ( int j = 0; j < d; ++j) {
                q_coords[ j]
                    += value * tco.access_coordinates_begin_d_object()(
                                                   q_points[ index-r])[ j];
            }
            q_support_indices.push_back( index-r);
        }
    }
    p_coords[ d] = q_coords[ d] = solver.variables_common_denominator();
@end

@! ----------------------------------------------------------------------------
\subsubsection{Choosing the Pricing Strategy}

@macro <Poly_dist_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_PARTIAL_EXACT_PRICING_H
    #  include <CGAL/_QP_solver/Partial_exact_pricing.h>
    #endif
    #ifndef CGAL_PARTIAL_FILTERED_PRICING_H
    #  include <CGAL/_QP_solver/Partial_filtered_pricing.h>
    #endif
@end

@macro <Poly_dist_d data members> += @begin
    
    typename Solver::Pricing_strategy*  // pricing strategy
                             strategyP; // of the QP solver
@end

@macro <Poly_dist_d QP-solver set-up> many = @begin
    set_pricing_strategy( NT());
@end

@macro <Poly_dist_d private member functions> += @begin
    
    template < class NT >
    void  set_pricing_strategy( NT)
        { strategyP = new CGAL::Partial_filtered_pricing<QP_rep>;
          solver.set_pricing_strategy( *strategyP); }

    #ifndef _MSC_VER
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

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

The function \ccc{test_Polytope_distance_d}, invoked with a set of points
and a traits class model, calls each function of \ccc{Polytope_distance_d}
at least once to ensure code coverage. If \ccc{verbose} is set to $-1$, the
function is ``silent'', otherwise some diagnosing output is written to the
standard error stream.

@macro <Poly_dist_d test function> = @begin
    #define COVER(text,code) \
                verr0.out().width( 26); verr0 << text << "..." << flush; \
                verrX.out().width(  0); verrX << "==> " << text << endl \
                  << "----------------------------------------" << endl; \
                { code } verr0 << "ok."; verr << endl;
    
    template < class ForwardIterator, class Traits >
    void
    test_Polytope_distance_d( ForwardIterator p_first, ForwardIterator p_last,
                              ForwardIterator q_first, ForwardIterator q_last,
                              const Traits& traits, int verbose)
    {
        using namespace std;
        
        typedef  CGAL::Polytope_distance_d< Traits >  Poly_dist;
        typedef  typename Traits::Point_d             Point;

        CGAL::Verbose_ostream verr ( verbose >= 0);
        CGAL::Verbose_ostream verr0( verbose == 0);
        CGAL::Verbose_ostream verrX( verbose >  0);
        CGAL::set_pretty_mode( verr.out());

        bool  is_valid_verbose = ( verbose > 0);

        // constructors
        COVER( "default constructor",
            Poly_dist  pd( traits, verbose, verr.out());
            assert( pd.is_valid( is_valid_verbose));
            assert( ! pd.is_finite());
        )

        COVER( "point set constructor",
            Poly_dist  pd( p_first, p_last, q_first, q_last,
                           traits, verbose, verr.out());
            assert( pd.is_valid( is_valid_verbose));
        )

        Poly_dist  poly_dist( p_first, p_last, q_first, q_last);
        COVER( "ambient dimension",
            Poly_dist  pd;
            assert( pd.ambient_dimension() == -1);
            verrX << poly_dist.ambient_dimension() << endl;
        )

        COVER( "(number of) points",
            verrX << poly_dist.number_of_points() << endl;
            verrX << endl << poly_dist.number_of_points_p() << endl;
            typename Poly_dist::Point_iterator
                point_p_it = poly_dist.points_p_begin();
            for ( ; point_p_it != poly_dist.points_p_end(); ++point_p_it) {
                verrX << *point_p_it << endl;
            }
            verrX << endl << poly_dist.number_of_points_q() << endl;
            typename Poly_dist::Point_iterator
                point_q_it = poly_dist.points_q_begin();
            for ( ; point_q_it != poly_dist.points_q_end(); ++point_q_it) {
                verrX << *point_q_it << endl;
            }
            assert( ( poly_dist.number_of_points_p()
                      + poly_dist.number_of_points_q())
                    == poly_dist.number_of_points());
            assert( ( poly_dist.points_p_end() - poly_dist.points_p_begin())
                    == poly_dist.number_of_points_p());
            assert( ( poly_dist.points_q_end() - poly_dist.points_q_begin())
                    == poly_dist.number_of_points_q());
        )

        COVER( "(number of) support points",
            verrX << poly_dist.number_of_support_points() << endl;
            verrX << endl << poly_dist.number_of_support_points_p() << endl;
            typename Poly_dist::Support_point_iterator
                point_p_it = poly_dist.support_points_p_begin();
            for ( ; point_p_it != poly_dist.support_points_p_end();
                  ++point_p_it) {
                verrX << *point_p_it << endl;
            }
            verrX << endl << poly_dist.number_of_support_points_q() << endl;
            typename Poly_dist::Support_point_iterator
                point_q_it = poly_dist.support_points_q_begin();
            for ( ; point_q_it != poly_dist.support_points_q_end();
                  ++point_q_it) {
                verrX << *point_q_it << endl;
            }
            assert( ( poly_dist.number_of_support_points_p()
                      + poly_dist.number_of_support_points_q())
                    == poly_dist.number_of_support_points());
            assert( ( poly_dist.support_points_p_end()
                      - poly_dist.support_points_p_begin())
                    == poly_dist.number_of_support_points_p());
            assert( ( poly_dist.support_points_q_end()
                      - poly_dist.support_points_q_begin())
                    == poly_dist.number_of_support_points_q());
        )

        COVER( "realizing points",
            typename Poly_dist::Coordinate_iterator  coord_it;
            verrX << "p:";
            for ( coord_it  = poly_dist.realizing_point_p_coordinates_begin();
                  coord_it != poly_dist.realizing_point_p_coordinates_end();
                  ++coord_it) {
                verrX << ' ' << *coord_it;
            }
            verrX << endl;
            verrX << "q:";
            for ( coord_it  = poly_dist.realizing_point_q_coordinates_begin();
                  coord_it != poly_dist.realizing_point_q_coordinates_end();
                  ++coord_it) {
                verrX << ' ' << *coord_it;
            }
            verrX << endl;
        )
 
        COVER( "squared distance",
            verrX << poly_dist.squared_distance_numerator()   << " / "
                  << poly_dist.squared_distance_denominator() << endl;
        )

        COVER( "clear",
            poly_dist.clear();
            verrX << "poly_dist is" << ( poly_dist.is_finite() ? "" : " not")
                  << " finite." << endl;
            assert( ! poly_dist.is_finite());
        )

        COVER( "insert (single point)",
            poly_dist.insert_p( *p_first);
            poly_dist.insert_q( *q_first);
            assert( poly_dist.is_valid( is_valid_verbose));
        )

        COVER( "insert (point set)",
            poly_dist.insert( p_first, p_last, q_first, q_last);
            assert( poly_dist.is_valid( is_valid_verbose));
            poly_dist.clear();
            poly_dist.insert_p( q_first, q_last);
            poly_dist.insert_q( p_first, p_last);
            assert( poly_dist.is_valid( is_valid_verbose));
        )

        COVER( "traits class access",
            poly_dist.traits();
        )

        COVER( "I/O",
            verrX << poly_dist;
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
$d$-dimensional \cgal~kernel. All three traits class models are used twice,
firstly with one exact number type (the ``default'' use) and secondly with
three different number types (the ``advanced'' use). Since the current
implementation of the underlying linear programming solver can only handle
input points with Cartesian representation, we use \cgal's Cartesian kernel
for testing.

Some of the following macros are parameterized with the dimension,
e.g.~with $2$, $3$, or $d$.

@macro <Poly_dist_d test: includes and typedefs>(1) many += @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Polytope_distance_d.h>
    #include <CGAL/Optimisation_d_traits_@1.h>
@end

We use the number type \ccc{leda_integer} from \leda{} for the first
variant.

@macro <Poly_dist_d test: includes and typedefs> += @begin

    // test variant 1 (needs LEDA)
    #ifdef CGAL_USE_LEDA
    # include <CGAL/leda_integer.h>
      typedef  CGAL::Cartesian<leda_integer>       K_1;
      typedef  CGAL::Optimisation_d_traits_@1<K_1>  Traits_1;
    # define TEST_VARIANT_1 \
        "Optimisation_d_traits_@1< Cartesian<leda_integer> >"
    #endif
@end

The second variant uses points with \ccc{int} coordinates. The exact number
type used by the underlying quadratic programming solver is
\ccc{GMP::Double}, i.e.~an arbitrary precise floating-point type based on
\textsc{Gmp}'s integers. To speed up the pricing, we use \ccc{double}
arithmetic.

@macro <Poly_dist_d test: includes and typedefs> += @begin

    // test variant 2 (needs GMP)
    #ifdef CGAL_USE_GMP
    # include <CGAL/_QP_solver/Double.h>
      typedef  CGAL::Cartesian< int >                                 K_2;
      typedef  CGAL::Optimisation_d_traits_@1<K_2,GMP::Double,double>  Traits_2;
    # define TEST_VARIANT_2 \
        "Optimisation_d_traits_@1< Cartesian<int>, GMP::Double, double >"
    #endif
@end

The test sets consist of $50$ points with $20$-bit random integer
coordinates uniformly distributed in a $d$-cube.

@macro <Poly_dist_d test: includes and typedefs> += @begin

    #include <CGAL/Random.h>
    #include <vector>
@end

@macro <Poly_dist_d test: generate point set (2D)>(1) = @begin
    std::vector<K_@1::Point_2>  p_points_@1, q_points_@1;
    p_points_@1.reserve( 50);
    q_points_@1.reserve( 50);
    {
        int  i;
        for ( i = 0; i < 50; ++i) {
            p_points_@1.push_back( K_@1::Point_2(
                CGAL::default_random( 0x100000),
                CGAL::default_random( 0x100000)));
        }
        for ( i = 0; i < 50; ++i) {
            q_points_@1.push_back( K_@1::Point_2(
                -CGAL::default_random( 0x100000),
                -CGAL::default_random( 0x100000)));
        }
    }
@end

@macro <Poly_dist_d test: generate point set (3D)>(1) = @begin
    std::vector<K_@1::Point_3>  p_points_@1, q_points_@1;
    p_points_@1.reserve( 50);
    q_points_@1.reserve( 50);
    {
        int  i;
        for ( i = 0; i < 50; ++i) {
            p_points_@1.push_back( K_@1::Point_3(
                CGAL::default_random( 0x100000),
                CGAL::default_random( 0x100000),
                CGAL::default_random( 0x100000)));
        }
        for ( i = 0; i < 50; ++i) {
            q_points_@1.push_back( K_@1::Point_3(
                -CGAL::default_random( 0x100000),
                -CGAL::default_random( 0x100000),
                -CGAL::default_random( 0x100000)));
        }
    }
@end

The traits class model with $d$-dimensional points is tested with $d = 5$
(variant 1) and $d = 10$ (variant 2).

@macro <Poly_dist_d test: generate point set (dD)>(1) = @begin
    std::vector<K_@1::Point_d>  p_points_@1, q_points_@1;
    p_points_@1.reserve( 50);
    q_points_@1.reserve( 50);
    {
        int d = 5*@1;
        std::vector<int>  coords( d);
        int  i, j;
        for ( i = 0; i < 50; ++i) {
            for ( j = 0; j < d; ++j)
                coords[ j] = CGAL::default_random( 0x100000);
            p_points_@1.push_back( K_@1::Point_d( d, coords.begin(),
                                                     coords.end()));
        }
        for ( i = 0; i < 50; ++i) {
            for ( j = 0; j < d; ++j)
                coords[ j] = -CGAL::default_random( 0x100000);
            q_points_@1.push_back( K_@1::Point_d( d, coords.begin(),
                                                     coords.end()));
        }
    }
@end

Finally we call the test function (described in the last section).

@macro <Poly_dist_d test: includes and typedefs> += @begin

    #include "test_Polytope_distance_d.h"
@end

@macro <Poly_dist_d test: call test function>(1) many = @begin
    CGAL::test_Polytope_distance_d( p_points_@1.begin(), p_points_@1.end(),
                                    q_points_@1.begin(), q_points_@1.end(),
                                    Traits_@1(), verbose);
@end

Each of the two test variants is compiled and executed only if the
respective number type is available.

@macro <Poly_dist_d test: test variant output>(1) many = @begin
    verr << endl
         << "==================================="
         << "===================================" << endl
         << "Testing `Polytope_distance_d' with traits class model" <<endl
         << "==> " << TEST_VARIANT_@1 << endl
         << "==================================="
         << "===================================" << endl
         << endl;
@end

@macro <Poly_dist_d test: test variant (2D)>(1) many = @begin
    // test variant @1
    // --------------
    #ifdef TEST_VARIANT_@1

        @<Poly_dist_d test: test variant output>(@1)

        // generate point set
        @<Poly_dist_d test: generate point set (2D)>(@1)

        // call test function
        @<Poly_dist_d test: call test function>(@1)

    #endif
@end

@macro <Poly_dist_d test: test variant (3D)>(1) many = @begin
    // test variant @1
    // --------------
    #ifdef TEST_VARIANT_@1

        @<Poly_dist_d test: test variant output>(@1)

        // generate point set
        @<Poly_dist_d test: generate point set (3D)>(@1)

        // call test function
        @<Poly_dist_d test: call test function>(@1)

    #endif
@end

@macro <Poly_dist_d test: test variant (dD)>(1) many = @begin
    // test variant @1
    // --------------
    #ifdef TEST_VARIANT_@1

        @<Poly_dist_d test: test variant output>(@1)

        // generate point set
        @<Poly_dist_d test: generate point set (dD)>(@1)

        // call test function
        @<Poly_dist_d test: call test function>(@1)

    #endif
@end

The complete bodies of the test programs look as follows. Verbose output
can be enabled by giving a number between 0 and 3 at the command line.

@macro <Poly_dist_d test: command line argument> many = @begin
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr ( verbose >= 0); verr  << "";
@end

@macro <Poly_dist_d test: main (2D)> = @begin
    using namespace std;
    
    @<Poly_dist_d test: command line argument>
    
    @<Poly_dist_d test: test variant (2D)>(1)

    @<Poly_dist_d test: test variant (2D)>(2)

    return 0;
@end

@macro <Poly_dist_d test: main (3D)> = @begin
    using namespace std;
    
    @<Poly_dist_d test: command line argument>
    
    @<Poly_dist_d test: test variant (3D)>(1)

    @<Poly_dist_d test: test variant (3D)>(2)

    return 0;
@end

@macro <Poly_dist_d test: main (dD)> = @begin
    using namespace std;
    
    @<Poly_dist_d test: command line argument>
    
    @<Poly_dist_d test: test variant (dD)>(1)

    @<Poly_dist_d test: test variant (dD)>(2)

    return 0;
@end


@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Polytope_distance_d.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Polytope\_distance\_d.h}

@file <include/CGAL/Polytope_distance_d.h> = @begin
    @<file header>(
        "include/CGAL/Polytope_distance_d.h",
        "Distance of convex polytopes in arbitrary dimension")

    #ifndef CGAL_POLYTOPE_DISTANCE_D_H
    #define CGAL_POLYTOPE_DISTANCE_D_H

    // includes
    // --------
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<Poly_dist_d CGAL includes>
    @<Poly_dist_d CGAL/QP_solver includes>
    @<Poly_dist_d standard includes>

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Poly_dist_d declarations>
    
    // Class interfaces
    // ================
    @<Poly_dist_d interface>

    @<Poly_dist_d signed-point iterator>

    @<Poly_dist_d signed-inner-product function class>

    @<Poly_dist_d row-of-D function class>

    @<Poly_dist_d QP representation>

    // Function declarations
    // =====================
    @<Poly_dist_d I/O operators declaration>

    @<dividing line>
    
    // Class implementation
    // ====================

    @<Poly_dist_d validity check>
    
    @<Poly_dist_d I/O operators>

    @<namespace end>("CGAL")
    
    #endif // CGAL_POLYTOPE_DISTANCE_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Polytope_distance_d.h
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d.h}

@file <test/Polytope_distance_d/test_Polytope_distance_d.h> = @begin
    @<file header>(
        "test/Polytope_distance_d/test_Polytope_distance_d.h",
        "test function for smallest enclosing sphere")

    #ifndef CGAL_TEST_POLYTOPE_DISTANCE_D_H
    #define CGAL_TEST_POLYTOPE_DISTANCE_D_H

    // includes
    #ifndef CGAL_IO_VERBOSE_OSTREAM_H
    #  include <CGAL/IO/Verbose_ostream.h>
    #endif
    #include <cassert>
    
    @<namespace begin>("CGAL")

    @<Poly_dist_d test function>
    
    @<namespace end>("CGAL")

    #endif // CGAL_TEST_POLYTOPE_DISTANCE_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Polytope_distance_d_2.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_2.C}

@file <test/Polytope_distance_d/test_Polytope_distance_d_2.C> = @begin
    @<file header>(
        "test/Polytope_distance_d/test_Polytope_distance_d_2.C",
        "test program for polytope distance (2D traits class)")

    // includes and typedefs
    // ---------------------
    @<Poly_dist_d test: includes and typedefs>(2)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Poly_dist_d test: main (2D)>
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Polytope_distance_d_3.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_3.C}

@file <test/Polytope_distance_d/test_Polytope_distance_d_3.C> = @begin
    @<file header>(
        "test/Polytope_distance_d/test_Polytope_distance_d_3.C",
        "test program for polytope distance (3D traits class)")

    // includes and typedefs
    // ---------------------
    @<Poly_dist_d test: includes and typedefs>(3)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Poly_dist_d test: main (3D)>
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Polytope_distance_d_d.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_sphere\_d/test\_Min\_sphere\_d\_d.C}

@file <test/Polytope_distance_d/test_Polytope_distance_d_d.C> = @begin
    @<file header>(
        "test/Polytope_distance_d/test_Polytope_distance_d_d.C",
        "test program for polytope distance (dD traits class")

    // includes and typedefs
    // ---------------------
    @<Poly_dist_d test: includes and typedefs>(d)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Poly_dist_d test: main (dD)>
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
        "Polytope_distance_d",
        "Geometric Optimisation",
        "Polytope_distance_d",
        "$Id$","$Date$",
        "Sven Schönherr <sven@@inf.ethz.ch>",
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
