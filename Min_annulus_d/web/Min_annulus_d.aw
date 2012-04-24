@! ============================================================================
@! The CGAL Library
@! Implementation: Smallest Enclosing Annulus in Arbitrary Dimension
@! ----------------------------------------------------------------------------
@! file  : web/Min_annulus_d.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Min_annulus_d WIP $
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
\newcommand{\cgalMinAnnulusLayout}{%
  \ccSetThreeColumns{}{min_annulus.center()\,}{returns
    \ccGlobalScope\ccc{ON_BOUNDED_SIDE},
    \ccGlobalScope\ccc{ON_BOUNDARY},}
  \ccPropagateThreeToTwoColumns}

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
@t title titlefont centre "Smallest Enclosing Annulus"
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
    \verb+$CGAL_Package: Min_annulus_d WIP+\cgalWIP\verb+$+ \\
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
  We provide an implementation for computing the smallest enclosing annulus
  (region between two concentric spheres that minimizes the difference
  between the squared radii) of a finite point set in arbitrary dimension.
  The problem is formulated as a linear program and a dedicated
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

We consider the problem of finding the annulus (region between two
concentric spheres) enclosing a finite set of points in $d$-dimensional
Euclidean space $\E_d$ that minimizes the difference between its squared
radii. It can be formulated as an optimization problem with linear
constraints and a linear objective function~\cite{gs-eegqp-00}.

@! ----------------------------------------------------------------------------
@! Smallest Enclosing Annulus as a Linear Programming Problem
@! ----------------------------------------------------------------------------

\subsection{Smallest Enclosing Annulus as a Linear Programming Problem}

Let $r$ and $R$ be the inner and outer radii of the annulus and $c$ its center.
If the point set is given as $P = \{p_1,\dots,p_n\}$, we want to minimize
%
\begin{equation} \label{eq:objective_function}
  R^2-r^2
\end{equation}
%
subject to the constraints
%
\begin{equation*}
  r \leq ||p_i-c|| \leq R, \quad \forall i \in \{1,\dots,n\}.
\end{equation*}
%
If $p_i = (p^i_1,\dots,p^i_d)$ and $c = (c_1,\dots,c_d)$ we can
equivalently write this as
%
\begin{equation*}
  r^2 \leq (p^i_1-c_1)^2 + \dots + (p^i_d-c_d)^2 \leq R^2,
  \quad \forall i \in \{1,\dots,n\}
\end{equation*}
%
or
%
\begin{equation} \label{eq:constraints}
  \begin{array}{@@{}r@@{\,}}
    (r^2 - c_1^2 - \dots - c_d^2) + 2p^i_1c_1 + \dots + 2p^i_dc_d
    \leq (p^i_1)^2 + \dots + (p^i_d)^2 \\[1ex]
    (R^2 - c_1^2 - \dots - c_d^2) + 2p^i_1c_1 + \dots + 2p^i_dc_d
    \geq (p^i_1)^2 + \dots + (p^i_d)^2 \\
  \end{array},
  \quad \forall i \in \{1,\dots,n\}.
\end{equation}
%
Defining
%
\begin{alignat}{2}
  \alpha &:= \mbox{} & r^2 - c_1^2 - \dots - c_d^2, & \label{eq:alpha} \\
  \beta  &:= \mbox{} & R^2 - c_1^2 - \dots - c_d^2, & \label{eq:beta}
\end{alignat}
%
then (\ref{eq:objective_function}) equals $\beta-\alpha$, while
(\ref{eq:constraints}) becomes
%
\begin{equation*}
  \begin{array}{@@{}r@@{\,}}
    \alpha + 2p^i_1c_1 + \dots + 2p^i_dc_d
    \leq (p^i_1)^2 + \dots + (p^i_d)^2 \\[1ex]
    \beta  + 2p^i_1c_1 + \dots + 2p^i_dc_d
    \geq (p^i_1)^2 + \dots + (p^i_d)^2 \\
  \end{array},
  \quad \forall i \in \{1,\dots,n\},
\end{equation*}
%
so we get a linear program (LP) with $2n$ constraints and $d\!+\!2$
variables $\alpha,\beta,c_1,\dots,c_d$. From a solution to this problem the
desired radii $r$ and $R$ can obviously be reconstructed via
(\ref{eq:alpha}) and (\ref{eq:beta}).

For efficient solvability by the dedicated solver described
in~\cite{s-qpego1-00}, the large number of constraints is not suitable.
Instead, we would need a large number of variables and only few
constraints. Fortunately, we may consider the dual version of this problem
which has exactly this feature.

Consider a general LP in the form
%
\begin{equation*}
  \begin{array}{ll}
    \text{minimize} & c^T x \\[0.5ex]
    \text{subject to} & Ax \leq b,
  \end{array}
\end{equation*}
%
then its dual is
%
\begin{equation*}
  \begin{array}{ll}
    \text{minimize} & -b^T y \\[0.5ex]
    \text{subject to} & A^Ty = c, \\
    & y \geq 0.
  \end{array}
\end{equation*}

Thus, the dual of the LP we want to solve has $2n$ variables, namely
$\lambda = (\lambda_1,\dots,\lambda_d)$ and $\mu = (\mu_1,\dots,\mu_d)$,
and $d\!+\!2$~constraints, given as
%
\begin{equation*}
  \begin{array}{ll}
    \text{minimize} & \sum_{i=1}^n ((p^i_1)^2 + \dots + (p^i_d)^2)\lambda_i
    - \sum_{i=1}^n ((p^i_1)^2 + \dots + (p^i_d)^2)\mu_i \\[1ex]
    \text{subject to}
    & \sum_{i=1}^n \lambda_i = 1, \\[0.5ex]
    & \sum_{i=1}^n -\mu_i = -1, \\[0.5ex]
    & \sum_{i=1}^n 2p^i_1\lambda_i - \sum_{i=1}^n 2p^i_1\mu_i = 0, \\[0.5ex]
    & \vdots \\[0.5ex]
    & \sum_{i=1}^n 2p^i_d\lambda_i - \sum_{i=1}^n 2p^i_d\mu_i = 0, \\[0.5ex]
    & \lambda,\mu \geq 0.
  \end{array}
\end{equation*}
%
In the following matrix notation we slightly changed the column order by
putting together the two columns corresponding to $p_i$, for every $i$.
%
\begin{equation} \label{eq:MA_as_LP}
  \begin{array}{lll}
    \text{(MA)} & \text{minimize}
    & \left( \begin{array}{@@{}ccccc@@{}}
      +p_1^T p_1 & -p_1^T p_1 & \dots & +p_n^T p_n & -p_n^T p_n
    \end{array} \right) \, y \\[1ex]
  & \text{subject to}
  & \left( \begin{array}{ccccc}
      \\
      +2p_1 & -2p_1 & \dots & +2p_n & -2p_n \\[0.5ex]
      \\
      +1    &  0    & \dots & +1    &  0    \\
       0    & -1    & \dots &  0    & -1
    \end{array} \right)
  y =
  \left( \begin{array}{c}
      0 \\[-1ex] \vdots \\ 0 \\ +1 \\ -1 \end{array} \right), \\
  & & y \geq 0,
  \end{array}
\end{equation}
%
where $y^T = (\lambda_1,\mu_1,\dots,\lambda_n,\mu_n)$. The optimal
values to the primal problem are determined by
%
\begin{equation*}
  (\alpha,\beta,c_1,\dots,c_d) = c_B^T A_B^{-1},
\end{equation*}
%
where $A_B^{-1}$ is the final basis inverse of the dual problem.


@! ============================================================================
@! Reference Pages
@! ============================================================================

\clearpage
\section{Reference Pages} \label{sec:reference_pages}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

@p maximum_input_line_length = 102

@! ----------------------------------------------------------------------------
@! Class: Min_annulus_d
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{CGAL::Min\_annulus\_d\texttt{<}Traits\texttt{>}}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/Min_annulus_d.tex}

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
@! The Class Template CGAL::Min_annulus_d<Traits>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Min\_annulus\_d\texttt{<}Traits\texttt{>}}

The class template \ccc{Min_annulus_d} expects a model of the concept
\ccc{OptimisationDTraits} (see
Section~\ref{ccRef_OptimisationDTraits}.2) as its template argument.

@macro <Min_annulus_d declarations> += @begin
    template < class Traits_ >
    class Min_annulus_d;
@end

The interface consists of the public types and member functions described
in Section~\ref{ccRef_CGAL::Min_annulus_d<Traits>}.1 and of some private
types, private member functions, and data members.

@macro <Min_annulus_d interface> = @begin
    template < class Traits_ >
    class Min_annulus_d {
      public:
        // self
        typedef  Traits_                    Traits;
        typedef  Min_annulus_d<Traits>      Self;

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
        @<Min_annulus_d Solver type>
        @<Min_annulus_d QP_solver types>
        @<Min_annulus_d private types>

      public:
        @<Min_annulus_d types>
        
        @<Min_annulus_d member functions>
        
      private:
        @<Min_annulus_d data members>

        @<Min_annulus_d private member functions>
    };
@end

@! ----------------------------------------------------------------------------
\subsubsection{Data Members}

Mainly, we have to store the given input points, the center and the squared
radii of the smallest enclosing annulus, and an instance of the linear
programming solver. Additional variables, that are used in the member
functions described below, are introduced when they appear for the first
time.

We start with the traits class object.

@macro <Min_annulus_d data members> += @begin

    Traits                   tco;       // traits class object
@end

The inputs points are kept in a vector to have random access to them.
Their dimension is stored separately.

@macro <Min_annulus_d standard includes> += @begin
    #include <vector>
@end

@macro <Min_annulus_d private types> += @begin
    // private types
    typedef  std::vector<Point>         Point_vector;
@end

@macro <Min_annulus_d data members> += @begin

    Point_vector             points;    // input points
    int                      d;         // dimension of input points
@end

The center and the squared radii of the smallest enclosing annulus are
stored with rational representation, i.e.~numerators and denominators are
kept separately. The vector \ccc{center_coords} contains $d+1$ entries,
the numerators of the $d$ coordinates and the common denominator.

@macro <Min_annulus_d private types> += @begin
    typedef  std::vector<ET>            ET_vector;
@end

@macro <Min_annulus_d data members> += @begin
    
    ET_vector                center_coords;     // center of small.encl.annulus

    ET                       sqr_i_rad_numer;   // squared inner radius of
    ET                       sqr_o_rad_numer;   // ---"--- outer ----"----
    ET                       sqr_rad_denom;     // smallest enclosing annulus
@end

We store an instance of the linear and quadratic programming solver
described in~\cite{s-qpego1-00}. The details are given in
Section~\ref{sec:using_qp_solver} below, here it suffice to know that there
is a variable \ccc{solver} of type \ccc{Solver}.

@macro <Min_annulus_d private types: linear programming solver> zero = @begin
    typedef  ...                        Solver;
@end

@macro <Min_annulus_d data members> += @begin

    Solver                   solver;    // linear programming solver
@end

@! ----------------------------------------------------------------------------
\subsubsection{Creation}

Two constructors are provided. If the user wants to get some verbose output
(of the underlying QP solver), he can override the default arguments of
\ccc{verbose} and \ccc{stream}.

@macro <Min_annulus_d standard includes> += @begin
    #include <iostream>
@end

@macro <Min_annulus_d member functions> += @begin
    // creation
    Min_annulus_d( const Traits&  traits  = Traits(),
                   int            verbose = 0,
                   std::ostream&  stream  = std::cout)
      : tco( traits), d( -1), solver( verbose, stream)
        {
            @<Min_annulus_d QP-solver set-up>
        }
@end

The second constructor expects a set of points given via an iterator range.
It calls the \ccc{set} member function described in
Subsection~\ref{sec:modifiers} to store the points and to compute the
smallest enclosing annulus of the given point set.

@macro <Min_annulus_d member functions> += @begin

    template < class InputIterator >
    Min_annulus_d( InputIterator  first,
                   InputIterator  last,
                   const Traits&  traits = Traits(),
                   int            verbose = 0,
                   std::ostream&  stream  = std::cout)
      : tco( traits), solver( verbose, stream)
        {
            @<Min_annulus_d QP-solver set-up>
            set( first, last);
        }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Access}

The following types and member functions give access to the set of points
contained in the smallest enclosing annulus.

@macro <Min_annulus_d types> += @begin
    // public types
    typedef  typename Point_vector::const_iterator
                                        Point_iterator;
@end

@macro <Min_annulus_d member functions> += @begin

    // access to point set
    int  ambient_dimension( ) const { return d; }

    int  number_of_points( ) const { return points.size(); }

    Point_iterator  points_begin( ) const { return points.begin(); }
    Point_iterator  points_end  ( ) const { return points.end  (); }
@end

To access the support points, we exploit the following fact. A point~$p_i$
is a support point, iff one of its corresponding variable $\lambda_i$ or
$\mu_i$ (of the QP solver) is basic. Thus the number of support points is
equal to the number of basic variables, if the smallest enclosing annulus
is not empty.

@macro <Min_annulus_d member functions> += @begin

    // access to support points
    int
    number_of_support_points( ) const
        { return number_of_points() < 2 ? number_of_points()
                                        : solver.number_of_basic_variables(); }
@end

If $i$ is the index of the $k$-th basic variable, then $p_{i/2}$ is the
$k$-th support point. To access a point given its index, we use the
following function class.

@macro <Min_annulus_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_ACCESS_BY_INDEX_H
    #  include <CGAL/_QP_solver/Access_by_index.h>
    #endif
@end

@macro <Min_annulus_d private types> += @begin

    typedef  CGAL::Access_by_index<typename std::vector<Point>::const_iterator>
                                        Point_by_index;
@end

Another function class is used to divide the index by 2.

@macro <Min_annulus_d standard includes> += @begin
    #include <functional>
@end

@macro<Min_annulus_d private types> += @begin

    typedef  std::binder2nd< std::divides<int> >
                                        Divide;
@end

The indices of the basic variables can be accessed with the following
iterator.

@macro <Min_annulus_d QP_solver types> += @begin
    // types from the QP solver
    typedef  typename Solver::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;
@end

Combining the function classes with the index iterator gives the support
point iterator.

@macro <Min_annulus_d CGAL includes> += @begin
    #ifndef CGAL_FUNCTION_OBJECTS_H
    #  include <CGAL/function_objects.h>
    #endif
@end

@macro <Min_annulus_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
    #  include <CGAL/_QP_solver/Join_random_access_iterator.h>
    #endif
@end

@macro <Min_annulus_d types> += @begin

    typedef  CGAL::Join_random_access_iterator_1<
                 Basic_variable_index_iterator,
                 CGAL::Unary_compose_1<Point_by_index,Divide> >
                                        Support_point_iterator;
@end

@macro <Min_annulus_d member functions> += @begin

    Support_point_iterator
    support_points_begin() const
        { return Support_point_iterator(
                     solver.basic_variables_index_begin(),
                     CGAL::compose1_1(
                         Point_by_index( points.begin()),
                         std::bind2nd( std::divides<int>(), 2)));}

    Support_point_iterator
    support_points_end() const
        { return Support_point_iterator( number_of_points() < 2
                     ? solver.basic_variables_index_begin()
                     : solver.basic_variables_index_end(),
                     CGAL::compose1_1(
                         Point_by_index( points.begin()),
                         std::bind2nd( std::divides<int>(), 2)));}
@end

Before we can access the inner and outer support points, we have to divide
the set of basic variables into two sets corresponding to the inner and
outer support points, respectively. The indices are stored in
\ccc{inner_indices} and \ccc{outer_indices}, while the actual split-up is
done in the private member function \ccc{compute_distance} described below
in Section~\ref{sec:using_qp_solver}.

@macro <Min_annulus_d private types> += @begin

    typedef  std::vector<int>           Index_vector;
@end

@macro <Min_annulus_d data members> += @begin

    Index_vector             inner_indices;
    Index_vector             outer_indices;
@end

@macro <Min_annulus_d member functions> += @begin

    int  number_of_inner_support_points() const { return inner_indices.size();}
    int  number_of_outer_support_points() const { return outer_indices.size();}
@end

@macro <Min_annulus_d types> += @begin

    typedef  typename Index_vector::const_iterator IVCI;
    typedef  CGAL::Join_random_access_iterator_1<
                 IVCI, Point_by_index >
                                        Inner_support_point_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 IVCI, Point_by_index >
                                        Outer_support_point_iterator;
@end

@macro <Min_annulus_d member functions> += @begin

    Inner_support_point_iterator
    inner_support_points_begin() const
        { return Inner_support_point_iterator(
                     inner_indices.begin(),
                     Point_by_index( points.begin())); }

    Inner_support_point_iterator
    inner_support_points_end() const
        { return Inner_support_point_iterator(
                     inner_indices.end(),
                     Point_by_index( points.begin())); }

    Outer_support_point_iterator
    outer_support_points_begin() const
        { return Outer_support_point_iterator(
                     outer_indices.begin(),
                     Point_by_index( points.begin())); }

    Outer_support_point_iterator
    outer_support_points_end() const
        { return Outer_support_point_iterator(
                     outer_indices.end(),
                     Point_by_index( points.begin())); }
@end


The following types and member functions give access to the center and the
squared radii of the smallest enclosing annulus.

@macro <Min_annulus_d types> += @begin

    typedef  typename ET_vector::const_iterator
                                        Coordinate_iterator;
@end

@macro <Min_annulus_d member functions> += @begin

    // access to center (rational representation)
    Coordinate_iterator
    center_coordinates_begin( ) const { return center_coords.begin(); }

    Coordinate_iterator
    center_coordinates_end  ( ) const { return center_coords.end  (); }

    // access to squared radii (rational representation)
    ET  squared_inner_radius_numerator( ) const { return sqr_i_rad_numer; }
    ET  squared_outer_radius_numerator( ) const { return sqr_o_rad_numer; }
    ET  squared_radii_denominator     ( ) const { return sqr_rad_denom; }
@end

For convinience, we also provide member functions for accessing the center
as a single point of type \ccc{Point} and the squared radii as single
numbers of type \ccc{FT}. Both functions only work, if an implicit
conversion from number type \ccc{ET} to number type \ccc{RT} is available,
e.g.~if both types are the same.

@macro <Min_annulus_d member functions> += @begin

    // access to center and squared radii
    // NOTE: an implicit conversion from ET to RT must be available!
    Point
    center( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return tco.construct_point_d_object()( ambient_dimension(),
                                                 center_coordinates_begin(),
                                                 center_coordinates_end()); }

    FT
    squared_inner_radius( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_inner_radius_numerator()) /
                 FT( squared_radii_denominator()); }

    FT
    squared_outer_radius( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_outer_radius_numerator()) /
                 FT( squared_radii_denominator()); }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

We use the private member function \ccc{sqr_dist} to compute the
squared distance of a given point to the center of the smallest
enclosing annulus.

@macro <Min_annulus_d CGAL includes> += @begin
    #ifndef CGAL_IDENTITY_H
    #  include <CGAL/_QP_solver/identity.h>
    #endif
@end

@macro <Min_annulus_d private member functions> += @begin

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

@macro <Min_annulus_d member functions> += @begin

    // predicates
    CGAL::Bounded_side
    bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return CGAL::Bounded_side(
                       CGAL_NTS sign( sqr_d - sqr_i_rad_numer)
                     * CGAL_NTS sign( sqr_o_rad_numer - sqr_d)); }

    bool
    has_on_bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
            is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return ( ( sqr_i_rad_numer < sqr_d) && ( sqr_d < sqr_o_rad_numer)); }

    bool
    has_on_boundary( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return (( sqr_d == sqr_i_rad_numer) || ( sqr_d == sqr_o_rad_numer));}

    bool
    has_on_unbounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return ( ( sqr_d < sqr_i_rad_numer) || ( sqr_o_rad_numer < sqr_d)); }
@end

The smallest enclosing annulus is \emph{empty}, if it contains no points,
and it is \emph{degenerate}, if it has less than two support points.

@macro <Min_annulus_d member functions> += @begin

    bool  is_empty     ( ) const { return number_of_points() == 0; }
    bool  is_degenerate( ) const
        { return ! CGAL_NTS is_positive( sqr_o_rad_numer); }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers} \label{sec:modifiers}

These private member functions are used by the following \ccc{set} and
\ccc{insert} member functions to set and check the dimension of the input
points, respectively.

@macro <Min_annulus_d private member functions> += @begin

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
\ccc{compute_min_annulus} (described in Section~\ref{sec:using_qp_solver})
to compute the smallest enclosing annulus.

@macro <Min_annulus_d member functions> += @begin

    // modifiers
    template < class InputIterator >
    void
    set( InputIterator first, InputIterator last)
        { if ( points.size() > 0) points.erase( points.begin(), points.end());
          std::copy( first, last, std::back_inserter( points));
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension(),
              "Not all points have the same dimension.");
          compute_min_annulus(); }
@end

The \ccc{insert} member functions append the given point(s) to the point
set and recompute the smallest enclosing annulus.

@macro <Min_annulus_d member functions> += @begin

    void
    insert( const Point& p)
        { if ( is_empty()) d = tco.access_dimension_d_object()( p);
          CGAL_optimisation_precondition( 
              tco.access_dimension_d_object()( p) == d);
          points.push_back( p);
          compute_min_annulus(); }

    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
        { CGAL_optimisation_precondition_code( int old_n = points.size());
          points.insert( points.end(), first, last);
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension( old_n),
              "Not all points have the same dimension.");
          compute_min_annulus(); }
@end

The \ccc{clear} member function deletes all points and resets the smallest
enclosing annulus to the empty annulus.

@macro <Min_annulus_d member functions> += @begin

    void
    clear( )
        { points.erase( points.begin(), points.end());
          compute_min_annulus(); }
@end


@! ----------------------------------------------------------------------------
\subsubsection{Validity Check}

A \ccc{Min_annulus_d<Traits>} object can be checked for validity. This
means, it is checked whether (a) the annulus contains all points of its
defining set $P$, and (b) the annulus is the smallest annulus spanned by its
support set $S$ and the support set is minimal, i.e.~no support point is
redundant. The function \ccc{is_valid} is mainly intended for debugging
user supplied traits classes but also for convincing the anxious user that
the traits class implementation is correct. If \ccc{verbose} is \ccc{true},
some messages concerning the performed checks are written to standard error
stream. The second parameter \ccc{level} is not used, we provide it only
for consistency with interfaces of other classes.

@macro <Min_annulus_d member functions> += @begin

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
@end

@macro <Min_annulus_d validity check> = @begin
    // validity check
    template < class Traits_ >
    bool
    Min_annulus_d<Traits_>::
    is_valid( bool verbose, int level) const
    {
        using namespace std;
        
        CGAL::Verbose_ostream verr( verbose);
        verr << "CGAL::Min_annulus_d<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;
        
        // containment check (a)
        // ---------------------
        @<Min_annulus_d validity check: containment>

        // support set check (b)
        // ---------------------
        @<Min_annulus_d validity check: support set>

        verr << "  object is valid!" << endl;
        return( true);
    }
@end

The containment check (a) is easy to perform, just a loop over all
points in $|P|$.

@macro <Min_annulus_d validity check: containment> = @begin
    verr << "  (a) containment check..." << flush;

    Point_iterator  point_it = points_begin();
    for ( ; point_it != points_end(); ++point_it) {
        if ( has_on_unbounded_side( *point_it)) 
            return CGAL::_optimisation_is_valid_fail( verr,
                       "annulus does not contain all points");
    }

    verr << "passed." << endl;
@end

To validate the support set, we check whether all inner and outer support
points lie on the inner and outer boundary of the smallest enclosing
annulus, respectively.

@macro <Min_annulus_d validity check: support set> = @begin
    verr << "  (b) support set check..." << flush;

    // all inner support points on inner boundary?
    Inner_support_point_iterator  i_pt_it = inner_support_points_begin();
    for ( ; i_pt_it != inner_support_points_end(); ++i_pt_it) {
        if ( sqr_dist( *i_pt_it) != sqr_i_rad_numer)
            return CGAL::_optimisation_is_valid_fail( verr,
       "annulus does not have all inner support points on its inner boundary");
    }

    // all outer support points on outer boundary?
    Outer_support_point_iterator  o_pt_it = outer_support_points_begin();
    for ( ; o_pt_it != outer_support_points_end(); ++o_pt_it) {
        if ( sqr_dist( *o_pt_it) != sqr_o_rad_numer)
            return CGAL::_optimisation_is_valid_fail( verr,
       "annulus does not have all outer support points on its outer boundary");
    }
    /*
    // center strictly in convex hull of support points?
    typename Solver::Basic_variable_numerator_iterator
        num_it = solver.basic_variables_numerator_begin();
    for ( ; num_it != solver.basic_variables_numerator_end(); ++num_it) {
        if ( ! (    CGAL_NTS is_positive( *num_it)
                 && *num_it <= solver.variables_common_denominator()))
            return CGAL::_optimisation_is_valid_fail( verr,
              "center does not lie strictly in convex hull of support points");
    }
    */
    
    verr << "passed." << endl;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Miscellaneous}

The member function \ccc{traits} returns a const reference to the
traits class object.

@macro <Min_annulus_d member functions> += @begin

    // traits class access
    const Traits&  traits( ) const { return tco; }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Min_annulus_d I/O operators declaration> = @begin
    // I/O operators
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os, const Min_annulus_d<Traits_>& min_annulus);

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is,       Min_annulus_d<Traits_>& min_annulus);
@end

@macro <Min_annulus_d I/O operators> = @begin
    // output operator
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Min_annulus_d<Traits_>& min_annulus)
    {
        using namespace std;

        typedef  Min_annulus_d<Traits_>::Point  Point;
        typedef  ostream_iterator<Point>       Os_it;
        typedef  typename Traits_::ET          ET;
        typedef  ostream_iterator<ET>          Et_it;

        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Min_annulus_d( |P| = "
               << min_annulus.number_of_points() << ", |S| = "
               << min_annulus.number_of_inner_support_points() << '+'
               << min_annulus.number_of_outer_support_points() << endl;
            os << "  P = {" << endl;
            os << "    ";
            copy( min_annulus.points_begin(), min_annulus.points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S_i = {" << endl;
            os << "    ";
            copy( min_annulus.inner_support_points_begin(),
                  min_annulus.inner_support_points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S_o = {" << endl;
            os << "    ";
            copy( min_annulus.outer_support_points_begin(),
                  min_annulus.outer_support_points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  center = ( ";
            copy( min_annulus.center_coordinates_begin(),
                  min_annulus.center_coordinates_end(),
                  Et_it( os, " "));
            os << ")" << endl;
            os << "  squared inner radius = "
               << min_annulus.squared_inner_radius_numerator() << " / "
               << min_annulus.squared_radii_denominator() << endl;
            os << "  squared outer radius = "
               << min_annulus.squared_outer_radius_numerator() << " / "
               << min_annulus.squared_radii_denominator() << endl;
            break;

          case CGAL::IO::ASCII:
            copy( min_annulus.points_begin(), min_annulus.points_end(),
                  Os_it( os, "\n"));
            break;

          case CGAL::IO::BINARY:
            copy( min_annulus.points_begin(), min_annulus.points_end(),
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
    operator >> ( std::istream& is, CGAL::Min_annulus_d<Traits_>& min_annulus)
    {
        using namespace std;
        
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            typedef  CGAL::Min_annulus_d<Traits_>::Point  Point;
            typedef  istream_iterator<Point>             Is_it;
            min_annulus.set( Is_it( is), Is_it());
            break;

          default:
            CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end


@! ----------------------------------------------------------------------------
@! Using the Linear Programming Solver
@! ----------------------------------------------------------------------------

\subsection{Using the Linear Programming Solver}
\label{sec:using_qp_solver}

We use the solver described in~\cite{s-qpego1-00} to determine the solution
of the linear programming problem~(\ref{eq:MA_as_LP}).

@macro <Min_annulus_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_QP_SOLVER_H
    #  include <CGAL/_QP_solver/QP_solver.h>
    #endif
@end


@! ----------------------------------------------------------------------------
\subsubsection{Representing the Linear Program}

We need a model of the concept \ccc{QP_representation}, which defines the
number types and iterators used by the QP solver.

@macro <Min_annulus_d declarations> += @begin
    
    template < class ET_, class NT_, class Point, class PointIterator,
               class Access_coord, class Access_dim >
    struct LP_rep_min_annulus_d;
@end

@macro <Min_annulus_d LP representation> = @begin
    template < class ET_, class NT_, class Point, class PointIterator,
               class Access_coord, class Access_dim >
    struct LP_rep_min_annulus_d {
        typedef  ET_                    ET;
        typedef  NT_                    NT;

        @<Min_annulus_d LP representation: iterator types>

        typedef  CGAL::Tag_true         Is_lp;
    };
@end

The matrix $A$ and the vectors $b$ and $c$ are stored in the data members
\ccc{a_matrix}, \ccc{b_vector}, and \ccc{c_vector}, respectively.

@macro <Min_annulus_d private types> += @begin

    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;
@end

@macro <Min_annulus_d data members> += @begin

    NT_matrix                a_matrix;  // matrix `A' of dual LP
    NT_vector                b_vector;  // vector `b' of dual LP
    NT_vector                c_vector;  // vector `c' of dual LP
@end

@macro <Min_annulus_d declarations> += @begin

    template < class NT >
    struct LP_rep_row_of_a {
        typedef  std::vector<NT>                         argument_type;
        typedef  typename argument_type::const_iterator  result_type;
      
        result_type
        operator ( ) ( const argument_type& v) const { return v.begin(); }
    };
@end
    
@macro <Min_annulus_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
    #  include <CGAL/_QP_solver/Join_random_access_iterator.h>
    #endif
@end

@macro <Min_annulus_d LP representation: iterator types> += @begin
    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;

    typedef  typename NT_matrix::const_iterator NTMCI;
    typedef  CGAL::Join_random_access_iterator_1<
                 NTMCI, LP_rep_row_of_a<NT> >  A_iterator;
    typedef  typename NT_vector::const_iterator
                                        B_iterator;
    typedef  typename NT_vector::const_iterator
                                        C_iterator;

    typedef  A_iterator                 D_iterator;     // dummy
@end

Now we are able to define the fully specialized type of the LP solver.

@macro <Min_annulus_d Solver type> = @begin
    // LP solver
    typedef  CGAL::LP_rep_min_annulus_d<
                 ET, NT, Point, typename std::vector<Point>::const_iterator,
                 Access_coordinates_begin_d, Access_dimension_d >
                                        LP_rep;
    typedef  CGAL::QP_solver< LP_rep >  Solver;
    typedef  typename Solver::Pricing_strategy
                                        Pricing_strategy;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Computing the Smallest Enclosing Annulus}

We set up the dual of the linear program, solve it, and compute center and
squared radii of the smallest enclosing annulus.

@macro <Min_annulus_d private member functions> += @begin
    
    // compute smallest enclosing annulus
    void
    compute_min_annulus( )
    {
        // clear inner and outer support points
        inner_indices.erase( inner_indices.begin(), inner_indices.end());
        outer_indices.erase( outer_indices.begin(), outer_indices.end());

        if ( is_empty()) {
            center_coords.resize( 1);
            sqr_i_rad_numer = -ET( 1);
            sqr_o_rad_numer = -ET( 1);
            return;
        }
        
        if ( number_of_points() == 1) {
            inner_indices.push_back( 0);
            outer_indices.push_back( 0);
            center_coords.resize( d+1);
            std::copy( tco.access_coordinates_begin_d_object()( points[ 0]),
                       tco.access_coordinates_begin_d_object()( points[ 0])+d,
                       center_coords.begin());
            center_coords[ d] = ET( 1);
            sqr_i_rad_numer = ET( 0);
            sqr_o_rad_numer = ET( 0);
            sqr_rad_denom   = ET( 1);
            return;
        }
        
        // set up and solve dual LP
        @<Min_annulus_d compute_min_annulus: set up and solve dual LP>

        // compute center and squared radius
        @<Min_annulus_d compute_min_annulus: compute center and ...>
    }
@end

@macro <Min_annulus_d compute_min_annulus: set up and solve dual LP> = @begin
    int i, j;
    NT  nt_0 = 0, nt_1 = 1, nt_2 = 2;
    NT  nt_minus_1 = -nt_1, nt_minus_2 = -nt_2;

    // vector b
    b_vector.resize( d+2);
    for ( j = 0; j < d; ++j) b_vector[ j] = nt_0;
    b_vector[ d  ] = nt_1;
    b_vector[ d+1] = nt_minus_1;
    
    // matrix A, vector c
    a_matrix.erase( a_matrix.begin(), a_matrix.end());
    a_matrix.insert( a_matrix.end(), 2*points.size(), NT_vector( d+2));
    c_vector.resize( 2*points.size());
    for ( i = 0; i < number_of_points(); ++i) {
        typename Traits::Access_coordinates_begin_d::Coordinate_iterator
            coord_it = tco.access_coordinates_begin_d_object()( points[i]);
        NT  sum = 0;
        for ( j = 0; j < d; ++j) {
            a_matrix[ 2*i  ][ j] = nt_2*coord_it[ j];
            a_matrix[ 2*i+1][ j] = nt_minus_2*coord_it[ j];
            sum += NT( coord_it[ j])*NT( coord_it[ j]);
        }
        a_matrix[ 2*i  ][ d  ] = nt_1;
        a_matrix[ 2*i+1][ d  ] = nt_0;
        a_matrix[ 2*i  ][ d+1] = nt_0;
        a_matrix[ 2*i+1][ d+1] = nt_minus_1;
        c_vector[ 2*i  ] =  sum;
        c_vector[ 2*i+1] = -sum;
    }
    typedef  typename LP_rep::A_iterator  A_it;
    typedef  typename LP_rep::D_iterator  D_it;
    solver.set( 2*points.size(), d+2, d+2,
                A_it( a_matrix.begin()), b_vector.begin(),
                c_vector.begin(), D_it());
    solver.init();
    solver.solve();
@end

@macro <Min_annulus_d compute_min_annulus: compute center and ...> = @begin
    ET sqr_sum = 0;
    center_coords.resize( ambient_dimension()+1);
    for ( i = 0; i < d; ++i) {
        center_coords[ i] = -solver.dual_variable( i);
        sqr_sum += center_coords[ i] * center_coords[ i];
    }
    center_coords[ d] = solver.variables_common_denominator();
    sqr_i_rad_numer = sqr_sum
                      - solver.dual_variable( d  )*center_coords[ d];
    sqr_o_rad_numer = sqr_sum
                      - solver.dual_variable( d+1)*center_coords[ d];
    sqr_rad_denom   = center_coords[ d] * center_coords[ d];

    // split up support points
    for ( i = 0; i < solver.number_of_basic_variables(); ++i) {
        int index = solver.basic_variables_index_begin()[ i];
        if ( index % 2 == 0) {
            inner_indices.push_back( index/2);
        } else {
            outer_indices.push_back( index/2);
        }
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Choosing the Pricing Strategy}

@macro <Min_annulus_d CGAL/QP_solver includes> += @begin
    #ifndef CGAL_PARTIAL_EXACT_PRICING_H
    #  include <CGAL/_QP_solver/Partial_exact_pricing.h>
    #endif
    #ifndef CGAL_PARTIAL_FILTERED_PRICING_H
    #  include <CGAL/_QP_solver/Partial_filtered_pricing.h>
    #endif
@end

@macro <Min_annulus_d data members> += @begin
    
    typename Solver::Pricing_strategy*  // pricing strategy
                             strategyP; // of the QP solver
@end

@macro <Min_annulus_d QP-solver set-up> many = @begin
    set_pricing_strategy( NT());
@end

@macro <Min_annulus_d private member functions> += @begin
    
    template < class NT >
    void  set_pricing_strategy( NT)
        { strategyP = new CGAL::Partial_filtered_pricing<LP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    
    #ifndef _MSC_VER
    void  set_pricing_strategy( ET)
        { strategyP = new CGAL::Partial_exact_pricing<LP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    #endif
@end


@! ============================================================================
@! Test Programs
@! ============================================================================

\clearpage
\section{Test Programs} \label{sec:test_program}

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

The function \ccc{test_Min_annulus_d}, invoked with a set of points and a
traits class model, calls each function of \ccc{Min_annulus_d} at least once
to ensure code coverage. If \ccc{verbose} is set to $-1$, the function is
``silent'', otherwise some diagnosing output is written to the standard
error stream.

@macro <Min_annulus_d test function> = @begin
    #define COVER(text,code) \
                verr0.out().width( 32); verr0 << text << "..." << flush; \
                verrX.out().width(  0); verrX << "==> " << text << endl \
                  << "----------------------------------------" << endl; \
                { code } verr0 << "ok."; verr << endl;
    
    template < class ForwardIterator, class Traits >
    void
    test_Min_annulus_d( ForwardIterator first, ForwardIterator last,
                        const Traits& traits, int verbose)
    {
        using namespace std;
        
        typedef  CGAL::Min_annulus_d< Traits >  Min_annulus;
        typedef  typename Traits::Point_d      Point;

        CGAL::Verbose_ostream verr ( verbose >= 0);
        CGAL::Verbose_ostream verr0( verbose == 0);
        CGAL::Verbose_ostream verrX( verbose >  0);
        CGAL::set_pretty_mode( verr.out());

        bool  is_valid_verbose = ( verbose > 0);

        // constructors
        COVER( "default constructor",
            Min_annulus  ms( traits, verbose, verr.out());
            assert( ms.is_valid( is_valid_verbose));
            assert( ms.is_empty());
        )

        COVER( "point set constructor",
            Min_annulus  ms( first, last, traits, verbose, verr.out());
            assert( ms.is_valid( is_valid_verbose));
        )

        Min_annulus  min_annulus( first, last);
        COVER( "ambient dimension",
            Min_annulus  ms;
            assert( ms.ambient_dimension() == -1);
            verrX << min_annulus.ambient_dimension() << endl;
        )

        COVER( "(number of) points",
            verrX << min_annulus.number_of_points() << endl;
            typename Min_annulus::Point_iterator
                point_it = min_annulus.points_begin();
            for ( ; point_it != min_annulus.points_end(); ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_annulus.points_end() - min_annulus.points_begin())
                    == min_annulus.number_of_points());
        )

        COVER( "(number of) support points",
            verrX << min_annulus.number_of_support_points() << endl;
            typename Min_annulus::Support_point_iterator
                point_it = min_annulus.support_points_begin();
            for ( ; point_it != min_annulus.support_points_end(); ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_annulus.support_points_end()
                      - min_annulus.support_points_begin())
                    == min_annulus.number_of_support_points());
        )

        COVER( "(number of) inner support points",
            verrX << min_annulus.number_of_inner_support_points() << endl;
            typename Min_annulus::Inner_support_point_iterator
                point_it = min_annulus.inner_support_points_begin();
            for ( ; point_it != min_annulus.inner_support_points_end();
                  ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_annulus.inner_support_points_end()
                      - min_annulus.inner_support_points_begin())
                    == min_annulus.number_of_inner_support_points());
        )

        COVER( "(number of) outer support points",
            verrX << min_annulus.number_of_outer_support_points() << endl;
            typename Min_annulus::Outer_support_point_iterator
                point_it = min_annulus.outer_support_points_begin();
            for ( ; point_it != min_annulus.outer_support_points_end();
                  ++point_it) {
                verrX << *point_it << endl;
            }
            assert( ( min_annulus.outer_support_points_end()
                      - min_annulus.outer_support_points_begin())
                    == min_annulus.number_of_outer_support_points());
        )

        COVER( "center and squared radii",
            verrX << "center:";
            typename Min_annulus::Coordinate_iterator  coord_it;
            for ( coord_it  = min_annulus.center_coordinates_begin();
                  coord_it != min_annulus.center_coordinates_end();
                  ++coord_it) {
                verrX << ' ' << *coord_it;
            }
            verrX << endl << "squared inner radius: "
                  << min_annulus.squared_inner_radius_numerator()   << " / "
                  << min_annulus.squared_radii_denominator() << endl;
            verrX << endl << "squared outer radius: "
                  << min_annulus.squared_outer_radius_numerator()   << " / "
                  << min_annulus.squared_radii_denominator() << endl;
        )

        COVER( "predicates",
            CGAL::Bounded_side  bounded_side;
            bool                has_on_bounded_side;
            bool                has_on_boundary;
            bool                has_on_unbounded_side;
            Point               p;
            typename Min_annulus::Point_iterator
                point_it = min_annulus.points_begin();
            for ( ; point_it != min_annulus.points_end(); ++point_it) {
                p = *point_it;
                bounded_side          = min_annulus.bounded_side( p);
                has_on_bounded_side   = min_annulus.has_on_bounded_side( p);
                has_on_boundary       = min_annulus.has_on_boundary( p);
                has_on_unbounded_side = min_annulus.has_on_unbounded_side( p);
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
            min_annulus.clear();
            verrX << "min_annulus is" << ( min_annulus.is_empty() ? "" : " not")
                  << " empty." << endl;
            assert( min_annulus.is_empty());
        )

        COVER( "insert (single point)",
            min_annulus.insert( *first);
            assert( min_annulus.is_valid( is_valid_verbose));
            assert( min_annulus.is_degenerate());
        )

        COVER( "insert (point set)",
            min_annulus.insert( first, last);
            assert( min_annulus.is_valid( is_valid_verbose));
        )

        COVER( "traits class access",
            min_annulus.traits();
        )

        COVER( "I/O",
            verrX << min_annulus;
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

@macro <Min_annulus_d test: includes and typedefs>(1) many += @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Min_annulus_d.h>
    #include <CGAL/Optimisation_d_traits_@1.h>
@end

We use the number type \ccc{leda_integer} from \leda{} for the first
variant.

@macro <Min_annulus_d test: includes and typedefs> += @begin

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
type used by the underlying linear programming solver is
\ccc{GMP::Double}, i.e.~an arbitrary precise floating-point type based on
\textsc{Gmp}'s integers. To speed up the pricing, we use \ccc{double}
arithmetic.

@macro <Min_annulus_d test: includes and typedefs> += @begin

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

@macro <Min_annulus_d test: includes and typedefs> += @begin

    #include <CGAL/Random.h>
    #include <vector>
@end

@macro <Min_annulus_d test: includes and typedefs (2/3D)>(1) many = @begin
    #include <CGAL/point_generators_@1.h>
    #include <CGAL/algorithm.h>
    #include <iterator>
@end

@macro <Min_annulus_d test: generate point set>(3) = @begin
    std::vector<K_@1::Point_@2>  points_@1;
    points_@1.reserve( 100);
    CGAL::cpp0x::copy_n( CGAL::Random_points_on_@3_@2<K_@1::Point_@2>( 0x100000),
                  100, std::back_inserter( points_@1));
@end

The traits class model with $d$-dimensional points is tested with $d = 5$
(variant 1) and $d = 10$ (variant 2). The points are distributed uniformly
in a $d$-cube.

@macro <Min_annulus_d test: generate point set (dD)>(1) = @begin
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

@macro <Min_annulus_d test: includes and typedefs> += @begin

    #include "test_Min_annulus_d.h"
@end

@macro <Min_annulus_d test: call test function>(1) many = @begin
    CGAL::test_Min_annulus_d( points_@1.begin(), points_@1.end(),
                              Traits_@1(), verbose);
@end

Each of the two test variants is compiled and executed only if the
respective number type is available.

@macro <Min_annulus_d test: test variant output>(1) many = @begin
    verr << endl
         << "==================================="
         << "===================================" << endl
         << "Testing `Min_annulus_d' with traits class model" << endl
         << "==> " << TEST_VARIANT_@1 << endl
         << "==================================="
         << "===================================" << endl
         << endl;
@end

@macro <Min_annulus_d test: test variant>(3) many = @begin
    // test variant @1
    // --------------
    #ifdef TEST_VARIANT_@1

        @<Min_annulus_d test: test variant output>(@1)

        // generate point set
        @<Min_annulus_d test: generate point set>(@1,@2,@3)

        // call test function
        @<Min_annulus_d test: call test function>(@1)

    #endif
@end

@macro <Min_annulus_d test: test variant (dD)>(1) many = @begin
    // test variant @1
    // --------------
    #ifdef TEST_VARIANT_@1

        @<Min_annulus_d test: test variant output>(@1)

        // generate point set
        @<Min_annulus_d test: generate point set (dD)>(@1)

        // call test function
        @<Min_annulus_d test: call test function>(@1)

    #endif
@end

The complete bodies of the test programs look as follows. Verbose output
can be enabled by giving a number between 0 and 3 at the command line.

@macro <Min_annulus_d test: command line argument> many = @begin
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr ( verbose >= 0); verr  << "";
@end

@macro <Min_annulus_d test: main>(2) many = @begin
    using namespace std;

    @<Min_annulus_d test: command line argument>
    
    @<Min_annulus_d test: test variant>(1,@1,@2)

    @<Min_annulus_d test: test variant>(2,@1,@2)

    return 0;
@end

@macro <Min_annulus_d test: main (dD)> = @begin
    using namespace std;

    @<Min_annulus_d test: command line argument>
    
    @<Min_annulus_d test: test variant (dD)>(1)

    @<Min_annulus_d test: test variant (dD)>(2)

    return 0;
@end

@! ============================================================================
@! Files
@! ============================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Min_annulus_d.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Min\_annulus\_d.h}

@file <include/CGAL/Min_annulus_d.h> = @begin
    @<file header>(
        "include/CGAL/Min_annulus_d.h",
        "Smallest enclosing annulus in arbitrary dimension")

    #ifndef CGAL_MIN_ANNULUS_D_H
    #define CGAL_MIN_ANNULUS_D_H

    // includes
    // --------
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<Min_annulus_d CGAL includes>
    @<Min_annulus_d CGAL/QP_solver includes>
    @<Min_annulus_d standard includes>

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_annulus_d declarations>
    
    // Class interfaces
    // ================
    @<Min_annulus_d interface>

    @<Min_annulus_d LP representation>

    // Function declarations
    // =====================
    @<Min_annulus_d I/O operators declaration>

    @<dividing line>
    
    // Class implementation
    // ====================

    @<Min_annulus_d validity check>
    
    @<Min_annulus_d I/O operators>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_ANNULUS_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_annulus_d.h
@! ----------------------------------------------------------------------------

\subsection{test/Min\_annulus\_d/test\_Min\_annulus\_d.h}

@file <test/Min_annulus_d/test_Min_annulus_d.h> = @begin
    @<file header>(
        "test/Min_annulus_d/test_Min_annulus_d.h",
        "test function for smallest enclosing annulus")

    #ifndef CGAL_TEST_MIN_ANNULUS_D_H
    #define CGAL_TEST_MIN_ANNULUS_D_H

    // includes
    #ifndef CGAL_IO_VERBOSE_OSTREAM_H
    #  include <CGAL/IO/Verbose_ostream.h>
    #endif
    #include <cassert>
    
    @<namespace begin>("CGAL")

    @<Min_annulus_d test function>
    
    @<namespace end>("CGAL")

    #endif // CGAL_TEST_MIN_ANNULUS_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_annulus_d_2.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_annulus\_d/test\_Min\_annulus\_d\_2.C}

@file <test/Min_annulus_d/test_Min_annulus_d_2.C> = @begin
    @<file header>(
        "test/Min_annulus_d/test_Min_annulus_d_2.C",
        "test program for smallest enclosing annulus (2D traits class)")

    // includes and typedefs
    // ---------------------
    @<Min_annulus_d test: includes and typedefs>(2)
    @<Min_annulus_d test: includes and typedefs (2/3D)>(2)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Min_annulus_d test: main>(2,circle)
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_annulus_d_3.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_annulus\_d/test\_Min\_annulus\_d\_3.C}

@file <test/Min_annulus_d/test_Min_annulus_d_3.C> = @begin
    @<file header>(
        "test/Min_annulus_d/test_Min_annulus_d_3.C",
        "test program for smallest enclosing annulus (3D traits class)")

    // includes and typedefs
    // ---------------------
    @<Min_annulus_d test: includes and typedefs>(3)
    @<Min_annulus_d test: includes and typedefs (2/3D)>(3)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Min_annulus_d test: main>(3,sphere)
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_annulus_d_d.C
@! ----------------------------------------------------------------------------

\subsection{test/Min\_annulus\_d/test\_Min\_annulus\_d\_d.C}

@file <test/Min_annulus_d/test_Min_annulus_d_d.C> = @begin
    @<file header>(
        "test/Min_annulus_d/test_Min_annulus_d_d.C",
        "test program for smallest enclosing annulus (dD traits class")

    // includes and typedefs
    // ---------------------
    @<Min_annulus_d test: includes and typedefs>(d)

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        @<Min_annulus_d test: main (dD)>
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
        "Min_annulus_d",
        "Geometric Optimisation",
        "Min_annulus_d",
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
