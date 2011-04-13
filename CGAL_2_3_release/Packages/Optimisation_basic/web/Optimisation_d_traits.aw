@! ============================================================================
@! The CGAL Library
@! Implementation: Traits Classes for dD Optimisation Algorithms
@! ----------------------------------------------------------------------------
@! file  : web/Optimisation_d_traits.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Optimisation_basic WIP $
@! $Revision$
@! $Date$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amsmath}
@usepackage{amssymb}
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

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 20 mm
@t title titlefont centre "Traits Class Models for"
@t vskip 0 mm
@t title titlefont centre "dD Optimisation Algorithms*"
@t vskip 10 mm
@t title smalltitlefont centre "Sven Schönherr"
\begin{center}
  \textbf{ETH Z{\"u}rich}
\end{center}
@t vskip 10 mm
{\small
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Package: Optimisation_basic WIP+\cgalWIP\verb+$+ \\
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
  We provide traits class models for $d$-dimensional optimisation
  algorithms using the two-, three-, and $d$-dimensional \cgal~kernel,
  respectively.
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

We provide traits class models for $d$-dimensional optimisation
algorithms using the two-, three-, and $d$-dimensional \cgal~kernel,
respectively.

\ldots


@! ============================================================================
@! Reference Pages
@! ============================================================================

\clearpage
\section{Reference Pages} \label{sec:reference_pages}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

@p maximum_input_line_length = 82

@! ----------------------------------------------------------------------------
@! Concept: OptimisationDTraits
@! ----------------------------------------------------------------------------

\subsectionRef{Concept}{OptimisationDTraits}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/OptimisationDTraits.tex}

@! ----------------------------------------------------------------------------
@! Class: Optimisation_d_traits_2
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{%
  CGAL::Optimisation\_d\_traits\_2\texttt{<}K,ET,NT\texttt{>}}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/Optimisation_d_traits_2.tex}

@! ----------------------------------------------------------------------------
@! Class: Optimisation_d_traits_3
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{%
  CGAL::Optimisation\_d\_traits\_3\texttt{<}K,ET,NT\texttt{>}}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/Optimisation_d_traits_3.tex}

@! ----------------------------------------------------------------------------
@! Class: Optimisation_d_traits_d
@! ----------------------------------------------------------------------------

\subsectionRef{Class}{%
  CGAL::Optimisation\_d\_traits\_d\texttt{<}K,ET,NT\texttt{>}}
\input{../doc_tex/basic/Optimisation/Optimisation_ref/Optimisation_d_traits_d.tex}

@p maximum_input_line_length = 80


@! ============================================================================
@! Traits Class Models
@! ============================================================================

\clearpage
\section{Traits Class Models} \label{sec:traits_class_models}

@! ----------------------------------------------------------------------------
@! The Class Template CGAL::Optimisation_d_traits_2<K,ET,NT>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Optimisation\_d\_traits\_2\texttt{<}K,ET,NT\texttt{>}}
\label{sec:Optimisation_d_traits_2}

The first template argument of \ccc{Optimisation_d_traits_2} is expected to
be a \cgal\ representation class. The second and third template argument
are expected to be number types fulfilling the requirements of a \cgal\ 
number type. They have default type \ccc{K::RT}.

@macro <Optimisation_d_traits_2 declaration> = @begin
    template < class K_, class ET_ = CGAL_TYPENAME_MSVC_NULL K_::RT,
                         class NT_ = CGAL_TYPENAME_MSVC_NULL K_::RT >
    class Optimisation_d_traits_2;
@end

The interface consists of the types and member functions described in
Section~\ref{ccRef_CGAL::Optimisation_d_traits_2<K,ET,NT>}.3.

@macro <Optimisation_d_traits_2 interface> = @begin
    template < class K_, class ET_, class NT_>
    class Optimisation_d_traits_2 {
      public:
        // self
        typedef  K_                         K;
        typedef  ET_                        ET;
        typedef  NT_                        NT;
        typedef  Optimisation_d_traits_2<K,ET,NT>
                                            Self;

        // types
        typedef  typename K::Point_2        Point_d;

        typedef  typename K::Rep_tag        Rep_tag;

        typedef  typename K::RT             RT;
        typedef  typename K::FT             FT;

        typedef  Access_dimension_2<K>      Access_dimension_d;
        typedef  Access_coordinates_begin_2<K>
                                            Access_coordinates_begin_d;

        typedef  Construct_point_2<K>       Construct_point_d;

        // creation
        Optimisation_d_traits_2( ) { }
        Optimisation_d_traits_2( const Optimisation_d_traits_2<K_,ET_,NT_>&) {}

        // operations
        Access_dimension_d
        access_dimension_d_object( ) const
            { return Access_dimension_d(); }

        Access_coordinates_begin_d
        access_coordinates_begin_d_object( ) const
            { return Access_coordinates_begin_d(); }

        Construct_point_d
        construct_point_d_object( ) const
            { return Construct_point_d(); }
    };
@end


@! ----------------------------------------------------------------------------
@! The Class Template CGAL::Optimisation_d_traits_3<K,ET,NT>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Optimisation\_d\_traits\_3\texttt{<}K,ET,NT\texttt{>}}
\label{sec:Optimisation_d_traits_3}

The first template argument of \ccc{Optimisation_d_traits_3} is expected to
be a \cgal\ representation class. The second and third template argument
are expected to be number types fulfilling the requirements of a \cgal\ 
number type. They have default type \ccc{K::RT}.

@macro <Optimisation_d_traits_3 declaration> = @begin
    template < class K_, class ET_ = CGAL_TYPENAME_MSVC_NULL K_::RT,
                         class NT_ = CGAL_TYPENAME_MSVC_NULL K_::RT >
    class Optimisation_d_traits_3;
@end

The interface consists of the types and member functions described in
Section~\ref{ccRef_CGAL::Optimisation_d_traits_3<K,ET,NT>}.4.

@macro <Optimisation_d_traits_3 interface> = @begin
    template < class K_, class ET_, class NT_>
    class Optimisation_d_traits_3 {
      public:
        // self
        typedef  K_                         K;
        typedef  ET_                        ET;
        typedef  NT_                        NT;
        typedef  Optimisation_d_traits_3<K,ET,NT>
                                            Self;

        // types
        typedef  typename K::Point_3        Point_d;

        typedef  typename K::Rep_tag        Rep_tag;

        typedef  typename K::RT             RT;
        typedef  typename K::FT             FT;

        typedef  Access_dimension_3<K>      Access_dimension_d;
        typedef  Access_coordinates_begin_3<K>
                                            Access_coordinates_begin_d;

        typedef  Construct_point_3<K>       Construct_point_d;

        // creation
        Optimisation_d_traits_3( ) { }
        Optimisation_d_traits_3( const Optimisation_d_traits_3<K_,ET_,NT_>&) {}

        // operations
        Access_dimension_d
        access_dimension_d_object( ) const
            { return Access_dimension_d(); }

        Access_coordinates_begin_d
        access_coordinates_begin_d_object( ) const
            { return Access_coordinates_begin_d(); }

        Construct_point_d
        construct_point_d_object( ) const
            { return Construct_point_d(); }
    };
@end


@! ----------------------------------------------------------------------------
@! The Class Template CGAL::Optimisation_d_traits_d<K,ET,NT>
@! ----------------------------------------------------------------------------

\subsection{The Class Template \ccFont
  CGAL::Optimisation\_d\_traits\_d\texttt{<}K,ET,NT\texttt{>}}
\label{sec:Optimisation_d_traits_d}

The first template argument of \ccc{Optimisation_d_traits_d} is expected to
be a \cgal\ representation class. The second and third template argument
are expected to be number types fulfilling the requirements of a \cgal\ 
number type. They have default type \ccc{K::RT}.

@macro <Optimisation_d_traits_d declaration> = @begin
    template < class K_, class ET_ = CGAL_TYPENAME_MSVC_NULL K_::RT,
                         class NT_ = CGAL_TYPENAME_MSVC_NULL K_::RT >
    class Optimisation_d_traits_d;
@end

The interface consists of the types and member functions described in
Section~\ref{ccRef_CGAL::Optimisation_d_traits_d<K,ET,NT>}.5.

@macro <Optimisation_d_traits_d interface> = @begin
    template < class K_, class ET_, class NT_>
    class Optimisation_d_traits_d {
      public:
        // self
        typedef  K_                         K;
        typedef  ET_                        ET;
        typedef  NT_                        NT;
        typedef  Optimisation_d_traits_d<K,ET,NT>
                                            Self;

        // types
        typedef  typename K::Point_d        Point_d;

        typedef  typename K::Rep_tag        Rep_tag;

        typedef  typename K::RT             RT;
        typedef  typename K::FT             FT;

        typedef  Access_dimension_d<K>      Access_dimension_d;
        typedef  Access_coordinates_begin_d<K>
                                            Access_coordinates_begin_d;

        typedef  Construct_point_d<K>       Construct_point_d;

        // creation
        Optimisation_d_traits_d( ) { }
        Optimisation_d_traits_d( const Optimisation_d_traits_d<K_,ET_,NT_>&) {}

        // operations
        Access_dimension_d
        access_dimension_d_object( ) const
            { return Access_dimension_d(); }

        Access_coordinates_begin_d
        access_coordinates_begin_d_object( ) const
            { return Access_coordinates_begin_d(); }

        Construct_point_d
        construct_point_d_object( ) const
            { return Construct_point_d(); }
    };
@end

        
@! ============================================================================
@! Files
@! ============================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Optimisation_d_traits_2.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation\_d\_traits\_2.h}

@file <include/CGAL/Optimisation_d_traits_2.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation_d_traits_2.h",
        "Traits class (2D) for dD optimisation algorithms")

    #ifndef CGAL_OPTIMISATION_D_TRAITS_2_H
    #define CGAL_OPTIMISATION_D_TRAITS_2_H

    // includes
    #ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_2_H
    #  include <CGAL/Optimisation/Access_dimension_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_2_H
    #  include <CGAL/Optimisation/Access_coordinates_begin_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_2_H
    #  include <CGAL/Optimisation/Construct_point_2.h>
    #endif

    @<namespace begin>("CGAL")

    // Class declaration
    // =================
    @<Optimisation_d_traits_2 declaration>
    
    // Class interface
    // ===============
    @<Optimisation_d_traits_2 interface>

    @<namespace end>("CGAL")
    
    #endif // CGAL_OPTIMISATION_D_TRAITS_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_d_traits_3.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation\_d\_traits\_3.h}

@file <include/CGAL/Optimisation_d_traits_3.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation_d_traits_3.h",
        "Traits class (3D) for dD optimisation algorithms")

    #ifndef CGAL_OPTIMISATION_D_TRAITS_3_H
    #define CGAL_OPTIMISATION_D_TRAITS_3_H

    // includes
    #ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_3_H
    #  include <CGAL/Optimisation/Access_dimension_3.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_3_H
    #  include <CGAL/Optimisation/Access_coordinates_begin_3.h>
    #endif
    #ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_3_H
    #  include <CGAL/Optimisation/Construct_point_3.h>
    #endif

    @<namespace begin>("CGAL")

    // Class declaration
    // =================
    @<Optimisation_d_traits_3 declaration>
    
    // Class interface
    // ===============
    @<Optimisation_d_traits_3 interface>

    @<namespace end>("CGAL")
    
    #endif // CGAL_OPTIMISATION_D_TRAITS_3_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_d_traits_d.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation\_d\_traits\_d.h}

@file <include/CGAL/Optimisation_d_traits_d.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation_d_traits_d.h",
        "Traits class (dD) for dD optimisation algorithms")

    #ifndef CGAL_OPTIMISATION_D_TRAITS_D_H
    #define CGAL_OPTIMISATION_D_TRAITS_D_H

    // includes
    #ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_D_H
    #  include <CGAL/Optimisation/Access_dimension_d.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_D_H
    #  include <CGAL/Optimisation/Access_coordinates_begin_d.h>
    #endif
    #ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H
    #  include <CGAL/Optimisation/Construct_point_d.h>
    #endif

    @<namespace begin>("CGAL")

    // Class declaration
    // =================
    @<Optimisation_d_traits_d declaration>
    
    // Class interface
    // ===============
    @<Optimisation_d_traits_d interface>

    @<namespace end>("CGAL")
    
    #endif // CGAL_OPTIMISATION_D_TRAITS_D_H

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
        "Optimisation_basic",
        "Geometric Optimisation",
        "Optimisation_d_traits",
        "$Revision$","$Date$",
        "Sven Schönherr <sven@@inf.ethz.ch>",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ===== EOF ==================================================================
