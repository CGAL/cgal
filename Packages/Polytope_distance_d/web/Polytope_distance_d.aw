@! ============================================================================
@! The CGAL Library
@! Implementation: Distance of Polytopes in Arbitrary Dimension
@! ----------------------------------------------------------------------------
@! file  : Polytope_distance_d/web/Polytope_distance_d.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Polytope_distance_d WIP $
@! $Revision$
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
\setlength{\skip\footins}{3ex}

@! LaTeX macros
\newcommand{\remark}[2]{[\textbf{#1:} \emph{#2}]}

\newcommand{\linebreakByHand}{\ccTexHtml{\linebreak[4]}{}}
\newcommand{  \newlineByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

@! settings for `cc_manual.sty'
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{\rm}
\newcommand{\cgalColumnLayout}{%
  \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
  \ccPropagateThreeToTwoColumns}
\newcommand{\cgalPolytopeDistanceLayout}{%
  \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
  \ccPropagateThreeToTwoColumns}


@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsRevision}{$Revision$}
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
  Education and Sciences for the latter project.}
\renewcommand{\thefootnote}{\arabic{footnote}}

@! --------
@! Abstract
@! --------

\begin{abstract}
  \ldots
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



@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Polytope_distance_d.h
@! ----------------------------------------------------------------------------

\subsection{Polytope\_distance\_d.h}

@file <include/CGAL/Polytope_distance_d.h> = @begin
    @<file header>(
        "include/CGAL/Polytope_distance_d.h",
        "Distance of convex polytopes in arbitrary dimension")

    #ifndef CGAL_POLYTOPE_DISTANCE_D_H
    #define CGAL_POLYTOPE_DISTANCE_D_H

    // includes
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<namespace begin>("CGAL")
    

    @<namespace end>("CGAL")
    
    #endif // CGAL_POLYTOPE_DISTANCE_D_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Polytope_distance_d.C
@! ----------------------------------------------------------------------------

\subsection{test\_Polytope\_distance\_d.C}

@file <test/Polytope_distance_d/test_Polytope_distance_d.C> = @begin
    @<file header>(
        "test/Polytope_distance_d/test_Polytope_distance_d.C",
        "test program for distance of convex polytopes")

    // main
    // ----
    int
    main( int argc, char* argv[])
    {

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
        "Geometric Optimisation",
        "Polytope_distance_d", "Polytope_distance_d",
        "$Revision$","$Date$",
        "Sven Schönherr",
        "Sven Schönherr <sven@@inf.fu-berlin.de>",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ===== EOF ==================================================================
