%%%%%%%%%%%%%%%%%%%%%% start of egpubl.cls %%%%%%%%%%%%%%%%%%%%%%
\def\fileversion{3.30}
\def\filedate{2004/01/30}
\def\docdate {2003/01/30}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{egpubl}[\filedate \space
                       Document Style for EG Publications v\fileversion ]
%
% EG Publications document style
% (based on EG Annual Event document style, which in turn was based on CFG Journal)
% Copyright (c) 1997-2002 D.Fellner
%
% based on the original LaTeX ARTICLE DOCUMENT STYLE
%    for LaTeX version 2.09 by Leslie Lamport
% ...modifications by Sue Brooks, BUSS Ltd
% ... and Graham White 8.ix.92
% ...redesign for LaTeX2e and adaption to Conf. Proc. by D.Fellner
%
% 30-01-04 DF  3.30 support both latex and pdflatex
% 11-08-03 DF  3.20 adapt to new citation style (via eg-alpha.bst)
%                   and fix refs to Blasckwell (address, name)
% 30-04-03 DF  3.15 replace use of \tt by \ttfamily
% 17-03-03 DF  3.14 add string constants to describe function of editor ...
% 09-02-03 DF  3.13 fix shape of page numbers in running heads (upshape)
% 22-04-02 DF  3.12 cleanup of some redundant structures; 
%                   commented-out sections deleted
% 22-02-02 DF  3.11 delete page range information on top line of first page;
% 22-11-01 DF  3.10 include keyword-support;
%                   increase area top-left at first page
% 23-06-01 DF  3.00 integrate support for EG WS Proc
%                   \WsSubmission, \WsPaper
% 15-06-01 DF  2.90 activate computation of page range through an
%                   .lbl file holding a \pEndPage command
% 07-06-01 DF  2.80 fix problem with empty last page 
% 27-03-00 DF  2.71 replace 'Short Paper' by 'Short Presentation'
% 23-03-00 DF  2.70 inlude support for ShortPaper
% 25-02-00 DF  2.60 use \volume to determine year on first page
%                   for all types of publications (not submissions);
%                   use \pStartPage to initialize page counter and pagerange;
% 25-02-00 DF  2.50 final polishing of (c) lines;
%                   introduce \issue, \pStartPage, \pEndPage to
%                   replace \pnumber and \pagerange;
%                   selection between printed and electronic version
% 10-10-99 DF  2.40 implement commands \JournalSubmission, \JournalPaper;
%                   define p@volume and compute as function of year;
%                   parametrize copyright text (to cater for different 
%                   def's of year);
%                   replaced \editors by \ConfEditors;
%                   customize for EG 00
% 20-11-98 DF  2.20 customize for EG 99 and introduce \ConfYear
% 26-02-98 DF  2.10 enhance level of customization
%                   \editors, \EGyear, 
% 11-04-97 DF  2.00 port to LaTeX 2e and customization for EG'97
%                   Introduced \STAR, \TUTORIAL
% 06-01-94 AJW 1.10 NFSS compatible
% 15-03-93 MR  1.00 Added check for \CUP@mtlplain@loaded and
%                   \CUP@cmtt@loaded flags for use at CUP.
%                   This allows the style file to check for the use
%                   of `mtclplain' instead of `mtlplain'. `mtclplain'
%                   has all it's mttt fonts changed to cmtt.
% 10-03-93 RM       \pnumber, \looseness in catchline address,
%                   spaces in running head text.
% 25-02-93 SFB 0.01.01 keywords,
%              0.01.02 catchline default
%              0.01.03 copyright in running foot
%              0.01.04 onecolumn style option
%              0.01.05 appendix titles...
%              0.01.06 section headings -- bold and italic
%
\newif\ifoldfss
\newif\ifnfssone
\newif\ifnfsstwo
\def\f@s@s{}%
%
\@ifundefined{DeclareFontShape}%
  {%
   \@ifundefined{selectfont}%
    {\global\oldfsstrue\let\reset@font=\relax}% old font selection scheme
    {\global\nfssonetrue\def\f@s@s{with NFSS release 1,\space}}%
  }%
  {\global\nfsstwotrue\def\f@s@s{with NFSS release 2,\space}}%
%
% \typeout{Document Style `CGF' \f@s@s v1.10 released 6th January 1994}
%
%  Journals use two-sided printing.
%
\@twosidetrue
\@twocolumntrue
\@mparswitchtrue
\def\ds@draft{\overfullrule 5pt}
\def\ds@twocolumn{\@twocolumntrue}
\def\ds@onecolumn{\@twocolumnfalse}
\newif\ifSFB@galley
\def\ds@galley{\SFB@galleytrue\@twosidefalse\@mparswitchfalse}
\@options
%
%  ****************************************
%  *               FONTS                  *
%  ****************************************
%
\lineskip 1pt
\normallineskip 1pt
\def\baselinestretch{1}
  \def\@normalsize{\@setsize\normalsize{11pt}\ixpt\@ixpt
   \abovedisplayskip 5.5pt plus 2pt minus .5pt
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip 5.5pt plus 2pt
   \belowdisplayshortskip \abovedisplayshortskip
   \let\@listi\@listI}
%% DF
\let\normalsize=\@normalsize
%%
\def\small{\@setsize\small{10pt}\viiipt\@viiipt
 \abovedisplayskip 5pt plus 1pt minus .5pt
 \belowdisplayskip \abovedisplayskip
 \abovedisplayshortskip 5pt plus 1pt
 \belowdisplayshortskip \abovedisplayshortskip
 \def\@listi{\leftmargin\leftmargini
  \topsep  \z@ plus 1pt minus .5pt
  \parsep  \z@
  \itemsep \z@}}
\def\footnotesize{\@setsize\footnotesize{10pt}\viiipt\@viiipt
 \abovedisplayskip 5pt plus 1pt minus .5pt
 \belowdisplayskip \abovedisplayskip
 \abovedisplayshortskip 5pt plus 1pt
 \belowdisplayshortskip \abovedisplayshortskip
 \def\@listi{\leftmargin\leftmargini
  \topsep  \z@ plus 1pt minus .5pt
  \parsep  \z@
  \itemsep \z@}}
\def\scriptsize{\@setsize\scriptsize{7pt}\vipt\@vipt} % ajw: was 7/8
\def\tiny{\@setsize\tiny{8pt}\vipt\@vipt}
\def\large{\@setsize\large{11pt}\xpt\@xpt}
\def\Large{\@setsize\Large{14pt}\xiipt\@xiipt}
\def\LARGE{\@setsize\LARGE{17pt}\xivpt\@xivpt}
\def\huge{\@setsize\huge{20pt}\xviipt\@xviipt}
\def\Huge{\@setsize\Huge{25pt}\xxpt\@xxpt}
\normalsize

%  ****************************************
%  *            PAGE LAYOUT               *
%  ****************************************
%
% SIDE MARGINS:
\oddsidemargin \z@ \evensidemargin \z@
\marginparwidth 4pc
\marginparsep 0.5pc

% VERTICAL SPACING:
% \topmargin -3pc
\topmargin -1pc
\headheight 11pt
\headsep    11pt
\topskip     9pt
%% DF 
%\footheight 11pt
%%
\footskip   22pt

% DIMENSION OF TEXT:
\textwidth 38pc
\textheight = 55\baselineskip \advance\textheight by \topskip
\columnsep 2pc
\columnseprule 0pt

% FOOTNOTES:
\footnotesep 6pt
\skip\footins 19.5pt plus 12pt minus 1pt

% FLOATS:
% for floats on a text page:
%  one-column mode or single-column floats in two-column mode:
\floatsep     11pt plus 5.5pt minus 2.25pt
\textfloatsep 22pt plus 11pt minus 5.5pt
\intextsep    11pt plus 5.5pt minus 2.25pt
%% DF
% \@maxsep      22pt
%%
%    two-column floats in two-column mode:
\dblfloatsep     11pt plus 5.5pt minus 2.25pt
\dbltextfloatsep 22pt plus 11pt minus 5.5pt
%% DF
% \@dblmaxsep      22pt
%%
% for floats on a separate float page or column:
% one-column mode or single-column floats in two-column mode:
\@fptop  0pt plus 0fil
\@fpsep 11pt plus 0.5fil
\@fpbot  0pt plus 3fil
% double-column floats in two-column mode.
\@dblfptop  0pt plus 0fil
\@dblfpsep 11pt plus 0.5fil
\@dblfpbot  0pt plus 3fil

\marginparpush 5.5pt

%  ****************************************
%  *           PARAGRAPHING               *
%  ****************************************
%
\parskip 5.5pt plus 1pt minus 1.5pt
\parindent 1em
\partopsep 0pt plus 1pt
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty -\@lowpenalty
\@itempenalty -\@lowpenalty
\clubpenalty=0
\widowpenalty=10000

%  ****************************************
%  *         FIGURES AND TABLES           *
%  ****************************************
%
\setcounter{topnumber}{2}
\def\topfraction{.9}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.9}
\setcounter{totalnumber}{3}
\def\textfraction{.1}
\def\floatpagefraction{.75}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.9}
\def\dblfloatpagefraction{.75}

\newcounter{table}
\def\thetable{\@arabic\c@table}
\ifSFB@galley \def\fps@table{p}
        \else \def\fps@table{tbp}
\fi
\def\ftype@table{1}
\def\ext@table{lot}
\def\fnum@table{Table \thetable}
\def\table{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}

\newcounter{figure}
\def\thefigure{\@arabic\c@figure}
\ifSFB@galley \def\fps@figure{p}
        \else \def\fps@figure{tbp}
\fi
\def\ftype@figure{2}
\def\ext@figure{lof}
\def\fnum@figure{Figure \thefigure}
\def\figure{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}

\long\def\@makecaption#1#2{\vskip 6pt
 \setbox\@tempboxa\hbox{\textbf{#1:} \reset@font\textit{#2}}
 \ifdim \wd\@tempboxa >\hsize
  \textbf{#1:} \reset@font\itshape #2\par
 \else
  \hbox to\hsize{\hfil\box\@tempboxa\hfil}
 \fi
 \vskip 6pt
}

%  ****************************************
%  *                PARTS                 *
%  ****************************************
%
\def\part{\par \addvspace{4ex}\@afterindentfalse
 \secdef\@part\@spart}
\def\@part[#1]#2{\ifnum \c@secnumdepth >\m@ne
  \refstepcounter{part}
  \addcontentsline{toc}{part}{Part \thepart: #1}
 \else \addcontentsline{toc}{part}{#1}
 \fi
 {\parindent 0pt \raggedright
  \ifnum \c@secnumdepth >\m@ne
   \reset@font\large\rmfamily PART
   \ifcase\thepart \or ONE \or TWO \or THREE \or FOUR \or FIVE
    \or SIX \or SEVEN \or EIGHT \or NINE \or TEN \else \fi
   \par \nobreak
  \fi
  \reset@font\LARGE \rmfamily #2 \markboth{}{}\par }
 \nobreak \vskip 3ex \@afterheading}
\def\@spart#1{{\parindent 0pt \raggedright
  \reset@font\LARGE \rmfamily #1\par}
 \nobreak \vskip 3ex \@afterheading}

%  ****************************************
%  *               SECTIONS               *
%  ****************************************
%
\def\section{\@startsection {section}{1}{\z@}
 {-11pt plus -5.5pt minus -1pt}
 {0.001pt}
 {\raggedright\reset@font\normalsize\bfseries}}
\def\subsection{\@startsection{subsection}{2}{\z@}
 {-11pt plus -5.5pt minus -1pt}
 {0.001pt}
 {\raggedright\reset@font\normalsize\bfseries}}
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}
 {-5.5pt plus -5.5pt minus -1pt}
 {0.001pt}
 {\raggedright\reset@font\normalsize\bfseries}}
\def\paragraph{\@startsection{paragraph}{4}{\z@}
 {5.5pt plus 2.25pt minus 1pt}
 {-0.5em}
 {\reset@font\normalsize\bfseries}}
\def\subparagraph{\@startsection{subparagraph}{4}{1em}
 {5.5pt plus 2.25pt minus 1pt}
 {-0.5em}
 {\reset@font\normalsize\bfseries}}
\setcounter{secnumdepth}{5}

\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
 \par \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
 \if@nobreak \everypar{}\else
  \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi
 \@ifstar{\SFB@ssect{#2}{#3}{#4}{#5}{#6}}%
         {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%
% \@sect{NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}[ARG1]{ARG2}
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
 \ifnum #2>\c@secnumdepth
  \def\@svsec{}%
 \else
  \refstepcounter{#1}
  \edef\@svsec{\csname the#1\endcsname.\hskip 0.5em}%
 \fi
 \@tempskipa #5\relax
 \ifdim \@tempskipa>\z@
  \begingroup #6\relax
   \ifnum #2=\@ne
    \ifSFB@appendix
     \ifSFB@appendixstar
      \@hangfrom{\hskip #3\relax}{\interlinepenalty \@M
       Appendix \csname the#1\endcsname\par}%
     \else
      \@hangfrom{\hskip #3\relax}{\interlinepenalty \@M
       Appendix \csname the#1\endcsname:\hskip 0.5em\reset@font\rmfamily #8\par}%
     \fi
    \else
     \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
    \fi
   \else
    \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
   \fi
  \endgroup
  \csname #1mark\endcsname{#7}%
  \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else
   \protect\numberline{\csname the#1\endcsname}\fi #7}%
 \else
  \def\@svsechd{#6\hskip #3\@svsec #8%\par
   \csname #1mark\endcsname{#7}
   \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else
    \protect\numberline{\csname the#1\endcsname}\fi#7}}%
 \fi
 \@xsect{#5}}
%
% \SFB@ssect{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}{ARG} ==
%
\def\SFB@ssect#1#2#3#4#5#6{\@tempskipa #4\relax
 \ifSFB@appendix
  \ifnum #1>\@ne             \edef\@svsec{#6}%
   \else \ifSFB@appendixstar \edef\@svsec{Appendix}%
                       \else \edef\@svsec{Appendix:\hskip 0.5em\reset@font\rmfamily #6}%
   \fi
  \fi
 \else
  \ifnum #1>\@ne             \edef\@svsec{#6}%
           \else             \edef\@svsec{#6}%
  \fi
 \fi
 \ifdim \@tempskipa>\z@
  \begingroup #5\relax
    \@hangfrom{\hskip #2\relax}{\interlinepenalty \@M \@svsec \par}%
  \endgroup
 \else
  \def\@svsechd{#5\hskip #2\relax \@svsec}%
 \fi
 \@xsect{#4}}

\newif\ifSFB@appendix
\newif\ifSFB@appendixstar
\def\APPENDIX{\appendix\SFB@appendixstartrue}
\def\appendix{\SFB@appendixstarfalse\par
 \SFB@appendixtrue
 \setcounter{secnumdepth}{1}
 \setcounter{section}{0}
 \setcounter{subsection}{0}
 \def\thesection{\Alph{section}}}

%  ****************************************
%  *       PARTS & SECTIONS NUMBERS       *
%  ****************************************
%
\newcounter{part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\def\thepart          {\arabic{part}}
\def\thesection       {\arabic{section}}
\def\thesubsection    {\thesection.\arabic{subsection}}
\def\thesubsubsection {\thesubsection .\arabic{subsubsection}}
\def\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}

%  ****************************************
%  *                LISTS                 *
%  ****************************************
%
\leftmargini   1em
\leftmarginii  1pc
\leftmarginiii 1pc
\leftmarginiv  1pc
\leftmarginv   1pc
\leftmarginvi  1pc
\leftmargin\leftmargini
\labelsep 0.5em
\labelwidth\leftmargini\advance\labelwidth-\labelsep
%
% label macros for Range-Left and Range-Right labels
\def\makeRLlabel#1{\rlap{#1}\hss}
\def\makeRRlabel#1{\hss\llap{#1}}
%
\def\@listI{\leftmargin\leftmargini
 \parsep  \z@
 \partopsep \z@
 \topsep  \z@ plus 1pt minus .5pt
 \itemsep \z@ plus .1pt
}
\let\@listi\@listI
\@listi
\def\@listii{\leftmargin\leftmarginii
 \labelwidth\leftmarginii\advance\labelwidth-\labelsep
 \topsep  5.5pt plus 1pt minus 1pt
 \parsep  \z@
 \partopsep \z@
 \itemsep \z@ plus .1pt
}
\def\@listiii{\leftmargin\leftmarginiii
 \labelwidth\leftmarginiii\advance\labelwidth-\labelsep
 \topsep 5.5pt plus 1pt minus 1pt
 \parsep \z@
 \partopsep \z@
 \itemsep \z@
}
\def\@listiv{\leftmargin\leftmarginiv
 \labelwidth\leftmarginiv\advance\labelwidth-\labelsep
}
\def\@listv{\leftmargin\leftmarginv
 \labelwidth\leftmarginv\advance\labelwidth-\labelsep
}
\def\@listvi{\leftmargin\leftmarginvi
 \labelwidth\leftmarginvi\advance\labelwidth-\labelsep
}
%
% ITEMIZE
\def\itemize{\ifnum \@itemdepth >3 \@toodeep
  \else \advance\@itemdepth \@ne
   \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
   \list{\csname\@itemitem\endcsname}%
    {\let\makelabel\makeRLlabel}%
  \fi}
\def\labelitemi{$\bullet$}
\def\labelitemii{--}
\def\labelitemiii{$\circ$}
\def\labelitemiv{$\diamond$}
%
% ENUMERATE
% SFB -- enumerate with optional argument to set left margin...
\def\enumerate{\ifnum \@enumdepth >3 \@toodeep \else
  \advance\@enumdepth \@ne
  \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
 \fi
 \@ifnextchar [{\@enumeratetwo}{\@enumerateone}%
}
\def\@enumeratetwo[#1]{\list{\csname label\@enumctr\endcsname}%
  {\settowidth\labelwidth{#1}
   \leftmargin\labelwidth \advance\leftmargin\labelsep
   \usecounter{\@enumctr}
   \let\makelabel\makeRRlabel}%
}
\def\@enumerateone{\list{\csname label\@enumctr\endcsname}%
  {\settowidth\labelwidth{9.}
   \leftmargin\labelwidth \advance\leftmargin\labelsep
   \usecounter{\@enumctr}
   \let\makelabel\makeRRlabel}%
}
\def\labelenumi{\theenumi}            \def\theenumi  {\arabic{enumi}.}
\def\labelenumii{\theenumii}          \def\theenumii {\alph{enumii}.}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii}        \def\theenumiii{\roman{enumiii}.}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv}          \def\theenumiv {\Alph{enumiv}.}
\def\p@enumiv{\p@enumiii\theenumiii}
%
% DESCRIPTION
\def\descriptionlabel#1{\hspace\labelsep \bfseries #1}
\def\description{\list{}{\labelwidth\z@ \itemindent-\leftmargin
  \let\makelabel\descriptionlabel}}
\let\enddescription\endlist

%  ****************************************
%  *            DISPLAYED TEXT            *
%  ****************************************
%
% VERSE
\def\verse{\let\\=\@centercr
 \list{}{\itemsep\z@
  \itemindent -1em
  \listparindent \itemindent
  \rightmargin\leftmargin
  \advance\leftmargin 1em}\item[]}
\let\endverse\endlist
%
% QUOTATION
\def\quotation{\list{}{%\listparindent 1em
%  \itemindent\listparindent
  \leftmargin 1.5pc\rightmargin\leftmargin
  \parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist
%
% QUOTE
\def\quote{\list{}{\leftmargin 1.5pc\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

%  ****************************************
%  *         OTHER ENVIRONMENTS           *
%  ****************************************
%
% THEOREM
\def\@begintheorem#1#2{\reset@font\rmfamily\list{}%
  {\leftmargin\z@\rightmargin\leftmargin}\item[]{\bfseries #1\ #2 }}
\def\@opargbegintheorem#1#2#3{\reset@font\rmfamily\list{}%
 {\leftmargin\z@\rightmargin\leftmargin}\item[]{\bfseries #1\ #2 (#3) }}
\def\@endtheorem{\endlist}
%
% PROOF
\def\proof{\reset@font\rmfamily\list{}%
  {\leftmargin\z@\rightmargin\leftmargin}\item[]{\reset@font\itshape Proof\/}}
\def\endproof{\hspace*{1em}{\begin{picture}(6.5,6.5)%
  \put(0,0){\framebox(6.5,6.5){}}\end{picture}}\endlist}
\@namedef{proof*}{\reset@font\rmfamily\list{}%
  {\leftmargin\z@\rightmargin\leftmargin}\item[]{\reset@font\itshape Proof\/}}
\@namedef{endproof*}{\endlist}
\def\proofbox{\begin{picture}(6.5,6.5)%
  \put(0,0){\framebox(6.5,6.5){}}\end{picture}}
%
% TITLEPAGE
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \else \newpage \fi \thispagestyle{empty}\c@page\z@}
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

% ARRAY AND TABULAR
%
\arraycolsep   5pt
\tabcolsep     6pt
\arrayrulewidth .5pt
\doublerulesep 1.5pt

\def\tabular{\def\@halignto{}
 \def\hline{\noalign{\ifnum0=`}\fi
  \vskip 3pt
  \hrule \@height \arrayrulewidth
  \vskip 3pt
  \futurelet \@tempa\@xhline}
 \def\fullhline{\noalign{\ifnum0=`}\fi
  \vskip 3pt
  \hrule \@height \arrayrulewidth
  \vskip 3pt
  \futurelet \@tempa\@xhline}
 \def\@xhline{\ifx\@tempa\hline
   \vskip -6pt
   \vskip \doublerulesep
  \fi
  \ifnum0=`{\fi}}
  \def\@arrayrule{\@addtopreamble{\hskip -.5\arrayrulewidth
%                                  \vrule \@width \arrayrulewidth
                                  \hskip .5\arrayrulewidth}}
\@tabular
}
% EQNARRAY
%
\def\eqnarray{\stepcounter{equation}%
 \let\@currentlabel=\theequation
 \global\@eqnswtrue
 \global\@eqcnt\z@
 \tabskip\@centering
 \let\\=\@eqncr
 $$\halign to \displaywidth\bgroup\@eqnsel\hskip\@centering
  $\displaystyle\tabskip\z@{##}$&\global\@eqcnt\@ne
  \hskip \arraycolsep
  \hfil${##}$\hfil
  &\global\@eqcnt\tw@
  \hskip \arraycolsep
  $\displaystyle\tabskip\z@{##}$\hfil
   \tabskip\@centering&\llap{##}\tabskip\z@\cr}

% TABBING
%
\tabbingsep \labelsep

% MINIPAGE
%
\skip\@mpfootins = \skip\footins

% FRAMEBOX
%
% Space left between box and text by \fbox and \framebox.
\fboxsep = 3pt
% Width of rules in box made by \fbox and \framebox.
\fboxrule = \arrayrulewidth

%  ****************************************
%  *         TITLE AND ABSTRACT           *
%  ****************************************
%
\def\maketitle{\par
 \begingroup
 \setcounter{footnote}{0}
  \def\thefootnote{\fnsymbol{footnote}}
  \def\@makefnmark{\mbox{$^{\@thefnmark}$}}
  \if@twocolumn
   \twocolumn[\vspace*{27pt}\@maketitle]
  \else
   \newpage
   \global\@topnum\z@
   \@maketitle
  \fi
  \thispagestyle{titlepage}
 \endgroup
 \setcounter{footnote}{0}
 \let\maketitle\relax
 \let\@maketitle\relax
 \gdef\@author{}
 \gdef\@title{}
 \let\thanks\relax
}
\def\and{\end{author@tabular}\vskip 6pt\par
 \begin{author@tabular}[t]{@{}c@{}}}
\def\@maketitle{\newpage
 \vspace*{10pt}
 {\centering \sloppy
  {\reset@font\huge \bfseries \@title \par}
  \vskip 35pt
  {\reset@font\normalsize
   \begin{author@tabular}[t]{@{}c@{}}\@author
   \end{author@tabular}\par}
  \vskip 46pt
 }%
}
%
\newif\ifSFB@abstract
\def\abstract{\if@twocolumn
  \SFB@abstracttrue\start@SFBbox\SFB@widerule\addvspace{2pt}\bgroup\@abstract
 \else
  \@abstract
 \fi}
\def\endabstract{\if@twocolumn
  \SFB@abstractfalse\endlist\egroup
  \SFB@widerule\addvspace{25pt}\finish@SFBbox
  \@thanks \gdef\@thanks{}%
 \else
  \endlist
 \fi}
\def\@abstract{\list{}{\leftmargin 2pc\rightmargin\leftmargin
  \parsep 0pt plus .1pt}\item[]{\textbf{Abstract}}\\\reset@font\itshape}
%
\def\classification{\vskip 5.5pt\par\reset@font\rmfamily%
Categories and Subject Descriptors {\small (according to ACM CCS)}: }
\def\endclassification{\relax}
\newcommand\CCScat[3]{#1 [#2]: #3}
%
\def\keywords{\vskip 5.5pt\par\reset@font\rmfamily\textbf{Keywords: }}
\def\endkeywords{\relax}
\newcommand\keyword[1]{[#1]}
\def\nokeywords{\relax}

\def\SFB@widerule{\bgroup\parskip 0pt\par\rule{\textwidth}{0.25pt}\par\egroup}

\def\author@tabular{\def\@halignto{}\@authortable}
\let\endauthor@tabular=\endtabular
\def\author@tabcrone{{\ifnum0=`}\fi\@xtabularcr[11pt]\reset@font\small\rmfamily
 \let\\=\author@tabcrtwo\ignorespaces}
\def\author@tabcrtwo{{\ifnum0=`}\fi\@xtabularcr[-1pt]\reset@font\small\rmfamily
 \let\\=\author@tabcrtwo\ignorespaces}
\def\@authortable{\leavevmode \hbox \bgroup $\let\@acol\@tabacol
 \let\@classz\@tabclassz \let\@classiv\@tabclassiv
 \let\\=\author@tabcrone \ignorespaces \@tabarray}

\def\start@SFBbox{\@next\@currbox\@freelist{}{}%
 \global\setbox\@currbox
 \vbox\bgroup
  \hsize \textwidth
  \@parboxrestore
}
\def\finish@SFBbox{\par\vskip -\dbltextfloatsep
  \egroup
  \global\count\@currbox\tw@
  \global\@dbltopnum\@ne
  \global\@topnum\z@
  \global\@dbltoproom\maxdimen\@addtodblcol
  \global\vsize\@colht
  \global\@colroom\@colht
}

%  ****************************************
%  *            PAGE STYLES               *
%  ****************************************
%
\mark{{}{}}
\gdef\@author{\mbox{}}
\def\author{\@ifnextchar [{\@authortwo}{\@authorone}}
\def\@authortwo[#1]#2{\gdef\@author{#2}\gdef\@shortauthor{#1}}
\def\@authorone#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\gdef\@shortauthor{}
\gdef\@title{\mbox{}}
\def\title{\@ifnextchar [{\@titletwo}{\@titleone}}
\def\@titletwo[#1]#2{\gdef\@title{#2}\gdef\@shorttitle{#1}}
\def\@titleone#1{\gdef\@title{#1}\gdef\@shorttitle{#1}}
\gdef\@shorttitle{}
%
%% DF add p@volume which is a counter automatically adjusting by the year
%% DF (the volumes volume can be computed as year-1981)
%
\gdef\j@year{1981}
\gdef\j@volume{0}
\gdef\j@issue{0}
%
\newcounter{p@volume}
% \setcounter{p@volume}{\number\year} % compute volume based on current year
% \addtocounter{p@volume}{-1981}%
%
\newcommand\volume[1]{%
 \gdef\j@volume{#1}%
 \setcounter{p@volume}{#1}%     compute year for the given volume
 \addtocounter{p@volume}{1981}%
 \gdef\j@year{\thep@volume}%
}
%
\newcommand\issue[1]{\gdef\j@issue{#1}}
%
\newcommand\pStartPage[1]{\gdef\j@pStartPage{#1}}
\gdef\j@pStartPage{1}
%
\newcommand\pEndPage[1]{\gdef\j@pEndPage{#1}}
\gdef\j@pEndPage{000}
%
%** \newcommand\pagerange[1]{\gdef\j@pagerange{#1}}
%** \gdef\j@pagerange{000--000}
%
\newcommand\journal[1]{\gdef\@journal{#1}}
\gdef\@journal{{\reset@font\large\sffamily COMPUTER GRAPHICS \reset@font\itshape forum}}
%
%% DF - add a few commands for customization
%
\newcommand\ConfEditors[1]{\gdef\p@ConfEditors{#1}}
\gdef\p@ConfEditors{N.N. and N.N.}
%
\newcommand\ConfEditorStrg[1]{\gdef\p@ConfEditorStrg{#1}}
\gdef\p@ConfEditorStrg{(Editors)}
%
\newcommand\ShortPresEditors[1]{\gdef\p@ShortPresEditors{#1}}
\gdef\p@ShortPresEditors{N.N. and N.N.}
%
\newcommand\EGyear[1]{\gdef\p@EGyear{#1}}
\gdef\p@EGyear{'0x}
%
\newcommand\ConfYear[1]{\gdef\p@ConfYear{#1}}
\gdef\p@ConfYear{200x}
%
\newcommand\ConfName[1]{\gdef\p@ConfName{#1}}
\gdef\p@ConfName{EUROGRAPHICS Workshop on ...}
%
% --- Handling of copyright text
%     start with cmds to set the text for title page and running page
%
\newcommand\copyrightTextTitPag[1]{%
 \renewcommand\p@copyrightTextTitPag{#1}%
}
\newcommand\p@copyrightTextTitPag{}%
%
\newcommand\copyrightTextRunPag[1]{%
 \renewcommand\p@copyrightTextRunPag{#1}%
}
\newcommand\p@copyrightTextRunPag{}%
%
% --- and here are the parameterized texts defined
%
\newcommand\p@copyrightTextLong[1]{%
\copyright\ The Eurographics Association 
and Blackwell Publishing #1. %
Published by Blackwell Publishing, 9600 Garsington Road, Oxford OX4 2DQ, UK
and 350 Main Street, Malden, MA 02148, USA.%
}%
\newcommand\p@copyrightTextShort[1]{%
\copyright\ The Eurographics Association and Blackwell Publishing #1.%
}%
\newcommand\p@copyrightTextEG[1]{%
\copyright\ The Eurographics Association #1.%
}%
\newcommand\p@journalSubmissionText[1]{%
submitted to COMPUTER GRAPHICS \textit{Forum} #1.%
}%
\newcommand\p@confSubmissionText[1]{%
submitted to EUROGRAPHICS #1.% 
}%
\newcommand\p@wsSubmissionText[1]{%
submitted to \textit{#1}% 
}%
% ----
%
% \ps@empty and \ps@plain defined in LATEX.TEX
%
\newcommand\EmptySectionmark[1]{}
\newcommand\EmptySubsectionmark[1]{}
%
\def\noEGpagenumber{\def\p@EGpagenumber{\mbox{}}}
\def\EGpagenumber{\def\p@EGpagenumber{\thepage}}
\EGpagenumber%
%
% ------ define Document Types
%          ConferenceSubmission, ConferencePaper, 
%          STAR, Tutorial, ShortPresentation,
%          JournalSubmission, JournalPaper,
%          WsSubmission, WsPaper
%
\def\ConferenceSubmission{%
\EGpagenumber%
\copyrightTextTitPag{\p@confSubmissionText{\p@ConfYear}}
\copyrightTextRunPag{\p@confSubmissionText{\p@ConfYear}}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   \parbox[t]{.7\textwidth}{\sloppy\raggedright
     EUROGRAPHICS \p@EGyear{} / \p@ConfEditors{}\\
     \p@ConfEditorStrg}% (Guest Editors)}
   \hfill
   \textit{Volume \j@volume} (\textit{\p@ConfYear\/}),
   \textit{Number \j@issue}}%
 }%
}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\ConferencePaper{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextLong{\number\year}}
\copyrightTextRunPag{\p@copyrightTextShort{\number\year}}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   \parbox[t]{.7\textwidth}{\sloppy\raggedright
     EUROGRAPHICS \p@EGyear{} / \p@ConfEditors{}\\
     \p@ConfEditorStrg}% (Guest Editors)}
   \hfill
%   \textit{Volume \j@volume} (\textit{\p@ConfYear\/}),
   \textit{Volume \j@volume} (\textit{\j@year\/}),
   \textit{Number \j@issue}}%
 }%
}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\JournalSubmission{%
\EGpagenumber%
\copyrightTextTitPag{\p@journalSubmissionText{(\number\month/\number\year)}}
\copyrightTextRunPag{\p@journalSubmissionText{(\number\month/\number\year)}}
%% \pagerange{000--000}%
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
%    \textit{Volume \j@volume} (\textit{\number\year}\/),
   \textit{Volume xx} (\textit{200y}\/),         % Vol + year yet unknown !
   \textit{Number z, \ pp.\ \j@pStartPage--\j@pEndPage}%
   \hfill
 }}}%
%
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\JournalPaper{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextLong{\number\year}}
\copyrightTextRunPag{\p@copyrightTextShort{\number\year}}
%% \pagerange{\mbox{}}%
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
%   \textit{Volume \j@volume} (\textit{\number\year}\/),
   \textit{Volume \j@volume} (\textit{\j@year}\/),
   \textit{Number \j@issue\ pp.\ \j@pStartPage--\j@pEndPage}%
   \hfill
 }}}%
%
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\Tutorial{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextEG{\p@ConfYear}}
\copyrightTextRunPag{\p@copyrightTextEG{\p@ConfYear}}
\gdef\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   EUROGRAPHICS \p@EGyear \hfill {\itshape Tutorial}}}}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
   \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\STAR{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextEG{\p@ConfYear}}
\copyrightTextRunPag{\p@copyrightTextEG{\p@ConfYear}}
\gdef\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   EUROGRAPHICS \p@EGyear \hfill {\itshape STAR -- State of The Art Report}}}}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
%
%---
\def\ShortPresentation{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextEG{\p@ConfYear}}
\copyrightTextRunPag{\p@copyrightTextEG{\p@ConfYear}}
\gdef\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   EUROGRAPHICS \p@EGyear{} / \p@ShortPresEditors{}
   \hfill {\itshape Short Presentations}}}}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
%
%---
\def\WsSubmission{%
\EGpagenumber%
\copyrightTextTitPag{\p@wsSubmissionText{\p@ConfName{} (\p@ConfYear})}
\copyrightTextRunPag{\p@wsSubmissionText{\p@ConfName{} (\p@ConfYear})}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   \parbox[t]{.7\textwidth}{\sloppy\raggedright
     \p@ConfName{} (\p@ConfYear), \ pp.\ \j@pStartPage--\j@pEndPage \\
     \p@ConfEditors{} \p@ConfEditorStrg}}%
 }%
}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
\def\WsPaper{%
\noEGpagenumber%
\copyrightTextTitPag{\p@copyrightTextEG{\p@ConfYear}}
\copyrightTextRunPag{\p@copyrightTextEG{\p@ConfYear}}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
 \def\@oddhead{\raisebox{\z@}[8pt][1pt]{\parbox{\textwidth}{\small
   \parbox[t]{.7\textwidth}{\sloppy\raggedright
%     \p@ConfName{} (\p@ConfYear), \ pp.\ \j@pStartPage--\j@pEndPage \\
     \p@ConfName{} (\p@ConfYear)\\
     \p@ConfEditors{} \p@ConfEditorStrg}}%
 }%
}
 \def\@oddfoot{{\tiny\raisebox{\z@}[8pt][1pt]{\parbox[t]{18pc}{\sloppy
  \p@copyrightTextTitPag}}}\hfill}
 \let\@evenhead=\@oddhead
 \let\@evenfoot=\@oddfoot
 \let\sectionmark=\EmptySectionmark
 \let\subsectionmark=\EmptySubsectionmark
}}
% -------
% Headings for Running Pages
%
\def\ps@headings{\let\@mkboth\markboth
%% DF
 \def\@oddhead{\small \hfill \itshape \@shortauthor\ / \@shorttitle
               \hfill \upshape \textrm{\p@EGpagenumber}}
%%- \def\@oddfoot{\mbox{}\hfill submitted to EUROGRAPHICS \EGyear}
 \def\@oddfoot{\tiny \p@copyrightTextRunPag \hfill\mbox{}}
%% DF
 \def\@evenhead{\small \upshape \textrm{\p@EGpagenumber} \hfill 
                \itshape \@shortauthor\ / \@shorttitle \hfill}
%%- \def\@evenfoot{submitted to EUROGRAPHICS '98\hfill\mbox{}}
 \def\@evenfoot{\mbox{}\hfill\tiny \p@copyrightTextRunPag}
 \def\sectionmark##1{}
 \def\subsectionmark##1{}}
%
\def\ps@myheadings{\let\@mkboth\@gobbletwo
%% DF
%%  \def\@oddhead{\small \it \rightmark \hfill \textrm{\thepage}}
 \def\@oddhead{\small \itshape \rightmark \hfill \mbox{}}
%%
 \def\@oddfoot{}
%% DF
%%  \def\@evenhead{\small \it \leftmark \hfill \textrm{\thepage}}
 \def\@evenhead{\small \itshape \leftmark \hfill \mbox{}}
%%
 \def\@evenfoot{}
 \def\sectionmark##1{}
 \def\subsectionmark##1{}}

%  ****************************************
%  *         TABLE OF CONTENTS, ETC.      *
%  ****************************************
%
\def\@pnumwidth{1.55em}
\def\@tocrmarg {2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{2}
\def\@undottedtocline#1#2#3#4#5{\ifnum #1>\c@tocdepth
 \else
  \vskip \z@ plus .2pt
  {\hangindent #2\relax
   \rightskip \@tocrmarg \parfillskip -\rightskip
   \parindent #2\relax \@afterindenttrue
   \interlinepenalty\@M \leavevmode
   \@tempdima #3\relax #4\nobreak \hfill \nobreak
%% DF
%    \hbox to\@pnumwidth{\hfil\reset@font\rm \ }\par}\fi}
   \hbox to\@pnumwidth{\hfil\reset@font\rmfamily \ }\par}\fi}
%%
\def\tableofcontents{\@restonecolfalse
 \if@twocolumn\@restonecoltrue\onecolumn\fi
 \section*{Contents} \@starttoc{toc}
 \if@restonecol\twocolumn\fi \par\vspace{11pt}}
\def\l@part#1#2{\addpenalty{-\@highpenalty}
 \addvspace{2.25em plus 1pt}
 \begingroup
  \parindent \z@ \rightskip \@pnumwidth
  \parfillskip -\@pnumwidth
%% DF
%   {\reset@font\normalsize\rm
  {\reset@font\normalsize\rmfamily
%%
   \leavevmode \hspace*{3pc}
   #1\hfil \hbox to\@pnumwidth{\hss \ }}\par
   \nobreak \global\@nobreaktrue
   \everypar{\global\@nobreakfalse\everypar{}}\endgroup}
\def\l@section#1#2{\addpenalty{\@secpenalty}
 \@tempdima 1.5em
 \begingroup
  \parindent \z@ \rightskip \@pnumwidth
%% DF
%   \parfillskip -\@pnumwidth \reset@font\rm \leavevmode
  \parfillskip -\@pnumwidth \reset@font\rmfamily \leavevmode
%%
  \advance\leftskip\@tempdima \hskip -\leftskip
  #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss \ }\par
 \endgroup}
\def\l@subsection{\@undottedtocline{2}{1.5em}{2.3em}}
\def\l@subsubsection{\@undottedtocline{3}{3.8em}{3.2em}}
\def\l@paragraph{\@undottedtocline{4}{7.0em}{4.1em}}
\def\l@subparagraph{\@undottedtocline{5}{10em}{5em}}
\def\listoffigures{\@restonecolfalse
 \if@twocolumn\@restonecoltrue\onecolumn\fi
 \section*{List of Figures\@mkboth{List of Figures}{List of Figures}}
 \@starttoc{lof} \if@restonecol\twocolumn\fi}
\def\l@figure{\@undottedtocline{1}{1.5em}{2.3em}}
\def\listoftables{\@restonecolfalse
 \if@twocolumn\@restonecoltrue\onecolumn\fi
 \section*{List of Tables\@mkboth{List of Tables}{List of Tables}}
 \@starttoc{lot} \if@restonecol\twocolumn\fi}
\let\l@table\l@figure

%  ****************************************
%  *             BIBLIOGRAPHY             *
%  ****************************************
%
% \def\thebibliography#1{\subsection*{References}
%  \addcontentsline{toc}{section}{References}
%  \list{\arabic{enumi}.}{\usecounter{enumi}%
%     \leftmargin 18pt
%     \labelwidth\leftmargin\advance\labelwidth-\labelsep
%     \let\makelabel\makeRLlabel
%     \itemsep 5.5pt plus 1pt minus 1pt%
% }
%  \parindent\z@
%  \parskip\z@ plus .1pt\relax
%  \def\newblock{\hskip .11em plus .33em minus .07em}
%  \sloppy\clubpenalty4000\widowpenalty4000
%  \sfcode`\.=1000\relax
% }
% 
\def\thebibliography#1{\subsection*{References}
 \addcontentsline{toc}{section}{References}
 \list{\@biblabel {\arabic{enumi}}}{\usecounter{enumi}%
    \settowidth \labelwidth {\@biblabel {#1}}
    \leftmargin \labelwidth \advance \leftmargin \labelsep 
    \let\makelabel\makeRLlabel
    \itemsep 5.5pt plus 1pt minus 1pt%
}
 \parindent\z@
 \parskip\z@ plus .1pt\relax
 \def\newblock{\hskip .11em plus .33em minus .07em}
 \sloppy\clubpenalty4000\widowpenalty4000
 \sfcode`\.=1000\relax
}
%
% \thebibliography=\long macro:
% #1->\section *{\refname \@mkboth {\MakeUppercase \refname }{\MakeUppercase \refname }}
% \list {\@biblabel {\@arabic \c@enumiv }}%
%    {\settowidth \labelwidth {\@biblabel {#1}}
%     \leftmargin \labelwidth \advance \leftmargin \labelsep 
%     \@openbib@code
%     \usecounter {enumiv}
%     \let \p@enumiv \@empty 
%     \renewcommand \theenumiv {\@arabic \c@enumiv }}
%  \sloppy \clubpenalty 4000 \@clubpenalty \clubpenalty \widowpenalty 
% 4000\sfcode `\.\@m .
% 
\let\endthebibliography=\endlist
% \def\@biblabel#1{\hspace*{-\labelsep}#1.}
\def\@biblabel#1{[#1]}
% \def\@cite#1{{\reset@font\scriptsize\mbox{$^{#1}$}}}
\def\@cite#1{\mbox{[#1]}}

%  ****************************************
%  *              THE INDEX               *
%  ****************************************
%
\newif\if@restonecol
\def\theindex{\section*{Index}
 \addcontentsline{toc}{section}{Index}
 \reset@font\footnotesize \parindent\z@ \parskip\z@ plus .1pt\relax
 \let\item\@idxitem}
\def\@idxitem{\par\hangindent 1em}
\def\subitem {\par\hangindent 1em \hspace*{1em}}
\def\subsubitem{\par\hangindent 1em\hspace*{1em}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par\vskip\baselineskip \relax}

%  ****************************************
%  *             FOOTNOTES                *
%  ****************************************
%
% \def\footnoterule{\kern-3\p@ \hrule width 12pc height \z@ \kern 3\p@}
%
% DF
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}

\renewcommand{\thefootnote}{\mbox{{$\fnsymbol{footnote}$}}}

% \def\@fnsymbol#1{\ifcase#1\or \dagger\or \ddagger\or
%  \S \or \P \or \|\or **\or \dagger\dagger \or
%  \ddagger\ddagger\or \S\S\or \P\P\or \|\|\else ***\fi\relax}
%
% DF 
\def\@fnsymbol#1{\ifcase#1\or \dagger\or \ddagger\or
 \mathchar "278\or \mathchar "27B\or \|\or **\or \dagger\dagger \or
 \ddagger\ddagger\or \mathchar "278\mathchar "278\or 
 \mathchar "27B\mathchar "27B\or \|\|\else ***\fi\relax}

\long\def\@makefntext#1{\parindent 1em\noindent 
 $^{\@thefnmark}$\hspace{4pt}#1}

% redefine footnotes to issue a warning...
\def\footnote{%
 \@warning{Footnotes are not permitted in CGF articles -- please remove}%
 \@ifnextchar[{\@xfootnote}{\stepcounter{\@mpfn}%
  \xdef\@thefnmark{\thempfn}\@footnotemark\@footnotetext}%
}
%    ****************************************
%    *               VERBATIM               *
%    ****************************************
% verbatim.sty
%%
%% This is file `verbatim.sty' generated
%% on <1991/11/8> with the docstrip utility (v1.1l test).
%%
%% The original source file was `verbatim.doc'.
%%
%%
%% Copyright (C) 1989,1990,1991 by Rainer Schoepf. All rights reserved.
%%
\def\addto@hook#1#2{#1\expandafter{\the#1#2}}
\newtoks\every@verbatim
\every@verbatim={}
\def\@makeother#1{\catcode`#112\relax}
\begingroup
 \catcode`\ =\active%
\gdef\@vobeyspaces{\catcode`\ \active\let \@xobeysp}%
\endgroup
\def\@xobeysp{\leavevmode\penalty\@M\ }
\newtoks\verbatim@line
\def\verbatim@startline{\verbatim@line{}}
\def\verbatim@addtoline#1{%
  \verbatim@line\expandafter{\the\verbatim@line#1}}
\def\verbatim@processline{\the\verbatim@line\par}
\def\verbatim@finish{\ifcat$\the\verbatim@line$\else
  \verbatim@processline\fi}
\begingroup
 \catcode`\`=\active
 \gdef\verbatim@font{\ttfamily \catcode96\active
   \def`{\leavevmode\kern\z@\char96 }}
\endgroup
\def\@verbatim{\the\every@verbatim
  \trivlist \item[]%
  \addvspace{5.5pt}
  \leftskip\@totalleftmargin\rightskip\z@
  \parindent\z@\parfillskip\@flushglue\parskip\z@
  \@@par
  \def\par{\leavevmode\null\@@par\penalty\interlinepenalty}%
  \obeylines
  \verbatim@font
  \let\do\@makeother \dospecials}
\def\verbatim{\@verbatim \frenchspacing\@vobeyspaces\verbatim@start}
\@namedef{verbatim*}{\@verbatim\verbatim@start}
\let\endverbatim=\endtrivlist
\expandafter\let\csname endverbatim*\endcsname =\endtrivlist
\def\comment{\@bsphack
             \let\do\@makeother\dospecials\catcode`\^^M\active
             \let\verbatim@startline\relax
             \let\verbatim@addtoline\@gobble
             \let\verbatim@processline\relax
             \let\verbatim@finish\relax
             \verbatim@}
\let\endcomment=\@esphack
\@ifundefined{vrb@catcodes}%
  {\def\vrb@catcodes{%
     \catcode`\!12\catcode`\[12\catcode`\]12}}{}
\begingroup
 \vrb@catcodes
 \lccode`\!=`\\ \lccode`\[=`\{ \lccode`\]=`\}
 \catcode`\~=\active \lccode`\~=`\^^M
 \lccode`\C=`\C
 \lowercase{%
    \gdef\verbatim@start#1{%
      \verbatim@startline
      \if\noexpand#1\noexpand~%
        \let\next\verbatim@
      \else \def\next{\verbatim@#1}\fi
      \next}%
    \gdef\verbatim@#1~{\verbatim@@#1!end\@nil}%
    \gdef\verbatim@@#1!end{%
       \verbatim@addtoline{#1}%
       \futurelet\next\verbatim@@@}%
    \gdef\verbatim@@@#1\@nil{%
       \ifx\next\@nil
         \verbatim@processline
         \verbatim@startline
         \let\next\verbatim@
       \else
         \def\@tempa##1!end\@nil{\toks@{##1}}%
           \@tempa#1\@nil
         \@temptokena{!end}%
         \edef\next{\noexpand\verbatim@test\the\toks@\noexpand~}%
       \fi \next}%
 \gdef\verbatim@test#1{%
           \let\next\verbatim@test
           \if\noexpand#1\noexpand~%
             \expandafter\verbatim@addtoline
               \expandafter{\the\@temptokena}%
             \verbatim@processline
             \verbatim@startline
             \let\next\verbatim@
           \else \if\noexpand#1
             \@temptokena\expandafter{\the\@temptokena#1}%
           \else \if\noexpand#1\noexpand[%
             \let\@tempc\@empty
             \let\next\verbatim@testend
           \else
             \expandafter\verbatim@addtoline
               \expandafter{\the\@temptokena}%
             \def\next{\verbatim@#1}%
           \fi\fi\fi
           \next}%
    \gdef\verbatim@testend#1{%
         \if\noexpand#1\noexpand~%
           \expandafter\verbatim@addtoline
             \expandafter{\the\@temptokena[}%
           \expandafter\verbatim@addtoline
             \expandafter{\@tempc}%
           \verbatim@processline
           \verbatim@startline
           \let\next\verbatim@
         \else\if\noexpand#1\noexpand]%
           \let\next\verbatim@@testend
         \else\if\noexpand#1\noexpand!%
           \expandafter\verbatim@addtoline
             \expandafter{\the\@temptokena[}%
           \expandafter\verbatim@addtoline
             \expandafter{\@tempc}%
           \def\next{\verbatim@!}%
         \else \expandafter\def\expandafter\@tempc\expandafter
           {\@tempc#1}\fi\fi\fi
         \next}%
    \gdef\verbatim@@testend{%
       \ifx\@tempc\@currenvir
         \verbatim@finish
         \edef\next{\noexpand\end{\@currenvir}%
                    \noexpand\verbatim@rescan{\@currenvir}}%
       \else
         \expandafter\verbatim@addtoline
           \expandafter{\the\@temptokena[}%
           \expandafter\verbatim@addtoline
             \expandafter{\@tempc]}%
         \let\next\verbatim@
       \fi
       \next}%
  \gdef\verbatim@rescan#1#2~{\if\noexpand~\noexpand#2~\else
        \@warning{Characters dropped after `\string\end{#1}'}\fi}}
\endgroup
\def\verbatiminput{\begingroup
  \@ifstar\sverbatim@input\verbatim@input}
\def\sverbatim@input#1{\@verbatim
  \@input{#1}\endtrivlist\endgroup\@doendpe}
\def\verbatim@input#1{\@verbatim
  \frenchspacing \@vobeyspaces
  \@input{#1}\endtrivlist\endgroup\@doendpe}
\begingroup
  \lccode`\~=`\^^M
  \lowercase{%
    \gdef\verb{\begingroup
      \verbatim@font
      \catcode`\^^M\active
      \def~{\endgroup\@latexerr{\string\verb\space command ended by
                                end of line.}\@ehc}%
    \let\do\@makeother \dospecials
    \@ifstar\@sverb{\@vobeyspaces \frenchspacing \@sverb}}}
\endgroup
\def\@sverb#1{%
  \catcode`#1\active
  \lccode`\~`#1%
  \lowercase{\let~\endgroup}%
  \leavevmode\null}
%\endinput
%%
%% End of file `verbatim.sty'.
%
% SFB modifications to the above
%
\begingroup
  \lccode`\~=`\^^M
  \lowercase{%
    \gdef\verb{\begingroup
%      \verbatim@font
      \ttfamily
      \catcode`\^^M\active
      \def~{\endgroup\@latexerr{\string\verb\space command ended by
                                end of line.}\@ehc}%
    \let\do\@makeother \dospecials
    \@ifstar\@sverb{\@vobeyspaces \frenchspacing \@sverb}}}
\endgroup
%
% end of myverb.sty
%
\def\verbatim@font{\reset@font\small\ttfamily}

%  ****************************************
%  *            MISCELLANEOUS             *
%  ****************************************
%
% DATE
\def\today{\number\day\ \ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December
 \fi \ \number\year}

%  ****************************************
%  *         HYPERLINK SUPPORT            *
%  ****************************************
%

\newif\ifpdf  % determine if we are running PDFLaTeX or not

\ifx\pdfoutput\undefined \pdffalse % we are not running PDFLaTeX 
\else 
  \pdfoutput=1 % we are running PDFLaTeX 
  \pdftrue 
\fi 

\newcommand{\pdf@Subject}{Vol \j@volume, No \j@issue}
\newcommand{\pdfSubject}[1]{\renewcommand{\pdf@Subject}{#1}}

\newif\ifproduce@print
\produce@printtrue

\newcommand{\printedVersion}{\global\produce@printtrue}

\newcommand{\printed@Version}{}

\newcommand{\electronicVersion}{\global\produce@printfalse}

\ifpdf 
  \newcommand{\electronic@Version}{%
   \usepackage[pdftex,
    pagebackref=true,
    colorlinks,linkcolor=blue,citecolor=blue,urlcolor=blue,
    bookmarks=false,
    pdfpagemode=UseNone,
    pdftitle={\@shorttitle},
    pdfauthor={\@shortauthor},
    pdfsubject={\pdf@Subject},
    pdfkeywords={Computer Graphics Forum, EUROGRAPHICS}]{hyperref}
  }
\else 
  \newcommand{\electronic@Version}{%
   \usepackage[dvips,
    pagebackref=true,
    colorlinks,linkcolor=blue,citecolor=blue,urlcolor=blue,
    bookmarks=false,
    pdfpagemode=UseNone,
    pdftitle={\@shorttitle},
    pdfauthor={\@shortauthor},
    pdfsubject={\pdf@Subject},
    pdfkeywords={Computer Graphics Forum, EUROGRAPHICS}]{hyperref}
  }
\fi

\newcommand{\PrintedOrElectronic}{%
\ifproduce@print \printed@Version \else \electronic@Version \fi
}

%  ****************************************
%  *           INITIALIZATION             *
%  ****************************************
%
\pagenumbering{arabic}
\if@twocolumn
 \twocolumn
 \sloppy
\fi
\ifSFB@galley
 \ps@empty
 \raggedbottom
\else
 \ps@headings
 \flushbottom
\fi
\frenchspacing

\AtBeginDocument{{\setcounter{page}{\j@pStartPage}}}

\gdef\jobname@lbl{\jobname.lbl}%
\InputIfFileExists{\jobname.lbl}{}{}%
\def\label@out{9}%
\immediate\openout\label@out=\jobname@lbl%

\AtEndDocument{%
\immediate\write\label@out{\string\pEndPage{\thepage}\string\endinput}%
}

\JournalSubmission

% ---------------------
%% DF add-on to set variable parameters for the publication run
%
% for journal papers - if known at compile time
% 
% \issue{3}	   % the issue to appear in (year + volume are computed)
% \pStartPage{x}   % set starting page

\endinput
%%%%%%%%%%%%%%%%%%%%%% end of egpubl.cls %%%%%%%%%%%%%%%%%%%%%%
