% EG:	bibliographies style for EUROGRAPHICS publications
%    derived from BI_NUM (which originates from ACM Transactions bibliography style)
%            and  ALPHA  (with \etalchar{*}).
%
%    modifications:
%            - labels for bibitem are computed according to alpha.bst
%            - names are included in {\sc ... :}
%            - no period after line listing the authors
%            - no comma within a name
%            - names are separated by comma
%            - new database item STANDARD; its like MANUAL but an optional
%              VERSION field is added after TITLE
%
%    History
%     7.08.90 (DF) above modifications
%    20.11.90 (DF) add title handling for GARTICLE: format.german.title
%    24.05.91 (DF) port to BibTeX 0.99 (using ACM.BST as the source to start from)
%    10.04.93 (DF) fix semicolon between editors and "(Ed[s].)": scapify;
%                   use \uppercase{longest.label} to reserve enough space
%    11.09.03 (DF) adapt to LateX 2e conventions
%
%
% ACM Transactions bibliography style (24-Jan-88 version)
%    A lot like abbrv.bst, but names come out "Last, initials", and in \sc.
%    Some dates are parenthesized.
%
%    History
%     2/ 6/86    (HWT)    Original version, by Howard Trickey.
%     3/ 5/86    (HWT)    Put in pp. everywhere but articles, as per ACM style.
%     1/24/88    (OP&HWT) Updated for BibTeX version 0.99a, Oren Patashnik;
%            corrected the abbreviations to "Mar." and "Sept.";
%            THIS `acm' VERSION DOES NOT WORK WITH BIBTEX 0.98i.
 
ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    version
    volume
    year
  }
  {}
  { label extra.label sort.label }
 
INTEGERS { output.state before.all mid.sentence after.sentence after.block }
 
FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}
 
STRINGS { s t }
 
FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}
 
FUNCTION {output.nonnull.noperiod}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}
 
FUNCTION {output.noperiod}
{ duplicate$ empty$
    'pop$
    'output.nonnull.noperiod
  if$
}
 
FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}
 
FUNCTION {output.check.noperiod}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull.noperiod
  if$
}
 
FUNCTION {output.bibitem}
{ newline$
  "\bibitem[" write$
  label write$
  "]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}
 
FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}
 
FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}
 
FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}
 
FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}
 
FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}
 
FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}
 
FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}
 
FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}
 
FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}
 
FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\emph{" swap$ * "}" * }
  if$
}
 
FUNCTION {emphasizeic}
{ duplicate$ empty$
    { pop$ "" }
    { "\emph{" swap$ * "}" * }
  if$
}
 
FUNCTION {scapify}
{ duplicate$ empty$
    { pop$ "" }
    { "\textsc{" swap$ * "}" * }
  if$
}
 
INTEGERS { nameptr namesleft numnames }
 
FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
    { s nameptr "{vv~}{ll}{ jj}{ f.}" format.name$ 't :=
    nameptr #1 >
      { namesleft #1 >
          { ", " * t * }
          { t "others" =
              { ", et~al." * }
              { ", " * t * }
            if$
          }
        if$
      }
      't
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 - 'namesleft :=
    }
  while$
}
 
% For names inside entries (e.g., editors of an "In ...");
% this is exactly ABBRV.BST's `format.names' function.
FUNCTION {format.innames}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
    { s nameptr "{vv~}{ll}{ jj}{ f.}" format.name$ 't :=
    nameptr #1 >
      { namesleft #1 >
          { ", " * t * }
          { numnames #2 >
              { "," * }
              'skip$
            if$
            t "others" =
              { " et~al." * }
              { ", " * t * }
            if$
          }
        if$
      }
      't
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 - 'namesleft :=
    }
  while$
}
 
FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names scapify ":" * }
  if$
}
 
FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names scapify
      editor num.names$ #1 >
        { " (Eds.):" * }
        { " (Ed.):" * }
      if$
    }
  if$
}
 
FUNCTION {format.ineditors}
{ editor empty$
    { "" }
    { editor format.innames
      editor num.names$ #1 >
        { ", (Eds.)" * }
        { ", (Ed.)" * }
      if$
    }
  if$
}
 
FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}
 
FUNCTION {format.german.title}
{ title field.or.null
}
 
FUNCTION {n.dashify}
{ 't :=
  ""
  { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            { { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}
 
FUNCTION {format.date}
{ year empty$
    { month empty$
        { "" }
        { "there's a month but no year in " cite$ * warning$
          month
        }
      if$
    }
    { month empty$
        'year
        { month " " * year * }
      if$
    }
  if$
}
 
FUNCTION {format.btitle}
{ title emphasize
}
 
FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}
 
FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}
 
FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "vol.~" volume *
      series empty$
        'skip$
        { " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}
 
FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { output.state mid.sentence =
            { "no.~" }
            { "No.~" }
          if$
          number *
          series empty$
            { "there's a number but no series in " cite$ * warning$ }
            { " in " * series * }
          if$
        }
      if$
    }
    { "" }
  if$
}
 
FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
        { edition "l" change.case$ "~ed." * }
        { edition "t" change.case$ "~ed." * }
      if$
    }
  if$
}
 
FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages n.dashify }
  if$
}
 
INTEGERS { multiresult }
 
FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
  { multiresult not  t empty$ not  and }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}
 
FUNCTION {format.pp.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
        { "pp.~" pages n.dashify * }
        { "p.~" pages * }
      if$
    }
  if$
}
 
FUNCTION {format.journal.vol.num.date}
{ journal empty$
    { "empty journal in " cite$ * warning$
      ""
    }
    { journal
      volume empty$
        'skip$
        { " " * volume * }
      if$
      number empty$
        'emphasizeic
        { emphasize ", " * number * }
      if$
      year empty$
        { "empty year in " cite$ * warning$ }
        { " (" * format.date * ")" * }
      if$
    }
  if$
}
 
FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pp.pages
    { type empty$
        { "ch.~" chapter * }
        { type "l" change.case$ chapter tie.or.space.connect }
      if$
      pages empty$
        'skip$
        { ", " * format.pp.pages * }
      if$
    }
  if$
}
 
FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { "In " booktitle emphasize *
      editor empty$
        'skip$
        { ", " * format.ineditors * }
      if$
    }
  if$
}
 
% The proceedings title (it's on the stack) gets an (address, date) appended
FUNCTION {format.proc.date}
{ duplicate$ empty$
    { pop$ "" }
    { year empty$
        { "empty year in " cite$ * warning$
          address empty$
            'emphasize
            { emphasizeic
              " (" * address * ")" *
            }
          if$
        }
        { emphasizeic
          " (" *
          address empty$
            'skip$
            { address * ", " * }
          if$
          format.date *
          ")" *
        }
      if$
    }
  if$
}
 
FUNCTION {format.in.proc.date}
{ booktitle empty$
    { "" }
    { "In " booktitle format.proc.date * }
  if$
}
 
FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}
 
FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}
 
FUNCTION {format.tr.number}
{ type empty$
    { "Tech. Rep." }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}
 
FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
        { "need key or journal for " cite$ * " to crossref " * crossref *
          warning$
          ""
        }
        { "In \emph{" journal * "}" * }
      if$
    }
    { "In " key * }
  if$
  " \cite{" * crossref * "}" *
}
 
FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
        'skip$
        { editor #2 "{vv~}{ll}{ jj}{ f.}" format.name$ "others" =
            { " et~al." * }
            { " and " * editor #2 "{vv~}{ll}" format.name$ * }
          if$
        }
      if$
    }
  if$
}
 
FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Vol.~" volume *
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { series empty$
            { "need editor, key, or series for " cite$ * " to crossref " *
              crossref * warning$
              "" *
            }
            { "\emph{" * series * "}" * }
          if$
        }
        { key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}
 
FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { booktitle empty$
            { "need editor, key, or booktitle for " cite$ * " to crossref " *
              crossref * warning$
              ""
            }
            { "In \emph{" booktitle * "}" * }
          if$
        }
        { "In " key * }
      if$
    }
    { "In " format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}
 
FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.title "title" output.check.noperiod
  new.block
  crossref missing$
    { format.journal.vol.num.date output
      format.pages output
    }
    { format.article.crossref output.nonnull
      format.pp.pages output
    }
  if$
  new.block
  note output
  fin.entry
}
 
FUNCTION {garticle}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.german.title "title" output.check.noperiod
  new.block
  crossref missing$
    { format.journal.vol.num.date output
      format.pages output
    }
    { format.article.crossref output.nonnull
      format.pp.pages output
    }
  if$
  new.block
  note output
  fin.entry
}
 
FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  new.block
  format.btitle "title" output.check.noperiod
  format.edition output
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.date "year" output.check
  new.block
  note output
  fin.entry
}
 
FUNCTION {booklet}
{ output.bibitem
  format.authors output
  new.block
  format.title "title" output.check.noperiod
  howpublished address new.block.checkb
  howpublished output
  address output
  format.date output
  new.block
  note output
  fin.entry
}
 
FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  new.block
  format.btitle "title" output.check.noperiod
  format.edition output
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.date "year" output.check
  format.chapter.pages "chapter and pages" output.check
  new.block
  note output
  fin.entry
}
 
FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.title "title" output.check.noperiod
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.edition output
      format.bvolume output
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
      format.date "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull }
  if$
  format.chapter.pages output
  new.block
  note output
  fin.entry
}
 
FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.title "title" output.check.noperiod
  new.block
  crossref missing$
    { format.in.proc.date "booktitle" output.check
      format.ineditors output
      format.bvolume output
      format.number.series output
      organization output
      publisher output
    }
    { format.incoll.inproc.crossref output.nonnull }
  if$
  format.pp.pages output
  new.block
  note output
  fin.entry
}
 
FUNCTION {conference} { inproceedings }
 
FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization scapify ":" * output }
    { format.authors output.nonnull }
  if$
  new.block
  format.btitle "title" output.check.noperiod
  format.edition output
  author empty$
    { address new.block.checka }
    { organization address new.block.checkb
      organization output
    }
  if$
  address output
  format.date output
  new.block
  note output
  fin.entry
}
 
FUNCTION {standard}
{ output.bibitem
  author empty$
    { organization scapify ":" * output }
    { format.authors output.nonnull }
  if$
  new.block
  format.btitle "title" output.check.noperiod
  version output
  format.edition output
  author empty$
    { address new.block.checka }
    { organization address new.block.checkb
      organization output
    }
  if$
  address output
  format.date output
  new.block
  note output
  fin.entry
}
 
FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.btitle "title" output.check.noperiod
  new.block
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}
 
FUNCTION {misc}
{ output.bibitem
  format.authors output
  title howpublished new.block.checkb
  format.title output.noperiod
  howpublished new.block.checka
  howpublished output
  format.date output
  new.block
  note output
  fin.entry
  empty.misc.check
}
 
FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.btitle "title" output.check.noperiod
  new.block
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}
 
FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization scapify ":" * output }
    { format.editors output.nonnull }
  if$
  new.block
  title format.proc.date "title" output.check
  format.bvolume output
  format.number.series output
  editor empty$
    'skip$
    { organization output }
  if$
  publisher output
  new.block
  note output
  fin.entry
}
 
FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.btitle "title" output.check.noperiod
  new.block
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}
 
FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.title "title" output.check.noperiod
  new.block
  note "note" output.check
  format.date output
  fin.entry
}
 
FUNCTION {default.type} { misc }
 
MACRO {jan} {"Jan."}
 
MACRO {feb} {"Feb."}
 
MACRO {mar} {"Mar."}
 
MACRO {apr} {"Apr."}
 
MACRO {may} {"May"}
 
MACRO {jun} {"June"}
 
MACRO {jul} {"July"}
 
MACRO {aug} {"Aug."}
 
MACRO {sep} {"Sept."}
 
MACRO {oct} {"Oct."}
 
MACRO {nov} {"Nov."}
 
MACRO {dec} {"Dec."}
 
MACRO {acmcg} {"Computer Graphics"}

MACRO {acmcs} {"ACM Comput. Surv."}
 
MACRO {acta} {"Acta Inf."}
 
MACRO {c&g} {"Computers \& Graphics"}

MACRO {cacm} {"Commun. ACM"}
 
MACRO {cgforum} {"Computer Graphics Forum"}

MACRO {computer} {"IEEE Computer"}

MACRO {ibmjrd} {"IBM J. Res. Dev."}
 
MACRO {ibmsj} {"IBM Syst.~J."}
 
MACRO {ieeecg} {"IEEE CG\&A"}

MACRO {ieeese} {"IEEE Trans. Softw. Eng."}
 
MACRO {ieeetc} {"IEEE Trans. Comput."}
 
MACRO {ieeetcad} {"IEEE Trans. Comput.-Aided Design Integrated Circuits"}
 
MACRO {ipl} {"Inf. Process. Lett."}
 
MACRO {jacm} {"J.~ACM"}
 
MACRO {jcss} {"J.~Comput. Syst. Sci."}
 
MACRO {scp} {"Sci. Comput. Programming"}
 
MACRO {sicomp} {"SIAM J. Comput."}
 
MACRO {tocs} {"ACM Trans. Comput. Syst."}
 
MACRO {tods} {"ACM Trans. Database Syst."}
 
MACRO {tog} {"ACM TOG"}
 
MACRO {tosem} {"ACM Trans. Softw. Eng. Method."}

MACRO {toms} {"ACM Trans. Math. Softw."}
 
MACRO {toois} {"ACM Trans. Office Inf. Syst."}
 
MACRO {toplas} {"ACM Trans. Program. Lang. Syst."}
 
MACRO {tcs} {"Theoretical Comput. Sci."}
 
READ
 
FUNCTION {sortify}
{ purify$
  "l" change.case$
}
 
INTEGERS { len }
 
FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}
 
%% DF 24.5.91  
%%   insert from ALPHA.bst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INTEGERS { et.al.char.used }

FUNCTION {initialize.et.al.char.used}
{ #0 'et.al.char.used :=
}

EXECUTE {initialize.et.al.char.used}

FUNCTION {format.lab.names}
{ 's :=
  s num.names$ 'numnames :=
  numnames #1 >
    { numnames #4 >
        { #3 'namesleft := }
        { numnames 'namesleft := }
      if$
      #1 'nameptr :=
      ""
      { namesleft #0 > }
        { nameptr numnames =
            { s nameptr "{ff }{vv }{ll}{ jj}" format.name$ "others" =
                { "{\etalchar{*}}" *
                  #1 'et.al.char.used :=
                }
                { s nameptr "{v{}}{l{}}" format.name$ * }
              if$
            }
            { s nameptr "{v{}}{l{}}" format.name$ * }
          if$
          nameptr #1 + 'nameptr :=
          namesleft #1 - 'namesleft :=
        }
      while$
      numnames #4 >
        { "{\etalchar{*}}" *
          #1 'et.al.char.used :=
        }
        'skip$
      if$
    }
    { s #1 "{v{}}{l{}}" format.name$
      duplicate$ text.length$ #2 <
        { pop$ s #1 "{ll}" format.name$ #3 text.prefix$ }
        'skip$
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        { key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            { key #3 text.prefix$ }
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.key.organization.label}
{ author empty$
    { key empty$
        { organization empty$
            { cite$ #1 #3 substring$ }
            { "The " #4 organization chop.word #3 text.prefix$ }
          if$
        }
        { key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.organization.label}
{ editor empty$
    { key empty$
        { organization empty$
            { cite$ #1 #3 substring$ }
            { "The " #4 organization chop.word #3 text.prefix$ }
          if$
        }
        { key #3 text.prefix$ }
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.organization.label
        { type$ "manual" =
            'author.key.organization.label
            'author.key.label
          if$
        }
      if$
    }
  if$
  duplicate$
  year field.or.null purify$ #-1 #2 substring$
  *
  'label :=
  year field.or.null purify$ #-1 #4 substring$
  *
  sortify 'sort.label :=
}

%%   end insert from ALPHA.bst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
    { nameptr #1 >
        { "   " * }
        'skip$
      if$
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
        { "et al" * }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}
 
FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}
 
FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}
 
FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}
 
FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
        { key empty$
            { "to sort, need author, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}
 
FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
        { key empty$
            { "to sort, need editor, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}
 
FUNCTION {presort}
{ calc.label		%% DF 24.5.91:  inserted from ALPHA
  sort.label		%% DF 24.5.91:  inserted from ALPHA
  "    "		%% DF 24.5.91:  inserted from ALPHA
  *			%% DF 24.5.91:  inserted from ALPHA
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.organization.sort
        { type$ "manual" =
          type$ "standard" =
          or
            'author.organization.sort
            'author.sort
          if$
        }
      if$
    }
  if$
  *			%% DF 24.5.91:  inserted from ALPHA
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}
 
ITERATE {presort}
 
SORT
 
STRINGS { longest.label last.sort.label next.extra }

INTEGERS { longest.label.width last.extra.num }
 
FUNCTION {initialize.longest.label}	%% DF 24.5.91:  whole function replaced from ALPHA
{ "" 'longest.label :=
  #0 int.to.chr$ 'last.sort.label :=
  "" 'next.extra :=
  #0 'longest.label.width :=
  #0 'last.extra.num :=
}

%% DF 24.5.91
%%   FUNCTION {longest.label.pass} ...
%% replaced by

FUNCTION {forward.pass}
{ last.sort.label sort.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      sort.label 'last.sort.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  label extra.label * 'label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
  extra.label 'next.extra :=
}

%% end replace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EXECUTE {initialize.longest.label}
 
%% DF 24.5.91: 
%%   replace  ITERATE {longest.label.pass}
%% by:

ITERATE {forward.pass}

REVERSE {reverse.pass}

%% end replacement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
FUNCTION {begin.bib}
{ et.al.char.used						%% insert from ALPHA
    { "\newcommand{\etalchar}[1]{$^{#1}$}" write$ newline$ }	%% insert from ALPHA
    'skip$							%% insert from ALPHA
  if$								%% insert from ALPHA
  preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{\uppercase{"  longest.label  * "}}" * write$ newline$
}
 
EXECUTE {begin.bib}
 
EXECUTE {init.state.consts}
 
ITERATE {call.type$}
 
FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
 
EXECUTE {end.bib}
