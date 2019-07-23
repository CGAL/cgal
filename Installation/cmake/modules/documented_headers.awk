#! /bin/awk

function rec_wrap(str)
{
    matches=""
    return rec_func(str)
}
function rec_func(str2)
{
    where=match(str2, /# *include *(&lt;|[<"])(CGAL\/[^>&"]*)([>"]|&gt;)/, arr)
    if(where!=0) {
        matches=(matches "\n" arr[2])
        rec_func(substr(str2, RSTART+RLENGTH, length(str2)))
    }
    return matches
}

    { matches=rec_wrap($0); if(length(matches) > 0) print matches}
