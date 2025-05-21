/* Author: Maxime Gimeno <maxime.gimeno@gmail.com> */

/***
  This is designed to treat lines from CGAL testsuites results files,
  of the form Package_name [r/y/w/n].

  Input : 2 arrays of string, sorted by alphabetical order.

  Output: the arrays, alphabetically sorted, of the same size, filled with
  empty strings.

  Short: Equalizes the sizes of the two input arrays by adding empty strings.

  Detailed: for each element of the smaller input,
  if base[i] != newtest[i] (not taking the last char into account),
  adds an empty string in the input that is missing this entry,
  based on the alphabetical order. Once it is done, will fill the smaller array
  with empty strings if necessary.


***/
function addMissingLines(base, newtext){
  for(i=0, j=0; i<Math.min(base.length, newtext.length); i++, j++)
  {
    if(base[i].substr(0, base[i].length-2) != newtext[i].substr(0,newtext[i].length-2) ) {
      var tmp = [base[i], newtext[i]];
      tmp.sort();
      if(tmp[0] == base[i]){
        newtext.splice(i, 0, "");
      } else {
        base.splice(i, 0, "");
      }
    }
  }
  var minsize = base.length;
  var maxsize = newtext.length;
  var to_fill = base;
  if(base.length>newtext.length){
    minsize = newtext.length;
    maxsize = base.length;
    to_fill = newtext;
  }
  for(i = minsize; i < maxsize; i++)
  {
    to_fill.push("");
  }
}
