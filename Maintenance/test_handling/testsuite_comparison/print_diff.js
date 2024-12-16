/* Author: Maxime Gimeno <maxime.gimeno@gmail.com> */
/***
Need difflib.js

Input : 2 strings

Output: res. Contains the difference between the two input strings.

According to the result of difflib, will output :
- replace : -base, +newtext (the two input differs, but none are empty)
- insert: +newtext (the first input is empty)
- delete: -base (the second input is empty)
***/

function print_diff(base, newtext){
  // create a SequenceMatcher instance that diffs the two sets of lines
  var sm = new difflib.SequenceMatcher(base, newtext);
  // get the opcodes from the SequenceMatcher instance
  // opcodes is a list of 3-tuples describing what changes should be made to the base text
  // in order to yield the new text
  var opcodes = sm.get_opcodes();
  var res=["",""];
  for (var idx = 0; idx < opcodes.length; idx++) {
    code = opcodes[idx];
    change = code[0];
    var b = code[1];
    var be = code[2];
    var n = code[3];
    var ne = code[4];
    if(newtext.charAt(newtext.length-1) !== "y"
       && change != "equal"){
        res[0]=base.charAt(base.length-1);
        res[1]=newtext.charAt(newtext.length-1);
        if (res[0]=="w" && res[1]=="t")
        {
          res[0]="";
          res[1]="";
        }
      }
      //else if(change == "insert") {
      //  res=newtext.charAt(newtext.length-1);
      //}
      //else if(change == "delete") {
      //  res="-";
      //}
    //}
  }
  return res;
}
