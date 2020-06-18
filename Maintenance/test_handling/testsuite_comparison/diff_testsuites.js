/* Author: Maxime Gimeno <maxime.gimeno@gmail.com> */

/***
  Needs difflib.js, fill_empty_lines.js and print_diff.js to work.

  Input : 2 strings of the form "CGAL-M.m.mm-I[c]-XXX"

  Output: the list of differences between the testsuites.

***/

//v1.7
function diff_testsuites(baseTest, newTest){
  var URL_suff='https://cgal.geometryfactory.com/~mgimeno/testsuite_comparison/list_of_suffixes.txt';
  var URL_testsuite='https://cgal.geometryfactory.com/CGAL/Members/testsuite/';
  //get the list of suffixes
  var xhr = new XMLHttpRequest();
  xhr.open('GET', URL_suff, false);
  xhr.send(null);
  var tmp=xhr.responseText;
  var suffixes=tmp.split("\n");
  suffixes.sort();
  suffixes.reverse();
  var myArray = new Array();
    
    //contains the diff of existing platforms.
    
  var diffArray = new Array();
    
  var init = false;
  for(s = 0; s < suffixes.length; s++) {
    var new_column = new Array();
    xhr.open('GET', URL_testsuite+baseTest+'/'+suffixes[s], false);
    xhr.send(null);
    var base_exists = (xhr.status === 200);
    var base=xhr.responseText;
    xhr.open('GET', URL_testsuite+newTest+'/'+suffixes[s], false);
    xhr.send(null);
    var new_exists = (xhr.status === 200)
      
    if(base_exists && !new_exists)
    {
        diffArray.push("-"+suffixes[s]);
        continue;
    }
    if(!base_exists && new_exists)
    {
        diffArray.push("+"+suffixes[s]);
        continue;
    }
    if(!base_exists && !new_exists)
    {
        continue;
    }
   
    var newtext=xhr.responseText;
    var sp_base=base.split("\n");
    sp_base.sort();
    var sp_newtext=newtext.split("\n");
    sp_newtext.sort();
    addMissingLines(sp_base, sp_newtext);
    if(!init)
    {
      var first_column = new Array();
        first_column.push("Platform")
      for(i=0; i< sp_base.length; i++){
          if(sp_base[i] !== ""){
            first_column.push(sp_base[i].substr(0, sp_base[i].length-2));
          } else {
            first_column.push(sp_newtext[i].substr(0, sp_newtext[i].length-2));
          }
      }
      myArray.push(first_column);
      init = true;
    }
    var fragments = suffixes[s].split("_");
    fragments.shift();
    var name = fragments.join("_");
    name = name.replace('.txt', '');
    if(name !== ""){
      new_column.push(name.replace(fragments[0]+"_", ''));
      for(i=0; i< sp_base.length; i++){
        var broken = false;
        var res = print_diff(sp_base[i], sp_newtext[i]);
          var compensator=0;
          if(sp_base[i] !== ""){
            while(sp_base[i].substr(0, sp_base[i].length-2) !== first_column[i+compensator]){
                if(compensator >10){
                  broken=true;
                  break;
                }
              compensator++;
            }
          }
          else{
            while(sp_newtext[i].substr(0, sp_newtext[i].length-2) !== first_column[i+compensator]){
                if(compensator >10){
                  broken=true;
                  break;
                }
              compensator++;
            }
          }
          if(broken)
          {
            continue;
          }
        var new_line=first_column[i+compensator]+"||"+"<td style='width: 25px; text-align: right;";
        var result="";
        if(res[0]===""){
          if(res[1] !==""){
            result+="+";
          }
        }
        else {
          result+="<a href=\""+URL_testsuite+baseTest+"/"+first_column[i+compensator]+"/TestReport_"+name+".gz\">"+res[0]+"</a>\/";
        }

        if(res[1]===""){
            if(res[0]!==""){
              result+="-";
            }
        }
        else {
            result+="<a href=\""+URL_testsuite+newTest+"/"+first_column[i+compensator]+"/TestReport_"+name+".gz\">"+res[1]+"</a>";
        }
        if(res[1] === "w"){
          new_line+=" background-color: rgb(100%,100%,50%)'> "+result;
        } else if(res[1]=== "r"){
          new_line+=" background-color: rgb(65%,65%,100%)'> "+result;
        } else if(res[1]=== "n"){
          new_line+=" background-color: rgb(100%,50%,50%)'> "+result;
        }
        else if(res[1]==="" && res[0]!==""){
          new_line+=" background-color: rgb(50%,25%,75%)'>"+result;
        }
        else{
          new_line+="'>";
        }
        new_line+="</td>";
        new_column.push(new_line);
      }
      myArray.push(new_column);
    }
  }
  
  return [myArray, diffArray];
}
