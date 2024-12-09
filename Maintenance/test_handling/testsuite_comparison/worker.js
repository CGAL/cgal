

onmessage=function(e){
    importScripts('difflib.js',  'diff_testsuites.js');
    importScripts('fill_empty_lines.js','print_diff.js');
  var res_arrays = diff_testsuites(e.data[0],e.data[1]);
  var myArray = res_arrays[0];
  var diffArray = res_arrays[1];
    for(var i=0; i<diffArray.length; i++){
      postMessage(["diffPlatforms","<tr><td>"+ diffArray[i]+"</td></tr>"]);
    }
    postMessage(["diffPlatforms", "</table>"]);
  //pass over My Array to equalize columns length

    var max_length=myArray[0].length;
    for( i = 1; i< myArray.length; ++i){
      var length = myArray[i].length;
        if(length === max_length){ continue; }
      for(var j=2; j<length;j++){
      var base_line = myArray[i][j];
      var max_line = myArray[0][j];
      var base_identifier=base_line.split("||")[0];
      var max_identifier=max_line.split("||")[0];
      if(base_identifier !== max_identifier) {
          myArray[i].splice(j, 0, "<td></td>");
      }
    }
    var minsize = myArray[i].length;

    for(k = minsize; k < myArray[0].length; k++)
    {
      myArray[i].push("<td></td>");
    }
  }

  for(i=1; i<myArray.length; i++){
    postMessage(["namesTable","<tr><td>"+i+"</td><td>"+ myArray[i][0]+"</td></tr>"]);
  }
  postMessage(["namesTable", "</table>"]);

  postMessage(["testTable", "<tr><td style='width: 25px;'>Platform:  </td>"]);
  for(i=1; i<myArray.length; i++)
  {
    postMessage(["testTable", "<td style='width: 25px; text-align: right;'>" + i + "</td>"]);
  }
  postMessage(["testTable", "</tr>"]);
  for (var j=2; j<myArray[0].length; j++) {
    var p =  myArray[0][j];
    postMessage(["testTable", "<tr><td style='width: 25px;'><a class=\"package_name\" href=\"#"+p+"\" name=\"" + p + "\">" + p + "</a></td>"]);
    for(i=1; i<myArray.length; i++)
    {
      var sp = myArray[i][j].split("||");
      if(sp.length === 1){
        postMessage(["testTable", myArray[i][j]]);
      } else{
        postMessage(["testTable", sp[1]]);
      }
    }
    postMessage(["testTable", "</tr>"]);
  }
  postMessage(["testTable", "</table></body>"]);
  postMessage(["finished"]);
}
