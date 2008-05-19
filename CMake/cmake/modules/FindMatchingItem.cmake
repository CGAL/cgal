macro( find_matching_item list regex output )

  set( ${output} "NOTFOUND" )
  
  foreach( ${list}_row__ ${${list}} )
  
    if ( ${list}_row__ MATCHES ${regex} )
    
      set( ${output} "${${list}_row__}" )
      
    endif( ${list}_row__ MATCHES ${regex} )
    
  endforeach( ${list}_row__ ${${list}} )
  
endmacro( find_matching_item list regex output )

