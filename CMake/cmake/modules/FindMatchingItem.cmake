macro( find_matching_item list regex output )

  set( ${output} "NOTFOUND" )
  
  foreach( ${list}_lfmi_idx ${${list}} )
  
    list( GET ${list}_lfmi_idx ${list} ${list}_lfmi_row )
    
    if ( ${list}_lfmi_row MATCHES ${regex} )
    
      set( ${output} "${${list}_lfmi_row}" )
      
    endif( ${list}_lfmi_row MATCHES ${regex} )  
    
  endforeach( ${list}_lfmi_idx ${${list}} )
  
endmacro( find_matching_item list regex output )
