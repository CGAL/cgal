macro( readlines file lines )

  if ( EXISTS ${file} )
  
    file(READ ${file} ${lines}_file_contents )
    
    string( REPLACE "\n" ";" "${lines}" "${${lines}_file_contents}" )
    
  else()
  
    set( ${lines} "NOTFOUND" )
    
  endif()
  
endmacro( readlines file lines )
 
