macro( readlines file lines )

  file(READ ${file} ${file}_contents )
  
  string( REPLACE "\n" ";" ${lines} ${${file}_contents} )
  
endmacro( readlines file lines )
 
