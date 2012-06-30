if ( NOT FIND_OPENGL_WRAPPER )

  set ( FIND_OPENGL_WRAPPER 1 )
  
  set(SAVED_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} )
  
  set(CMAKE_MODULE_PATH ${ORIG_CMAKE_MODULE_PATH} )
  
  find_package(OpenGL)
  
  if ( OPENGL_FOUND AND CGAL_APPLE_LEOPARD )
  
    if ( BUILD_SHARED_LIBS )
    
      uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS
                        "-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib"
                        )
                        
    else()
    
      uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS
                        "-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib"
                        )
    endif()
                      
                      
  endif()
  
  set(CMAKE_MODULE_PATH ${SAVED_CMAKE_MODULE_PATH} )

endif()
