set CGAL_LIB=%1:\uni\libs\CGAL-3.7\bin\lib
set CGAL_INCLUDE=%1:\uni\libs\CGAL-3.7\include
set CGAL_COMPILER_CONFIG=%1:\uni\libs\CGAL-3.7\bin\include
set CGAL_GMP=%1:\uni\libs\CGAL-3.7\auxiliary\gmp\include
set CGAL_GMP_LIB=%1:\uni\libs\CGAL-3.7\auxiliary\gmp\lib
set BOOST_INCLUDE=%1:\uni\libs\boost_1_45_0
set BOOST_LIB=%1:\uni\libs\boost_1_45_0\stage\lib
set QT_INCLUDE=%1:\uni\libs\qt\include
set QT_LIB=%1:\uni\libs\qt\lib\
set PATH=%1:\uni\libs\CGAL-3.7\auxiliary\gmp\lib;%1:\uni\libs\qt\bin;%PATH%
"C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\devenv.exe"