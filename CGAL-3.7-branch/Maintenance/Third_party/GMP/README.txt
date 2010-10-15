This README contains instructions to create compiled binaries for gmp and mpfr for Visual Studio
with the names expected by the auto-linking feature of CGAL (as those provided with the
Windows Installer)

STEP 1.

  Download the Brian Gladman's ports for gmp and mpfr for Visual Studio from:

    http://fp.gladman.plus.com/computing/gmp4win.htm

  You would choose a version of gmp and mpfr according to the ports available.
  
STEP 2.

  Download the corresponding versions of gmp and mpfr depending on the Gladman's port downloaded
  
STEP 3.

  Follow the instructions supplied by Gladman's port and build one of the static versions
  of the gmp and mpfr libs.
    
  Typically, in this step you would
  
     untar the gmp, mpfr and ports archives under a common root folder
     launch gmp.sln
     build the preliminary "gen_*" targets 
     build one of the static library targets, such as "gmp_lib_gc"
     launch mpfr.sln
     renamed "mparam_h.in" to "mparam.h" 
     build the target lib_mpfr.
  
  The details may vary depending on the actual port you are using. Just make sure to follow
  the port instructions and build ONE static version of gmp and mpfr. You can choose
  the generic version (_gc) or one of the optimized (_pN) version but make sure to build only one.
  
  You would end up with a structure similar to the following:
    
        <base directory>
          gmp-4.2.4
            build.vc9
              lib_gmp_gc
                Win32
                  Debug
                    gmp.lib
                    gmp.pdb
                  Release
                    gmp.lib
                    gmp.pdb
          mpfr-2.3.2
            build.vc9 
              lib_mpfr
                Win32
                  Debug
                    mpfr.lib
                    mpfr.pdb
                  Release
                    mpfr.lib
                    mpfr.pdb


STEP 4

  copy the batch script "create_mangled_versions.bat" into the base directory
  containing gmp and mpfr.
  
  Run it passing the directory containing gmp and mpfr. For example:
  
    C:\>create_mangled_versions gmp-4.2.4
    C:\>create_mangled_versions mpfr-2.3.2
    
  This script scans the given directory in search for gmp and mpfr binaries.
  Each time one is found, it creates a copy of it using a "mangled" named that
  conforms to the autolinking requirements. It constructs the mangled name based
  on the relative downpath to the library found. For example, if it finds
    
    build.vc9\lib_gmp_gc\Win32\Debug\gmp.lib
    
  it creates a copy, in the same directory, named:
  
    gmp_vc90_mt_gd.lib
    
  All mangled binaries are also copied to a folder named "mangled_binaries" that
  the script creates in the directory where it is run.
  

      
