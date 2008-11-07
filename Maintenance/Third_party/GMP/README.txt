This directory contains the visual studio project files used to generate the precompiled binaries for gmp and mpfr.

To install it:

1.Download gmp-4.2.2 and mpfr-2.3.1 from 

2.Download the Gladman's ports to Visual Studio from

http://fp.gladman.plus.com/computing/gmp4win.htm

and follow the instructions to properly unzip the archive, but do not build any libraries.

3.Unzip "mangled_gmp_VC.zip" at the root directory of your gmp-4.2.2 installation.
4.Unzip "mangled_mpfr_VC.zip" at the rot directory of your mpfr-2.3.1 installation.

5.To build the gmp libraries open <gmp-4.2.2 root>/build.vc<8|9>/lib_gmp_gc/mangled_lib_gmp_gc.vcproj in Visual Studio and build it.

6.TO build the mpfr libs, *after* having built gmp, use <mpfr-2.3.1 root>/build.vc<8|9>/lib_mpfr/mangled_lib_mpfr.vcproj

