Instructions for creating the MacOS X Installer

The instructions below assume that you will be creating CGAL-X.Y
Make the appropriate modifications for building other releases.
The instructions below work for versions of CGAL from 3.2 and up.

1. Put the CGAL root directory under Package_root. Make sure you
   have stripped all directories and files that do not go to the
   release. The directory inside Package_root should be named CGAL-X.Y
2. Modify the scripts and files in the Resources directory to reflect
   the welcoming message (Welcome.rtf), the installations instructions
   to the user (ReadMe.rtf) and the license (License.rtf) to reflect the
   status of the release you are distributing.
3. This step requires the PackageMaker application that comes with Xcode
   developer tools. Open CGAL.pmproj with PackageMaker in order to make
   the appropriate modifications for version X.Y. In particular
   a. in the Configuration section make sure that
       i. the "Authentication" is set to "None"
      ii. the "Relocatable" and "Follow Symbolic Links" flags are checked
          and that all other flags are unchecked.
   b. Modify the package version and info string to reflect the current
      version and copyright year.
   The rest should be okay are they are.
4. Save your new project and type Apple+B to build the package.
   Save it as CGAL.pkg. Eventually commit your new project to the
   svn server.
5. Open the DiskUtility and create a read-write disk image (named
   CGAL-X.Y-RW.dmg) with the appropriate size. Copy CGAL.pkg in the disk
   image. Eject the disk image. Using the DiskUtility application convert
   the CGAL-X.Y-RW.dmg disk image to a read-only disk image. Save the new
   disk image as CGAL-X.Y.dmg. You all set!!!