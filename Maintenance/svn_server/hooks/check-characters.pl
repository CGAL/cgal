#!/usr/bin/perl -w

# ====================================================================
# check-characters.pl
# Laurent Saboret, INRIA, 2006
#
# This script can be called from a pre-commit hook on either Windows or a Unix
# like operating system.  It implements the checks required to ensure that the
# repository accepts only file/folder names compatible with Windows/Unix/Mac
# file systems (plus SVN own requirements).
#
# When a file is added this script checks the file name. If it
# fails one of the requirements below, the commit is rejected.
# SVN requirements:
# - Control characters 0x0 to 0x1f and 0x7f are not supported
#   (see http://subversion.tigris.org/issues/show_bug.cgi?id=1954).
# Unix file system (including Linux and MacOSX):
# - Characters \0 and / are not supported
#   (see http://svn.haxx.se/dev/archive-2002-04/0144.shtml).
# Windows file systems FAT and NTFS:
# - Control characters 0x0 to 0x1f are not supported.
# - Characters / \ : * ? " < > | are not supported.
# - Filenames CON, PRN, AUX, NUL, COM1 to COM9, LPT1 to LPT9 and CLOCK$ are reserved,
#   even with a path or an extension.
# - A filename cannot end by . or SPACE characters.
# (see http://msdn.microsoft.com/library/default.asp?url=/library/en-us/fileio/fs/naming_a_file.asp)
#
# On a Unix system put this script in the hooks directory and add this to the
# pre-commit script:
#  $REPOS/hooks/check-characters.pl "$REPOS" "$TXN" || exit 1
#
# On a windows machine add this to pre-commit.bat:
#  perl <path-to-script>\check-characters.pl %1 %2
#  if errorlevel 1 goto :ERROR
#  exit 0
#  :ERROR
#  echo Error found in commit 1>&2
#  exit 1
#
# You may need to change the setting of $svnlook to the path to the
# executable on your system.
#
# Turn on debug by adding a -debug option as the first option in the list.
#
# Note: this script was created from Collab Net check-case-insensitive.pl
# ====================================================================

use strict;
require 5.004; # This is when locale support was added.
use Encode;
$ENV{'LANG'} = 'en_US.UTF-8';
$ENV{'LC_CTYPE'} = 'en_US.UTF-8';
binmode STDERR, ":utf8";

# Please check the path to svnlook is correct...
my $svnlook;
if ($^O eq 'MSWin32') {
    $svnlook = '"c:\Program Files\subversion\bin\svnlook.exe"';
} else {
    $svnlook = '/usr/bin/svnlook';
}

# Shift off any debug options.
my $debug = 0;
while (@ARGV and $ARGV[0] =~ /^-d(ebug)?$/) {
  $debug++;
  shift;
}

# Usage
unless (@ARGV > 1) {
print <<"END_USAGE";
usage: $0 [-d] repos txn [--revision]
This script can be called from a pre-commit hook on either Windows or a Unix
like operating system.  It implements the checks required to ensure that the
repository accepts only file/folder names compatible with Windows/Unix/Mac
file systems (plus SVN own requirements).
END_USAGE
exit(0);
}

# Fetch the command line arguments.
my $repos = shift;
my $txn = shift;

# Jeremy Bettis <jeremy@deadbeef.com> wrote the $flag code and has this to
# say about it:
#
# The reason I did that was so that I could test the hook without actually
# doing a commit.  Whenever I had a commit that succeeded in making a bad file
# or directory, or when a commit took too long I just did a sequence of
# operations like this:
#
# svnlook youngest path
# (it tells me that HEAD is 987 or whatever)
# check-characters.pl -debug path 987 -r
# and then the check-characters.pl passes -r to svnlook instead of
# --transaction.
my $flag = '--transaction';
$flag = shift if @ARGV;

# Each added path put here.
my @added;

# Command being executed.
my $cmd;

print STDERR "LANG=", $ENV{'LANG'}, "\n" if ($debug and defined($ENV{'LANG'}));
print STDERR "LC_CTYPE=", $ENV{'LC_CTYPE'}, "\n" if ($debug and defined($ENV{'LC_CTYPE'}));

# Get a list of added files.
local *SVNLOOK;
$cmd = "$svnlook changed \"$repos\" $flag $txn";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, '-|:utf8', $cmd)
    or die("$0: cannot open '$cmd' pipe for reading: $!\n");
while (<SVNLOOK>) {
  chomp;
  if (/^A\s+(\S.*)/) {
    push @added, $1;
  }
}
close SVNLOOK;

# Trace
if ($debug) {
  print STDERR "Added " . ($#added + 1) . " item(s):\n";
  foreach my $itm (@added) {
    print STDERR " $itm\n";
  }
}

# Stop if no added files
unless (@added) {
  print STDERR "No files added\n" if ($debug);
  # No added files so no problem.
  exit(0);
}

# Check each added file
my $failmsg;
foreach my $newfile (@added) {
    # Get basename
    $newfile =~ s/\/$//;                      # Remove trailing slash from folders
    $newfile =~ s/.*\///;                     # Remove path
    if ($debug) {
        print STDERR "Checking $newfile\n";

        print STDERR "    as unicode: ";
        my @unicode_chars = split //, $newfile;
        map { print STDERR "$_ (" . ord($_) . ") " } @unicode_chars;
        print STDERR "\n";

        print STDERR "    as UTF-8:   ";
        my $utf8_raw_string = Encode::encode("UTF-8", $newfile);
        my @utf8_octets = split //, $utf8_raw_string;
        map { print STDERR "$_ (" . ord($_) . ") " } @utf8_octets;
        print STDERR "\n";
    }

    # SVN requirements:
    # - Control characters 0x0 to 0x1f and 0x7f are not supported
    #   (see http://subversion.tigris.org/issues/show_bug.cgi?id=1954).
    if ($newfile =~ m/([\x00-\x1F\x7f])/) {
        $failmsg .= "\n  $newfile contains an invalid control character (0x00-0x1f 0x7f).\n";
        next;
    }

    # Unix file system (including Linux and MacOSX):
    # - Characters \0 and / are not supported
    #   (see http://svn.haxx.se/dev/archive-2002-04/0144.shtml).
    if ($newfile =~ m/[\x00\/]/) {
        $failmsg .= "\n  $newfile contains a character invalid on Unix (\\0 \/).\n";
        next;
    }

    # Windows file systems FAT and NTFS:
    # - Control characters 0x0 to 0x1f are not supported.
    # - Characters / \ : * ? " < > | are not supported.
    # - Filenames CON, PRN, AUX, NUL, COM1 to COM9, LPT1 to LPT9 and CLOCK$ are reserved,
    #   even with a path or an extension.
    # - A filename cannot end by . or SPACE characters.
    # (see http://msdn.microsoft.com/library/default.asp?url=/library/en-us/fileio/fs/naming_a_file.asp)
    if ($newfile =~ m/[\x00-\x1F\/\\\:\*\?\"\<\>\|]/) {
        $failmsg .= "\n  $newfile contains a character invalid on Windows (0x00-0x1f \/ \\ \: \* \? \" \< \> \|).\n";
        next;
    }
    if ($newfile =~ m/^(CON|PRN|AUX|NUL|COM\d|LPT\d|CLOCK\$)(\..*)?$/i) {
        $failmsg .= "\n  $newfile is a reserved filename on Windows.\n";
        next;
    }
    if ($newfile =~ m/[\. ]$/) {
        $failmsg .= "\n  $newfile is illegal. Windows filenames cannot end by . or SPACE characters.\n";
        next;
    }
}
if (defined($failmsg)) {
    print STDERR "\nInvalid file name(s) found:\n" . $failmsg . "\n";
    exit 1;
}

exit 0;
