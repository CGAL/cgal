#!/usr/bin/perl -w

# ====================================================================
# check-symlinks.pl
# Laurent Saboret, INRIA, 2006
#
# This script can be called from a pre-commit hook on either Windows or a Unix
# like operating system.  It implements the checks required to reject Unix
# symbolic links, which are not valid (in fact: ignored) on Windows.
#
# When a file is added this script checks if it is a symbolic link. If it
# is, the commit is rejected.
#
# On a Unix system put this script in the hooks directory and add this to the
# pre-commit script:
#  $REPOS/hooks/check-symlinks.pl "$REPOS" "$TXN" || exit 1
#
# On a windows machine add this to pre-commit.bat:
#  perl <path-to-script>\check-symlinks.pl %1 %2
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
$ENV{'LANG'} = 'en_US.UTF-8';

# Please check the path to svnlook is correct...
my $svnlook;
if ($^O eq 'MSWin32') {
    $svnlook = '"c:\Program Files\subversion\bin\svnlook.exe"';
} else {
    $svnlook = '/usr/bin/svnlook';
}

my $openstr = '-|';
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
like operating system. It implements the checks required to reject Unix
symbolic links, which are not valid (in fact: ignored) on Windows.
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
# check-symlinks.pl -debug path 987 -r
# and then the check-symlinks.pl passes -r to svnlook instead of
# --transaction.
my $flag = '--transaction';
$flag = shift if @ARGV;

# Each added path put here.
my @added;

# Command being executed.
my $cmd;

print STDERR "LANG=", $ENV{'LANG'}, "\n" if ($debug and defined($ENV{'LANG'}));

# Get a list of added files.
local *SVNLOOK;
$cmd = "$svnlook changed \"$repos\" $flag $txn";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, $openstr, $cmd)
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
  $cmd = "$svnlook proplist \"$repos\" \"$newfile\" $flag $txn";
  print STDERR "$cmd\n" if ($debug);
  open(SVNLOOK, $openstr, $cmd)
    or die("$0: cannot open '$cmd' pipe for reading: $!\n");
  while (<SVNLOOK>) {
    chomp;
    if (/\bsvn:special\b/) {
      $failmsg .= "\n  $newfile is a symbolic link.\n";
      next;
    }
  }
  close SVNLOOK;
}
if (defined($failmsg)) {
    print STDERR "\nInvalid file(s) found:\n" . $failmsg . "\n";
    exit 1;
}

exit 0;
