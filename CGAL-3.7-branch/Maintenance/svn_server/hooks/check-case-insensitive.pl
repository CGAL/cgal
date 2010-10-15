#!/usr/bin/perl -w

# ====================================================================
# Copyright (c) 2000-2004 Collab Net.  All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at http://subversion.tigris.org/license-1.html.
# If newer versions of this license are posted there, you may use a
# newer version instead, at your option.
#
# This software consists of voluntary contributions made by many
# individuals.  For exact contribution history, see the revision
# history and logs, available at http://subversion.tigris.org/.
# ====================================================================

use strict;
require 5.004; # This is when locale support was added.

# This 'use encoding' and setting the LANG environment variable has the
# desired effect of handling the comparison of extended characters and
# preventing a commit.
#use encoding "utf8"; # LS 06/2006: commented out because it displays incorrectly
                      #             filenames with diacritic characters on InriaGForge
$ENV{'LANG'} = 'en_US.UTF-8'; # LS 06/2006: 'en_GB.UTF-8' breaks svnlook on InriaGForge

# Please check the path to svn and svnlook is correct...
# LS 08/2006: request also path to svn.
my $svnlook;
my $svn;
if ($^O eq 'MSWin32') {
  $svnlook = '"c:\Program Files\subversion\bin\svnlook.exe"';
  $svn     = '"c:\Program Files\subversion\bin\svn.exe"';
} else {
  $svnlook = '/usr/bin/svnlook';
  $svn     = '/usr/bin/svn';
}

# This script can be called from a pre-commit hook on either Windows or a Unix
# like operating system.  It implements the checks required to ensure that the
# repository acts in a way which is compatible with a case preserving but
# case insensitive file system.
#
# When a file is added this script checks the file tree in the repository for
# files which would be the same name on a case insensitive file system and
# rejects the commit.
#
# On a Unix system put this script in the hooks directory and add this to the
# pre-commit script:
#
#  $REPOS/hooks/check-case-insensitive.pl "$REPOS" "$TXN" || exit 1
#
# On a windows machine add this to pre-commit.bat:
#
#  perl <path-to-script>\check-case-insensitive.pl %1 %2
#  if errorlevel 1 goto :ERROR
#  exit 0
#  :ERROR
#  echo Error found in commit 1>&2
#  exit 1
#
# You may need to change the setting of $svn and $svnlook to the path to the
# executable on your system.
#
# Turn on debug by adding up to three -debug options as the first options in
# the list.  The more -debug options the more output.  If you specify more
# than one the output goes into a file.
#
# If you have any problems with this script feel free to contact
# Martin Tomes <martin@tomes.org.uk>

# Bugfixes and some debug code added by Jeremy Bettis <jeremy@deadbeef.com>

my $openstr = '-|';
# Shift off any debug options.
my $debug = 0;
while (@ARGV and $ARGV[0] =~ /^-d(ebug)?$/) {
  $debug++;
  shift;
}

# If there is too much debug output to STDERR subversion doesn't like it, so,
# if a lot of output is expected send it to a file instead.
if ($debug > 1) { # LS 08/2006: compare to 1 to match comments
  if ($^O eq 'MSWin32') {
    open(STDERR, ">c:/svnlog.txt")
      or die "$0: cannot open 'c:/svnlog.txt' for writing: $!\n";
  } else {
    open(STDERR, ">/tmp/svnlog.txt")
      or die "$0: cannot open '/tmp/svnlog.txt' for writing: $!\n";
  }
}

# Usage
unless (@ARGV > 1) {
  die "usage: $0 [-d [-d [-d]]] repos txn [--revision]\n";
}

# Fetch the command line arguments.
my $repos = shift;
$repos =~ s/\/$//; # LS 08/2006: remove trailing slash
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
# check-case-insensitive.pl -debug path 987 -r
# and then the check-case-insensitive.pl passes -r to svn/svnlook instead of
# --transaction.
my $flag = '--transaction';
$flag = shift if @ARGV;

# Each added path put here.
my @added;

# LS 08/2006: Each removed path put here.
my @removed;

# LS 08/2006: Each modified directory put here.
my @modified_dirs;

# The file tree as a hash, index lower cased name, value actual name.
my %tree;

# Command being executed.
my $cmd;

print STDERR "LANG=", $ENV{'LANG'}, "\n" if ($debug and defined($ENV{'LANG'}));

# Get lists of added and removed files.
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
  # LS 08/2006: Get also the list of removed files.
  elsif (/^D\s+(\S.*)/) {
    push @removed, $1;
  }
}
close SVNLOOK;

# Trace
if ($debug) {
  print STDERR "Added " . ($#added + 1) . " item(s):\n";
  foreach my $itm (@added) {
    print STDERR " $itm\n";
  }

  print STDERR "Removed " . ($#removed + 1) . " item(s):\n";
  foreach my $itm (@removed) {
    print STDERR " $itm\n";
  }
}

unless (@added) {
  print STDERR "No files added\n" if ($debug);
  # No added files so no problem.
  exit(0);
}

# LS 08/2006: Get the list of *all* modified directories.
$cmd = "$svnlook dirs-changed \"$repos\" $flag $txn";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, $openstr, $cmd)
  or die("$0: cannot open '$cmd' pipe for reading: $!\n");
while (<SVNLOOK>) {
  chomp;
  print STDERR "  ", $_, "\n" if ($debug > 2);
  # There isn't a leading slash on $changed but there is a trailing one.  When
  # it is the root of the repository the / is a pain, so always remove the
  # trailing slash and put it back in where needed.
  $_ =~ s/\/$//;
  push @modified_dirs, $_;
}
close SVNLOOK;

# LS 08/2006: We will check transations against HEAD and revisions against the previous commit.
my $lastrev = "";
unless ($flag eq '--transaction') {
    $lastrev = "--revision " . ($txn-1);
}

# Get the modified directories' content at the previous revision and turn the output into
# complete paths for each file.
# LS 08/2006: large rewrite with 'svn ls' instead of 'svnlook tree' (which is too slow)
local *SVN;
foreach my $changed (@modified_dirs) {
  $cmd = "$svn ls \"file://$repos/$changed\" $lastrev";
  print STDERR "$cmd\n" if ($debug);
  open(SVN, $openstr, $cmd)
    or die("$0: cannot open '$cmd' pipe for reading: $!\n");
  while (<SVN>) {
    chomp;
    if (/^(.*)\/$/) { # Is a directory.
      my $name = $1 . '/';
      my $index;
      if ($changed eq '') {
        $index = $name;
      } else {
        $index = $changed . '/' . $name;
      }
      $index =~ s/\/$//; # LS 06/2006: Remove trailing slash to compare file and folders
      $index = lc($index);
      $tree{$index} = $name; # Index the hash with case folded name.
      print STDERR "\$tree{$index}=$name (dir)\n" if ($debug > 1);
    } else {  # This is a real file name, not a directory.
      my $name = $_;
      my $index;
      if ($changed eq '') {
        $index = $name;
      } else {
        $index = $changed . '/' . $name;
      }
      $index = lc($index);
      $tree{$index} = $name; # Index the hash with case folded name.
      print STDERR "\$tree{$index}=$name\n" if ($debug > 1);
    }
  }
  close SVN;
}

# LS 08/2006: remove deleted paths from %tree
foreach my $newfile (@removed) {
  $newfile =~ s/\/$//;
  my $lcnewfile = lc($newfile);
  print STDERR "Delete \$tree{$lcnewfile}\n" if ($debug > 1);
  delete($tree{$lcnewfile});
}

# Check each added file
my $failmsg;
my %newtree;
foreach my $newfile (@added) {
  $newfile =~ s/\/$//; # LS 06/2006: Remove trailing slash to compare file and folders
  my $lcnewfile = lc($newfile);
  print STDERR "Checking \$tree{$lcnewfile}\n" if ($debug > 1);
  if (exists($tree{$lcnewfile})) {
    $failmsg .= "\n  $newfile already exists as " . $tree{$lcnewfile};
  }
  elsif (exists($newtree{$lcnewfile})) {
    $failmsg .= "\n  $newfile also added as " . $newtree{$lcnewfile};
  }
  $newtree{$lcnewfile} = $newfile;
}
if (defined($failmsg)) {
  print STDERR "\nFile name case conflict(s) found:\n" . $failmsg . "\n";
  exit 1;
}

exit 0;
