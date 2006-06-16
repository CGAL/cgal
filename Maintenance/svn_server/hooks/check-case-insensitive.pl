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
# preventing a commit.  However, if any of the files in conflict have
# extended characters in them this is the error displayed by the client:
#
#   Commit failed (details follow):
#   svn: MERGE request failed on '/svn/play/martinto/trunk'
#   svn: General svn error from server
#
# It should list the file names which are in conflict.  But it does stop the
# commit.
use encoding "utf8";
#$ENV{'LANG'} = 'en_GB.UTF-8';  # breaks svnlook on InriaGForge
$ENV{'LANG'} = 'fr_FR.UTF-8';   # ok on InriaGForge

# Please check the path to svnlook is correct...
my $svnlook;
if ($^O eq 'MSWin32') {
  $svnlook = '"c:\Program Files\subversion\bin\svnlook.exe"';
} else {
  $svnlook = '/usr/bin/svnlook';
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
# You may need to change the setting of $svnlook to the path to the
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
if ($debug > 0) {
  if ($^O eq 'MSWin32') {
    open(STDERR, ">c:/svnlog.txt")
      or die "$0: cannot open 'c:/svnlog.txt' for writing: $!\n";
  } else {
    open(STDERR, ">/tmp/svnlog.txt")
      or die "$0: cannot open '/tmp/svnlog.txt' for writing: $!\n";
  }
}

# Fetch the command line arguments.
unless (@ARGV > 1) {
  die "usage: $0 [-d [-d [-d]]] repos txn [--revision]\n";
}

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
# check-case-insensitive.pl -debug path 987 -r
# and then the check-case-insensitive.pl passes -r to svnlook instead of
# --transaction.
#
# Of course when it gets down to # Get the file tree at the previous revision,
# then it doesn't work, but most of my problems were found before that point.
my $flag = '--transaction';
$flag = shift if @ARGV;

# Each added path put here.
my @added;

# The file tree as a hash, index lower cased name, value actual name.
my %tree;

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

if ($debug) {
  print STDERR "Added " . ($#added + 1) . " items:\n";
  foreach my $itm (@added) {
    print STDERR " $itm\n";
  }
}

unless (@added) {
  print STDERR "No files added\n" if ($debug);
  # No added files so no problem.
  exit(0);
}

# Get the shortest directory name which has changed, this will be the path
# into the repository to use to get the history.
$cmd = "$svnlook dirs-changed \"$repos\" $flag $txn";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, $openstr, $cmd)
  or die("$0: cannot open '$cmd' pipe for reading: $!\n");
my $shortest=999999;
my $changed;
while (<SVNLOOK>) {
  chomp;
  print STDERR "  ", $_, "\n" if ($debug > 2);
  if (length($_) < $shortest) {
    $changed = $_;
    $shortest = length($_);
  }
}
close SVNLOOK;
# There isn't a leading slash on $changed but there is a trailing one.  When
# it is the root of the repository the / is a pain, so always remove the
# trailing slash and put it back in where needed.
$changed =~ s/\/$//;

# Use the history of $changed path to find the revision of the previous commit.
$cmd = "$svnlook history \"$repos\" \"$changed/\"";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, $openstr, $cmd)
  or die("$0: cannot open '$cmd' pipe for reading: $!\n");
my $lastrev;
while (<SVNLOOK>) {
  chomp;
  if (/(\d+)/) {
    $lastrev = $1;
    last;
  }
}
close SVNLOOK;

# Get the file tree at the previous revision and turn the output into
# complete paths for each file.
my @path;
$cmd = "$svnlook tree \"$repos\" \"$changed/\" --revision $lastrev";
print STDERR "$cmd\n" if ($debug);
open(SVNLOOK, $openstr, $cmd)
  or die("$0: cannot open '$cmd' pipe for reading: $!\n");
while (<SVNLOOK>) {
  chomp;
  print STDERR "tree: '", $_, "'\n" if ($debug > 2);
  next if (/^\/{1,2}$/); # Ignore the root node.  Two /'s at root of the repos.
  if (/^(\s+)(.*)\/$/) { # Is a directory.
    $#path = length($1)-2; # Number of spaces at start of line is nest level.
    push @path, $2;
    my $name = join('/', @path) . '/';
    my $index;
    if ($changed eq '') {
      $index = $name;
    } else {
      $index = $changed . '/' . $name;
    }
    $index =~ s/\/$//; # Remove trailing slash to compare file and folders
    $tree{lc($index)} = $name; # Index the hash with case folded name.
    print STDERR "\$tree{lc($index)}=$name (dir)\n" if ($debug > 1);
  } elsif (/^(\s+)(.*)$/) {  # This is a real file name, not a directory.
    $#path = length($1)-2; # Number of spaces at start of line is nest level.
    my $name;
    if ($#path eq -1) {
      $name = $2;
    } else {
      $name = join('/', @path) . '/' . $2;
    }
    my $index;
    if ($changed eq '') {
      $index = $name;
    } else {
      $index = $changed . '/' . $name;
    }
    $tree{lc($index)} = $name; # Index the hash with case folded name.
    print STDERR "\$tree{lc($index)}=$name\n" if ($debug > 1);
  }
}
close SVNLOOK;

my $failmsg;

my %newtree;
foreach my $newfile (@added) {
  $newfile =~ s/\/$//; # Remove trailing slash to compare file and folders
  print STDERR "Checking \$tree{lc($newfile)}\n" if ($debug > 1);
  # Without the following line it gets the lc() wrong.
  my $junk = "x$newfile";
  my $lcnewfile = lc($newfile);
  if (exists($tree{$lcnewfile})) {
    $failmsg .= "\n  $newfile already exists as " . $tree{lc($newfile)};
  }
  elsif (exists($newtree{$lcnewfile})) {
    $failmsg .= "\n  $newfile also added as " . $newtree{lc($newfile)};
  }
  $newtree{$lcnewfile} = $newfile;
}
if (defined($failmsg)) {
  print STDERR "\nFile name case conflict found:\n" . $failmsg . "\n";
  exit 1;
}

exit 0;
