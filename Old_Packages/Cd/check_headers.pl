#!/usr/local/bin/perl -w

use strict;
#use Cwd;
use File::Copy;
#use File::Basename;
#use File::Find;
use Getopt::Std;

my ($authors,$revision,$revision_date);

#
# This script Checks if a file header complies to the CGAL rules.

sub usage()
{
    print STDERR<<'EOF';
usage: $0 (-h|-u|-c) file1 ...

Exactly one of the options -h, -u and -c must be present.
    -h show this message and quit
    -u update headers
    -c check but do no updates

EOF
}

sub gjmove($$)
{
    return 1 if rename($_[0], $_[1] );
    return (system('mv', "$_[0]", "$_[1]") == 0);
}

#----------------------------------------------------------------#
#                     initialization                             #
#----------------------------------------------------------------#

my $TEMPFILE;

sub print_license()
{
    print TEMPFILE <<'END_OF_LICENSE';
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
END_OF_LICENSE
}

#----------------------------------------------------------------#
#                     add_header                                 #
#
#
#----------------------------------------------------------------#

sub add_header($$$$$$$)
{
    my ($version, $date, $package, $full_filename, $revision, $revision_date, $authors) = @_;
    $version = "" if (!defined($version));
    $package = "" if (!defined($package));
    my $year = `date '+%Y'`;
    chomp $year;
    print TEMPFILE "// ", '=' x 70,"\n";
    print TEMPFILE "//\n";
    print TEMPFILE "// Copyright \(c\) $year The CGAL Consortium\n";
    print_license;
    print TEMPFILE "// ", '-' x 70,"\n";
    print TEMPFILE <<"END_OF_HEADER";
//
// release       : $version
// release_date  : $date//
// file          : $full_filename
// package       : $package
$revision$revision_date$authors// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann\@sophia.inria.fr)
//
END_OF_HEADER
    print TEMPFILE "// ", '=' x 70,"\n";
    while( <SOURCE_FILE> ) {
	print TEMPFILE $_ if !m<^\s*//\s*revision_date\s*:>
	                  && !m<^\s*//\s*revision\s*:>
			  && !m<^\s*//\s*author(s|\(s\))?\s*:>;
    }
}

sub remove_header()
{
    my $no_header = 0;
    while( <SOURCE_FILE> ) {
        $no_header=1 if !m<^\s*//>;
	print TEMPFILE $_ if $no_header
	                  || m<^\s*//\s*revision_date\s*:>
	                  || m<^\s*//\s*revision\s*:>
			  || m<^\s*//\s*author(s|\(s\))?\s*:>;
    }
}

#----------------------------------------------------------------#
#                     check_and_update_2                         #
#----------------------------------------------------------------#

sub check_and_update_2($$$$)
{
    my ($version, $date, $package, $filename) = @_;
    my $status=0;
    # 0: no header 1: required header 2: header

    my ($authors_seen, $revision_seen, $revision_date_seen) = (0,0,0);
    while ( <SOURCE_FILE> ) {
	++$authors_seen if m<^\s*//\s*author(s|\(s\))?\s*:>;
	$authors = $_ if m<^\s*//\s*author(s|\(s\))?\s*:>;
	++$revision_seen if m<^\s*//\s*revision\s*:>;
	$revision = $_ if m<^\s*//\s*revision\s*:>;
	++$revision_date_seen if m<^\s*//\s*revision_date\s*:>;
	$revision_date = $_ if m<^\s*//\s*revision_date\s*:>;
	$status = 2 if m<\s*//\s==========>;
	$status = 2 if m<\s*//\s*Copyright>;
	$status = 2 if m<\s*//\s*release\s*:>;
	$status = 2 if m<\s*//\s*release\s*:>;
	$status = 2 if m<\s*//\s*release_date\s*:>;
	$status = 2 if m<\s*//\s*file\s*:>;
	$status = 2 if m<\s*//\s*package\s*:>;
	$status = 2 if m<\s*//\s*coordinator\s*:>;
	print TEMPFILE $_;
    }

    if ($revision_seen == 0) {
	print FILE_CHECKS
	    "$filename: revision field missing in header.\n";
	print TEMPFILE "// revision      : ?\n";
    } elsif ($revision_seen > 1) {
	print FILE_CHECKS
	    "$filename: Multiple revision fields in header.\n";
    }
    if ($revision_date_seen == 0) {
	print FILE_CHECKS
	    "$filename: revision_date field missing in header.\n";
	print TEMPFILE "// revision_date : ?\n";
    } elsif ($revision_date_seen > 1) {
	print FILE_CHECKS
	    "$filename: Multiple revision_date fields in header.\n";
    }
    if ($authors_seen == 0) {
	print FILE_CHECKS
	    "$filename: authors field missing in header.\n";
	print TEMPFILE "// author(s)     : ?\n";
    } elsif ($authors_seen > 1) {
	print FILE_CHECKS
	    "$filename: Multiple authors fields in header.\n";
    }
    if ($status==0 && $authors_seen==1 && $revision_seen==1 && $revision_date_seen==1) {
        $status=1;
    }

    my ($lines_exceeding_length, $has_line_directives) = (0, 0);
    while ( <SOURCE_FILE> ) {
	$lines_exceeding_length +=1 if length $_ > 80;
	$has_line_directives = 1 if m|^\s*#\s*line\s|;
	print TEMPFILE $_;
    }
    if ($lines_exceeding_length) {
	print FILE_CHECKS
	    "$filename has $lines_exceeding_length",
	    " lines over 80 characters.\n";
    }
    if ($has_line_directives) {
	print FILE_CHECKS "$filename has line directives.\n";
    }
    return $status;
}

sub check_and_update_file($$$$)
{
    my ($filename, $version, $date, $package) = @_;
    my $check_status;
    open SOURCE_FILE, "<$filename" || die "Error opening $filename: $!\n";
    open TEMPFILE, ">$TEMPFILE" || die;
    $check_status =check_and_update_2($version, $date, $package, $filename);
    close SOURCE_FILE || die "Error closing $filename: $!";
    close TEMPFILE || die "Error closing temporary file: $!\n";
    if ($check_status == 1) {
	print FILE_CHECKS "Completing header for $filename.\n";
	open SOURCE_FILE, "<$filename"
	    || die "Error opening $filename: $!\n";
	open TEMPFILE, ">$TEMPFILE" || die;
	add_header($version, $date, $package, $filename, $revision, $revision_date, $authors);
	close SOURCE_FILE || die "Error closing $filename: $!";
	close TEMPFILE || die "Error closing temporary file: $!\n";
    } elsif ($check_status == 2) {
	print FILE_CHECKS "Stripping header for $filename.\n";
	open SOURCE_FILE, "<$filename"
	    || die "Error opening $filename: $!\n";
	open TEMPFILE, ">$TEMPFILE" || die;
	remove_header();
	close SOURCE_FILE || die "Error closing $filename: $!";
	close TEMPFILE || die "Error closing temporary file: $!\n";
    }
    if ($::opt_u) {
	gjmove($TEMPFILE, $filename )
	|| warn "Could not update file $filename\n";
    }
}


#$::PARENT_DIR=cwd();

sub main()
{
    umask(002);
    
    my $version = `cat version`;
    my $date=`date +"%-e %b %Y"`;
    $version=~s/ (.*)\s*$//;
    getopts('hucv:d:p:');
    if ($::opt_h ) {
	usage;
	return;
    }
    $::opt_h = 0;
    if ($::opt_u and $::opt_c) {
	usage;
	die "Both -c and -u option present\n";
    }
    if ($::opt_u ) {
	$TEMPFILE="tmp.$$";
    } else {  # no updates, only checking
	die if !$::opt_c; # mainly put here for shutting up warnings
	$TEMPFILE="/dev/null";
    }
    
    open FILE_CHECKS, ">-";
    my $filename;
    foreach $filename (@ARGV) {
	check_and_update_file($filename, $version, $date, "Cd ($version)");
    }
    close  FILE_CHECKS;
}

main;

