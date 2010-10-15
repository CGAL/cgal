#!/usr/bin/perl -w

use strict;
use lib "/home/groups/cgal/hooks/";
use Mail::Sender;

# Shift off any debug options.
my $debug = 0;
while (@ARGV and $ARGV[0] =~ /^-d(ebug)?$/) {
  $debug++;
  shift;
}
# Usage
unless (@ARGV > 2) {
    print <<"EOF";
usage: $0 [-d] repos user capabilities
This script can be called from a pre-commit hook on either Windows or a Unix
like operating system. It implements the checks required to reject Unix
symbolic links, which are not valid (in fact: ignored) on Windows.
EOF
    exit(0);
}

my $repo = shift @ARGV;
my $user = shift @ARGV;
my $capabilities_line = shift @ARGV;

my @capabilities = split(':', $capabilities_line);

if( not grep "mergeinfo",@capabilities ) {
}

my $subject = "CGAL svn: check-capabilities.pl";

my $header;
$header->{'to'} = "Laurent.Rineau\@geometryfactory.com";
$header->{'from'} = "$user\@users.gforge.inria.fr";
$header->{'subject'} = $subject;
$header->{'smtp'} = 'localhost';
$header->{'multipart'} = 'alternative';
#$header->{'debug'} = *main::STDERR;

my $mail = new Mail::Sender($header) or die "$0: cannot open mail connection: $!\n";
$mail->Open($header);
$mail->SendLineEnc("REPO: $repo\n",
                   "USER: $user\n",
                   "CAPS: $capabilities_line\n");
$mail->Close();

exit 0
