#! /usr/bin/env perl

use Git;

#print STDERR "$ARGV[0]\n";

my $repo = Git->repository();

$last_date=$repo->command_oneline("log", "--pretty=format:\"%ad\"",  "-1",  "--",  $ARGV[0]);

while(<STDIN>) {
    s/\$Date[^\$]*\$/\$Date: $last_date\$/;
    print;
  }
