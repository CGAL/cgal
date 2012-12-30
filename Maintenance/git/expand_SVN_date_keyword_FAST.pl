#! /usr/bin/env perl


# This filter is not correct.
# It sets the same modification date for all files.

use Git;

#print STDERR "$ARGV[0]\n";

my $repo = Git->repository();

$last_date=$repo->command_oneline("log", "--pretty=format:\"%ad\"",  "-1");

while(<STDIN>) {
    s/\$Date[^\$]*\$/\$Date: $last_date\$/;
    print;
  }
