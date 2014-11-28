#!/usr/bin/env perl

## Run that script on a file (files one by one) that is under QPL, to turn
## its license into GPLv3+.

use warnings;

sub print_gplv3
{
    my $mark = shift;
    my $software = shift;
    print <<"EOF"
$mark This file is part of $software.
$mark You can redistribute it and/or modify it under the terms of the GNU
$mark General Public License as published by the Free Software Foundation,
$mark either version 3 of the License, or (at your option) any later version.
EOF
}

unshift(@ARGV, '-') unless @ARGV;
while ($ARGV = shift) {
    my $first_line=1;

    my $previous_line;
    my $current_line;

    my $backup = $ARGV . ".bak";
    rename($ARGV, $backup);
    open(my $fh, "<", $backup);
    open(my $out, ">", $ARGV);
    select($out);
    while ($current_line = <$fh>) {
        if($first_line) {
            $first_line=0;
            $previous_line = $current_line;        
        } else {
          if( ((my $mark, $software) = $previous_line =~ /(\/\/|\#| \*|%+) *This file is part of (.*);/) &&
                ($current_line =~ /the terms of the Q Public License /) )
            {
                print_gplv3($mark, $software);
                $previous_line = <$fh>;
                $previous_line = <$fh>;
            }
            else {
                print $previous_line;
                $previous_line = $current_line;        
            }
        }
    }
    print $previous_line;
    close($fh);
    select(STDOUT);
    close($out);
    unlink($backup);
}
