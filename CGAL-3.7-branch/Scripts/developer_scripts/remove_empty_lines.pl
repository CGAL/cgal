#!/usr/bin/perl

# This script-let removes empty lines at the beginning and at the end of the file.
#
# As complement cleanup, one can also remove trailing white spaces at the end of lines:
# perl -pi.bak -e 's/\s+$/\n/' */examples/*/*.cpp

local($/) = undef;
my $text = <>;

$text =~ s/\A\n+//mg;
$text =~ s/\n+\Z/\n/mg;

print "$text";
