#! /sw/bin/perl5

$prog = substr($0,rindex($0,'/')+1) ;
$Usage = <<USAGE ;
Usage: $prog args
USAGE
sub Usage { die "$_[0]$Usage" ; }
sub Error { die "$prog: $_[0]\n" ; }

# usage: &NGetOpt(ARG,ARG,ARG,..)
# ARG = 'ID' | 'ID=SPC' | 'ID:SPC' for no-arg, required-arg or optional-arg
# ID  = perl identifier
# SPC = i|f|s for integer, fixedpoint real or string argument
# defines $opt_ID as 1 or user specified value

require 'newgetopt.pl' ;
$newgetopt'ignorecase = 0 ; # options are case-sensitive
&Usage() unless &NGetOpt() ;
&Usage("Arg count\n") unless @ARGV >= 0 ;

@ARGV = ( '/users/www/CGAL/GET/get-info' ) unless @ARGV ;

while ( <> )
  { chomp ;
    $uml=0;
    $email_address="";
    for $pair ( split(/&/, $_))
      { ($name, $value) = split(/=/, $pair) ;
        $value =~ tr/+/ / ;
        $name =~ s/%([0-9A-Za-z][0-9A-Za-z])/pack(c,hex($1))/eg ;
        $value =~ s/%([0-9A-Za-z][0-9A-Za-z])/pack(c,hex($1))/eg ;

	$uml = 1 if ($name eq 'uml' and $value eq 'yes');
	$email_address = $value if $name eq 'email';

      }
    if ($uml and $email_address) {
	print "$email_address\n" ;
    }
  }
