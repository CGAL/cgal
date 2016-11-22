#!/usr/bin/perl

use strict;

my $dir = @ARGV[0];
my $htmldir = @ARGV[1];

if ( ! -d $htmldir ) {
	`mkdir $htmldir`;
}
if ( ! -d "$htmldir/html" ) {
	`mkdir $htmldir/html`;
}
if ( ! -d "$htmldir/html/picts" ) {
	`mkdir $htmldir/html/picts`;
}

my $html = "$htmldir/html";

my @meshes = glob($dir.'/ProgramOutput.*');
my @corenames;
foreach (@meshes) {
	#get name
	(my $corename) = ( $_ =~ m/ProgramOutput\.(.*)\.txt/ );
	push(@corenames,$corename);

	# generate histogram svg
	my $svg_out = "picts/$corename.svg";
	`create_histogram.pl $dir/$corename-histo.txt > $html/$svg_out`;
}


my $out_filename = "$html/results.html";
open(my $out, ">", $out_filename) || die("Impossible d'ouvrir $out_filename: $!");

print $out "<html>
<style>.results{border:1px solid black; margin:10px;} .result{padding-top:1ex; background-color:#EEEEEE; margin-left:220px; padding-left:1em;} .spacer{clear:both;}
.mesh{color:#884416;} .error{color:red;} .addinfo{color:darkgreen; font-style:italic} .time{color:darkblue;}
h1,h2{text-align:center;} h1{border: 3px solid black; background-color:#AAAAAA; }
h3{font-size:100%; margin:0 1em 1.1ex -.5em; border-bottom:1px solid #888888; padding-left:.5em;}</style>
<body>\n";

foreach (@corenames) {
	if ( (my $mesh_file) = ( $_ =~ m/(^[^_]*)_0$/ ) ) {
		print $out "<h1>$mesh_file</h1>\n";
	}
	if ( (my $mesh_file) = ( $_ =~ m/(^[^_]*_[0-9]+)$/ ) ) {
		print $out "<h2>$mesh_file</h2>\n";
	}
	print $out "<div class=\"results\">\n";

	my $svg_file = svg_name($_);
	if ( -f "$html/$svg_file" ) {
		print $out "<object type=\"image/svg+xml\" data=\"$svg_file\" style=\"{float:left;}\"></object>\n";
	}
	else {
		print $out "<img src=\"dummy\">\n";
	}
	print $out "<div class=result>\n";
	print $out "<h3>$_: ";
	my $btitle = 1;

	my $out_file = out_name($_);
	open(my $fic, $out_file) || die("Impossible d'ouvrir $out_file: $!");
	
	while( <$fic> ) {
		chomp $_;
		$_ =~ s/^(--.*)/<span class=mesh>$1<\/span>/;
		$_ =~ s/^(ERROR.*|  ! .*)/<span class=error>$1<\/span>/;
		$_ =~ s/(\[.*\])$/<span class=addinfo>$1<\/span>/;
		$_ =~ s/\(([0-9]{3}[0-9]*|[0-9\.]{4}).*s\)/<span class=time>\($1s\)<\/span>/;
		
		if ( $btitle==1 ) { $_ = $_."</h3>"; $btitle=0;}
		else { $_ = $_."<br>"; }

		print $out "$_";
	}
	print $out "\n<div class=\"spacer\"> </div></div></div>\n";
}

print $out "</body></html>";




sub svg_name() {
	my $arg = shift;
	return "picts/$arg.svg";
}

sub out_name() {
	my $arg = shift;
	return "$dir/ProgramOutput.$arg.txt";
}
