# Mail::Sender.pm version 0.8.13
#
# Copyright (c) 2001 Jan Krynicky <Jenda@Krynicky.cz>. All rights reserved.
# This program is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.

package Mail::Sender; local $^W;
require 'Exporter.pm';
use vars qw(@ISA @EXPORT @EXPORT_OK);
@ISA = (Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(@error_str GuessCType);

$Mail::Sender::VERSION = '0.8.16';
$Mail::Sender::ver=$Mail::Sender::VERSION;

BEGIN {
	if ($] >= 5.008) {
		# fsck, the 5.8 is broken. The binmode() doesn't work for socket()s.
		require 'open.pm';
		'open'->import(OUT=>":raw :perlio");
	}
}

use strict;
use warnings;
no warnings 'uninitialized';
use Carp;
use FileHandle;
use IO::Socket::INET;
use File::Basename;

use MIME::Base64;
use MIME::QuotedPrint;
                    # if you do not use MailFile or SendFile and only send 7BIT or 8BIT "encoded"
					# messages you may comment out these lines.
                    #MIME::Base64 and MIME::QuotedPrint may be found at CPAN.

# include config file and libraries when packaging the script
if (0) {
	require 'Mail/Sender.config'; 	# local configuration
	require 'Symbol.pm'; 			# for debuging and GetHandle() method
	require 'Tie/Handle.pm';	# for debuging and GetHandle() method
	require 'IO/Handle.pm';	# for debuging and GetHandle() method
	require 'Digest/HMAC_MD5.pm'; # for CRAM-MD5 authentication only
	require 'Authen/NTLM.pm'; # for NTLM authentication only
} # this block above is there to let PAR, PerlApp, PerlCtrl, PerlSvc and Perl2Exe know I may need those files.

BEGIN {
    my $config = $INC{'Mail/Sender.pm'};
    die "Wrong case in use statement or Mail::Sender module renamed. Perl is case sensitive!!!\n" unless $config;
	my $compiled = !(-e $config); # if the module was not read from disk => the script has been "compiled"
    $config =~ s/\.pm$/.config/;
	if ($compiled or -e $config) {
		# in a Perl2Exe or PerlApp created executable or PerlCtrl generated COM object
		# or the config is known to exist
		eval {require $config};
		if ($@ and $@ !~ /Can't locate /) {
			print STDERR "Error in Mail::Sender.config : $@" ;
		}
	}
}

#local IP address and name
my $local_name =  $ENV{HOSTNAME} || $ENV{HTTP_HOST} || (gethostbyname 'localhost')[0];
$local_name =~ s/:.*$//; # the HTTP_HOST may be set to something like "foo.bar.com:1000"
my $local_IP =  join('.',unpack('CCCC',(gethostbyname $local_name)[4]));

#time diference to GMT - Windows will not set $ENV{'TZ'}, if you know a better way ...
my $GMTdiff;

use Time::Local;
sub ResetGMTdiff {
	my $local = time;
	my $gm = timelocal( gmtime $local );
	my $sign = qw( + + - )[ $local <=> $gm ];
	$GMTdiff = sprintf "%s%02d%02d", $sign, (gmtime abs( $local - $gm ))[2,1];
}
ResetGMTdiff();

#
my @priority = ('','1 (Highest)','2 (High)', '3 (Normal)','4 (Low)','5 (Lowest)');

#data encoding
my $chunksize=1024*4;
my $chunksize64=71*57; # must be divisible by 57 !

sub enc_base64 {my $s = encode_base64($_[0]); $s =~ s/\x0A/\x0D\x0A/sg; return $s;}
my $enc_base64_chunk = 57;

sub enc_qp {my $s = $_[0];$s =~ s/\x0D\x0A/\n/g;$s = encode_qp($s); $s=~s/^\./../gm; $s =~ s/\x0A/\x0D\x0A/sg; return $s}

sub enc_plain {my $s = shift; $s=~s/^\./../gm; $s =~ s/(?:\x0D\x0A?|\x0A)/\x0D\x0A/sg; return $s}

{ my $username;
sub getusername () {
	return $username if defined($username);
	return $username=eval{getlogin || getpwuid($<)} || $ENV{USERNAME};
}
}

#IO
use vars qw($debug);
$debug = 0;

#reads the whole SMTP response
# converts
#	nnn-very
#	nnn-long
#	nnn message
# to
#	nnn very
#	long
#	message
sub get_response ($) {
	my $s = shift;
	my $res = <$s>;
	if ($res =~ s/^(\d\d\d)-/$1 /) {
		my $nextline = <$s>;
		while ($nextline =~ s/^\d\d\d-//) {
			$res .= $nextline;
			$nextline = <$s>;
		}
		$nextline =~ s/^\d\d\d //;
		$res .= $nextline;
	}
	$Mail::Sender::LastResponse = $res;
	return $res;
}

sub send_cmd ($$) {
	my ($s, $cmd) = @_;
	chomp $cmd;
	if ($s->opened()) {
		print $s "$cmd\x0D\x0A";
		get_response($s);
	} else {
		return '400 connection lost';
	}
}

sub print_hdr {
	my ($s, $hdr, $str, $charset) = @_;
	return if !defined $str or $str eq '';
	$str =~ s/[\x0D\x0A\s]+$//;

	if ($charset && $str =~ /[^[:ascii:]]/) {
		$str = encode_qp($str);
		$str =~ s/=\r?\n$//;
		$str = "=?$charset?Q?" . $str . "?=";
	}

	$str =~ s/(?:\x0D\x0A?|\x0A)/\x0D\x0A/sg; # \n or \r => \r\n
	$str =~ s/\x0D\x0A([^\t])/\x0D\x0A\t$1/sg;
	if (length($str)+length($hdr) > 997) { # header too long, max 1000 chars
		$str =~ s/(.{1,997}[;,])\s+/$1\x0D\x0A\t/g;
	}
	print $s "$hdr: $str\x0D\x0A";
}


sub say_helo {
	my ($self, $s) = @_;
	my $res = send_cmd $s, "EHLO $self->{'client'}";
	if ($res !~  /^[123]/) {
		$res = send_cmd $s, "HELO $self->{'client'}";
		if ($res !~ /^[123]/) { return $self->Error(COMMERROR($_));}
		return;
	}

	$res =~ s/^.*\n//;
	$self->{'supports'} = {map {split /\s+/, $_, 2} split /\n/, $res};

	if (exists $self->{'supports'}{AUTH}) {
		my @auth = split /\s+/, uc($self->{'supports'}{AUTH});
		$self->{'auth_protocols'} = {map {$_, 1} @auth};
			# create a hash with accepted authentication protocols
	}

	$self->{esmtp}{_MAIL_FROM} = '';
	$self->{esmtp}{_RCPT_TO} = '';
	if (exists $self->{'supports'}{DSN} and exists $self->{esmtp}) {
		for (qw(RET ENVID)) {
			$self->{esmtp}{_MAIL_FROM} .= " $_=$self->{esmtp}{$_}"
				if $self->{esmtp}{$_} ne '';
		}
		for (qw(NOTIFY ORCPT)) {
			$self->{esmtp}{_RCPT_TO} .= " $_=$self->{esmtp}{$_}"
				if $self->{esmtp}{$_} ne '';
		}
	}
	return;
}

sub login {
	my $self = shift();
	my $auth = uc( $self->{'auth'}) || 'LOGIN';
	if (! $self->{'auth_protocols'}->{$auth}) {
		return $self->Error(INVALIDAUTH($auth));
	}

	$self->{'authid'} = $self->{'username'}
		if (exists $self->{'username'} and !exists $self->{'authid'});

	$self->{'authpwd'} = $self->{'password'}
		if (exists $self->{'password'} and !exists $self->{'authpwd'});

	$auth =~ tr/a-zA-Z0-9_/_/c; # change all characters except letters, numbers and underscores to underscores
	no strict qw'subs refs';
	&{"Mail::Sender::Auth::".$auth}($self);
}

# authentication code stolen from http://support.zeitform.de/techinfo/e-mail_prot.html
sub Mail::Sender::Auth::LOGIN {
	my $self = shift();
	my $s = $self->{'socket'};

	$_ = send_cmd $s, 'AUTH LOGIN';
	if (!/^[123]/) { return $self->Error(INVALIDAUTH('LOGIN', $_)); }

	if ($self->{auth_encoded}) {
		# I assume the username and password had been base64 encoded already!
		$_ = send_cmd $s, $self->{'authid'};
		if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }

		$_ = send_cmd $s, $self->{'authpwd'};
		if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }
	} else {
		$_ = send_cmd $s, &encode_base64($self->{'authid'});
		if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }

		$_ = send_cmd $s, &encode_base64($self->{'authpwd'});
		if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }
	}
	return;
}

use vars qw($MD5_loaded);
$MD5_loaded = 0;
sub Mail::Sender::Auth::CRAM_MD5 {
	my $self = shift();
	my $s = $self->{'socket'};

	$_ = send_cmd $s, "AUTH CRAM-MD5";
	if (!/^[123]/) { return $self->Error(INVALIDAUTH('CRAM-MD5', $_)); }
	my $stamp = $1 if /^\d{3}\s+(.*)$/;

	unless ($MD5_loaded) {
		eval 'use Digest::HMAC_MD5 qw(hmac_md5_hex)';
		die "$@\n" if $@;
		$MD5_loaded = 1;
	}

	my $user = $self->{'authid'};
	my $secret = $self->{'authpwd'};

	my $decoded_stamp = decode_base64($stamp);
	my $hmac = hmac_md5_hex($decoded_stamp, $secret);
	my $answer = encode_base64($user . ' ' . $hmac);
	$_ = send_cmd $s, $answer;
	if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }
	return;
}

sub Mail::Sender::Auth::PLAIN {
	my $self = shift();
	my $s = $self->{'socket'};

	$_ = send_cmd $s, "AUTH PLAIN";
	if (!/^[123]/) { return $self->Error(INVALIDAUTH('PLAIN', $_)); }

	$_ = send_cmd $s, encode_base64("\000" . $self->{'authid'} . "\000" . $self->{'authpwd'});
	if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }
	return;
}

{
my $NTLM_loaded=0;
sub Mail::Sender::Auth::NTLM {
	unless ($NTLM_loaded) {
		eval "use Authen::NTLM qw();";
		die "$@\n" if $@;
		$NTLM_loaded = 1;
	}
	my $self = shift();
	my $s = $self->{'socket'};

	$_ = send_cmd $s, "AUTH NTLM";
	if (!/^[123]/) { return $self->Error(INVALIDAUTH('NTLM', $_)); }

	Authen::NTLM::ntlm_user($self->{'authid'});
	Authen::NTLM::ntlm_password($self->{'authpwd'});
	Authen::NTLM::ntlm_domain($self->{'authdomain'})
		if defined $self->{'authdomain'};

	$_ = send_cmd $s, Authen::NTLM::ntlm();
	if (!/^3\d\d (.*)$/s) { return $self->Error(LOGINERROR($_)); }
	my $response = $1;
	$_ = send_cmd $s, Authen::NTLM::ntlm($response);
	if (!/^[123]/) { return $self->Error(LOGINERROR($_)); }
	return;
}}

sub Mail::Sender::Auth::AUTOLOAD {
    (my $auth = $Mail::Sender::Auth::AUTOLOAD) =~ s/.*:://;
	my $self = shift();
	my $s = $self->{'socket'};
	send_cmd $s, "QUIT";
	close $s;
	delete $self->{'socket'};
	return $self->Error( UNKNOWNAUTH($auth));
}

my $debug_code;
sub __Debug {
	my ($socket, $file) = @_;
	if (defined $file) {
		unless (defined @Mail::Sender::DBIO::ISA) {
			eval "use Symbol;";
			eval $debug_code;
			die $@ if $@;
		}
		my $handle = gensym();
		*$handle = \$socket;
		if (! ref $file) {
			my $DEBUG = new FileHandle;
			open $DEBUG, "> $file" or die "Cannot open the debug file '$file': $^E\n";
			binmode $DEBUG;
			$DEBUG->autoflush();
			tie *$handle, 'Mail::Sender::DBIO', $socket, $DEBUG, 1;
		} else {
			my $DEBUG = $file;
			tie *$handle, 'Mail::Sender::DBIO', $socket, $DEBUG, 0;
		}
		bless $handle, 'Mail::Sender::DBIO';
		return $handle;
	} else {
		return $socket;
	}
}

#internale

sub HOSTNOTFOUND {
	$!=2;
	$Mail::Sender::Error="The SMTP server $_[0] was not found";
	return -1, $Mail::Sender::Error;
}

sub SOCKFAILED {
	$Mail::Sender::Error='socket() failed: $^E';
	$!=5;
	return -2, $Mail::Sender::Error;
}

sub CONNFAILED {
	$Mail::Sender::Error="connect() failed: $^E";
	$!=5;
	return -3, $Mail::Sender::Error;
}

sub SERVNOTAVAIL {
	$!=40;
	$Mail::Sender::Error="Service not available. Reply: $_[0]";
	return -4, $Mail::Sender::Error;
}

sub COMMERROR {
	$!=5;
	if ($_[0] eq '') {
		$Mail::Sender::Error="No response from server";
	} else {
		$Mail::Sender::Error="Server error: $_[0]";
	}
	return -5, $Mail::Sender::Error;
}

sub USERUNKNOWN {
	$!=2;
	if ($_[2] and $_[2] !~ /Local user/i) {
		my $err= $_[2];
		$err =~ s/^\d+\s*//;
		$err =~ s/\s*$//s;
		$err ||= "Error";
		$Mail::Sender::Error="$err for \"$_[0]\" on host \"$_[1]\"";
	} else {
		$Mail::Sender::Error="Local user \"$_[0]\" unknown on host \"$_[1]\"";
	}
	return -6, $Mail::Sender::Error;
}

sub TRANSFAILED {
	$!=5;
	$Mail::Sender::Error="Transmission of message failed ($_[0])";
	return -7, $Mail::Sender::Error;
}

sub TOEMPTY {
	$!=14;
	$Mail::Sender::Error="Argument \$to empty";
	return -8, $Mail::Sender::Error;
}

sub NOMSG {
	$!=22;
	$Mail::Sender::Error="No message specified";
	return -9, $Mail::Sender::Error;
}

sub NOFILE {
	$!=22;
	$Mail::Sender::Error="No file name specified";
	return -10, $Mail::Sender::Error;
}

sub FILENOTFOUND {
	$!=2;
	$Mail::Sender::Error="File \"$_[0]\" not found";
	return -11, $Mail::Sender::Error;
}

sub NOTMULTIPART {
	$!=40;
	$Mail::Sender::Error="$_[0] not available in singlepart mode";
	return -12, $Mail::Sender::Error;
}

sub SITEERROR {
	$!=15;
	$Mail::Sender::Error="Site specific error";
	return -13, $Mail::Sender::Error;
}

sub NOTCONNECTED {
	$!=1;
	$Mail::Sender::Error="Connection not established";
	return -14, $Mail::Sender::Error;
}

sub NOSERVER {
	$!=22;
	$Mail::Sender::Error="No SMTP server specified";
	return -15, $Mail::Sender::Error;
}

sub NOFROMSPECIFIED {
	$!=22;
	$Mail::Sender::Error="No From: address specified";
	return -16, $Mail::Sender::Error;
}

sub INVALIDAUTH {
	$!=22;
	$Mail::Sender::Error="Authentication protocol $_[0] is not accepted by the server";
	$Mail::Sender::Error.=",\nresponse: $_[1]" if defined $_[1];
	return -17, $Mail::Sender::Error;
}

sub LOGINERROR {
	$!=22;
	$Mail::Sender::Error="Login not accepted";
	return -18, $Mail::Sender::Error;
}

sub UNKNOWNAUTH {
	$!=22;
	$Mail::Sender::Error="Authentication protocol $_[0] is not implemented by Mail::Sender";
	return -19, $Mail::Sender::Error;
}

sub ALLRECIPIENTSBAD {
	$!=2;
	return -20, $Mail::Sender::Error;
}

sub FILECANTREAD {
	$Mail::Sender::Error="File \"$_[0]\" cannot be read: $^E";
	return -21, $Mail::Sender::Error;
}

sub DEBUGFILE {
	$Mail::Sender::Error=$_[0];
	return -22, $Mail::Sender::Error;
}

@Mail::Sender::Errors = (
	'OK',
	'debug file cannot be opened',
	'file cannot be read',
	'all recipients have been rejected',
	'authentication protocol is not implemented',
	'login not accepted',
	'authentication protocol not accepted by the server',
	'no From: address specified',
	'no SMTP server specified',
	'connection not established. Did you mean MailFile instead of SendFile?',
	'site specific error',
	'not available in singlepart mode',
	'file not found',
	'no file name specified in call to MailFile or SendFile',
	'no message specified in call to MailMsg or MailFile',
	'argument $to empty',
	'transmission of message failed',
	'local user $to unknown on host $smtp',
	'unspecified communication error',
	'service not available',
	'connect() failed',
	'socket() failed',
	'$smtphost unknown'
);

=head1 NAME

Mail::Sender - module for sending mails with attachments through an SMTP server

Version 0.8.16

=head1 SYNOPSIS

 use Mail::Sender;
 $sender = new Mail::Sender
  {smtp => 'mail.yourdomain.com', from => 'your@address.com'};
 $sender->MailFile({to => 'some@address.com',
  subject => 'Here is the file',
  msg => "I'm sending you the list you wanted.",
  file => 'filename.txt'});

=head1 DESCRIPTION

C<Mail::Sender> provides an object oriented interface to sending mails.
It doesn't need any outer program. It connects to a mail server
directly from Perl, using Socket.

Sends mails directly from Perl through a socket connection.

=head1 new Mail::Sender

 new Mail::Sender ([from [,replyto [,to [,smtp [,subject [,headers [,boundary]]]]]]])
 new Mail::Sender {[from => 'somebody@somewhere.com'] , [to => 'else@nowhere.com'] [...]}

Prepares a sender. This doesn't start any connection to the server. You
have to use C<$Sender->Open> or C<$Sender->OpenMultipart> to start
talking to the server.

The parameters are used in subsequent calls to C<$Sender->Open> and
C<$Sender->OpenMultipart>. Each such call changes the saved variables.
You can set C<smtp>, C<from> and other options here and then use the info
in all messages.

=head2 Parameters

=over 4

=item from

C<>=> the sender's e-mail address

=item fake_from

C<>=> the address that will be shown in headers.

If not specified we use the value of C<from>.

=item replyto

C<>=> the reply-to address

=item to

C<>=> the recipient's address(es)

This parameter may be either a comma separated list of email addresses
or a reference to a list of addresses.

=item fake_to

C<>=> the recipient's address that will be shown in headers.
If not specified we use the value of "to".

If the list of addresses you want to send your message to is long or if you do not want
the recipients to see each other's address set the C<fake_to> parameter to some informative,
yet bogus, address or to the address of your mailing/distribution list.

=item cc

C<>=> address(es) to send a copy (CC:) to

=item fake_cc

C<>=> the address that will be shown in headers.

If not specified we use the value of "cc".

=item bcc

C<>=> address(es) to send a copy (BCC: or blind carbon copy).
these addresses will not be visible in the mail!

=item smtp

C<>=> the IP or domain address of your SMTP (mail) server

This is the name of your LOCAL mail server, do NOT try
to contact directly the adressee's mailserver! That would be slow and buggy,
your script should only pass the messages to the nearest mail server and leave
the rest to it. Keep in mind that the recipient's server may be down temporarily.

=item port

C<>=> the TCP/IP port used form the connection. By default getservbyname('smtp', 'tcp')||25.
You should only need to use this option if your mail server waits on a nonstandard port.

=item subject

C<>=> the subject of the message

=item headers

C<>=> the additional headers

You may use this parameter to add custon headers into the message. The parameter may
be either a string containing the headers in the right format or a hash containing the headers
and their values.

=item boundary

C<>=> the message boundary

You usualy do not have to change this, it might only come in handy if you need
to attach a multipart mail created by Mail::Sender to your message as a single part.
Even in that case any problems are unlikely.

=item multipart

C<>=> the MIME subtype for the whole message (Mixed/Related/Alternative)

You may need to change this setting if you want to send a HTML body with some
inline images, or if you want to post the message in plain text as well as
HTML (alternative). See the examples at the end of the docs.
You may also use the nickname "subtype".

Please keep in mind though that it's not currently possible to create nested parts with Mail::Sender.
If you need that level of control you should try MIME::Lite.

=item ctype

C<>=> the content type of a single part message

Please do not confuse these two. The 'multipart' parameter is used to specify
the overall content type of a multipart message (for example a HTML document
with inlined images) while ctype is an ordinary content type for a single
part message. For example a HTML mail message without any inlines.

=item encoding

C<>=> encoding of a single part message or the body of a multipart message.

If the text of the message contains some extended characters or
very long lines you should use 'encoding => "Quoted-printable"' in the
call to Open(), OpenMultipart(), MailMsg() or MailFile().

Keep in mind that if you use some encoding you should either use SendEnc()
or encode the data yourself !

=item charset

C<>=> the charset of the message

=item client

C<>=> the name of the client computer.

During the connection you send
the mailserver your computer's name. By default Mail::Sender sends
C<(gethostbyname 'localhost')[0]>.
If that is not the address you need, you can specify a different one.

=item priority

C<>=> the message priority number

1 = highest, 2 = high, 3 = normal, 4 = low, 5 = lowest

=item confirm

C<>=> whether you request reading or delivery confirmations and to what addresses:

	"delivery" - only delivery, to the C<from> address
	"reading" - only reading, to the C<from> address
	"delivery, reading" - both confirmations, to the C<from> address
	"delivery: my.other@address.com" - only delivery, to my.other@address.com
	...

Keep in mind though that neither of those is guaranteed to work. Some servers/mail clients do not support
this feature and some users/admins may have disabled it. So it's possible that your mail was delivered and read,
but you wount get any confirmation!

=item ESMPT

	ESMTP => {
		NOTIFY => 'SUCCESS,FAILURE,DELAY',
		RET => 'HDRS',
		ORCPT => 'rfc822;my.other@address.com',
		ENVID => 'iuhsdfobwoe8t237',
	}

This option contains data for SMTP extensions, for example it allows you to request delivery
status notifications according to RFC1891.

NOTIFY - to specify the conditions under which a delivery status notification should be generated.
Should be either "NEVER" or a comma separated list of "SUCCESS", "FAILURE"  and "DELAY".

ORCPT - used to convey the "original" (sender-specified) recipient address

RET - to request that Delivery Status Notifications containing an indication of delivery
failure either return the entire contents of a message or only the message headers. Must be either
FULL or HDRS

ENVID - used to propagate an identifier for this message transmission envelope, which is also
known to the sender and will, if present, be returned in any Delivery Status Notifications  issued
for this transmission

You do not need to worry about encoding the ORCPT or ENVID parameters.

If the SMTP server you connect to doesn't support this extension, the options will be ignored.

=item debug

C<>=> C<"/path/to/debug/file.txt">

or

C<>=>  \*FILEHANDLE

or

C<>=> $FH

All the conversation with the server will be logged to that file or handle.
All lines in the file should end with CRLF (the Windows and Internet format).
If even a single one of them does not, please let me know!

If you pass the path to the log file, Mail::Sender will overwrite it. If you want to append to the file,
you have to open it yourself and pass the filehandle:

	open my $DEBUG, ">> /path/to/debug/file.txt"
		or die "Can't open the debug file: $!\n"
	$sender = new Mail::Sender ({
		...
		debug => $DEBUG,
	});

=item debug_level

Only taken into account if the C<debug> option is specified.

	1 - only log the conversation with the server, skip all message data
	2 - log the conversation and message headers
	3 - log the conversation and the message and part headers
	4 - log everything (default)

=item auth

the SMTP authentication protocol to use to login to the server
currently the only ones supported are LOGIN, PLAIN, CRAM-MD5 and NTLM.

Some protocols have module dependencies. CRAM-MD5 depends on
Digest::HMAC_MD5 and NTLM on Authen::NTLM.

You may add support for other authentication protocols yourself. See below.

=item authid

the username used to login to the server

=item authpwd

the password used to login to the server

=item authdomain

the domain name. Used optionaly by the NTLM authentication.

Other authentication protocols may use other options as well.
They should all start with "auth" though.

Please see the authentication section bellow.

=item auth_encoded

If set to a true value the LOGIN authentication assumes the authid and authpwd
is already base64 encoded.

=item keepconnection

If set to a true value causes the Mail::Sender to keep the connection open for several messages.
The connection will be closed if you call the Close() method with a true value or if you call Open,
OpenMultipart, MailMsg or MailFile with the "smtp" parameter.
This means that if you want the object to keep the connection you should pass the "smtp" either to "new Mail::Sender"
or only to the first Open, OpenMultipart, MailMsg or MailFile!

=item skip_bad_recipients

If this option is set to false or not specified then Mail::Sender stops trying to send a message as soon as
the first recipient's address fails. If it is set to a true value Mail::Sender skips the bad addresses and tries
to send the message at least to the good ones. If all addresses are rejected by the server it reports an
"All recipients were rejected" message.

If any addresses were skipped the C<$sender-E<gt>{'skipped_recipients'}> will be a reference to a hash
containing the failed address and the server's response.

=item createmessageid

This option allows you to overwrite the function that generates the message IDs for the emails.
The function gets the "pure" sender's address as it's only parameter and is supposed to return a string.
See the MessageID subroutine in Mail::Sender.pm.

If you want to specify a message id you can also use the "messageid" parameter for the Open, OpenMultipart,
MailMsg or MailFile methods.

=item	on_errors

This option allows you to affect the way Mail::Sender reports errors.

	=> 'die' - raise an exception
	=> 'code' - return the negative error code (default)
	=> 'undef' - return an undef

$Mail::Sender::Error, $sender->{'error'} and $sender->{'error_msg'} are set in all the cases.

All methods return the $sender object if they succeed.

P.S.: The die_on_errors option is deprecated. You may still use it, but it may be removed in future versions!

=back

=head2 Return codes

  ref to a Mail::Sender object =  success

  -1 = $smtphost unknown
  -2 = socket() failed
  -3 = connect() failed
  -4 = service not available
  -5 = unspecified communication error
  -6 = local user $to unknown on host $smtp
  -7 = transmission of message failed
  -8 = argument $to empty
  -9 = no message specified in call to MailMsg or MailFile
  -10 = no file name specified in call to SendFile or MailFile
  -11 = file not found
  -12 = not available in singlepart mode
  -13 = site specific error
  -14 = connection not established. Did you mean MailFile instead of SendFile?
  -15 = no SMTP server specified
  -16 = no From: address specified
  -17 = authentication protocol not accepted by the server
  -18 = login not accepted
  -19 = authentication protocol is not implemented

$Mail::Sender::Error contains a textual description of last error.

=cut

sub new {
	my $this = shift;
	my $self = {};
	my $class;
	if (ref($this)) {
		$class = ref($this);
		%$self = %$this;
	} else {
		$class = $this;
	}
	bless $self, $class;
	return $self->initialize(@_);
}

sub initialize {
		undef $Mail::Sender::Error;
	my $self = shift;

	delete $self->{'_buffer'};
	$self->{'debug'} = 0;
	$self->{'proto'} = (getprotobyname('tcp'))[2];
	$self->{'port'} = getservbyname('smtp', 'tcp')||25 if not defined $self->{'port'};

	$self->{'boundary'} = 'Message-Boundary-by-Mail-Sender-'.time();
	$self->{'multipart'} = 'mixed'; # default is multipart/mixed

	$self->{'client'} = $local_name;

	# Copy defaults from %Mail::Sender::default
	my $key;
	foreach $key (keys %Mail::Sender::default) {
		$self->{lc $key}=$Mail::Sender::default{$key};
	}

	if (@_ != 0) {
		if (ref $_[0] eq 'HASH') {
			my $hash=$_[0];
			foreach $key (keys %$hash) {
				$self->{lc $key}=$hash->{$key};
			}
			$self->{'reply'} = $self->{'replyto'} if (defined $self->{'replyto'} and !defined $self->{'reply'});
		} else {
			($self->{'from'}, $self->{'reply'}, $self->{'to'}, $self->{'smtp'},
			$self->{'subject'}, $self->{'headers'}, $self->{'boundary'}
			) = @_;
		}
	}

	$self->{'fromaddr'} = $self->{'from'};
	$self->{'replyaddr'} = $self->{'reply'};

	$self->_prepare_addresses('to') if defined $self->{'to'};
	$self->_prepare_addresses('cc') if defined $self->{'cc'};
	$self->_prepare_addresses('bcc') if defined $self->{'bcc'};

	$self->_prepare_ESMTP() if defined $self->{'esmtp'};

	$self->{'fromaddr'} =~ s/.*<([^\s]*?)>/$1/ if ($self->{'fromaddr'}); # get from email address
	if (defined $self->{'replyaddr'} and $self->{'replyaddr'}) {
		$self->{'replyaddr'} =~ s/.*<([^\s]*?)>/$1/; # get reply email address
		$self->{'replyaddr'} =~ s/^([^\s]+).*/$1/; # use first address
	}

	if (defined $self->{'smtp'}) {
		$self->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
		$self->{'smtp'} =~ s/\s+$//g;

		$self->{'smtpaddr'} = inet_aton($self->{'smtp'});
		if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }
		$self->{'smtpaddr'} = $1 if ($self->{'smtpaddr'} =~ /(.*)/s); # Untaint
	}

	$self->{'boundary'} =~ tr/=/-/ if defined $self->{'boundary'};

	$self->_prepare_headers() if (exists $self->{'headers'});

	return $self;
}

use vars qw(%CTypes);
%CTypes = (
	GIF => 'image/gif',
	JPE => 'image/jpeg',
	JPEG => 'image/jpeg',
	SHTML => 'text/html',
	SHTM => 'text/html',
	HTML => 'text/html',
	HTM => 'text/html',
	TXT => 'text/plain',
	INI => 'text/plain',
	DOC => 'application/x-msword',
	EML => 'message/rfc822',
);

sub GuessCType {
	my $ext = shift;
	$ext =~ s/^.*\.//;
	return $CTypes{uc $ext} || 'application/octet-stream';
}

sub Connect {
	my $self = shift();

	my $s = IO::Socket::INET->new(
		PeerHost    => $self->{'smtp'},
		PeerPort    => $self->{'port'},
		Proto       => "tcp",
		Timeout     => ($self->{'timeout'} || 120),
	) or return $self->Error(CONNFAILED);

	binmode($s)
		unless ($] >= 5.008);

### <???> Test only!!!
#binmode($s, ":utf8");
###

	$s->autoflush(1);

	if ($self->{'debug'}) {
		eval {
			$s = __Debug( $s, $self->{'debug'});
		}
		or return $self->Error(DEBUGFILE($@));
		$self->{'debug_level'} = 4 unless defined $self->{'debug_level'};
	}

	$_ = get_response($s); if (not $_ or !/^[123]/) { return $self->Error(SERVNOTAVAIL($_)); }
	$self->{'server'} = substr $_, 4;

	{	my $res = $self->say_helo($s);
		return $res if $res;
	}

	if ($self->{'auth'} or $self->{'username'}) {
		$self->{'socket'} = $s;
		my $res = $self->login();
		return $res if $res;
		delete $self->{'socket'}; # it's supposed to be added later
	}

	return $s;
}

sub Error {
	my $self = shift();
	if (@_) {
		if (defined $self->{'socket'}) {
			my $s = $self->{'socket'};
			print $s "quit\x0D\x0A";
			close $s;
			delete $self->{'socket'};
		}
		delete $self->{'_data'};
		($self->{'error'},$self->{'error_msg'}) = @_;
	}
	if ($self->{'die_on_errors'} or $self->{'on_errors'} eq 'die') {
		die $self->{'error_msg'}."\n" ;
	} elsif (exists $self->{'on_errors'} and (!defined($self->{'on_errors'}) or $self->{'on_errors'} eq 'undef')) {
		return
	} else {
		return $self->{'error'}
	}
}

sub ClearErrors {
	my $self = shift();
	delete $self->{'error'};
	delete $self->{'error_msg'};
	undef $Mail::Sender::Error;
}

sub _prepare_addresses {
	my ($self, $type) = @_;
	if (ref $self->{$type}) {
		$self->{$type.'_list'} = $self->{$type};
		$self->{$type} = join ', ', @{$self->{$type.'_list'}};
	} else {
		$self->{$type} =~ s/\s+/ /g;
		$self->{$type} =~ s/, ?,/,/g;
		$self->{$type.'_list'} = [map {s/\s+$//;$_} $self->{$type} =~ /((?:[^",]+|"[^"]*")+)(?:,\s*|\s*$)/g];
	}
}

sub _prepare_ESMTP {
	my $self = shift;
	$self->{esmtp} = {%{$self->{esmtp}}}; # make a copy of the hash. Just in case

	$self->{esmtp}{ORCPT} = 'rfc822;' . $self->{esmtp}{ORCPT} if $self->{esmtp}{ORCPT} ne '' and $self->{esmtp}{ORCPT} !~ /;/;
	for (qw(ENVID ORCPT)) {
		$self->{esmtp}{$_} = encode_qp($self->{esmtp}{$_})
#			if $self->{esmtp}{$_} =~ /[\x00-\x20+=\x7E-\xFFFF]/;
#			if $self->{esmtp}{$_} =~ /[^\x21-\x2A\x2C-\x3C\x3E-\x7E]/; # anything between ! (\x21) and ~ (\x7E) except + and =
			if $self->{esmtp}{$_} =~ /[^!-~]/; # anything between ! (\x21) and ~ (\x7E)
	}
}

sub _prepare_headers {
	my $self = shift;
	return unless exists $self->{'headers'};
	if ($self->{'headers'} eq '') {
		delete $self->{'headers'};
		return;
	}
	for ($self->{'headers'}) {
		if (ref($self->{'headers'}) eq 'HASH') {
			my $headers = '';
			while ( my ($hdr, $value) = each %{$self->{'headers'}}) {
				for ($hdr, $value) {
					s/(?:\x0D\x0A?|\x0A)/\x0D\x0A/sg; # convert all end-of-lines to CRLF
					s/^(?:\x0D\x0A)+//; # strip leading
					s/(?:\x0D\x0A)+$//;	# and trailing end-of-lines
					s/\x0D\x0A(\S)/\x0D\x0A\t$1/sg;
				}
				$headers .= "$hdr: $value\x0D\x0A";
			}
			$headers =~ s/(?:\x0D\x0A)+$//;	# and trailing end-of-lines
			$self->{'headers'} = $headers;
		} elsif (ref($self->{'headers'})) {
		} else {
			s/(?:\x0D\x0A?|\x0A)/\x0D\x0A/sg; # convert all end-of-lines to CRLF
			s/^(?:\x0D\x0A)+//; # strip leading
			s/(?:\x0D\x0A)+$//;	# and trailing end-of-lines
		}
	}
}
=head1 METHODS


=head2 Open

 Open([from [, replyto [, to [, smtp [, subject [, headers]]]]]])
 Open({[from => "somebody@somewhere.com"] , [to => "else@nowhere.com"] [...]})

Opens a new message. If some parameters are unspecified or empty, it uses
the parameters passed to the "C<$Sender=new Mail::Sender(...)>";

See C<new Mail::Sender> for info about the parameters.

The only additional parameter that may not be specified directly in the C<new Mail::Sender>
is messageid. If you set this option then the message will be sent with this Message-ID,
otherwise a new Message ID will be generated out of the sender's address, current date+time
and a random number (or by the function you specified in the C<createmessageid> option).

After the message is sent C<$sender-E<lt>{messageid}> will contain the Message-ID with
which the message was sent.

Returns ref to the Mail::Sender object if successfull.

=cut

sub Open {
		undef $Mail::Sender::Error;
	my $self = shift;
	local $_;
	if (!$self->{'keepconnection'} and $self->{'_data'}) { # the user did not Close() or Cancel() the previous mail
		if ($self->{'error'}) {
			$self->Cancel;
		} else {
			$self->Close;
		}
	}

	delete $self->{'error'};
	delete $self->{'encoding'};
	delete $self->{'messageid'};
	my %changed;
	$self->{'multipart'} = 0;
	$self->{'_had_newline'} = 1;

	if (ref $_[0] eq 'HASH') {
		my $key;
		my $hash=$_[0];
		$hash->{'reply'} = $hash->{'replyto'} if (defined $hash->{'replyto'} and !defined $hash->{'reply'});
		foreach $key (keys %$hash) {
			if (ref($hash->{$key}) eq 'HASH' and exists $self->{lc $key}) {
				$self->{lc $key} = { %{$self->{lc $key}}, %{$hash->{$key}} };
			} else {
				$self->{lc $key} = $hash->{$key};
			}
			$changed{lc $key}=1;
		}
	} else {
		my ($from, $reply, $to, $smtp, $subject, $headers) = @_;

		if ($from) {$self->{'from'}=$from;$changed{'from'}=1;}
		if ($reply) {$self->{'reply'}=$reply;$changed{'reply'}=1;}
		if ($to) {$self->{'to'}=$to;$changed{'to'}=1;}
		if ($smtp) {$self->{'smtp'}=$smtp;$changed{'smtp'}=1;}
		if ($subject) {$self->{'subject'}=$subject;$changed{'subject'}=1;}
		if ($headers) {$self->{'headers'}=$headers;$changed{'headers'}=1;}
	}

	$self->_prepare_addresses('to') if $changed{'to'};
	$self->_prepare_addresses('cc') if $changed{'cc'};
	$self->_prepare_addresses('bcc') if $changed{'bcc'};

	$self->_prepare_ESMTP() if defined $changed{'esmtp'};

	$self->{'boundary'} =~ tr/=/-/ if defined $changed{'boundary'};

	return $self->Error( NOFROMSPECIFIED) unless defined $self->{'from'};

	if ($changed{'from'}) {
		$self->{'fromaddr'} = $self->{'from'};
		$self->{'fromaddr'} =~ s/.*<([^\s]*?)>/$1/; # get from email address
	}

	if ($changed{'reply'}) {
		$self->{'replyaddr'} = $self->{'reply'};
		$self->{'replyaddr'} =~ s/.*<([^\s]*?)>/$1/; # get reply email address
		$self->{'replyaddr'} =~ s/^([^\s]+).*/$1/; # use first address
	}

	if ($changed{'smtp'}) {
		$self->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
		$self->{'smtp'} =~ s/\s+$//g;
		$self->{'smtpaddr'} = inet_aton($self->{'smtp'});
		if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }
		$self->{'smtpaddr'} = $1 if ($self->{'smtpaddr'} =~ /(.*)/s); # Untaint
		if (exists $self->{'socket'}) {
			my $s = $self->{'socket'};
			close $s;
			delete $self->{'socket'};
		}
	}

	$self->_prepare_headers() if ($changed{'headers'});

	if (!$self->{'to'}) { return $self->Error(TOEMPTY); }

	return $self->Error(NOSERVER) unless defined $self->{'smtp'};
#	if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }

	if ($Mail::Sender::{'SiteHook'} and !$self->SiteHook()) {
		return defined $self->{'error'} ? $self->{'error'} : $self->{'error'}=&SITEERROR;
	}

	my $s = $self->{'socket'} || $self->Connect();
	return $s unless ref $s; # return the error number if we did not get a socket
	$self->{'socket'} = $s;

	$_ = send_cmd $s, "MAIL FROM:<$self->{'fromaddr'}>$self->{esmtp}{_MAIL_FROM}";
	if (!/^[123]/) { return $self->Error(COMMERROR($_)); }

	{ local $^W;
		if ($self->{'skip_bad_recipients'}) {
			my $good_count = 0;
			my %failed;
			foreach my $addr ( @{$self->{'to_list'}}, @{$self->{'cc_list'}}, @{$self->{'bcc_list'}}) {
				if ($addr =~ /<(.*)>/) {
					$_ = send_cmd $s, "RCPT TO:<$1>$self->{esmtp}{_RCPT_TO}";
				} else {
					$_ = send_cmd $s, "RCPT TO:<$addr>$self->{esmtp}{_RCPT_TO}";
				}
				if (!/^[123]/) {
					chomp;
					s/^\d{3} //;
					$failed{$addr} = $_;
				} else {
					$good_count++
				}
			}
			$self->{'skipped_recipients'} = \%failed
				if %failed;
			if ($good_count == 0) {
				return $self->Error(ALLRECIPIENTSBAD);
			}
		} else {
			foreach my $addr ( @{$self->{'to_list'}}, @{$self->{'cc_list'}}, @{$self->{'bcc_list'}}) {
				if ($addr =~ /<(.*)>/) {
					$_ = send_cmd $s, "RCPT TO:<$1>$self->{esmtp}{_RCPT_TO}";
				} else {
					$_ = send_cmd $s, "RCPT TO:<$addr>$self->{esmtp}{_RCPT_TO}";
				}
				if (!/^[123]/) {
					return $self->Error(USERUNKNOWN($addr, $self->{'smtp'}, $_)); }
			}
		}
	}

	$_ = send_cmd $s, "DATA";
	if (!/^[123]/) { return $self->Error(COMMERROR($_)); }

	$self->{'socket'}->stop_logging("\x0D\x0A... message headers and data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} <= 1);
	$self->{'_data'} = 1;

	$self->{'ctype'} = 'text/plain' if (defined $self->{'charset'} and !defined $self->{'ctype'});

	my $headers;
	if (defined $self->{'encoding'} or defined $self->{'ctype'}) {
		$headers = 'MIME-Version: 1.0';
		$headers .= "\r\nContent-type: $self->{'ctype'}" if defined $self->{'ctype'};
		$headers .= "; charset=$self->{'charset'}" if defined $self->{'charset'};

		undef $self->{'chunk_size'};
		if (defined $self->{'encoding'}) {
			$headers .= "\r\nContent-transfer-encoding: $self->{'encoding'}";
			if ($self->{'encoding'} =~ /Base64/i) {
				$self->{'code'}=\&enc_base64;
				$self->{'chunk_size'} = $enc_base64_chunk;
			} elsif ($self->{'encoding'} =~ /Quoted[_\-]print/i) {
				$self->{'code'}=\&enc_qp;
			}
		}
	}
	$self->{'code'}=\&enc_plain unless $self->{'code'};

	print_hdr $s, "To" => (defined $self->{'fake_to'} ? $self->{'fake_to'} : $self->{'to'}), $self->{'charset'};
	print_hdr $s, "From" => (defined $self->{'fake_from'} ? $self->{'fake_from'} : $self->{'from'}), $self->{'charset'};
	if (defined $self->{'fake_cc'} and $self->{'fake_cc'}) {
		print_hdr $s, "Cc" => $self->{'fake_cc'}, $self->{'charset'};
	} elsif (defined $self->{'cc'} and $self->{'cc'}) {
		print_hdr $s, "Cc" => $self->{'cc'}, $self->{'charset'};
	}
	print_hdr $s, "Reply-to", $self->{'reply'},$self->{'charset'} if defined $self->{'reply'};

	$self->{'subject'} = "<No subject>" unless defined $self->{'subject'};
	print_hdr $s, "Subject" => $self->{'subject'}, $self->{'charset'};

	unless (defined $Mail::Sender::NO_DATE and $Mail::Sender::NO_DATE) {
		my $date = localtime(); $date =~ s/^(\w+)\s+(\w+)\s+(\d+)\s+(\d+:\d+:\d+)\s+(\d+)$/$1, $3 $2 $5 $4/;
		print_hdr $s, "Date" => "$date $GMTdiff";
	}

	if ($self->{'priority'}) {
		$self->{'priority'} = $priority[$self->{'priority'}]
			if ($self->{'priority'}+0 eq $self->{'priority'});
		print_hdr $s, "X-Priority" => $self->{'priority'};
	}

	if ($self->{'confirm'}) {
		for my $confirm (split /\s*,\s*/, $self->{'confirm'}) {
			if ($confirm =~ /^\s*reading\s*(?:\:\s*(.*))?/i) {
				print_hdr $s, "X-Confirm-Reading-To" => ($1 || $self->{'from'}), $self->{'charset'};
			} elsif ($confirm =~ /^\s*delivery\s*(?:\:\s*(.*))?/i) {
				print_hdr $s, "Return-receipt-to" => ($1 || $self->{'fromaddr'}), $self->{'charset'};
			}
		}
	}

	unless (defined $Mail::Sender::NO_X_MAILER) {
		my $script = basename($0);
		print_hdr $s, "X-Mailer" => qq{Perl script "$script"\r\n\tusing Mail::Sender $Mail::Sender::ver by Jenda Krynicky, Czechlands\r\n\trunning on $local_name ($local_IP)\r\n\tunder account "}.getusername().qq{"\r\n}
	}

	unless (defined $Mail::Sender::NO_MESSAGE_ID and $Mail::Sender::NO_MESSAGE_ID) {
		if (!defined $self->{'messageid'} or $self->{'messageid'} eq '') {
			if (defined $self->{'createmessageid'} and ref $self->{'createmessageid'} eq 'CODE') {
				$self->{'messageid'} = $self->{'createmessageid'}->($self->{'fromaddr'});
			} else {
				$self->{'messageid'} = MessageID($self->{'fromaddr'});
			}
		}
		print_hdr $s, "Message-ID" => $self->{'messageid'};
	}

	print $s $Mail::Sender::SITE_HEADERS,"\x0D\x0A" #<???> should handle \r\n at the end of the headers
		if (defined $Mail::Sender::SITE_HEADERS);

	print $s $self->{'headers'},"\x0D\x0A" if defined $self->{'headers'} and $self->{'headers'};
	print $s $headers,"\r\n" if defined $headers;

	print $s "\r\n";

	$self->{'socket'}->stop_logging("... message data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} <= 2);

	return $self;
}

=head2 OpenMultipart

 OpenMultipart([from [, replyto [, to [, smtp [, subject [, headers [, boundary]]]]]]])
 OpenMultipart({[from => "somebody@somewhere.com"] , [to => "else@nowhere.com"] [...]})

Opens a multipart message. If some parameters are unspecified or empty, it uses
the parameters passed to the C<$Sender=new Mail::Sender(...)>.

See C<new Mail::Sender> for info about the parameters.

Returns ref to the Mail::Sender object if successfull.

=cut

sub OpenMultipart {
	undef $Mail::Sender::Error;
	my $self = shift;

	local $_;
	if (!$self->{'keepconnection'} and $self->{'_data'}) { # the user did not Close() or Cancel() the previous mail
		if ($self->{'error'}) {
			$self->Cancel;
		} else {
			$self->Close;
		}
	}

	delete $self->{'error'};
	delete $self->{'encoding'};
	delete $self->{'messageid'};
	$self->{'_part'} = 0;

	my %changed;
	if (defined $self->{'type'} and $self->{'type'}) {
		$self->{'multipart'} = $1
			if $self->{'type'} =~ m{^multipart/(.*)}i;
	}
	$self->{'multipart'} ='Mixed' unless $self->{'multipart'};
	$self->{'idcounter'} = 0;

	if (ref $_[0] eq 'HASH') {
		my $key;
		my $hash=$_[0];
		$hash->{'multipart'} = $hash->{'subtype'} if defined $hash->{'subtype'};
		$hash->{'reply'} = $hash->{'replyto'} if (defined $hash->{'replyto'} and !defined $hash->{'reply'});
		foreach $key (keys %$hash) {
			if ((ref($hash->{$key}) eq 'HASH') and exists($self->{lc $key})) {
				$self->{lc $key} = { %{$self->{lc $key}}, %{$hash->{$key}} };
			} else {
				$self->{lc $key} = $hash->{$key};
			}
			$changed{lc $key}=1;
		}
	} else {
		my ($from, $reply, $to, $smtp, $subject, $headers, $boundary) = @_;

		if ($from) {$self->{'from'}=$from;$changed{'from'}=1;}
		if ($reply) {$self->{'reply'}=$reply;$changed{'reply'}=1;}
		if ($to) {$self->{'to'}=$to;$changed{'to'}=1;}
		if ($smtp) {$self->{'smtp'}=$smtp;$changed{'smtp'}=1;}
		if ($subject) {$self->{'subject'}=$subject;$changed{'subject'}=1;}
		if ($headers) {$self->{'headers'}=$headers;$changed{'headers'}=1;}
		if ($boundary) {$self->{'boundary'}=$boundary;}
	}

	$self->_prepare_addresses('to') if $changed{'to'};
	$self->_prepare_addresses('cc') if $changed{'cc'};
	$self->_prepare_addresses('bcc') if $changed{'bcc'};

	$self->_prepare_ESMTP() if defined $changed{'esmtp'};

	$self->{'boundary'} =~ tr/=/-/ if $changed{'boundary'};

	$self->_prepare_headers() if ($changed{'headers'});

	return $self->Error( NOFROMSPECIFIED) unless defined $self->{'from'};
	if ($changed{'from'}) {
		$self->{'fromaddr'} = $self->{'from'};
		$self->{'fromaddr'} =~ s/.*<([^\s]*?)>/$1/; # get from email address
	}

	if ($changed{'reply'}) {
		$self->{'replyaddr'} = $self->{'reply'};
		$self->{'replyaddr'} =~ s/.*<([^\s]*?)>/$1/; # get reply email address
		$self->{'replyaddr'} =~ s/^([^\s]+).*/$1/; # use first address
	}

	if ($changed{'smtp'}) {
		$self->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
		$self->{'smtp'} =~ s/\s+$//g;
		$self->{'smtpaddr'} = inet_aton($self->{'smtp'});
		if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }
		$self->{'smtpaddr'} = $1 if ($self->{'smtpaddr'} =~ /(.*)/s); # Untaint
		if (exists $self->{'socket'}) {
			my $s = $self->{'socket'};
			close $s;
			delete $self->{'socket'};
		}
	}

	if (!$self->{'to'}) { return $self->Error(TOEMPTY); }

	return $self->Error(NOSERVER) unless defined $self->{'smtp'};
#	if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }

	if ($Mail::Sender::{'SiteHook'} and !$self->SiteHook()) {
		return defined $self->{'error'} ? $self->{'error'} : $self->{'error'}=&SITEERROR;
	}

	my $s = $self->{'socket'} || $self->Connect();
	return $s unless ref $s; # return the error number if we did not get a socket
	$self->{'socket'} = $s;

	$_ = send_cmd $s, "MAIL FROM:<$self->{'fromaddr'}>$self->{esmtp}{_MAIL_FROM}";
	if (!/^[123]/) { return $self->Error(COMMERROR($_)); }

	{ local $^W;
		if ($self->{'skip_bad_recipients'}) {
			my $good_count = 0;
			my %failed;
			foreach my $addr ( @{$self->{'to_list'}}, @{$self->{'cc_list'}}, @{$self->{'bcc_list'}}) {
				if ($addr =~ /<(.*)>/) {
					$_ = send_cmd $s, "RCPT TO:<$1>$self->{esmtp}{_RCPT_TO}";
				} else {
					$_ = send_cmd $s, "RCPT TO:<$addr>$self->{esmtp}{_RCPT_TO}";
				}
				if (!/^[123]/) {
					s/^\d{3} //;
					$failed{$addr} = $_;
				} else {
					$good_count++
				}
			}
			$self->{'skipped_recipients'} = \%failed
				if %failed;
			if ($good_count == 0) {
				return $self->Error(ALLRECIPIENTSBAD);
			}
		} else {
			foreach my $addr ( @{$self->{'to_list'}}, @{$self->{'cc_list'}}, @{$self->{'bcc_list'}}) {
				if ($addr =~ /<(.*)>/) {
					$_ = send_cmd $s, "RCPT TO:<$1>$self->{esmtp}{_RCPT_TO}";
				} else {
					$_ = send_cmd $s, "RCPT TO:<$addr>$self->{esmtp}{_RCPT_TO}";
				}
				if (!/^[123]/) {
					return $self->Error(USERUNKNOWN($addr, $self->{'smtp'}, $_));
				}
			}
		}
	}

	$_ = send_cmd $s, "DATA";
	if (!/^[123]/) { return $self->Error(COMMERROR($_)); }

	$self->{'socket'}->stop_logging("\x0D\x0A... message headers and data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} <= 1);
	$self->{'_data'} = 1;

	print_hdr $s, "To" => (defined $self->{'fake_to'} ? $self->{'fake_to'} : $self->{'to'}), $self->{'charset'};
	print_hdr $s, "From" => (defined $self->{'fake_from'} ? $self->{'fake_from'} : $self->{'from'}), $self->{'charset'};
	if (defined $self->{'fake_cc'} and $self->{'fake_cc'}) {
		print_hdr $s, "Cc" => $self->{'fake_cc'}, $self->{'charset'};
	} elsif (defined $self->{'cc'} and $self->{'cc'}) {
		print_hdr $s, "Cc" => $self->{'cc'}, $self->{'charset'};
	}
	print_hdr $s, "Reply-to" => $self->{'reply'}, $self->{'charset'} if defined $self->{'reply'};

	$self->{'subject'} = "<No subject>" unless defined $self->{'subject'};
	print_hdr $s, "Subject" => $self->{'subject'}, $self->{'charset'};

	unless (defined $Mail::Sender::NO_DATE and $Mail::Sender::NO_DATE) {
		my $date = localtime(); $date =~ s/^(\w+)\s+(\w+)\s+(\d+)\s+(\d+:\d+:\d+)\s+(\d+)$/$1, $3 $2 $5 $4/;
		print_hdr $s, "Date" => "$date $GMTdiff";
	}

	if ($self->{'priority'}) {
		$self->{'priority'} = $priority[$self->{'priority'}]
			if ($self->{'priority'}+0 eq $self->{'priority'});
		print_hdr $s, "X-Priority" => $self->{'priority'};
	}

	if ($self->{'confirm'}) {
		for my $confirm (split /\s*,\s*/, $self->{'confirm'}) {
			if ($confirm =~ /^\s*reading\s*(?:\:\s*(.*))?/i) {
				print_hdr $s, "X-Confirm-Reading-To" => ($1 || $self->{'from'}), $self->{'charset'};
			} elsif ($confirm =~ /^\s*delivery\s*(?:\:\s*(.*))?/i) {
				print_hdr $s, "Return-receipt-to" => ($1 || $self->{'fromaddr'}), $self->{'charset'};
			}
		}
	}

	unless (defined $Mail::Sender::NO_X_MAILER and $Mail::Sender::NO_X_MAILER) {
		my $script = basename($0);
		print_hdr $s, "X-Mailer" => qq{Perl script "$script"\r\n\tusing Mail::Sender $Mail::Sender::ver by Jenda Krynicky, Czechlands\r\n\trunning on $local_name ($local_IP)\r\n\tunder account "}.getusername().qq{"\r\n}
	}

	print $s $Mail::Sender::SITE_HEADERS,"\r\n"
		if (defined $Mail::Sender::SITE_HEADERS);

	unless (defined $Mail::Sender::NO_MESSAGE_ID and $Mail::Sender::NO_MESSAGE_ID) {
		if (!defined $self->{'messageid'} or $self->{'messageid'} eq '') {
			if (defined $self->{'createmessageid'} and ref $self->{'createmessageid'} eq 'CODE') {
				$self->{'messageid'} = $self->{'createmessageid'}->($self->{'fromaddr'});
			} else {
				$self->{'messageid'} = MessageID($self->{'fromaddr'});
			}
		}
		print_hdr $s, "Message-ID" => $self->{'messageid'};
	}

	print $s $self->{'headers'},"\r\n" if defined $self->{'headers'} and $self->{'headers'};
	print $s "MIME-Version: 1.0\r\n";
	print_hdr $s, "Content-type", qq{multipart/$self->{'multipart'};\r\n\tboundary="$self->{'boundary'}"};

	print $s "\r\n";
	$self->{'socket'}->stop_logging("... message data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} <= 2);

	print $s "This message is in MIME format. Since your mail reader does not understand\r\n"
		. "this format, some or all of this message may not be legible.\r\n"
		. "\r\n--$self->{'boundary'}\r\n";

	return $self;
}

sub Connected {
	my $self = shift();
	return unless exists $self->{'socket'};
	my $s = $self->{'socket'};
	return $s->opened();
}



=head2 MailMsg

 MailMsg([from [, replyto [, to [, smtp [, subject [, headers]]]]]], message)
 MailMsg({[from => "somebody@somewhere.com"]
          [, to => "else@nowhere.com"] [...], msg => "Message"})

Sends a message. If a mail in $sender is opened it gets closed
and a new mail is created and sent. $sender is then closed.
If some parameters are unspecified or empty, it uses
the parameters passed to the "C<$Sender=new Mail::Sender(...)>";

See C<new Mail::Sender> for info about the parameters.

The module was made so that you could create an object initialized with
all the necesary options and then send several messages without need to
specify the SMTP server and others each time. If you need to send only
one mail using MailMsg() or MailFile() you do not have to create a named
object and then call the method. You may do it like this :

 (new Mail::Sender)->MailMsg({smtp => 'mail.company.com', ...});

Returns ref to the Mail::Sender object if successfull.

=cut

sub MailMsg {
	my $self = shift;
	my $msg;
	local $_;
	if (ref $_[0] eq 'HASH') {
		my $hash=$_[0];
		$msg=$hash->{'msg'};
	} else {
		$msg = pop;
	}
	return $self->Error(NOMSG) unless $msg;

	if (ref $self->Open(@_)
		and
		ref $self->SendEnc($msg)
		and
		ref $self->Close()
	) {
		return $self
	} else {
		return $self->{'error'}
	}
}


=head2 MailFile

 MailFile([from [, replyto [, to [, smtp [, subject [, headers]]]]]], message, file(s))
 MailFile({[from => "somebody@somewhere.com"]
           [, to => "else@nowhere.com"] [...],
           msg => "Message", file => "File"})

Sends one or more files by mail. If a mail in $sender is opened it gets closed
and a new mail is created and sent. $sender is then closed.
If some parameters are unspecified or empty, it uses
the parameters passed to the "C<$Sender=new Mail::Sender(...)>";

The C<file> parameter may be a "filename", a "list, of, file, names" or a \@list_of_file_names.

see C<new Mail::Sender> for info about the parameters.

Just keep in mind that parameters like ctype, charset and encoding
will be used for the attached file, not the body of the message.
If you want to specify those parameters for the body you have to use
b_ctype, b_charset and b_encoding. Sorry.

Returns ref to the Mail::Sender object if successfull.

=cut

sub MailFile {
	my $self = shift;
	my $msg;
	local $_;
	my ($file, $desc, $haddesc,$ctype,$charset,$encoding);
	my @files;
	if (ref $_[0] eq 'HASH') {
		my $hash = $_[0];
		$msg = $hash->{'msg'};
#		delete $hash->{'msg'};

		$file=$hash->{'file'};
#		delete $hash->{'file'};

		$desc=$hash->{'description'}; $haddesc = 1 if defined $desc;
#		delete $hash->{'description'};

		$ctype=$hash->{'ctype'};
#		delete $hash->{'ctype'};

		$charset=$hash->{'charset'};
#		delete $hash->{'charset'};

		$encoding=$hash->{'encoding'};
#		delete $hash->{'encoding'};

	} else {
		$desc=pop if ($#_ >=2); $haddesc = 1 if defined $desc;
		$file = pop;
		$msg = pop;
	}
	return $self->Error(NOMSG) unless $msg;
	return $self->Error(NOFILE) unless $file;

	if (ref $file eq 'ARRAY') {
		@files=@$file;
	} elsif ($file =~ /,/) {
		@files=split / *, */,$file;
	} else {
		@files = ($file);
	}
	foreach $file (@files) {
		return $self->Error(FILENOTFOUND($file)) unless ($file =~ /^&/ or -e $file);
	}

	ref $self->OpenMultipart(@_)
	and
	ref $self->Body(
		$self->{'b_charset'}||$self->{'charset'},
		$self->{'b_encoding'},
		$self->{'b_ctype'}
	)
	and
	$self->SendEnc($msg)
	or return $self->{'error'};

	$Mail::Sender::Error = '';
	foreach $file (@files) {
		my $cnt;
		my $filename = basename $file;
		my $ctype = $ctype || GuessCType $filename, $file;
		my $encoding = $encoding || ($ctype =~ m#^text/#i ? 'Quoted-printable' : 'Base64');

		$desc = $filename unless (defined $haddesc);

		$self->Part({encoding => $encoding,
				   disposition => (defined $self->{'disposition'} ? $self->{'disposition'} : "attachment; filename=\"$filename\""),
				   ctype => "$ctype; name=\"$filename\"" . (defined $charset ? "; charset=$charset" : ''),
				   description => $desc});

		my $code = $self->{'code'};

		my $FH = new FileHandle;
		open $FH, "<", $file
			or return $self->Error(FILECANTREAD($file));
		binmode $FH unless $ctype =~ m#^text/#i and $encoding =~ /Quoted[_\-]print|Base64/i;
		my $s;
		$s = $self->{'socket'};
		my $mychunksize = $chunksize;
		$mychunksize = $chunksize64 if defined $self->{'chunk_size'};
		while (read $FH, $cnt, $mychunksize) {
			$cnt = &$code($cnt);
			$cnt =~ s/^\.\././ unless $self->{'_had_newline'};
			print $s $cnt;
			$self->{'_had_newline'} = ($cnt =~ /[\n\r]$/);
		}
		close $FH;
	}

	if ($Mail::Sender::Error eq '') {
		undef $Mail::Sender::Error;
	} else {
		chomp $Mail::Sender::Error;
	}
	return $self->Close;
}



=head2 Send

 Send(@strings)

Prints the strings to the socket. Doesn't add any end-of-line characters.
Doesn't encode the data! You should use C<\r\n> as the end-of-line!

UNLESS YOU ARE ABSOLUTELY SURE YOU KNOW WHAT YOU ARE DOING
YOU SHOULD USE SendEnc() INSTEAD!

Returns the object if successfull.

=cut

sub Send {
	my $self = shift;
	my $s;
	$s = $self->{'socket'};
	print $s @_;
	return $self;
}

=head2 SendLine

 SendLine(@strings)

Prints the strings to the socket. Adds the end-of-line character at the end.
Doesn't encode the data! You should use C<\r\n> as the end-of-line!

UNLESS YOU ARE ABSOLUTELY SURE YOU KNOW WHAT YOU ARE DOING
YOU SHOULD USE SendLineEnc() INSTEAD!

Returns the object if successfull.

=cut

sub SendLine {
	my $self = shift;
	my $s = $self->{'socket'};
	print $s (@_,"\x0D\x0A");
	return $self;
}

=head2 print

Alias to SendEnc().

Keep in mind that you can't write :

	print $sender "...";

you have to use

	$sender->print("...");

If you want to be able to print into the message as if it was a normal file handle take a look at C<GetHandle>()

=head2 SendEnc

 SendEnc(@strings)

Prints the strings to the socket. Doesn't add any end-of-line characters.

Encodes the text using the selected encoding (none/Base64/Quoted-printable)

Returns the object if successfull.

=cut

sub SendEnc {
	my $self = shift;
	local $_;
	my $code = $self->{'code'};
	$self->{'code'}= $code = \&enc_plain
		unless defined $code;
	my $s;
	$s = $self->{'socket'}
		or return $self->Error(NOTCONNECTED);
	if (defined $self->{'chunk_size'}) {
		my $str;
		my $chunk = $self->{'chunk_size'};
		if (defined $self->{'_buffer'}) {
			$str=(join '',($self->{'_buffer'},@_));
		} else {
			$str=join '',@_;
		}
		my ($len,$blen);
		$len = length $str;
		if (($blen=($len % $chunk)) >0) {
			$self->{'_buffer'} = substr($str,($len-$blen));
			print $s (&$code(substr( $str,0,$len-$blen)));
		} else {
			delete $self->{'_buffer'};
			print $s (&$code($str));
		}
	} else {
		my $encoded = &$code(join('',@_));
		$encoded =~ s/^\.\././ unless $self->{'_had_newline'};
		print $s $encoded;
		$self->{'_had_newline'} = ($_[-1] =~ /[\n\r]$/);
	}
	return $self;
}

sub print;*print = \&SendEnc;

=head2 SendLineEnc

 SendLineEnc(@strings)

Prints the strings to the socket and adds the end-of-line character at the end.
Encodes the text using the selected encoding (none/Base64/Quoted-printable).

Do NOT mix up /Send(Line)?(Ex)?/ and /Send(Line)?Enc/! SendEnc does some buffering
necessary for correct Base64 encoding, and /Send(Ex)?/ is not aware of that!

Usage of /Send(Line)?(Ex)?/ in non xBIT parts not recommended.
Using C<Send(encode_base64($string))> may work, but more likely it will not!
In particular if you use several such to create one part,
the data is very likely to get crippled.

Returns the object if successfull.

=cut

sub SendLineEnc {
	push @_, "\r\n";
	goto &SendEnc;
}

=head2 SendEx

 SendEx(@strings)

Prints the strings to the socket. Doesn't add any end-of-line characters.
Changes all end-of-lines to C<\r\n>. Doesn't encode the data!

UNLESS YOU ARE ABSOLUTELY SURE YOU KNOW WHAT YOU ARE DOING
YOU SHOULD USE SendEnc() INSTEAD!

Returns the object if successfull.

=cut

sub SendEx {
	my $self = shift;
	my $s;
	$s = $self->{'socket'}
		or return $self->Error(NOTCONNECTED);
	my $str;my @data = @_;
	foreach $str (@data) {
		$str =~ s/(?:\x0D\x0A?|\x0A)/\x0D\x0A/sg;
		$str =~ s/^\./../mg;
	}
	print $s @data;
	return $self;
}

=head2 SendLineEx

 SendLineEx(@strings)

Prints the strings to the socket. Adds an end-of-line character at the end.
Changes all end-of-lines to C<\r\n>. Doesn't encode the data!

UNLESS YOU ARE ABSOLUTELY SURE YOU KNOW WHAT YOU ARE DOING
YOU SHOULD USE SendEnc() INSTEAD!

Returns the object if successfull.

=cut

sub SendLineEx {
	push @_, "\r\n";
	goto &SendEx;
}


=head2 Part

 Part( I<description>, I<ctype>, I<encoding>, I<disposition> [, I<content_id> [, I<msg>]]);
 Part( {[description => "desc"], [ctype => "content-type"], [encoding => "..."],
     [disposition => "..."], [content_id => "..."], [msg => ...]});

Prints a part header for the multipart message and (if specified) the contents.
The undefined or empty variables are ignored.

=over 2

=item description

The title for this part.

=item ctype

the content type (MIME type) of this part. May contain some other
parameters, such as B<charset> or B<name>.

Defaults to "application/octet-stream".

Since 0.8.00 you may use even "multipart/..." types. Such a multipart part should be
closed by a call to $sender->EndPart($ctype).

	...
	$sender->Part({ctype => "multipart/related", ...});
		$sender->Part({ctype => 'text/html', ...});
		$sender->Attach({file => 'some_image.gif', content_id => 'foo', ...});
	$sender->EndPart("multipart/related");
	...

Please see the examples below.

=item encoding

the encoding used for this part of message. Eg. Base64, Uuencode, 7BIT
...

Defaults to "7BIT".

=item disposition

This parts disposition. Eg: 'attachment; filename="send.pl"'.

Defaults to "attachment". If you specify "none" or "", the
Content-disposition: line will not be included in the headers.

=item content_id

The content id of the part, used in multipart/related.
If not specified, the header is not included.

=item msg

The content of the part. You do not have to specify the content here, you may use SendEnc()
to add content to the part.

=item charset

The charset of the part.

=back

Returns the Mail::Sender object if successfull, negative error code if not.

=cut

sub Part {
	my $self = shift;
	local $_;
	if (! $self->{'multipart'}) { return $self->Error(NOTMULTIPART("\$sender->Part()")); }
	$self->EndPart();

	my ($description, $ctype, $encoding, $disposition, $content_id, $msg, $charset);
	if (ref $_[0] eq 'HASH') {
		my $hash=$_[0];
		$description=$hash->{'description'};
		$ctype=$hash->{'ctype'};
		$encoding=$hash->{'encoding'};
		$disposition=$hash->{'disposition'};
		$content_id = $hash->{'content_id'};
		$msg = $hash->{'msg'};
		$charset = $hash->{'charset'};
	} else {
		($description, $ctype, $encoding, $disposition, $content_id, $msg) = @_;
	}

	$ctype = "application/octet-stream" unless defined $ctype;
	$disposition = "attachment" unless defined $disposition;
	$encoding="7BIT" unless defined $encoding;
	$self->{'encoding'} = $encoding;
	if (defined $charset and $charset and $ctype !~ /charset=/i) {
		$ctype .= qq{; charset="$charset"}
	}

	my $s;
	$s = $self->{'socket'}
		or return $self->Error(NOTCONNECTED);

	undef $self->{'chunk_size'};
	if ($encoding =~ /Base64/i) {
		$self->{'code'}=\&enc_base64;
		$self->{'chunk_size'} = $enc_base64_chunk;
	} elsif ($encoding =~ /Quoted[_\-]print/i) {
		$self->{'code'}=\&enc_qp;
	} else {
		$self->{'code'}=\&enc_plain;
	}

	$self->{'socket'}->start_logging() if ($self->{'debug'} and $self->{'debug_level'} == 3);

	if ($ctype =~ m{^multipart/}i) {
		$self->{'_part'}+=2;
		print $s "Content-Type: $ctype; boundary=\"Part-$self->{'boundary'}_$self->{'_part'}\"\r\n\r\n";
	} else {
		$self->{'_part'}++;
		print $s "Content-type: $ctype\r\n";
		if ($description) {print $s "Content-description: $description\r\n";}
		print $s "Content-transfer-encoding: $encoding\r\n";
		print $s "Content-disposition: $disposition\r\n" unless $disposition eq '' or uc($disposition) eq 'NONE';
		print $s "Content-ID: <$content_id>\r\n" if (defined $content_id);
		print $s "\r\n";

		$self->{'socket'}->stop_logging("... data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} == 3);
		$self->SendEnc($msg) if defined $msg;
	}

	#$self->{'_had_newline'} = 1;
	return $self;
}


=head2 Body

 Body([charset [, encoding [, content-type]]]);
 Body({charset => '...', encoding => '...', ctype => '...', msg => '...');

Sends the head of the multipart message body. You can specify the
charset and the encoding. Default is "US-ASCII","7BIT",'text/plain'.

If you pass undef or zero as the parameter, this function uses the default
value:

    Body(0,0,'text/html');

Returns the Mail::Sender object if successfull, negative error code if not.
You should NOT use this method in single part messages, that is, it works after OpenMultipart(),
but has no meaning after Open()!

=cut

sub Body {
	my $self = shift;
	if (! $self->{'multipart'}) {
		# ->Body() has no meanin in singlepart messages
		if (@_) {
			# they called it with some parameters? Too late for them, let's scream.
			return $self->Error(NOTMULTIPART("\$sender->Body()"));
		} else {
			# $sender->Body() ... OK, let's ignore it.
			return $self;
		}
	}
	my $hash;
	$hash = shift() if (ref $_[0] eq 'HASH');
	my $charset = shift || $hash->{'charset'} || 'US-ASCII';
	my $encoding = shift || $hash->{'encoding'} || $self->{'encoding'} || '7BIT';
	my $ctype = shift || $hash->{'ctype'} || $self->{'ctype'} || 'text/plain';

	$ctype .= qq{; charset="$charset"}
		unless $ctype =~ /charset=/i;

	$self->{'encoding'} = $encoding;
	$self->{'ctype'} = $ctype;

	$self->Part("Mail message body", $ctype,
		$encoding, 'inline', undef, $hash->{'msg'});
	return $self;
}

=head2 SendFile

Alias to Attach()

=head2 Attach

 Attach( I<description>, I<ctype>, I<encoding>, I<disposition>, I<file>);
 Attach( { [description => "desc"] , [ctype => "ctype"], [encoding => "encoding"],
             [disposition => "disposition"], file => "file"});

 Sends a file as a separate part of the mail message. Only in multipart mode.

=over 2

=item description

The title for this part.

=item ctype

the content type (MIME type) of this part. May contain some other
parameters, such as B<charset> or B<name>.

Defaults to "application/octet-stream".

=item encoding

the encoding used for this part of message. Eg. Base64, Uuencode, 7BIT
...

Defaults to "Base64".

=item disposition

This parts disposition. Eg: 'attachment; filename="send.pl"'. If you use
'attachment; filename=*' the * will be replaced by the respective names
of the sent files.

Defaults to "attachment; filename=*". If you do not want to include this header use
"" as the value.

=item file

The name of the file to send or a 'list, of, names' or a
['reference','to','a','list','of','filenames']. Each file will be sent as
a separate part.

Please keep in mind that if you pass a string as this parameter the module
will split it on commas! If your filenames may contain commas and you
want to be sure they are sent correctly you have to use the reference to array
format:

	file => [ $filename],

=item content_id

The content id of the message part. Used in multipart/related.

 Special values:
  "*" => the name of the file
  "#" => autoincremented number (starting from 0)

=back

Returns the Mail::Sender object if successfull, negative error code if not.

=cut

sub SendFile {
	my $self = shift;
	local $_;
	if (! $self->{'multipart'}) { return $self->Error(NOTMULTIPART("\$sender->SendFile()")); }
	if (! $self->{'socket'}) { return $self->Error(NOTCONNECTED); }

	my ($description, $ctype, $encoding, $disposition, $file, $content_id, @files);
	if (ref $_[0] eq 'HASH') {
		my $hash=$_[0];
		$description=$hash->{'description'};
		$ctype=$hash->{'ctype'};
		$encoding=$hash->{'encoding'};
		$disposition=$hash->{'disposition'};
		$file=$hash->{'file'};
		$content_id=$hash->{'content_id'};
	} else {
		($description, $ctype, $encoding, $disposition, $file, $content_id) = @_;
	}
	return ($self->{'error'}=NOFILE) unless $file;

	if (ref $file eq 'ARRAY') {
		@files=@$file;
	} elsif ($file =~ /,/) {
		@files=split / *, */,$file;
	} else {
		@files = ($file);
	}
	foreach $file (@files) {
		return $self->Error(FILENOTFOUND($file)) unless ($file =~ /^&/ or -e $file);
	}

	$disposition = "attachment; filename=*" unless defined $disposition;
	$encoding='Base64' unless $encoding;

	my $s=$self->{'socket'};

	if ($self->{'_buffer'}) {
		my $code = $self->{'code'};
		print $s (&$code($self->{'_buffer'}));
		delete $self->{'_buffer'};
	}

	my $code;
	if ($encoding =~ /Base64/i) {
		$code=\&enc_base64;
	} elsif ($encoding =~ /Quoted[_\-]print/i) {
		$code=\&enc_qp;
	} else {
		$code=\&enc_plain;
	}
	$self->{'code'}=$code;

	foreach $file (@files) {
		$self->EndPart();$self->{'_part'}++;
		$self->{'encoding'} = $encoding;
		my $cnt='';
		my $name =  basename $file;
		my $fctype = $ctype ? $ctype : GuessCType $name, $file;
		$self->{'ctype'} = $fctype;

		$self->{'socket'}->start_logging() if ($self->{'debug'} and $self->{'debug_level'} == 3);

		if ($fctype =~ /;\s*name=/) {
			print $s ("Content-type: $fctype\r\n");
		} else {
			print $s ("Content-type: $fctype; name=\"$name\"\r\n");
		}

		if ($description) {print $s ("Content-description: $description\r\n");}
		print $s ("Content-transfer-encoding: $encoding\r\n");

		if ($disposition =~ /^(.*)filename=\*(.*)$/i) {
			print $s ("Content-disposition: ${1}filename=\"$name\"$2\r\n");
		} elsif ($disposition and uc($disposition) ne 'NONE') {
			print $s ("Content-disposition: $disposition\r\n");
		}

		if ($content_id) {
			if ($content_id eq '*') {
				print $s ("Content-ID: <$name>\r\n");
			} elsif ($content_id eq '#') {
				print $s ("Content-ID: <id".$self->{'idcounter'}++.">\r\n");
			} else {
				print $s ("Content-ID: <$content_id>\r\n");
			}
		}
		print $s "\r\n";

		$self->{'socket'}->stop_logging("... data skipped ...") if ($self->{'debug'} and $self->{'debug_level'} == 3);

		my $FH = new FileHandle;
		open $FH, "<", $file
			or return $self->Error(FILECANTREAD($file));
		binmode $FH unless $fctype =~ m#^text/#i and $encoding =~ /Quoted[_\-]print|Base64/i;

		my $mychunksize = $chunksize;
		$mychunksize = $chunksize64 if lc($encoding) eq "base64";
		my $s;
		$s = $self->{'socket'}
			or return $self->Error(NOTCONNECTED);
		while (read $FH, $cnt, $mychunksize) {
			print $s (&$code($cnt));
		}
		close $FH;
	}

	return $self;
}

sub Attach; *Attach = \&SendFile;

=head2 EndPart

 $sender->EndPart($ctype);

Closes a multipart part.

If the $ctype is not present or evaluates to false, only the current SIMPLE part is closed!
Don't do that unless you are really sure you know what you are doing.

It's best to always pass to the ->EndPart() the content type of the corresponding ->Part().

=cut

sub EndPart {
	my $self = shift;
	return unless $self->{'_part'};
	my $end = shift();
	my $s;
	my $LN = "\x0D\x0A";
	$s = $self->{'socket'}
		or return $self->Error(NOTCONNECTED);
	# flush the buffer (if it contains anything)
	if ($self->{'_buffer'}) { # used only for base64
		my $code = $self->{'code'};
		if (defined $code) {
			print $s (&$code($self->{'_buffer'}));
		} else {
			print $s ($self->{'_buffer'});
		}
		delete $self->{'_buffer'};
	}
	if ($self->{'_had_newline'}) {
		$LN = '';
	} else {
		print $s "=" if !$self->{'bypass_outlook_bug'} and $self->{'encoding'} =~ /Quoted[_\-]print/i; # make sure we do not add a newline
	}

	$self->{'socket'}->start_logging() if ($self->{'debug'} and $self->{'debug_level'} == 3);

	if ($self->{'_part'}>1) { # end of a subpart
		print $s "$LN--Part-$self->{'boundary'}_$self->{'_part'}",
			($end ? "--" : ()),
			"\r\n";
	} else {
		print $s "$LN--$self->{'boundary'}",
			($end ? "--" : ()),
			"\r\n";
	}

	$self->{'_part'}--;
	$self->{'code'}=\&enc_plain;
	$self->{'encoding'} = '';
	return $self;
}

=head2 Close

 $sender->Close;
 $sender->Close(1);

Close and send the email message. If you pass a true value to the method the connection will be closed even
if the "keepconnection" was specified. You should only keep the connection open if you plan to send another
message immediately. And you should not keep it open for hundreds of emails even if you do send them all in a row.

This method should be called automatically when destructing the object, but you should not rely on it. If you want to be sure
your message WAS processed by the SMTP server you SHOULD call Close() explicitely.

Returns the Mail::Sender object if successfull, negative error code if not, zero if $sender was not connected at all.
The zero usualy means that the Open/OpenMultipart failed and you did not test its return value.

=cut

sub Close {
	my $self = shift;
	local $_;
	my $s = $self->{'socket'};
	return 0 unless $s;

	if ($self->{'_data'}) {
		# flush the buffer (if it contains anything)
		if ($self->{'_buffer'}) {
			my $code = $self->{'code'};
			if (defined $code) {
				print $s (&$code($self->{'_buffer'}));
			} else {
				print $s ($self->{'_buffer'});
			}
			delete $self->{'_buffer'};
		}

		if ($self->{'_part'}) {
			while ($self->{'_part'}) {
				$self->EndPart(1);
			}
		}

		$self->{'socket'}->start_logging() if ($self->{'debug'});
		print $s "\r\n.\r\n" ;
		$self->{'_data'} = 0;
		$_ = get_response($s); if (/^[45]\d* (.*)$/) { return $self->Error(TRANSFAILED($1)); }
		$self->{message_response} = $_;
	}

	delete $self->{'encoding'};
	delete $self->{'ctype'};

	if ($_[0] or !$self->{'keepconnection'}) {
		$_ = send_cmd $s, "QUIT";
		if (!/^[123]/) { return $self->Error(COMMERROR($_)); }
		close $s;
		delete $self->{'socket'};
		delete $self->{'debug'};
	}
	return $self;
}

=head2 Cancel

 $sender->Cancel;

Cancel an opened message.

SendFile and other methods may set $sender->{'error'}.
In that case "undef $sender" calls C<$sender->>Cancel not C<$sender->>Close!!!

Returns the Mail::Sender object if successfull, negative error code if not.

=cut

sub Cancel {
	my $self = shift;
	my $s;
	$s = $self->{'socket'}
		or return $self->Error(NOTCONNECTED);
	close $s;
	delete $self->{'socket'};
	delete $self->{'error'};
	return $self;
}

sub DESTROY {
	return if ref($_[0]) ne 'Mail::Sender';
	my $self = shift;
	if (defined $self->{'socket'}) {
		delete $self->{'keepconnection'};
		$self->Close;
	}
}

sub MessageID {
	my $from = shift;
	my ($sec,$min,$hour,$mday,$mon,$year)
		= gmtime(time);
	$mon++;$year+=1900;

	return sprintf "<%04d%02d%02d_%02d%02d%02d_%06d.%s>",
	$year,$mon,$mday,$hour,$min,$sec,rand(100000),
	$from;
}

=head2 QueryAuthProtocols

	@protocols = $sender->QueryAuthProtocols();
	@protocols = $sender->QueryAuthProtocols( $smtpserver);


Queryies the server (specified either in the default options for Mail::Sender,
the "new Mail::Sender" command or as a parameter to this method for
the authentication protocols it supports.

=cut

sub QueryAuthProtocols {
	my $self = shift;
	local $_;
	if (!defined $self) {
		croak "Mail::Sender::QueryAuthProtocols() called without any parameter!";
	} elsif (ref $self) { # $sender->QueryAuthProtocols() or $sender->QueryAuthProtocols('the.server.com)
		if ($self->{'socket'}) { # the user did not Close() or Cancel() the previous mail
			die "You forgot to close the mail before calling QueryAuthProtocols!\n"
		}
		if (@_) {
			$self->{'smtp'} = shift();
			$self->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
			$self->{'smtp'} =~ s/\s+$//g;
			$self->{'smtpaddr'} = inet_aton($self->{'smtp'});
			if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }
			$self->{'smtpaddr'} = $1 if ($self->{'smtpaddr'} =~ /(.*)/s); # Untaint
		}
	} elsif ($self =~ /::/) { # Mail::Sender->QueryAuthProtocols('the.server.com')
		croak "Mail::Sender->QueryAuthProtocols() called without any parameter!"
			if ! @_;
		$self = new Mail::Sender {smtp => $_[0]};
		return unless ref $self;
	} else { # Mail::Sender::QueryAuthProtocols('the.server.com')
		$self = new Mail::Sender {smtp => $self};
		return unless ref $self;
	}

	return $self->Error(NOSERVER) unless defined $self->{'smtp'};

	my $s = IO::Socket::INET->new(
		PeerHost    => $self->{'smtp'},
		PeerPort    => $self->{'port'},
		Proto       => "tcp",
		Timeout     => $self->{'timeout'} || 120,
	) or return $self->Error(CONNFAILED);

	binmode($s)
		unless ($] >= 5.008);

	$s->autoflush(1);

	$_ = get_response($s); if (not $_ or !/^[123]/) { return $self->Error(SERVNOTAVAIL($_)); }
	$self->{'server'} = substr $_, 4;

	{	my $res = $self->say_helo($s);
		return $res if $res;
	}

	$_ = send_cmd $s, "QUIT";
	close $s;
	delete $self->{'socket'};

	if (wantarray) {
		return keys %{$self->{'auth_protocols'}};
	} else {
		my $key = each %{$self->{'auth_protocols'}};
		return $key;
	}
}

sub printAuthProtocols {
	print "$_[1] supports: ",join(", ", Mail::Sender->QueryAuthProtocols($_[1] || 'localhost')),"\n";
}

sub TestServer {
	my $self = shift;
	local $_;
	if (!defined $self) {
		croak "Mail::Sender::TestServer() called without any parameter!";
	} elsif (ref $self) { # $sender->TestServer() or $sender->TestServer('the.server.com)
		if ($self->{'socket'}) { # the user did not Close() or Cancel() the previous mail
			die "You forgot to close the mail before calling TestServer!\n"
		}
		if (@_) {
			$self->{'smtp'} = shift();
			$self->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
			$self->{'smtp'} =~ s/\s+$//g;
			$self->{'smtpaddr'} = inet_aton($self->{'smtp'});
			if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }
			$self->{'smtpaddr'} = $1 if ($self->{'smtpaddr'} =~ /(.*)/s); # Untaint
		}
		$self->{'on_errors'} = 'die';
	} elsif ($self =~ /::/) { # Mail::Sender->TestServer('the.server.com')
		croak "Mail::Sender->TestServer() called without any parameter!"
			if ! @_;
		$self = new Mail::Sender {smtp => $_[0], on_errors => 'die'};
		return unless ref $self;
	} else { # Mail::Sender::QueryAuthProtocols('the.server.com')
		$self = new Mail::Sender {smtp => $self, on_errors => 'die'};
		return unless ref $self;
	}

	return $self->Error(NOSERVER) unless defined $self->{'smtp'};
#	if (!defined($self->{'smtpaddr'})) { return $self->Error(HOSTNOTFOUND($self->{'smtp'})); }

	if (exists $self->{'on_errors'} and (!defined($self->{'on_errors'}) or $self->{'on_errors'} eq 'undef')) {
		return $self->Connect() and $self->Close() and 1;
	} elsif (exists $self->{'on_errors'} and $self->{'on_errors'} eq 'die') {
		$self->Connect();
		$self->Close();
		return 1;
	} else {
		my $res = $self->Connect();
		return $res unless ref $res;
		$res = $self->Close();
		return $res unless ref $res;
		return $self;
	}
}

#====== Debuging bazmecks

$debug_code = <<'*END*';
package Mail::Sender::DBIO;
use IO::Handle;
use Tie::Handle;
@Mail::Sender::DBIO::ISA = qw(Tie::Handle);

sub SOCKET () {0}
sub LOG () {1}
sub ENDLINE () {2}
sub CLOSELOG () {3}
sub OFF () {4}

sub TIEHANDLE {
	my ($pkg,$socket,$debughandle, $mayCloseLog) = @_;
	return bless [$socket,$debughandle,1, $mayCloseLog,0], $pkg;
}

sub PRINT {
	my $self = shift;
	my $text = join(($\ || ''), @_);
	$self->[SOCKET]->print($text);
	return if $self->[OFF];
	$text =~ s/\x0D\x0A(?=.)/\x0D\x0A<< /g;
	$text = "<< ".$text if $self->[ENDLINE];
	$self->[ENDLINE] = ($text =~ /\x0D\x0A$/);
	$self->[LOG]->print($text);
}

sub READLINE {
	my $self = shift();
	my $socket = $self->[SOCKET];
	my $line = <$socket>;
	$self->[LOG]->print(">> $line") if defined $line and !$self->[OFF];
	return $line;
}

sub CLOSE {
	my $self = shift();
	$self->[SOCKET]->close();
	$self->[LOG]->close() if $self->[CLOSELOG];
	return $self->[SOCKET];
}

sub opened {
	our $SOCKET;
	local *SOCKET = $_[SOCKET];
	$SOCKET->opened();
}

use Data::Dumper;
sub stop_logging {
	my $self = tied(${$_[0]});

#print "stop_logging( ".$self." )\n";

	return if $self->[OFF];
	$self->[OFF] = 1;

	my $text = join(($\ || ''), $_[1])
		or return;
	$text .= "\x0D\x0A";
	$text =~ s/\x0D\x0A(?=.)/\x0D\x0A<< /g;
	$text = "<< ".$text if $self->[ENDLINE];
	$self->[ENDLINE] = ($text =~ /\x0D\x0A$/);
	$self->[LOG]->print($text);
}

sub start_logging {
	my $self = tied(${$_[0]});
	$self->[OFF] = 0;
}

*END*

my $pseudo_handle_code = <<'*END*';
package Mail::Sender::IO;
use IO::Handle;
use Tie::Handle;
@Mail::Sender::IO::ISA = qw(Tie::Handle);

sub TIEHANDLE {
	my ($pkg,$sender) = @_;
	return bless [$sender, $sender->{'_part'}], $pkg;
}

sub PRINT {
	my $self = shift;
	$self->[0]->SendEnc(@_);
}

sub PRINTF {
	my $self = shift;
	my $format = shift;
	$self->[0]->SendEnc( sprintf $format, @_);
}

sub CLOSE {
	my $self = shift();
	if ($self->[1]) {
		$self->[1]->EndPart();
	} else {
		$self->[0]->Close();
	}
}
*END*

=head2 GetHandle

Returns a "filehandle" to which you can print the message or file to attach or whatever.
The data you print to this handle will be encoded as necessary. Closing this handle closes
either the message (for single part messages) or the part.

	$sender->Open({...});
	my $handle = $sender->GetHandle();
	print $handle "Hello world.\n"
	my ($mday,$mon,$year) = (localtime())[3,4,5];
	printf $handle "Today is %04d/%02d/%02d.", $year+1900, $mon+1, $mday;
	close $handle;

P.S.: There is a big difference between the handle stored in $sender->{'socket'} and the handle
returned by this function ! If you print something to $sender->{'socket'} it will be sent to the server
without any modifications, encoding, escaping, ...
You should NOT touch the $sender->{'socket'} unless you really really know what you are doing.

=cut

package Mail::Sender;
sub GetHandle {
	my $self = shift();
	unless (defined @Mail::Sender::IO::ISA) {
		eval "use Symbol;";
		eval $pseudo_handle_code;
	}
	my $handle = gensym();
	tie *$handle, 'Mail::Sender::IO', $self;
	return $handle;
}

=head1 FUNCTIONS

=head2 GuessCType

	$ctype = GuessCType $filename, $filepath;

Guesses the content type based on the filename or the file contents.
This function is used when you attach a file and do not specify the content type.
It is not exported by default!

The builtin version uses the filename extension to guess the type.
Currently there are only a few extensions defined, you may add other extensions this way:

	$Mail::Sender::CTypes{'EXT'} = 'content/type';
	...

The extension has to be in UPPERCASE and will be matched case sensitively.

The package now includes two addins improving the guesswork. If you "use" one of them in your script,
it replaces the builtin GuessCType() subroutine with a better one:

	Mail::Sender::CType::Win32
		Win32 only, the content type is read from the registry
	Mail::Sender::CType::Ext
		any OS, a longer list of extensions from A. Guillaume

=head2 ResetGMTdiff

	ResetGMTdiff()

The module computes the local vs. GMT time difference to include in the timestamps
added into the message headers. As the time difference may change due to summer
savings time changes you may want to reset the time difference ocassionaly
in long running programs.

=head1 CONFIG

If you create a file named Sender.config in the same directory where
Sender.pm resides, this file will be "require"d as soon as you "use
Mail::Sender" in your script. Of course the Sender.config MUST "return a
true value", that is it has to be succesfully compiled and the last
statement must return a true value. You may use this to forbide the use
of Mail::Sender to some users.

You may define the default settings for new Mail::Sender objects and do
a few more things.

The default options are stored in hash %Mail::Sender::default. You may
use all the options you'd use in C<new>, C<Open>, C<OpenMultipart>,
C<MailMsg> or C<MailFile>.

 Eg.
  %default = (
    smtp => 'mail.yourhost.cz',
    from => getlogin.'yourhost.cz',
    client => getlogin.'.yourhost.cz'
  );
  # of course you will use your own mail server here !

The other options you may set here (or later of course) are
$Mail::Sender::SITE_HEADERS, $Mail::Sender::NO_X_MAILER and
$Mail::Sender::NO_DATE. (These are plain old scalar variables, there is no
function or method for modifying them. Just set them to anything you need.)

The $Mail::Sender::SITE_HEADERS may contain headers that will be added
to each mail message sent by this script, the $Mail::Sender::NO_X_MAILER
disables the header item specifying that the message was sent by
Mail::Sender and $Mail::Sender::NO_DATE turns off the Date: header generation.

!!! $Mail::Sender::SITE_HEADERS may NEVER end with \r\n !!!

If you want to set the $Mail::Sender::SITE_HEADERS for every script sent
from your server without your users being able to change it you may use
this hack:

 $loginname = something_that_identifies_the_user();
 *Mail::Sender::SITE_HEADERS = \"X-Sender: $loginname via $0";
 $Mail::Sender::NO_X_MAILER = 1;

You may even "install" your custom function that will be evaluated for
each message just before contacting the server. You may change all the
options from within as well as stop sending the message.

All you have to do is to create a function named SiteHook in
Mail::Sender package. This function will get the Mail::Sender object as
its first argument. If it returns a TRUE value the message is sent,
if it returns FALSE the sending is canceled and the user gets
"Site specific error" error message.

If you want to give some better error message you may do it like this :

 sub SiteHook {
  my $self = shift;
  if (whatever($self)) {
    $self->Error( SITEERROR);
    $Mail::Sender::Error = "I don't like this mail";
    return 0
  } else {
    return 1;
  }
 }


This example will ensure the from address is the users real address :

 sub SiteHook {
  my $self = shift;
  $self->{'fromaddr'} = getlogin.'@yoursite.com';
  $self->{'from'} = getlogin.'@yoursite.com';
  1;
 }

Please note that at this stage the from address is in two different
object properties.

$self->{'from'} is the address as it will appear in the mail, that is
it may include the full name of the user or any other comment
( "Jan Krynicky <jenda@krynicky.cz>" for example), while the
$self->{'fromaddr'} is realy just the email address per se and it will
be used in conversation with the SMTP server. It must be without
comments ("jenda@krynicky.cz" for example)!


Without write access to .../lib/Mail/Sender.pm or
.../lib/Mail/Sender.config your users will then be unable to get rid of
this header. Well ... everything is doable, if they are cheeky enough ... :-(

So if you take care of some site with virtual servers for several
clients and implement some policy via SiteHook() or
$Mail::Sender::SITE_HEADERS search the clients' scripts for "SiteHook"
and "SITE_HEADERS" from time to time. To see who's cheating.

=head1 AUTHENTICATION

If you get a "Local user "xxx@yyy.com" unknown on host "zzz"" message it usualy means that
your mail server is set up to forbid mail relay. That is it only accepts messages to or from a local user.
If you need to be able to send a message with both the sender's and recipient's address remote, you
need to somehow authenticate to the server. You may need the help of the mail server's administrator
to find out what username and password and/or what authentication protocol are you supposed to use.

There are many authentication protocols defined for ESTMP, Mail::Sender natively supports
only PLAIN, LOGIN, CRAM-MD5 and NTLM (please see the docs for C<new Mail::Sender>).

If you want to know what protocols are supported by your server you may get the list by this:

	/tmp# perl -MMail::Sender -e 'Mail::Sender->printAuthProtocols("the.server.com")'
  or
	c:\> perl -MMail::Sender -e "Mail::Sender->printAuthProtocols('the.server.com')"


There is one more way to authenticate. Some servers want you to login by POP3 before you
can send a message. You have to use Net::POP3 or Mail::POP3Client to do this.

=head2 Other protocols

It is possible to add new authentication protocols to Mail::Sender. All you have to do is
to define a function Mail::Sender::Auth::PROTOCOL_NAME that will implement
the login. The function gets one parameter ... the Mail::Sender object.
It can access these properties:

	$obj->{'socket'} : the socket to print to and read from
		you may use the send_cmd() function to send a request
		and read a response from the server
	$obj->{'authid'} : the username specified in the new Mail::Sender,
		Open or OpenMultipart call
	$obj->{'authpwd'} : the password
	$obj->{auth...} : all unknown parameters passed to the constructor or the mail
		opening/creation methods are preserved in the object. If the protocol requires
		any other options, please use names starting with "auth". Eg. "authdomain", ...
	$obj->{'error'} : this should be set to a negative error number. Please use numbers
		below -1000 for custom errors.
	$obj->{'error_msg'} : this should be set to the error message

	If the login fails you should
		1) Set $Mail::Sender::Error to the error message
		2) Set $obj->{'error_msg'} to the error message
		2) Set $obj->{'error'} to a negative number
		3) return a negative number
	If it succeeds, please return "nothing" :
		return;

Please use the protocols defined within Sender.pm as examples.

=head1 EXAMPLES

=head2 Object creation

 ref ($sender = new Mail::Sender { from => 'somebody@somewhere.com',
       smtp => 'mail.yourISP.com', boundary => 'This-is-a-mail-boundary-435427'})
 or die "Error in mailing : $Mail::Sender::Error\n";

or

 my $sender = new Mail::Sender { ... };
 die "Error in mailing : $Mail::Sender::Error\n" unless ref $sender;

or

 my $sender = new Mail::Sender { ..., on_errors => 'undef' }
   or die "Error in mailing : $Mail::Sender::Error\n";

You may specify the options either when creating the Mail::Sender object
or later when you open a message. You may also set the default options when
installing the module (See C<CONFIG> section). This way the admin may set
the SMTP server and even the authentication options and the users do not have
to specify it again.

You should keep in mind that the way Mail::Sender reports failures depends on the 'on_errors'=>
option. If you set it to 'die' it throws an exception, if you set it to C<undef> or C<'undef'> it returns
undef and otherwise it returns a negative error code!

=head2 Simple single part message

	$sender = new Mail::Sender {
		smtp => 'mail.yourISP.com',
		from => 'somebody@somewhere.com',
		on_errors => undef,
	}
		or die "Can't create the Mail::Sender object: $Mail::Sender::Error\n";
	$sender->Open({
		to => 'mama@home.org, papa@work.com',
		cc => 'somebody@somewhere.com',
		subject => 'Sorry, I\'ll come later.'
	})
		or die "Can't open the message: $sender->{'error_msg'}\n";
	$sender->SendLineEnc("I'm sorry, but thanks to the lusers,
		I'll come at 10pm at best.");
	$sender->SendLineEnc("\nHi, Jenda");
	$sender->Close()
		or die "Failed to send the message: $sender->{'error_msg'}\n";

or

	eval {
		$sender = new Mail::Sender {
			smtp => 'mail.yourISP.com',
			from => 'somebody@somewhere.com',
			on_errors => 'die',
		};
		$sender->Open({
			to => 'mama@home.org, papa@work.com',
			cc => 'somebody@somewhere.com',
			subject => 'Sorry, I\'ll come later.'
		});
		$sender->SendLineEnc("I'm sorry, but thanks to the lusers,
			I'll come at 10pm at best.");
		$sender->SendLineEnc("\nHi, Jenda");
		$sender->Close();
	};
	if ($@) {
		die "Failed to send the message: $@\n";
	}

or

	$sender = new Mail::Sender {
		smtp => 'mail.yourISP.com',
		from => 'somebody@somewhere.com',
		on_errors => 'code',
	};
	die "Can't create the Mail::Sender object: $Mail::Sender::Error\n"
		unless ref $sender;
	ref $sender->Open({
		to => 'mama@home.org, papa@work.com',
		cc => 'somebody@somewhere.com',
		subject => 'Sorry, I\'ll come later.'
	})
		or die "Can't open the message: $sender->{'error_msg'}\n";
	$sender->SendLineEnc("I'm sorry, but thanks to the lusers,
		I'll come at 10pm at best.");
	$sender->SendLineEnc("\nHi, Jenda");
	ref $sender->Close
		or die "Failed to send the message: $sender->{'error_msg'}\n";

=head2 Using GetHandle()

  ref $sender->Open({to => 'friend@other.com', subject => 'Hello dear friend'})
	 or die "Error: $Mail::Sender::Error\n";
  my $FH = $sender->GetHandle();
  print $FH "How are you?\n\n";
  print $FH <<'*END*';
  I've found these jokes.

   Doctor, I feel like a pack of cards.
   Sit down and I'll deal with you later.

   Doctor, I keep thinking I'm a dustbin.
   Don't talk rubbish.

  Hope you like'em. Jenda
  *END*

  $sender->Close;
  # or
  # close $FH;

or

  eval {
    $sender->Open({ on_errors => 'die',
			 to => 'mama@home.org, papa@work.com',
                cc => 'somebody@somewhere.com',
                subject => 'Sorry, I\'ll come later.'});
    $sender->SendLineEnc("I'm sorry, but due to a big load of work,
  I'll come at 10pm at best.");
    $sender->SendLineEnc("\nHi, Jenda");
    $sender->Close;
  };
  if ($@) {
    print "Error sending the email: $@\n";
  } else {
    print "The mail was sent.\n";
  }

=head2 Multipart message with attachment

 $sender->OpenMultipart({to => 'Perl-Win32-Users@activeware.foo',
                         subject => 'Mail::Sender.pm - new module'});
 $sender->Body;
 $sender->SendEnc(<<'*END*');
 Here is a new module Mail::Sender.
 It provides an object based interface to sending SMTP mails.
 It uses a direct socket connection, so it doesn't need any
 additional program.

 Enjoy, Jenda
 *END*
 $sender->Attach(
  {description => 'Perl module Mail::Sender.pm',
   ctype => 'application/x-zip-encoded',
   encoding => 'Base64',
   disposition => 'attachment; filename="Sender.zip"; type="ZIP archive"',
   file => 'sender.zip'
  });
 $sender->Close;

or

 $sender->OpenMultipart({to => 'Perl-Win32-Users@activeware.foo',
                         subject => 'Mail::Sender.pm - new version'});
 $sender->Body({ msg => <<'*END*' });
 Here is a new module Mail::Sender.
 It provides an object based interface to sending SMTP mails.
 It uses a direct socket connection, so it doesn't need any
 additional program.

 Enjoy, Jenda
 *END*
 $sender->Attach(
  {description => 'Perl module Mail::Sender.pm',
   ctype => 'application/x-zip-encoded',
   encoding => 'Base64',
   disposition => 'attachment; filename="Sender.zip"; type="ZIP archive"',
   file => 'sender.zip'
  });
 $sender->Close;

or (in case you have the file contents in a scalar)

 $sender->OpenMultipart({to => 'Perl-Win32-Users@activeware.foo',
                         subject => 'Mail::Sender.pm - new version'});
 $sender->Body({ msg => <<'*END*' });
 Here is a new module Mail::Sender.
 It provides an object based interface to sending SMTP mails.
 It uses a direct socket connection, so it doesn't need any
 additional program.

 Enjoy, Jenda
 *END*
 $sender->Part(
  {description => 'Perl module Mail::Sender.pm',
   ctype => 'application/x-zip-encoded',
   encoding => 'Base64',
   disposition => 'attachment; filename="Sender.zip"; type="ZIP archive"',
   msg => $sender_zip_contents,
  });
 $sender->Close;


=head2 Using exceptions (no need to test return values after each function)

 use Mail::Sender;
 eval {
 (new Mail::Sender {on_errors => 'die'})
 	->OpenMultipart({smtp=> 'jenda.krynicky.cz', to => 'jenda@krynicky.cz',subject => 'Mail::Sender.pm - new version'})
 	->Body({ msg => <<'*END*' })
 Here is a new module Mail::Sender.
 It provides an object based interface to sending SMTP mails.
 It uses a direct socket connection, so it doesn't need any
 additional program.

 Enjoy, Jenda
 *END*
 	->Attach({
 		description => 'Perl module Mail::Sender.pm',
 		ctype => 'application/x-zip-encoded',
 		encoding => 'Base64',
 		disposition => 'attachment; filename="Sender.zip"; type="ZIP archive"',
 		file => 'W:\jenda\packages\Mail\Sender\Mail-Sender-0.7.14.3.tar.gz'
 	})
 	->Close();
 } or print "Error sending mail: $@\n";

=head2 Using MailMsg() shortcut to send simple messages

If everything you need is to send a simple message you may use:

 if (ref ($sender->MailMsg({to =>'Jenda@Krynicky.czX', subject => 'this is a test',
                         msg => "Hi Johnie.\nHow are you?"}))) {
  print "Mail sent OK."
 } else {
  die "$Mail::Sender::Error\n";
 }

or

 if ($sender->MailMsg({
   smtp => 'mail.yourISP.com',
   from => 'somebody@somewhere.com',
   to =>'Jenda@Krynicky.czX',
   subject => 'this is a test',
   msg => "Hi Johnie.\nHow are you?"
 }) < 0) {
  die "$Mail::Sender::Error\n";
 }
 print "Mail sent OK."

=head2 Using MailMsg and authentication

 if ($sender->MailMsg({
   smtp => 'mail.yourISP.com',
   from => 'somebody@somewhere.com',
   to =>'Jenda@Krynicky.czX',
   subject => 'this is a test',
   msg => "Hi Johnie.\nHow are you?"
   auth => 'NTLM',
   authid => 'jenda',
   authpwd => 'benda',
 }) < 0) {
  die "$Mail::Sender::Error\n";
 }
 print "Mail sent OK."

=head2 Using MailFile() shortcut to send an attachment

If you want to attach some files:

 (ref ($sender->MailFile(
  {to =>'you@address.com', subject => 'this is a test',
   msg => "Hi Johnie.\nI'm sending you the pictures you wanted.",
   file => 'image1.jpg,image2.jpg'
  }))
  and print "Mail sent OK."
 )
 or die "$Mail::Sender::Error\n";

=head2 Sending HTML messages

If you are sure the HTML doesn't contain any accentuated characters (with codes above 127).

 open IN, $htmlfile or die "Cannot open $htmlfile : $!\n";
 $sender->Open({ from => 'your@address.com', to => 'other@address.com',
        subject => 'HTML test',
        ctype => "text/html",
        encoding => "7bit"
 }) or die $Mail::Sender::Error,"\n";

 while (<IN>) { $sender->SendEx($_) };
 close IN;
 $sender->Close();

Otherwise use SendEnc() instead of SendEx() and "quoted-printable" instead of "7bit".

Another ... quicker way ... would be:

 open IN, $htmlfile or die "Cannot open $htmlfile : $!\n";
 $sender->Open({ from => 'your@address.com', to => 'other@address.com',
        subject => 'HTML test',
        ctype => "text/html",
        encoding => "quoted-printable"
 }) or die $Mail::Sender::Error,"\n";

 while (read IN, $buff, 4096) { $sender->SendEnc($buff) };
 close IN;
 $sender->Close();

=head2 Sending HTML messages with inline images

	if (ref $sender->OpenMultipart({
		from => 'someone@somewhere.net', to => $recipients,
		subject => 'Embedded Image Test',
		boundary => 'boundary-test-1',
		multipart => 'related'})) {
		$sender->Attach(
			 {description => 'html body',
			 ctype => 'text/html; charset=us-ascii',
			 encoding => '7bit',
			 disposition => 'NONE',
			 file => 'test.html'
		});
		$sender->Attach({
			description => 'ed\'s gif',
			ctype => 'image/gif',
			encoding => 'base64',
			disposition => "inline; filename=\"apache_pb.gif\";\r\nContent-ID: <img1>",
			file => 'apache_pb.gif'
		});
		$sender->Close() or die "Close failed! $Mail::Sender::Error\n";
	} else {
		die "Cannot send mail: $Mail::Sender::Error\n";
	}

And in the HTML you'll have this :
 ... <IMG src="cid:img1"> ...
on the place where you want the inlined image.

Please keep in mind that the image name is unimportant, it's the Content-ID what counts!

# or using the eval{ $obj->Method()->Method()->...->Close()} trick ...

	use Mail::Sender;
	eval {
	(new Mail::Sender)
		->OpenMultipart({
			to => 'someone@somewhere.com',
			subject => 'Embedded Image Test',
			boundary => 'boundary-test-1',
			type => 'multipart/related'
		})
		->Attach({
			description => 'html body',
			ctype => 'text/html; charset=us-ascii',
			encoding => '7bit',
			disposition => 'NONE',
			file => 'c:\temp\zk\HTMLTest.htm'
		})
		->Attach({
			description => 'Test gif',
			ctype => 'image/gif',
			encoding => 'base64',
			disposition => "inline; filename=\"test.gif\";\r\nContent-ID: <img1>",
			file => 'test.gif'
		})
		->Close()
	}
	or die "Cannot send mail: $Mail::Sender::Error\n";

=head2 Sending message with plaintext and HTML alternatives

	use Mail::Sender;

	eval {
		(new Mail::Sender)
		->OpenMultipart({
			to => 'someone@somewhere.com',
			subject => 'Alternatives',
	#		debug => 'c:\temp\zkMailFlow.log',
			multipart => 'mixed',
		})
			->Part({ctype => 'multipart/alternative'})
				->Part({ ctype => 'text/plain', disposition => 'NONE', msg => <<'*END*' })
	A long
	mail
	message.
	*END*
				->Part({ctype => 'text/html', disposition => 'NONE', msg => <<'*END*'})
	<html><body><h1>A long</h1><p align=center>
	mail
	message.
	</p></body></html>
	*END*
			->EndPart("multipart/alternative")
		->Close();
	} or print "Error sending mail: $Mail::Sender::Error\n";

=head2 Sending message with plaintext and HTML alternatives with inline images

	use Mail::Sender;

	eval {
		(new Mail::Sender)
		->OpenMultipart({
			to => 'someone@somewhere.com',
			subject => 'Alternatives with images',
	#		debug => 'c:\temp\zkMailFlow.log',
			multipart => 'related',
		})
			->Part({ctype => 'multipart/alternative'})
				->Part({ ctype => 'text/plain', disposition => 'NONE', msg => <<'*END*' })
	A long
	mail
	message.
	*END*
				->Part({ctype => 'text/html', disposition => 'NONE', msg => <<'*END*'})
	<html><body><h1>A long</h1><p align=center>
	mail
	message.
	<img src="cid:img1">
	</p></body></html>
	*END*
			->EndPart("multipart/alternative")
			->Attach({
				description => 'ed\'s jpg',
				ctype => 'image/jpeg',
				encoding => 'base64',
				disposition => "inline; filename=\"0518m_b.jpg\";\r\nContent-ID: <img1>",
				file => 'E:\pix\humor\0518m_b.jpg'
			})
		->Close();
	} or print "Error sending mail: $Mail::Sender::Error\n";

Keep in mind please that different mail clients display messages differently. You may
need to try several ways to create messages so that they appear the way you need.
These two examples looked like I expected in Pegasus Email and MS Outlook.

If this doesn't work with your mail client, please let me know and we might find a way.


=head2 Sending a file that was just uploaded from an HTML form

 use CGI;
 use Mail::Sender;

 $query = new CGI;

 # uploading the file...
 $filename = $query->param('mailformFile');
 if ($filename ne ""){
  $tmp_file = $query->tmpFileName($filename);
 }

 $sender = new Mail::Sender {from => 'script@krynicky.cz',smtp => 'mail.krynicky.czX'};
 $sender->OpenMultipart({to=> 'jenda@krynicky.czX',subject=> 'test CGI attach'});
 $sender->Body();
 $sender->Send(<<"*END*");
 This is just a test of mail with an uploaded file.

 Jenda
 *END*
 $sender->Attach({
    encoding => 'Base64',
    description => $filename,
    ctype => $query->uploadInfo($filename)->{'Content-Type'},
    disposition => "attachment; filename = $filename",
    file => $tmp_file
 });
 $sender->Close();

 print "Content-type: text/plain\n\nYes, it's sent\n\n";

=head2 Listing the authentication protocols supported by the server

 use Mail::Sender;
 my $sender = new Mail::Sender {smtp => 'localhost'};
 die "Error: $Mail::Sender::Error\n" unless ref $sender;
 print join(', ', $sender->QueryAuthProtocols()),"\n";

or (if you have Mail::Sender 0.8.05 or newer)

 use Mail::Sender;
 print join(', ', Mail::Sender->QueryAuthProtocols('localhost')),"\n";

or

 use Mail::Sender;
 print join(', ', Mail::Sender::QueryAuthProtocols('localhost')),"\n";

=head2 FAQ

=head3 Forwarding the messages created by Mail::Sender removes accents. Why?

The most likely colprit is missing or incorrect charset specified for the body or
a part of the email. You should add something like

	charset => 'iso-8859-1',
	encoding => 'quoted-printable',

to the parameters passed to Open(), OpenMultipart(), MailMsg(), Body() or Part() or

	b_charset => 'iso-8859-1',
	b_encoding => 'quoted-printable',

to the parameters for MailFile().

If you use a different charset ('iso-8859-2', 'win-1250', ...) you will of course need
to specify that charset. If you are not sure, try to send a mail with some other mail client
and then look at the message/part headers.

=head2 Sometimes there is an equals sign at the end of an attached file when
I open the email in Outlook. What's wrong?

Outlook is. It has (had) a bug in its quoted printable decoding routines.
This problem happens only in quoted-printable encoded parts on multipart messages.
And only if the data in that part do not end with a newline. (This is new in 0.8.08, in older versions
it happened in all QP encoded parts.)

The problem is that an equals sign at the end of a line in a quoted printable encoded text means
"ignore the newline". That is

	fooo sdfg sdfg sdfh dfh =
	dfsgdsfg

should be decoded as

	fooo sdfg sdfg sdfh dfh dfsgdsfg

The problem is at the very end of a file. The part boundary (text separating different
parts of a multipart message) has to start on a new line, if the attached file ends by a newline everything is cool.
If it doesn't I need to add a newline and to denote that the newline is not part of the original file I add an equals:

	dfgd dsfgh dfh dfh dfhdfhdfhdfgh
	this is the last line.=
	--message-boundary-146464--

Otherwise I'd add a newline at the end of the file.
If you do not care about the newline and want to be sure Outlook doesn't add the equals to the file add

	bypass_outlook_bug => 1

parameter to C<new Mail::Sender> or C<Open>/C<OpenMultipart>.

=head2 WARNING

DO NOT mix Open(Multipart)|Send(Line)(Ex)|Close with MailMsg or MailFile.
Both Mail(Msg/File) close any Open-ed mail.
Do not try this:

 $sender = new Mail::Sender ...;
 $sender->OpenMultipart...;
 $sender->Body;
 $sender->Send("...");
 $sender->MailFile({file => 'something.ext');
 $sender->Close;

This WON'T work!!!

=head2 GOTCHAS

=head3 Local user "someone@somewhere.com" doesn't exist

"Thanks" to spammers mail servers usualy do not allow just anyone to post a message through them.
Most often they require that either the sender or the recipient is local to the server

=head3 Mail::Sendmail works, Mail::Sender doesn't

If you are able to connect to the mail server and scripts using Mail::Sendmail work, but Mail::Sender fails with
"connect() failed", please review the settings in /etc/services. The port for SMTP should be 25.

=head3 $/ and $\

If you change the $/ ($RS, $INPUT_RECORD_SEPARATOR) or $\ ($ORS, $OUTPUT_RECORD_SEPARATOR)
or $, ($OFS, $OUTPUT_FIELD_SEPARATOR) Mail::Sender may stop working! Keep in mind that those variables are global
and therefore they change the behaviour of <> and print everywhere.
And since the SMTP is a plain text protocol if you change the notion of lines you can break it.

If you have to fiddle with $/, $\ or $, do it in the smallest possible block of code and local()ize the change!

	open my $IN, '<', $filename or die "Can't open $filename: $!\n";
	my $data = do {local $/; <$IN>};
	close $IN;

=head1 BUGS

I'm sure there are many. Please let me know if you find any.

The problem with multiline responses from some SMTP servers (namely qmail) is solved. At last.

=head1 SEE ALSO

MIME::Lite, MIME::Entity, Mail::Sendmail, Mail::Mailer, ...

There are lots of mail related modules on CPAN, with different capabilities and interfaces. You
have to find the right one yourself :-)

=head1 DISCLAIMER

This module is based on SendMail.pm Version : 1.21 that appeared in
Perl-Win32-Users@activeware.com mailing list. I don't remember the name
of the poster and it's not mentioned in the script. Thank you mr. C<undef>.

=head1 AUTHOR

Jan Krynicky <Jenda@Krynicky.cz>
http://Jenda.Krynicky.cz

With help of Rodrigo Siqueira <rodrigo@insite.com.br>,
Ed McGuigan <itstech1@gate.net>,
John Sanche <john@quadrant.net>,
Brian Blakley <bblakley@mp5.net>,
and others.

=head1 COPYRIGHT

Copyright (c) 1997-2006 Jan Krynicky <Jenda@Krynicky.cz>. All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. There is only one aditional condition, you may
NOT use this module for SPAMing! NEVER! (see http://spam.abuse.net/ for definition)

=cut
