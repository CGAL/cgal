#!/usr/bin/env perl

# ====================================================================
# commit-email.pl: send a commit email for commit REVISION in
# repository REPOS to some email addresses.
#
# For usage, see the usage subroutine or run the script with no
# command line arguments.
#
# ====================================================================
# Copyright (c) 2000-2004 CollabNet.  All rights reserved.
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
#
# Changes for CGAL project:
# - Email content modified to be closer to CVS commit email format
# - Email sent in both HTML + text

use lib "/home/groups/cgal/hooks/";
use Mail::Sender;
require 5.004; # This is when locale support was added.

# Use UTF-8 to display properly diacritic characters
$ENV{'LANG'} = 'en_US.UTF-8';


# Keep beginning of $data_ref up to $max_size size
sub cut_array
{
    my $data_ref = shift @_;
    my $max_size = shift @_;

    my @new_data = ();
    my $line;
    my $total_size = 0;

    for $line (@{$data_ref}) {
        $total_size += length ($line);
        if ($total_size > $max_size) {
            push @new_data, "\n(snip)\n";
            return @new_data;
        }
        push @new_data, $line;
    }
    return @new_data;
}

# Encode raw text as HTML text
sub html_encode
{
    my $data_ref = shift @_;

    my @new_data = ();
    my $line;

    for $line (@{$data_ref}) {
        # Encode special characters ; < > SPACE
        $line =~ s/&/&amp;/g;
        $line =~ s/</&lt;/g;
        $line =~ s/>/&gt;/g;
        $line =~ s/ /&nbsp;/g;

        # Encode new line
        $line =~ s/\n/<br>\n/;

        push @new_data, $line;
    }

    return @new_data;
}

# Encode raw text for HTML <pre> tag
sub encode_for_pre
{
    my $data_ref = shift @_;

    my @new_data = ();
    my $line;

    for $line (@{$data_ref}) {
        # Encode special characters < >
        $line =~ s/</&lt;/g;
        $line =~ s/>/&gt;/g;

        push @new_data, $line;
    }

    return @new_data;
}


# Turn on warnings the best way depending on the Perl version.
BEGIN {
  if ( $] >= 5.006_000)
    { require warnings; import warnings; }
  else
    { $^W = 1; }
}

use strict;
use Carp;

######################################################################
# Configuration section.

# Svnlook path.
my $svnlook = "/usr/bin/svnlook";

# By default, when a file is deleted from the repository, svnlook diff
# prints the entire contents of the file.  If you want to save space
# in the log and email messages by not printing the file, then set
# $no_diff_deleted to 1.
my $no_diff_deleted = 1;

# ViewVC URL (specific to InriaGForge)
my $viewvc_url = "https://gforge.inria.fr/scm/viewvc.php";

# Since the path to svnlook depends upon the local installation
# preferences, check that the required programs exist to insure that
# the administrator has set up the script properly.
{
  my $ok = 1;
  foreach my $program ($svnlook)
    {
      if (-e $program)
        {
          unless (-x $program)
            {
              warn "$0: required program `$program' is not executable, ",
                   "edit $0.\n";
              $ok = 0;
            }
        }
      else
        {
          warn "$0: required program `$program' does not exist, edit $0.\n";
          $ok = 0;
        }
    }
  exit 1 unless $ok;
}


######################################################################
# Initial setup/command-line handling.

# Each value in this array holds a hash reference which contains the
# associated email information for one project.  Start with an
# implicit rule that matches all paths.
my @project_settings_list = (&new_project);

# Process the command line arguments till there are none left.  The
# first two arguments that are not used by a command line option are
# the repository path and the revision number.
my $repos;
my $rev;

# Use the reference to the first project to populate.
my $current_project = $project_settings_list[0];

# This hash matches the command line option to the hash key in the
# project.  If a key exists but has a false value (''), then the
# command line option is allowed but requires special handling.
my %opt_to_hash_key = ('--from' => 'from_address',
                       '-h'     => 'hostname',
                       '-l'     => 'log_file',
                       '-m'     => '',
                       '-r'     => 'reply_to',
                       '-s'     => 'subject_prefix');

while (@ARGV)
  {
    my $arg = shift @ARGV;
    if ($arg =~ /^-/)
      {
        my $hash_key = $opt_to_hash_key{$arg};
        unless (defined $hash_key)
          {
            die "$0: command line option `$arg' is not recognized.\n";
          }

        unless (@ARGV)
          {
            die "$0: command line option `$arg' is missing a value.\n";
          }
        my $value = shift @ARGV;

        if ($hash_key)
          {
            $current_project->{$hash_key} = $value;
          }
        else
          {
            # Here handle -m.
            unless ($arg eq '-m')
              {
                die "$0: internal error: should only handle -m here.\n";
              }
            $current_project                = &new_project;
            $current_project->{match_regex} = $value;
            push(@project_settings_list, $current_project);
          }
      }
    elsif ($arg =~ /^-/)
      {
        die "$0: command line option `$arg' is not recognized.\n";
      }
    else
      {
        if (! defined $repos)
          {
            $repos = $arg;
          }
        elsif (! defined $rev)
          {
            $rev = $arg;
          }
        else
          {
            push(@{$current_project->{email_addresses}}, $arg);
          }
      }
  }

# If the revision number is undefined, then there were not enough
# command line arguments.
&usage("$0: too few arguments.") unless defined $rev;

# Check the validity of the command line arguments.  Check that the
# revision is an integer greater than 0 and that the repository
# directory exists.
unless ($rev =~ /^\d+/ and $rev > 0)
  {
    &usage("$0: revision number `$rev' must be an integer > 0.");
  }
unless (-e $repos)
  {
    &usage("$0: repos directory `$repos' does not exist.");
  }
unless (-d _)
  {
    &usage("$0: repos directory `$repos' is not a directory.");
  }

# Check that all of the regular expressions can be compiled and
# compile them.
{
  my $ok = 1;
  for (my $i=0; $i<@project_settings_list; ++$i)
    {
      my $match_regex = $project_settings_list[$i]->{match_regex};

      # To help users that automatically write regular expressions
      # that match the root directory using ^/, remove the / character
      # because subversion paths, while they start at the root level,
      # do not begin with a /.
      $match_regex =~ s#^\^/#^#;

      my $match_re;
      eval { $match_re = qr/$match_regex/ };
      if ($@)
        {
          warn "$0: -m regex #$i `$match_regex' does not compile:\n$@\n";
          $ok = 0;
          next;
        }
      $project_settings_list[$i]->{match_re} = $match_re;
    }
  exit 1 unless $ok;
}

######################################################################
# Harvest data using svnlook.

# Change into /tmp so that svnlook diff can create its .svnlook
# directory.
my $tmp_dir = '/tmp';
chdir($tmp_dir)
  or die "$0: cannot chdir `$tmp_dir': $!\n";

# Get the author, date, and log from svnlook.
my @svnlooklines = &read_from_process($svnlook, 'info', $repos, '-r', $rev);
my $author = shift @svnlooklines;
my $date = shift @svnlooklines;
shift @svnlooklines;
my @log = map { "$_\n" } @svnlooklines;

# Figure out what directories have changed using svnlook.
my @dirschanged = &read_from_process($svnlook, 'dirs-changed', $repos,
                                     '-r', $rev);

# Lose the trailing slash in the directory names if one exists, except
# in the case of '/'.
my $rootchanged = 0;
for (my $i=0; $i<@dirschanged; ++$i)
  {
    if ($dirschanged[$i] eq '/')
      {
        $rootchanged = 1;
      }
    else
      {
        $dirschanged[$i] =~ s#^(.+)[/\\]$#$1#;
      }
  }

# Figure out what files have changed using svnlook.
@svnlooklines = &read_from_process($svnlook, 'changed', $repos, '-r', $rev);

# Parse the changed nodes.
my @adds;
my @dels;
my @mods;
foreach my $line (@svnlooklines)
  {
    my $path = '';
    my $code = '';

    # Split the line up into the modification code and path, ignoring
    # property modifications.
    if ($line =~ /^(.).  (.*)$/)
      {
        $code = $1;
        $path = $2;
      }

    if ($code eq 'A')
      {
        push(@adds, $path);
      }
    elsif ($code eq 'D')
      {
        push(@dels, $path);
      }
    else
      {
        push(@mods, $path);
      }
  }

# Get the diff from svnlook.
my @no_diff_deleted = $no_diff_deleted ? ('--no-diff-deleted') : ();
my @difflines = &read_from_process($svnlook, 'diff', $repos,
                                   '-r', $rev, @no_diff_deleted);

######################################################################
# Modified directory name collapsing.

# Collapse the list of changed directories only if the root directory
# was not modified, because otherwise everything is under root and
# there's no point in collapsing the directories, and only if more
# than one directory was modified.
my $commondir = '';
if (!$rootchanged and @dirschanged > 1)
  {
    my $firstline    = shift @dirschanged;
    my @commonpieces = split('/', $firstline);
    foreach my $line (@dirschanged)
      {
        my @pieces = split('/', $line);
        my $i = 0;
        while ($i < @pieces and $i < @commonpieces)
          {
            if ($pieces[$i] ne $commonpieces[$i])
              {
                splice(@commonpieces, $i, @commonpieces - $i);
                last;
              }
            $i++;
          }
      }
    unshift(@dirschanged, $firstline);

    if (@commonpieces)
      {
        $commondir = join('/', @commonpieces);
        my @new_dirschanged;
        foreach my $dir (@dirschanged)
          {
            if ($dir eq $commondir)
              {
                $dir = '.';
              }
            else
              {
                $dir =~ s#^$commondir/##;
              }
            push(@new_dirschanged, $dir);
          }
        @dirschanged = @new_dirschanged;
      }
  }
my $dirlist = join(' ', @dirschanged);

######################################################################
# Assembly of log message.

# Get project name (specific to InriaGForge)
my $project_name = $repos;
$project_name =~ s/\/svn\///;
$project_name =~ s/\/svnroot\///;

# Put together the body of the log message AS TEXT
my @body;
{
    my $prev_rev = $rev - 1;

    # Write summary
    push(@body, "*Summary*\n");
    push(@body, "\n");
    push(@body, "Revision in ViewVC <$viewvc_url?revision=$rev&root=$project_name&view=rev>\n");
    push(@body, "\n");
    push(@body, "New Revision: $rev\n");
    push(@body, "Author: $author\n");
    push(@body, "Date: $date\n");
    push(@body, "\n");
    push(@body, "Log message:\n");
    push(@body, @log);
    push(@body, "\n");

    # Added files list
    if (@adds)
    {
        @adds = sort @adds;
        push(@body, "*Added files*\n");
        push(@body, "\n");
        push(@body, map { /\/$/ ? "$_\n" : "$_ <$viewvc_url/$_?root=$project_name&revision=$rev&view=markup>\n" } @adds);
        push(@body, "\n");
    }

    # Deleted files list
    if (@dels)
    {
        @dels = sort @dels;
        push(@body, "*Removed files*\n");
        push(@body, "\n");
        push(@body, map { "$_\n" } @dels);
        push(@body, "\n");
    }

    # Modified files list
    if (@mods)
    {
        @mods = sort @mods;
        push(@body, "*Modified files*\n");
        push(@body, "\n");
        push(@body, map { "$_ <$viewvc_url/$_?root=$project_name&revision=$rev&r1=$prev_rev&r2=$rev>\n" } @mods);
        push(@body, "\n");
    }
    push(@body, "\n");

    # Write svn log
    push(@body, "*Differences as text*\n");
    push(@body, "\n");
    @difflines = map { /[\r\n]+$/ ? $_ : "$_\n" } @difflines;
    push(@body, @difflines);

    # Truncate very long body
    @body = cut_array(\@body, 30000);

    # Add empty line to have a prettier output
    # as Mailman will concatenate an "HTML attachment was scrubbed" comment
    push(@body, "\n");
}

# Put together the body of the log message AS HTML
my @body_html;
{
    my $prev_rev = $rev - 1;

    # Write HTML header
    push(@body_html, "<HTML>\n<HEAD>\n</HEAD>\n<BODY>\n");

    # Write summary
    push(@body_html, "<H3>Summary</H3>\n");
    push(@body_html, "<a href=\"$viewvc_url?revision=$rev&amp;root=$project_name&amp;view=rev\">Revision in ViewVC</a><br>\n");
    push(@body_html, "<br>\n");
    push(@body_html, "New Revision: $rev<br>\n");
    push(@body_html, "Author: $author<br>\n");
    push(@body_html, "Date: $date<br>\n");
    push(@body_html, "<br>\n");
    push(@body_html, "Log message:<br>\n");
    push(@body_html, "<pre>\n");
    push(@body_html, encode_for_pre(\@log));
    push(@body_html, "</pre>\n");

    # Added files list
    if (@adds)
    {
        @adds = sort @adds;
        push(@body_html, "<H3>Added files</H3>\n");
        push(@body_html, map { /\/$/ ? "$_<br>\n" : "<a href=\"$viewvc_url/$_?root=$project_name&amp;revision=$rev&amp;view=markup\">$_</a><br>\n" } @adds);
    }

    # Deleted files list
    if (@dels)
    {
        @dels = sort @dels;
        push(@body_html, "<H3>Removed files</H3>\n");
        push(@body_html, map { "$_<br>\n" } @dels);
    }

    # Modified files list
    if (@mods)
    {
        @mods = sort @mods;
        push(@body_html, "<H3>Modified files</H3>\n");
        push(@body_html, map { "<a href=\"$viewvc_url/$_?root=$project_name&amp;revision=$rev&amp;r1=$prev_rev&amp;r2=$rev\">$_</a><br>\n" } @mods);
    }

    # Write svn log (as plain text)
    push(@body_html, "<H3>Differences as text</H3>\n");
    push(@body_html, "<pre>\n");
    @difflines = map { /[\r\n]+$/ ? $_ : "$_\n" } @difflines;
    push(@body_html, encode_for_pre(\@difflines));
    push(@body_html, "</pre>\n");

    # Truncate very long body
    @body_html = cut_array(\@body_html, 30000);

    # Write HTML footer
    push(@body_html, "</BODY>\n</HTML>\n");
}

# Go through each project and see if there are any matches for this
# project.  If so, send the log out.
foreach my $project (@project_settings_list)
  {
    my $match_re = $project->{match_re};
    my $match    = 0;
    foreach my $path (@dirschanged, @adds, @dels, @mods)
      {
        if ($path =~ $match_re)
          {
            $match = 1;
            last;
          }
      }

    next unless $match;

    my @email_addresses = @{$project->{email_addresses}};
    my $userlist        = join(' ', @email_addresses);
    my $to              = join(', ', @email_addresses);
    my $from_address    = $project->{from_address};
    my $hostname        = $project->{hostname};
    my $log_file        = $project->{log_file};
    my $reply_to        = $project->{reply_to};
    my $subject_prefix  = $project->{subject_prefix};
    my $subject;

    if ($commondir ne '')
      {
        $subject = "r$rev - in $commondir: $dirlist";
      }
    else
      {
        $subject = "r$rev - $dirlist";
      }
    if ($subject_prefix =~ /\w/)
      {
        $subject = "$subject_prefix $subject";
      }
    my $mail_from = $author;

    if ($from_address =~ /\w/)
      {
        $mail_from = $from_address;
      }
    elsif ($hostname =~ /\w/)
      {
        $mail_from = "$mail_from\@$hostname";
      }

    my $header;
    #$header->{'debug'} = "/tmp/sendmail.log";
    $header->{'to'} = $to;
    $header->{'from'} = $mail_from;
    $header->{'subject'} = $subject;
    $header->{'replyto'} = $reply_to if $reply_to;
    $header->{'smtp'} = 'localhost';
    $header->{'multipart'} = 'alternative';

    ### Below, we set the content-type etc, but see these comments
    ### from Greg Stein on why this is not a full solution.
    #
    # From: Greg Stein <gstein@lyra.org>
    # Subject: Re: svn commit: rev 2599 - trunk/tools/cgi
    # To: dev@subversion.tigris.org
    # Date: Fri, 19 Jul 2002 23:42:32 -0700
    #
    # Well... that isn't strictly true. The contents of the files
    # might not be UTF-8, so the "diff" portion will be hosed.
    #
    # If you want a truly "proper" commit message, then you'd use
    # multipart MIME messages, with each file going into its own part,
    # and labeled with an appropriate MIME type and charset. Of
    # course, we haven't defined a charset property yet, but no biggy.
    #
    # Going with multipart will surely throw out the notion of "cut
    # out the patch from the email and apply." But then again: the
    # commit emailer could see that all portions are in the same
    # charset and skip the multipart thang.
    #
    # etc etc
    #
    # Basically: adding/tweaking the content-type is nice, but don't
    # think that is the proper solution.

    # Send mail via Mail::Sender
    if (@email_addresses)
      {
        my $mail;
	# Header
	$mail = new Mail::Sender($header) or die "$0: cannot open mail connection: $!\n";
        $mail->OpenMultipart($header) or die "$0: cannot open mail connection: $!\n";
        # Mail as text
        $mail->Part({ctype => 'text/plain; charset=UTF-8',
                    encoding => "quoted-printable",
                    disposition => 'NONE'
                    }) or warn "$0: cannot send mail: $!\n";
        $mail->SendEnc(@body) or warn "$0: cannot send mail: $!\n";
        # Mail as HTML
        $mail->Part({ctype => 'text/html; charset=UTF-8',
                    encoding => "quoted-printable",
                    disposition => 'NONE'
                    }) or warn "$0: cannot send mail: $!\n";
        $mail->SendEnc(@body_html) or warn "$0: cannot send mail: $!\n";
        $mail->Close() or warn "$0: cannot close mail connection: $!\n";
      }

    # Dump the output to logfile (if its name is not empty).
    if ($log_file =~ /\w/)
      {
        if (open(LOGFILE, ">> $log_file"))
          {
            print LOGFILE $header, @body;
            close LOGFILE
              or warn "$0: error in closing `$log_file' for appending: $!\n";
          }
        else
          {
            warn "$0: cannot open `$log_file' for appending: $!\n";
          }
      }
  }

exit 0;

sub usage
{
  warn "@_\n" if @_;
  die "usage: $0 REPOS REVNUM [[-m regex] [options] [email_addr ...]] ...\n",
      "options are\n",
      "  --from email_address  Email address for 'From:' (overrides -h)\n",
      "  -h hostname           Hostname to append to author for 'From:'\n",
      "  -l logfile            Append mail contents to this log file\n",
      "  -m regex              Regular expression to match committed path\n",
      "  -r email_address      Email address for 'Reply-To:'\n",
      "  -s subject_prefix     Subject line prefix\n",
      "\n",
      "This script supports a single repository with multiple projects,\n",
      "where each project receives email only for commits that modify that\n",
      "project.  A project is identified by using the -m command line\n",
      "with a regular expression argument.  If a commit has a path that\n",
      "matches the regular expression, then the entire commit matches.\n",
      "Any of the following -h, -l, -r and -s command line options and\n",
      "following email addresses are associated with this project.  The\n",
      "next -m resets the -h, -l, -r and -s command line options and the\n",
      "list of email addresses.\n",
      "\n",
      "To support a single project conveniently, the script initializes\n",
      "itself with an implicit -m . rule that matches any modifications\n",
      "to the repository.  Therefore, to use the script for a single\n",
      "project repository, just use the other command line options and\n",
      "a list of email addresses on the command line.  If you do not want\n",
      "a project that matches the entire repository, then use a -m with a\n",
      "regular expression before any other command line options or email\n",
      "addresses.\n",
      "\n",
      "Changes for CGAL project:\n",
      "- Email content modified to be closer to CVS commit email format\n",
      "- Email sent in both HTML + text\n";
}

# Return a new hash data structure for a new empty project that
# matches any modifications to the repository.
sub new_project
{
  return {email_addresses => [],
          from_address    => '',
          hostname        => '',
          log_file        => '',
          match_regex     => '.',
          reply_to        => '',
          subject_prefix  => ''};
}

# Start a child process safely without using /bin/sh.
sub safe_read_from_pipe
{
  unless (@_)
    {
      croak "$0: safe_read_from_pipe passed no arguments.\n";
    }

  my $pid = open(SAFE_READ, '-|');
  unless (defined $pid)
    {
      die "$0: cannot fork: $!\n";
    }
  unless ($pid)
    {
      open(STDERR, ">&STDOUT")
        or die "$0: cannot dup STDOUT: $!\n";
      exec(@_)
        or die "$0: cannot exec `@_': $!\n";
    }
  my @output;
  while (<SAFE_READ>)
    {
      s/[\r\n]+$//;
      push(@output, $_);
    }
  close(SAFE_READ);
  my $result = $?;
  my $exit   = $result >> 8;
  my $signal = $result & 127;
  my $cd     = $result & 128 ? "with core dump" : "";
  if ($signal or $cd)
    {
      warn "$0: pipe from `@_' failed $cd: exit=$exit signal=$signal\n";
    }
  if (wantarray)
    {
      return ($result, @output);
    }
  else
    {
      return $result;
    }
}

# Use safe_read_from_pipe to start a child process safely and return
# the output if it succeeded or an error message followed by the output
# if it failed.
sub read_from_process
{
  unless (@_)
    {
      croak "$0: read_from_process passed no arguments.\n";
    }
  my ($status, @output) = &safe_read_from_pipe(@_);
  if ($status)
    {
      return ("$0: `@_' failed with this output:", @output);
    }
  else
    {
      return @output;
    }
}
