#! /bin/env python3

from __future__ import print_function
from builtins import (bytes)

import github, sys, os, configparser, argparse, textwrap
from github import Github, Label, Issue, PullRequest
from xdg.BaseDirectory import load_first_config, xdg_config_home

from dulwich.repo import Repo
from dulwich.walk import Walker


config = configparser.ConfigParser()
if load_first_config('CGAL'):
    config_file = os.path.join(load_first_config('CGAL'), 'test-pygithub_rc')
else:
    config_file = os.path.join(xdg_config_home, 'test-pygithub_rc')
config.read(config_file)

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
Print the URLs or CGAL pull-request with a given label.''',\
                                 epilog=textwrap.dedent('''\
  The file {} must contain, in a `[main]` section, two keys:
    - user,
    - token
  to define Github credentials. For example:

      [main]
      user=lrineau
      token=7f3f81294bd965232ca3c657a23c16729dac663f

  '''.format(config_file)))
parser.add_argument('label', metavar='label', type=str,
                    help='the label used to select pull requests.\
                    If it starts with "4", then "Under testing in CGAL-"\
                    is added as a prefix.')
parser.add_argument('--verbose', action='store_const', const=True,
                    help='display verbose outputs to stderr')
args = parser.parse_args()

def print_verbose(str, *argv, **kwargs):
    if(verbose):
        print("... {}".format(str),\
              file=sys.stderr, *argv, **kwargs)

u=config['main']['user']
s=config['main']['token']
verbose=args.verbose

print_verbose("Getting commits in current branch", end='')
sys.stderr.flush()
repo = Repo.discover()
print_verbose(" branch: {}".format(repo.refs.read_ref('HEAD')),\
              end='')
sys.stderr.flush()
walk = repo.get_walker(include=repo.head(),\
                       exclude=[repo.get_refs()[b'refs/remotes/cgal/master']])
commits=list(w.commit.id for w in walk)
print_verbose(" {} commits".format(len(commits)))

print_verbose('Signing into Github using login "{}"...'.format(u))
g = Github(login_or_token=u,password=s);
CGAL = g.get_organization("CGAL")
cgal = CGAL.get_repo("cgal")

if(args.label[0]=='4'):
    label_name = "Under testing in CGAL-{}".format(args.label)
else:
    label_name = args.label

try:
    Ic = cgal.get_label(label_name)
except github.UnknownObjectException:
    print('ERROR: Unknown label "{}"'.format(label_name), file=sys.stderr)
    print('Know labels are: {}'.format([label.name\
                                        for label in cgal.get_labels()])\
          , file=sys.stderr)
    exit()

issues = {}

print_verbose("Gettings pulls...")
pulls = cgal.get_pulls()
[pr for pr in pulls]
print_verbose('Gettings issues with label "{}"...'.format(Ic.name))
for issue in cgal.get_issues(labels=[Ic]):
    try:
        pr = [pr for pr in pulls if issue.number == pr.number][0]
        issues[issue.number] = {"url": issue.html_url,\
                                         "sha": bytes(pr.head.sha, 'ascii')}
        print_verbose("PR {}".format(issue.html_url))
        print_verbose("  head:   {}".format(bytes(pr.head.sha, 'ascii')))
        print_verbose("  labels: {}".format([label.name for\
                                                   label in issue.labels]))
    except IndexError:
        print_verbose('Warning: issue #{} with label "{}"\
                      is not a pull-request'\
                      .format(issue.number, label_name))

for nb in issues.keys():
    print(issues[nb]["url"])
    if issues[nb]["sha"] not in commits:
        print_verbose("Warning: {} head {} is not in the current branch".\
                      format(issues[nb]["url"], issues[nb]["sha"]))
