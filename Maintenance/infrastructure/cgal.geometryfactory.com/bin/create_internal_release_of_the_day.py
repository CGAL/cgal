#! /bin/python3

'''This script is called by a cron job every day.
   It creates and publish a release tarball.
'''

import sys
import os
import datetime
import locale
from pathlib import Path
sys.path.append(Path(__file__))
from cgal_release import release, integration, master

# Define a dictionary that maps day of the week to an action
actions = {
    "Monday": integration,
    "Tuesday": integration,
    "Wednesday": integration,
    "Thursday": integration,
    "Friday": release("5.5"),
    "Saturday": release("5.4"),
    "Sunday": master
}

if __name__ == '__main__':
    # Get the current day of the week, or get it from the command line
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    try:
        DAY_OF_THE_WEEK = sys.argv[1]
    except IndexError:
        DAY_OF_THE_WEEK = datetime.datetime.now().strftime("%A")

    # Look up the action for the current day of the week in the dictionary
    create_release = actions[DAY_OF_THE_WEEK]

    # Then create the release tarball
    if os.system(create_release.command()) != 0 :
        raise RuntimeError(f"ERROR while creating release tarball")
