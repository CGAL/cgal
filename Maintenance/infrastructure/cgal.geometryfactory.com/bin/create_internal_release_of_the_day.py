#! /bin/python3

"""This script is called by a cron job every day.
   It creates and publish a release tarball.
"""

import os
import datetime
import locale
import argparse
from cgal_release import release, integration, master, beta_release, master, beta_release_from_master

# Define a dictionary that maps day of the week to an action
actions = {
    "Monday": integration,
    "Tuesday": integration,
    "Wednesday": integration,
    "Thursday": integration,
    "Friday": release("5.5"),
    "Saturday": release("5.6"),
    "Sunday": master,
}


def main():
    """Entry point of the script."""
    # Get the current day of the week, or get it from the command line
    locale.setlocale(locale.LC_ALL, "en_US.UTF-8")

    day_help = f"Day of the week (default: {datetime.datetime.now().strftime('%A')})"
    day_help += f" possible values: ({', '.join(actions.keys())})"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "day", help=day_help, nargs="?", default=datetime.datetime.now().strftime("%A")
    )
    parser.add_argument("-n", "--dry-run", action="store_true")
    args = parser.parse_args()

    # Look up the action for the current day of the week in the dictionary
    create_release = actions[args.day]

    if args.dry_run:
        print(create_release)
        return

    # Then create the release tarball
    if os.system(create_release.command()) != 0:
        raise RuntimeError("ERROR while creating release tarball")


if __name__ == "__main__":
    main()
