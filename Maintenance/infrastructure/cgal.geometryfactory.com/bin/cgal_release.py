'''Python module to create and publish CGAL releases from a branch'''

import os

class Release:
    '''class to create a CGAL release from a branch
       optionally, the release can be internal
    '''

    def __init__(self, branch, internal=False) :
        self.branch = branch
        self.internal = internal
        self.cwd = f'$HOME/CGAL/create_internal_release-{self.branch}-branch'
        self.repo = f'$HOME/CGAL/branches/CGAL-{self.branch}-branch.git'
        self.extra_options = ' --public'

    def command(self):
        '''return the command to create and publish the release'''
        return f"PATH=/home/lrineau/bin-cmake3:/bin:/usr/bin:/home/lrineau/bin; cd {self.cwd} && /usr/bin/time scl enable rh-git29 -- $HOME/bin/create_release {self.repo}{self.extra_options} --do-it"

    def __str__(self) :
        msg  = f"{'internal ' if self.internal else ''}release from {self.branch}\n" \
               f"cwd: {self.cwd}\nrepo: {self.repo}\n" \
               f"command:\n{self.command()}"
        return msg

    def __call__(self) :
        if os.system(self.command()) != 0 :
            raise RuntimeError(f"Error while creating {'internal ' if self.internal else ''}release from {self.branch}")

    INTERNAL = True

class InternalRelease(Release) :
    '''class to create an internal CGAL release from a branch'''
    def __init__(self, branch) :
        super().__init__(branch, Release.INTERNAL)
        self.extra_options = ' --integration'

integration = InternalRelease("integration")
integration.repo = '$HOME/CGAL/branches/integration.git'
integration.cwd = '$HOME/CGAL/create_internal_release'
master = Release("master")
master.repo = '$HOME/CGAL/branches/master.git'
master.cwd = '$HOME/CGAL/create_internal_release'

def release(branch) :
    '''Convenience function to create a release from a branch'''
    return Release(branch)

if __name__ == '__main__':
    print("This file is a Python module. Use create_internal_release_of_the_day.py instead.")
