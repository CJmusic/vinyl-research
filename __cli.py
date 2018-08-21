'''
This is the CLI to access 
'''

import cmd, sys
# import __audio
from __audio import *



class VinylShell(cmd.Cmd):
    intro = 'Welcome to the Vinyl shell.   Type help or ? to list commands.\n'
    prompt = '(((o))):'
    file = None

    def do_loadfile(self, arg): 
        'Loads a single audio file for analysis'
        # return __audio.loadfile(*parse(arg))
        # file_name = 
        if arg == '':
            return loadfile(raw_input("Please enter the location of the file: "))
        else: 
            return loadfile(arg)


def parse(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    return tuple(map(int, arg.split()))

if __name__ == '__main__':
    VinylShell().cmdloop()