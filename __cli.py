'''
This is the CLI to access 
'''

import cmd, sys, os
import __audio as audio
# from __audio import *



class VinylShell(cmd.Cmd):
    intro = 'Welcome to the Vinyl shell.   Type help or ? to list commands.\n'
    prompt = '(((o))):'
    current_file = None

    file_cache = {}

    def do_loadfile(self, arg): 
        'Loads a single audio file for analysis'
        # return __audio.loadfile(*parse(arg))
        # file_name =
        if arg == '':
            location = raw_input("Please enter the location of the file: ")
        else: 
            location = arg

        file_name = str(os.path.splitext(location)[0])
        file = audio.loadfile(location)

        self.file_cache[str(file_name)] = file
        self.update_current_file(file_name)
        

    def update_current_file(self, file_name):
        self.current_file = self.file_cache[file_name]
        print 'Current file has been updated to: ', file_name


def parse(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    return tuple(map(int, arg.split()))

if __name__ == '__main__':
    VinylShell().cmdloop()