'''
This is the CLI to access 
'''

import cmd, sys

class VinylShell(cmd.Cmd):
    intro = 'Welcome to the Vinyl shell.   Type help or ? to list commands.\n'
    prompt = '(((o))):'
    file = None

    # ----- basic turtle commands -----
    # def do_forward(self, arg):
    #     'Move the turtle forward by the specified distance:  FORWARD 10'
    #     forward(*parse(arg))
    # ----- record and playback -----
  

def parse(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    return tuple(map(int, arg.split()))

if __name__ == '__main__':
    VinylShell().cmdloop()