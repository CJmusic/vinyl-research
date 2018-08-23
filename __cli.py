'''
This is the CLI to access 
'''

import cmd, sys, os
import __audio as audio
import __database as database
import __plotting as plotting
# from __audio import *



class VinylShell(cmd.Cmd):
    intro = 'Welcome to the Vinyl shell.   Type help or ? to list commands.\n'
    prompt = '(((o))):'
    global current_file
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
        print 'file: ', file
        self.file_cache[str(file_name)] = file
        self.update_current_file(file_name)
        
    def do_plotwave(self, arg):
        'Plots the waveform of an audio file, args include the start point in either time or sample number and the number of points plotted (usually a power of two)'
        args = parse(arg)
        plotting.plot_wave(self.current_file, int(arg[0]), int(arg[1]))
    
    def do_plotfft(self, arg):
        plotting.plot_fft(self.current_file)
    

    def update_current_file(self, arg):
        'Updates the current audio file that will be acted upon'
        # file_name = parse(arg)[0]
        file_name = arg
        self.current_file = self.file_cache[file_name]
        print 'Current file has been updated to: ', file_name
        print 'Current file: ', self.current_file


def parse(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    return tuple(map(int, arg.split()))

if __name__ == '__main__':
    VinylShell().cmdloop()