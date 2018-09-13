'''
This is the CLI to access 
'''

import cmd, sys, os
import __audio as audio
import __database as database
import __plotting as plotting
# from __audio import *

import numpy as np

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
        location = fixstring(arg)
        file_name = str(os.path.splitext(location)[0])
        file = audio.loadfile(location)
        print 'file: ', file
        self.file_cache[str(file_name)] = file
        self.update_current_file(file_name)
    
    ##~~ __audio.py COMMANDS ~~##
    def do_RMS_level(self,arg):
        'Returns the RMS level over a section of audio'
        if arg == '': 
            start = int(raw_input('Enter the sample to start from: '))
            npoints = int(raw_input('Enter the number of samples to measure: '))
        else: 
            start = parse(args)[0]
            npoints = parse(args)[1]

        RMS_level = audio.RMS_level(self.current_file, start, npoints)
        print 'The RMS level of the audio is: ', 20.0*np.log10(RMS_level), 'dB FS'

    ##~~ __plotting.py COMMANDS ~~##
    def do_plotwave(self, arg):
        'Plots the waveform of an audio file, args include the start point in either time or sample number and the number of points plotted'
        args = parse(arg)
        if arg == '':
            start = int(raw_input("Please enter the sample to start plotting: "))
            Npoints = int(raw_input("Please enter the number of points to plot plotting: "))
        else: 
            start = int(arg[0])
            Npoints = int(arg[1])
        plotting.plot_wave(self.current_file, start, Npoints)
    
    def do_plotfft(self, arg):
        'Plots the fft of an audio file, args include the start point in either time or sample number and the number of points to be analyzed (usually a power of two)'
        plotting.plot_fft(self.current_file, int(args[0]), int(args[1]))
    

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

def fixstring(arg):
    for char in "\\":
        arg = arg.replace(char,'')
    return arg

if __name__ == '__main__':
    VinylShell().cmdloop()