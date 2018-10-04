'''
This is the CLI to access 
'''

import cmd, sys, os
import argparse

import __audio as audio
import __database as database
import __plotting as plotting
# from __audio import *

import numpy as np

# import warnings

# if not sys.warnoptions:
    # warnings.simplefilter("ignore")


class VinylShell(cmd.Cmd):
    intro = 'Welcome to the Vinyl shell.   Type help or ? to list commands.\n'
    prompt = '(((o))):'
    global current_files
    current_files = None

    file_cache = {}

    def do_loadfile(self, arg): 
        'Loads a single audio file for analysis'
        # return __audio.loadfile(*parse(arg))
        # file_name =
        if arg == '':
            arg = raw_input("Please enter the location of the file: ")
        # else: 
            # location = arg
        location = fixstring(arg)
        file_name = str(os.path.splitext(location)[0])
        file = audio.loadfile(location)
        print 'file: ', file
        self.file_cache[str(file_name)] = file
        self.update_current_files(file_name)

    def do_loadfolder(self, arg):
        'Loads a whole folder of .wav files for analysis'
        if arg == '':
            arg = raw_input("Please enter the path of the folder: ")
        location = fixstring(arg)
        for filename in os.listdir(location):
            if filename.endswith(".wav"):  
                file_name = str(os.path.splitext(location)[0])
                file = audio.loadfile(location+filename)
                print 'file: ', file
                print 'filename: ', filename
                self.file_cache[str(filename)] = file
                # self.update_current_files(file_name)
        self.current_files = self.file_cache

    def do_viewfiles(self,arg):
        print 'The current files are: '
        print self.current_files
        print '-----------------------'
        print 'The total file cache is: '
        print self.current_files
        print '-----------------------'




    ##~~ __audio.py COMMANDS ~~##
    def do_RMS_level(self,arg):
        'Returns the RMS level over a section of audio'
        if arg == '': 
            start = int(raw_input('Enter the sample to start from: '))
            npoints = int(raw_input('Enter the number of samples to measure: '))
        else: 
            # print sys.argv
            arg = self.parse(arg)
            start = int(arg[0])
            npoints = int(arg[1])

        # for file in self.current_files:
        for filename, file in self.current_files.iteritems():
            # print filename
            RMS_level = audio.RMS_level(file, start, npoints)
            print filename, ' : ', 20.0*np.log10(RMS_level), 'dB FS'
            # print 'The RMS level of the audio is: ', 20.0*np.log10(RMS_level), 'dB FS'


    def do_detect_leadin(self, arg):
        'Returns the samples between the time the needle drops and the signal starts'
        for filename, file in self.current_files.iteritems():
            leadin_start, signal_start = file.detect_leadin()
            print 'The leadin starts at: ', leadin_start
            print 'The signal starts at: ', signal_start

    def do_leadin_noise(self, arg): 
        'Automatically measures the noise in the leadin of the record'
        for filename, file in self.current_files.iteritems():
            print '%s : leadin noise: %f dB leadin start: %i leadin len: %i' % (filename, 20.0*np.log10(file.leadin_noise()), file.leadin_start, (file.signal_start - file.leadin_start))

    ##~~ __plotting.py COMMANDS ~~##
    def do_plotwave(self, arg):
        'Plots the waveform of an audio file, args include the start point in either time or sample number and the number of points plotted'
        args = self.parse(arg)
        if arg == '':
            start = int(raw_input("Please enter the sample to start plotting: "))
            Npoints = int(raw_input("Please enter the number of points to plot plotting: "))
        else: 
            start = int(args[0])
            Npoints = int(args[1])
            # print 'start do_plot: ', start
            # print 'Npoints do_plot: ', Npoints

        for filename, file in self.current_files.iteritems():
            plotting.plot_wave(file, start, Npoints)
    
    def do_plotfft(self, arg):
        'Plots the fft of an audio file, args include the start point in either time or sample number and the number of points to be analyzed (usually a power of two)'
        plotting.plot_fft(self.current_files, int(args[0]), int(args[1]))
    
    def do_plotleadin(self,arg):
        for filename, file in self.current_files.iteritems():
            plotting.plot_wave(file, file.leadin_start, (file.signal_start-file.leadin_start) )

    def do_scipyspectrograph(self,arg):
        args = self.parse(arg)
        if arg == '':
            start = int(raw_input("Please enter the sample to start plotting: "))
            Npoints = int(raw_input("Please enter the number of points to plot plotting: "))
        else: 
            start = int(args[0])
            Npoints = int(args[1])

        for filename, file in self.current_files.iteritems():
            if file.scipy_freq == False:
                file.scipy_spectograph(start, Npoints)

            plotting.scipy_plotspectograph(file)


    ##~~ generic COMMANDS ~~##

    def do_argtest(self,arg):
        # print sys.argv
        # print arg
        # print type(arg)
        # print len(arg)
        self.parse(arg)
        print arg
        # print sys.argv[0]
        # print sys.argv[1]
    


    def update_current_files(self, arg):
        'Updates the current audio file that will be acted upon'
        # file_name = parse(arg)[0]
        file_name = arg
        self.current_files = {file_name : self.file_cache[file_name]}
        print 'Current file has been updated to: ', file_name
        print 'Current file: ', self.current_files


    def parse(self, arg):
        # parser = argparse.ArgumentParser()
        for char in "\\":
            arg = arg.replace(char,'')

        arg = arg.split()
        # print 'parse(arg): ', (arg)
        # print type(arg)
        # print arg[0]
        if '-m' in arg: 
            print arg.index('-m')
            current_files = self.file_cache
        if '-a' in arg: 
            print "ALL FILES ENABLED"
            current_files = self.file_cache
            print current_files
        
        return arg


# def parse(arg):
#     for char in "\\":
#         arg = arg.replace(char,'')
#     return arg

def fixstring(arg):
    for char in "\\":
        arg = arg.replace(char,'')
    return arg

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-a', '--all', action='store_true', help='-a, --all runs the vinyl shell command on all audio files currently loaded') 

    VinylShell().cmdloop()