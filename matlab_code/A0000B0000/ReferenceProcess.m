%~~~~~~~~~~~~~~~~~ LOAD REFERENCE ~~~~~~~~~~~~~~~~~%
try 
    [ref, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
catch
    [ref, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\03141_A0000B0000r28a.wav');
end 
%~~~~~~~~~~~~~~~~~~ Reference info ~~~~~~~~~~~~~~~~%
folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\';


timestamps_ref = [0, 60, 90, 122, 158, 180, 246, 266, 304, 324, 362, 382, 417.5];
% this is how many seconds each signal is according to Chris Muths track listing
lengths = [60, 30, 31, 36, 21, 66, 20, 37, 19, 37, 19, 37, 19]; %starts with 1kHz
signal_names = {'leadin',    % 1
                '1kHz',      % 2
                '10kHz',     % 3
                '100Hz',     % 4
                'sweep',     % 5
                'quiet',     % 6
                '3150Hz',    % 7 
                '1kHzL',     % 8
                'sweepL',    % 9
                '1kHzR',     % 10
                'sweepR',    % 11
                '1kHzV',     % 12
                'sweepV',    % 13
                'transition',% 14
                '1kHz2',     % 15
                '10kHz2',    % 16
                '100Hz2',    % 17
                'freqsweep2',% 18
                'quiet2',    % 19
                '3150Hz2',   % 20
                '1kHzL2',    % 21
                'sweepL2',   % 22
                '1kHzR2',    % 23
                'sweepR2',   % 24
                '1kHzV2',    % 25
                'sweepV2',   % 26
                'leadout'    % 27
};
%% Reference 02072019_A0000B000r27a.wav 
offset = 10.625; 
transition = 517.375; 
% lockoutClipped = 953.746;
% lagdiff = []

% 031418_A0000B0000r27a.wav as reference timestamps
offset = 10.625; 
timestamps =       [[0, 61],    % 1. 1 kHz
                    [61,91],    % 2. 10 kHz
                    [91,121],   % 3. 100 Hz
                    [121,159],  % 4. sweep
                    [159,180],  % 5. quiet
                    [180,245],  % 6. 3150 Hz
                    [245,267],  % 7. 1 kHz left
                    [267, 302], % 8. sweep left
                    [302, 325], % 9. 1 kHz right
                    [325, 361], % 10. sweep right
                    [361, 383], % 11. 1 kHz vertical
                    [383, 418], % 12. sweep vertical
                    [418, 515], % 13. transition
                    [515, 578], % 14. 1 kHz
                    [578, 608], % 15. 10 kHz
                    [608, 639], % 16. 100 Hz
                    [639, 676], % 17. sweep
                    [676, 698], % 18. quiet
                    [698, 760], % 19. 3150 Hz
                    [760, 785], % 20. 1 kHz left
                    [785, 820], % 21. sweep left
                    [820, 842], % 22. 1 kHz right
                    [842, 878], % 23. sweep right
                    [878, 900], % 24. 1 kHz vertical
                    [900, 938]];% 25. sweep vertical  
                    % [938, 950]];               
                    %% dont forget lead in and leadout
timestamps = timestamps + offset;

timeref = (0:length(ref))/fs;
%~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%
t = 1;
sigtime = timeref(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
sig = ref(floor(timestamps(t,1)*fs) : floor(timestamps(t,2)*fs),:);
sigRMS=rms(sig);
normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
ref(:,1)=ref(:,1)/normalization(1);% now normalized to 40cm/s peak    
ref(:,2)=ref(:,2)/normalization(2);% now normalized to 40cm/s peak 
normalization_L = normalization(1);
normalization_R = normalization(2);

%save ref to a mat file

%run click detect on ref and save as a mat file 

[cref(:,1), REFS_L] = ClickDetect(ref(:,1),200,20);
[cref(:,2), REFS_R] = ClickDetect(ref(:,2),200,20);


% csvwrite(strcat(folder,'cref.txt'),cref)
% csvwrite(strcat(folder,'ref.txt'),ref)
csvwrite(strcat(folder,'03141_A0000B0000r28a-REFS_L.txt'),REFS_L)
csvwrite(strcat(folder,'03141_A0000B0000r28a-REFS_R.txt'),REFS_R)

