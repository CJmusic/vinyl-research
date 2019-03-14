%{
Lining up clicks: 

clicks = [ sample numbers ] 

% a click has a sustained high derivative 
% the energy of a click is between 2-4 kHz

%} 


function clicks = audio_clickdetect(data, fs);
   % implement a high pass filter 
   %fc = 1000.0 
   %[b, a] = butter(6, fc/(fs/2),'high'); % this should 
   %data = filter(b, a, data);

   time = (0:length(data)-1)/fs; 
   d_data = diff(data,1)*fs;% since delta_t = 1/fs:w;
   clicks = [];
   threshold = rms(d_data(1:1025));
   for i = (1:length(data)-1);
      if i > 512 && i + 512 < length(d_data);
        threshold = threshold - abs(d_data(i-511)) + abs(d_data(i+511));
      end
      if d_data(i) > threshold;
        click = i;
        %click_timestamp = i/fs;   
        clicks = [clicks, click];   
     end
   end
   disp('Number of clicks');
   size(clicks) 
   
   %plotting below
   
   %clf(figure(1));clf(figure(2)); 
   %figure(1); grid on; hold on;
   %plot(time, data); 
   % for xi = 1:length(clicks);
   %      x1 = time(clicks(xi));
   %      line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
   % end
   % title('Click detection: Amplitude vs. Time'); 
   % xlabel('Time [s]');
   % ylabel('Amplitude');

   %figure(2); grid on; hold on;
   %t_diff = length(time) - length(d_data)
   %plot(time(1:end-t_diff),d_data); grid on; hold on;
   %for xi = 1:length(clicks);
   %     x1 = time(clicks(xi));
   %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
   % end

   % title('Click detection: 1st Derivative vs. Time'); 
   % xlabel('Time [s]');
   % ylabel('dA/dt (s^-1)');
end 


