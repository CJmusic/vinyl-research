%{
Lining up clicks: 

clicks = [ sample numbers ] 


%} 


function clicks = audio_clickdetect(data, fs);
   time = (0:length(data)-1)/fs; 
   d_data = diff(data,1)*fs;% since delta_t = 1/fs:w;
   clicks = [];
   for i = (1:length(data));
      if i + 1024 < length(data);
          threshold = 17*rms(d_data(i:i+1024));
          if d_data(i) > threshold;
                click = i;
                %click_timestamp = i/fs;   
                clicks = [clicks, click];   
          end
      end
   end
   disp('Number of clicks');
   size(clicks) 
   
   %plotting below
   
   clf(figure(1));clf(figure(2)); 
   figure(1); grid on; hold on;
   plot(time, data); 
    for xi = 1:length(clicks);
         x1 = time(clicks(xi));
         line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end
    title('Click detection: Amplitude vs. Time'); 
    xlabel('Time [s]');
    ylabel('Amplitude');

   figure(2); grid on; hold on;
   t_diff = length(time) - length(d_data)
   plot(time(1:end-t_diff),d_data); grid on; hold on;
   for xi = 1:length(clicks);
        x1 = time(clicks(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end

    title('Click detection: 1st Derivative vs. Time'); 
    xlabel('Time [s]');
    ylabel('dA/dt (s^-1)');
end 

