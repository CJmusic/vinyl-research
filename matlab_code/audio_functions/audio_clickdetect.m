%{
Lining up clicks: 

clicks = [ sample numbers ] 


%} 


function clicks = audio_clickdetect(data, fs);
   time = (0:length(data)-1)/fs; 
   d_data = diff(data);%/diff(time);
   dd_data = diff(d_data);%/diff(time);
   clicks = [];
   for i = (1:length(data));
      if i + 1024 < length(data);
          threshold = rms(data(i:i+1024));
          if d_data(i) > threshold;
                click = i;
                %click_timestamp = i/fs;   
                clicks = [clicks, click];   
          end
      end
   end 
   
   %plotting below
   
  % x = zeros(length(clicks));
  % figure(1);
  % plot(time, data); grid on; hold on;
  % plot(clicks/fs,x, 'r.', 'MarkerSize', 20);
  % figure(2); 
  % subplot(2,1,1);
  % plot(time(1:end-1),d_data); grid on; hold on;
  % plot(clicks/fs,x,'r.', 'MarkerSize', 20);
  % subplot(2,1,2); 
  % plot(time(1:end-2),dd_data); grid on; hold on;
  % plot(clicks/fs,x,'r.', 'MarkerSize', 20);
end 


