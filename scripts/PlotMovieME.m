pth = 'C:\Users\munib\Desktop';
fn = 'vlc-record-2024-08-13-18h48m49s-JPV13_2023-10-03_behav_cam_0_date_2023_10_03_time_11_21_00_v001.mp4';

% Create a VideoReader object
video = VideoReader(fullfile(pth,fn));

% Initialize variables to store frames
frame1 = readFrame(video); % Read the first frame
frame1Gray = rgb2gray(frame1); % Convert to grayscale
frameDiffs = []; % To store differences
frameframe = [];

% Loop through each frame in the video
iii = 1;
while hasFrame(video)
    disp(num2str(iii))

    frame2 = readFrame(video); % Read the next frame
    frame2Gray = rgb2gray(frame2); % Convert to grayscale
    frameframe = cat(3,frameframe,frame2Gray);
    
    % Compute the difference between consecutive frames
    frameDiff = abs(double(frame2Gray) - double(frame1Gray));
    
    % Store the difference
    frameDiffs = cat(3, frameDiffs, frameDiff); % Append along the third dimension
    
    % Update frame1 for the next iteration
    frame1Gray = frame2Gray;

    iii = iii + 1;
end

%% save movie gif

clear s
close all
f = figure;

gifFile = 'concatMov_gif.gif';

for iframe = 1:size(frameframe,3)
    % Display the frame
    imshow(frameframe(:,:,iframe));
    
    s(iframe) = getframe(f);
    exportgraphics(f,gifFile,'Append',true,'Resolution',30);
end

%% save me gif

clear s
close all
f = figure('visible','on');
hold on;

gifFile = 'concatME_gif.gif';

cmap = inferno;
colormap(cmap)

for iframe = 1:size(frameDiffs,3)
    colormap(cmap)
    % Display the frame
    thisframe = frameDiffs(:,:,iframe);
    thisframe = imgaussfilt(thisframe);
    thisframe(1:350,:) = thisframe(1:350,:)./30;
    thisframe(351:end,:) = thisframe(351:end,:)./25;
    imshow(thisframe);
    colormap(cmap)
    
    s(iframe) = getframe(f);
    % exportgraphics(f,gifFile,'Append',true,'Resolution',30);
end


%% SAVE

outpth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
vOut = VideoWriter(fullfile(outpth,'concatVideo'),'Uncompressed AVI');
vOut.FrameRate = video.FrameRate;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)

