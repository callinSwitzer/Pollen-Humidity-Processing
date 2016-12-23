%% Callin Switzer
% Experiments from 08/04/2014

%% setup from 03/13/2016
% Aperture = 8 
% frequency of vibration = 280, sec = 0.125
% shutter 1/3000
% fps = 500

%% Last updated from experiments 08.02.2014
% note: make sure videos are stored w/ correct number of characters
%%


clear all
close all

% Part 1 Setup
cd('/Users/callinswitzer/Desktop/flowerShake03122016')
videos = dir2cell('*.avi') % function was downloaded 
% videos(1) = [] % gets rid of ampcheck

charNum = 15; % number of characters to extract from title
%%  Gets bottom of stamen and width of calibration
% object from the beginning of each video
for ii= length(videos)
vidnum = ii; 
% Define a handle to the video 

mm = VideoReader(char(videos(vidnum)));
 
% Load in an image, display, and select bottom of stamen

im = read(mm, 1);

% info on the image 
% whos im

% remove the color channel 
im = im(:,:,1);
im = im < 50; % we need a black & white image to find the sphere,
              % so we're setting a threshold of 50


% change the colormap
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(im);
colormap jet;

% Select an x, y range of interest
display('Select just below the bottom of the stamen')
bottomOfStamen = round(ginput(1)); 
bottomOfStamen = bottomOfStamen(2); % Want to cut out everything above this


% Select an x, y range of interest to exclude (i.e. the calibration stuff)
display('Select left and under the calibration object')
tapeExclude = round(ginput(1));  % cut everything above and left of this

% get sphere width for calibration -- all we need is the 
% radius of the circle
% we're only looking at the right-hand side of the image...
% remove the parentheses after im to look at whole image
[centers, radii] = imfindcircles(im(1:400, 500:1000),[30 45])
% imagesc(im(1:400, 500:1000));
% viscircles(centers, radii,'EdgeColor','b');


%display('Select the two sides of the sphere -- perpendicular to the wire')
%wireWidth = (ginput(2));
%ww = sqrt( (wireWidth(2)-wireWidth(1))^2 + (wireWidth(4)-wireWidth(3))^2 )

%display('select leftmost part of sphere')
%xCut = ginput(1); 

% ww = (abs(wireWidth(2) - wireWidth(1)));

ww = radii * 2; 
pxToMmConversion = ww / 1.39 % px/mm -- sphere is 1.39mm wide

% Do all of the above for each file in the folder

vidName = char(videos(vidnum));
WidthOfSphere = ww; 
stamEnd = table(WidthOfSphere, pxToMmConversion, bottomOfStamen,tapeExclude)
calFileName = [vidName(1:charNum) 'EndOFStamen' '.csv']

writetable(stamEnd, calFileName);

close all % close the figure after you're done

end
%% 
close all

%%  Look at 12 mm section of each video

for kk = length(videos)
    nums = kk; 
    vidName = char(videos(nums));
    fn = [vidName(1:charNum) 'EndOFStamen' '.csv']
    AA = importdata(fn, ',')
    ww = AA.data(1)
    pxToMmConversion = AA.data(2)
    bottomOfStamen = AA.data(3)
    tapeLeft = AA.data(4)
    tapeBottom = AA.data(5)
    xCut = 1; % round(AA.data(4))
    
    
    % Define a region that is 12 mm long and as wide as possible
    vertDist = 12; % millimeters
    
    xx = [xCut; 1024];
    yy = [bottomOfStamen; bottomOfStamen + vertDist*pxToMmConversion];
    xx = xx(1):xx(2);
    yy = yy(1):yy(2);
    
    %  Auto calibration  -- go along top of screen until you hit black wire, then
    % find the other side of the black wire -- do this later
    % Show the range of the selected image
    
    mm = VideoReader(vidName);
    
    % Load in an image, display, and select bottom of stamen
    
    im = read(mm, 1);
    
    % info on the image
    % whos im
    
    % remove the color channel
    im = im(:,:,1);
    figure(1)
    imagesc(im(yy, xx)) % images in matlab are y, x
    
end
close all
%% read in each video one at a time preprocess them all
% make black and white and save as .mp4
for ii =  length(videos)
    nums = ii; 
    vidName = char(videos(nums));
    fn = [vidName(1:charNum) 'EndOFStamen' '.csv'];
    AA = importdata(fn, ',');
    ww = AA.data(1);
    pxToMmConversion = AA.data(2);
    bottomOfStamen = AA.data(3);
    tapeLeft = AA.data(4);
    tapeBottom = AA.data(5);
    xCut = 1; %round(AA.data(4))
    
    % Define a region that is 12 mm long and as wide as possible
    vertDist = 12; % millimeters
    
    xx = [xCut; 1024];
    yy = [bottomOfStamen; bottomOfStamen + vertDist*pxToMmConversion];
    xx = xx(1):xx(2);
    yy = yy(1):yy(2);
    
    %  Auto calibration  -- go along top of screen until you hit black wire, then
    % find the other side of the black wire -- do this later
    % Show the range of the selected image
    
    mm = VideoReader(vidName);
    
    % Load in an image, display, and select bottom of stamen
    
    im = read(mm, 1);
    
    % info on the image
    % whos im
    
    % remove the color channel
    im = im(:,:,1);
    figure
    imagesc(im(yy, xx)) % images in matlab are y, x
    
    % Get some info on the video
    NumberOfFrames  = mm.NumberOfFrames;
    Width           = mm.Width;
    Height          = mm.Height;
    % closing figure seems to help
    % For speed, preallocate the array
    frames = uint8(zeros(length(yy), length(xx), NumberOfFrames));
    
    % Load in the frames into the preallocated array
    
    for kk=1:NumberOfFrames
        tmp = read(mm,kk); % load in the kk'th image
        %imagesc(tmp);
        tmp(1:tapeBottom, tapeLeft:1024) = 0; % remove tape from video
        frames(:,:,kk) = tmp(yy,xx);    % save the reduced x,y portion of the image
        %imagesc(tmp);
        kk
    end
    display('frames loaded -- subtracting background')
    
    % Lets play a video to find first pollen out
    % close all
    % f = figure(1)
    % display('Press spacebar to stop when pollen leaves the screen')
    % for kk=1:NumberOfFrames
    %
    %     imagesc(frames(:,:,kk));
    %     drawnow;
    %
    %     if get(f,'currentkey') ~= ' '
    %      disp('stopped');
    %      break;
    %     end
    %     %kk
    % end
    %
    % firstPollenOutOfScreen = kk
    
    % Make a background image from median filter
    % I could probably do this with less than the whole video -- maybe just the
    % first few frames
    bkg = median(frames, 3);
    
    %imagesc(bkg);
    imwrite(bkg, 'bkg.tif');
    
    display('backgound subtracted')
    % thresholded + mediand + mean light intensity Pollen video
    lightintensity = squeeze(sum(sum(frames)))/(size(frames,1)*size(frames,2));
    % colormap gray
    % figure(1)
    thresh = 0.9;
    
    display('finding first pollen out bottom of screen')
    % find the first pollen out of the bottom of the frame
    foo = zeros(mm.NumberOfFrames);
    for kk=1:(mm.NumberOfFrames)
        normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
        normalizedbkg = double(bkg)/mean(bkg(:));
        bwimage = ((normalizedimage./normalizedbkg < thresh));
        pollencount3(kk) = nnz(bwimage);
        imagesc((normalizedimage./normalizedbkg < thresh));
        bottomRow = bwimage((length(bwimage(:,1)) - 1), :);
        foo(kk) = nnz(bottomRow);
    end
    % find the frame number where pollen left the frame
    for jj = 1:length(foo)
        if foo(jj) > 1;
            firstPollenOutOfScreen = jj
            break
        end
        if jj == mm.NumberOfFrames
            firstPollenOutOfScreen = jj
            break
        end
    end
    
    
    
    % playback to check FPO is correct
    
    display('checking that first pollen out of frame is correct')
    colormap gray
    for kk = firstPollenOutOfScreen - 20: firstPollenOutOfScreen
        normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
        normalizedbkg = double(bkg)/mean(bkg(:));
        bwimage = ((normalizedimage./normalizedbkg < thresh));
        pollencount3(kk) = nnz(bwimage);
        imagesc((normalizedimage./normalizedbkg < thresh));
        bottomRow = bwimage((length(bwimage(:,1)) - 1), :);
        foo(kk) = nnz(bottomRow);
        drawnow;
    end

    % plot(foo, 'b')
    %
    % Convert pixels of pollen to mm of pollen
    % pxToMmConversion is about 40 px per mm -- 40^2 px per sq mm
    pollen = pollencount3 / pxToMmConversion^2;
    % figure
    % plot(pollen, 'r') % gives area of pollen in mm^2
    
    
    % Find first frame where pollen enters and blank frames
    
    for kk=1:mm.NumberOfFrames
        if pollencount3(kk) > 2
            firstPollen = kk % starts one frame before pollen enters
            break
        end
    end
    
    % remove zeros at the end
    for kk=firstPollen + 50:mm.NumberOfFrames
        if pollencount3(kk) == 0
            blankFrames = kk; % starts one frame before pollen enters
            break
        else blankFrames = mm.NumberOfFrames;
        end
    end
    
    blankFrames;
    
    hh = figure(2);
    plot([firstPollen:blankFrames]/1000,pollen(firstPollen:blankFrames), 'r')
    xlabel('Times (s)','FontSize',20, 'FontName', 'Times');
    ylabel('Area of pollen (mm^2)','FontSize',20, 'FontName', 'Times');
    set(gca,'FontSize',20, 'FontName', 'Times');
    axis([0 inf -inf inf]);
    
    set(gcf,'Color','w');
    
    % saveas(hh, fyle, 'pdf') % saves the figure
    
    % save the median light intensity Pollen video from just before pollen comes out
    % until after all the pollen is gone
    % mkdir('PPPollen')
    fyle = [vidName(1:charNum)];
    
    video = VideoWriter(fyle,'MPEG-4');
    video.FrameRate = 150;
    video.Quality = 50;
    
    
    open(video);
    
    cmap = gray(2);
    display('saving black and white video')
    for kk=firstPollen:blankFrames
        normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
        normalizedbkg = double(bkg)/mean(bkg(:));
        bwimage = ((normalizedimage./normalizedbkg < thresh));
        % pollencount3(kk) = nnz(bwimage);
        tmp = cmap(bwimage + 1, :);
        tmp = reshape(tmp, size(frames,1), size(frames,2), 3);
        writeVideo(video, tmp);
        kk
    end
    
    close(video)
    
    
    % write a csv file with titled with the filename and pollencount
    % mkdir('PollenCSVs')
    csvFile = [vidName(1:charNum) '.csv'];
    pollen1 = pollen(firstPollen:blankFrames);
    %
    outp = table(pxToMmConversion, firstPollen, firstPollenOutOfScreen, blankFrames);
    metaFileName = [vidName(1:charNum) 'Meta' '.csv'];
    csvwrite(csvFile, pollen1(:));
    writetable(outp, metaFileName);
    
    display('file saved')
    nums
end

close all

%%
