clc
clear

%Phase 2 -------------------------------------------------------------
%Input sound from phase 1
[y,Fs] = audioread("./252Crowd.wav");

% filter_first_use = filter_design(90,100,144,154);
% 
% newAudio = filter_first_use.filter(y);
% sound(newAudio, Fs);

filtered_audio = []; % this will contain all the audio that undergoes filteration
root = nthroot((8000/100), 12); % this will allow us to create a log space using the space we have

% this for loop can loop through the channels to create the filters and
% filtered audios accordingly
for channel = 1:12
    fpass1 = 100 * root^(channel-1); % setting the first passband
    fstop1 = fpass1 - 10; % setting the first stopband to allow increase to first passband
    fpass2 = 100 * root^(channel); % setting the second pass band to halt any frequencies beyond it
    fstop2 = fpass2 + 10; % where the decrease should happen to 

    % had this as an error check in case
    if channel == 12
        fpass2 = 8000;
        fstop2 = 8010;
    end

    % Creating a new filter using the specified parameters
    bandpass_filters(channel) = filter_design(fstop1,fpass1, fpass2, fstop2);

    %Using the filter just created to filter the audio and store it in an
    %array accordingly
    filtered_audio(:, channel) = bandpass_filters(channel).filter(y);
end

% Plotting the filtered audio from the lowest channel
figure;
plot(filtered_audio(:, 1));

% Plotting the filtered audio from the highest channel
figure;
plot(filtered_audio(:, 12));

% rectifying the filtered audio to pass to the low pass filter
filtered_audio = abs(filtered_audio);

% collecting the the audios in this array to store after low pass filtering
lp_filtered_audio=[];

for channel=1:12
    lp_filtered_audio(:, channel) = filterTask8().filter(filtered_audio(:, channel));
end

%Phase 3 ---------------------------------------------------------------
% We will need to access the center frequencies of each of these
bandpass_filters;

% This is the envelope of each channel that is going to be multiplied
lp_filtered_audio; 

% This array will hold each of the cosines with the central frequency of
% the bandpasses
cosinesArray = [];

% This will store the carrier x envelope 
carrierEnvelope = [];

% just to make the cosine function
[newInputSize,newChannels] = size(y);
t = linspace(0,(newInputSize/16000),newInputSize);

% This is for task 12 where we need to play the audio
carrierEnvelopeTotal = 0;

for channel=1:12
    lp_filtered_audio(:, channel) = filterTask8().filter(filtered_audio(:, channel));
end

%Phase 3 ---------------------------------------------------------------
% We will need to access the center frequencies of each of these
bandpass_filters;

% This array will hold each of the cosines with the central frequency of
% the bandpasses
cosinesArray = [];

% This will store the carrier x envelope 
carrierEnvelope = [];

% just to make the cosine function
[newInputSize,newChannels] = size(y);
t = linspace(0,(newInputSize/16000),newInputSize);

% This is for task 12 where we need to play the audio
carrierEnvelopeTotal = 0;

for channel=1:12
%   Find two cutoff frequencies
    fpass1 = 100 * root^(channel-1); % setting the first passband
    fpass2 = 100 * root^(channel); % setting the second pass band to halt any frequencies beyond it

%   Calculate the center frequency of each bandpass filter
    central_frequency = sqrt(fpass1*fpass2);
%   Create a cosine with each center frequency and store them in an array
    cosinesArray(:,channel) = cos(2*pi*central_frequency*t);
end

for channel=1:12
    %   Amplitude modulate each cosine signal using the low pass filter
    %   rectified wave
    carrierEnvelope(:,channel) = cosinesArray(:,channel).*lp_filtered_audio(:,channel);
end

for channel = 1:12
    carrierEnvelopeTotal = carrierEnvelopeTotal + carrierEnvelope(:, channel);
end

normalizeCarrier = max(abs(carrierEnvelopeTotal));
carrierEnvelopeTotal = carrierEnvelopeTotal/normalizeCarrier;

sound(carrierEnvelopeTotal,Fs);
