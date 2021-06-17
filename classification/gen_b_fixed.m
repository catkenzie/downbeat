
% PARAMETERS TO EVALUATE:
%   FFT_length (DONE)
%   Hop_length (DONE)
%   Window type for STFT (DONE)
%   Frequency band partitioning (arrayOfFrequencies)
%   Matrix Factorization Approaches (DONE)
%   Onset thresholding
%       Rising thresholds
%       Falling thresholds
%       Decay times
%   Additive noise
%
% PROGRESS SO FAR:
%   The resulting amplitude envelopes of the pieces of the kit played in
%   the test beat are, for all intents and purposes, the same regardless of
%   the approach used for determining the basis matricies of each piece
%   (4A averaging vs 4B-1 NNMF fixed # iters vs 4A-2 NNMF thresholding).
%   This makes sense because, as we discovered, the factorization reaches a
%   close approximation after 2 iterations
%
%   The resulting amplitude envelopes of the pieces of the kit played in
%   the test beat are, for all intents and purposes, the same regardless of
%   the approach used for running the NNMF (5A NNMF fixed # iters vs 5B NNMF
%   thresholding) when num_nnmf_iters = 100 and nnmf_cost_threshold =
%   0.001. One difference is the total number of iterations that are
%   performed is higher when using 5A.  On the other hand, 5B requires
%   calculating the divergence cost and comparing it to the threshold
%
%   Using FFT_length = 2048 and Hop_length=441 produces a false hi-hat onset
%   at the third snare hit. The amplitude of this hit is higher than that
%   of the previously detected correct hit, which poses a challenge for
%   setting an appropriate thresholding value.  Using FFT_length = 4096 and
%   Hop_length=441 greatly decreases the time resolution and may be
%   computationally too heavy for our use, but it does reduce the amplitude
%   of the false hit to the point where an appropriately-set threshold
%   could ignore it.  Using FFT_length = 1024 and Hop_length between 441
%   and 661 may be a better approach.  This produces more false onsets with
%   higher amplitudes for each instrument, but these onsets have very fast
%   decays compared to the real hits.  This could be accounted for in our
%   onset thresholding approach (ex: see how long it takes for the
%   amplitude to fall below a certain level after it rises above another
%   threshold, if it falls too quickly then it's likely a false onset and
%   can be ignored).
%
%   UPDATE (10/24)
%   Evaluate effect of window type on performance. No significant
%   difference, as reported in the thesis. The Hamming seems ever so
%   marginally better than a rectangular window, but should computational
%   abilities run tight than using a rectangular window will suffice.
%
%   Onset detection thresholding functionality added (Step #6). The process
%   checks the amplitude envelope of each piece of the kit at each frame to
%   see if it is above a certain threshold. If it is, then it finds when
%   that threshold reaches its peak and then evaluates the amplitude of the
%   envelope against a 2nd threshold after a certain number of frames
%   have passed (3 frames for hi-hat and kick, 4 for snare). If the
%   amplitude is below the 2nd threshold, then it is dismissed as a false
%   onset. Otherwise, it is treated as a correctly detected onset and
%   correctly formatted data is printed to the transcription.csv file. This
%   method currently fixes our issues with false onsets occuring with
%   higher amplitude than correct onsets and with faster decays than
%   correct onsets.  That said, the 2 thresholds and the decay time (number
%   of frames that pass before comparing to threshold #2) have been set
%   to work specifically in this case based on empirical observations.
%   Extra consideration needs to be given in order to find a way to set
%   these parameters for more general performance.
%
%   UPDATE (10/26)
%   SNR Evaluation added. Assesses the % of correct onsets detected and the
%   number of false onsets detected for each piece of the kit.  No
%   significant change in results for different NNMF approaches (step #4
%   and #5).  Avoid adding noise to training data.  Thresholds tuned for
%   improved SNR performance.  SNR Evaluation slightly naive since it
%   compares the number of onsets detected, not when they were detected
%   (this is due to noise causing the peaks to occasionally shift by 1
%   frame, which is inconsequential with regards to time resolution).  For
%   some reason the false Hi-Hat onset and % correct Snare onsets do not
%   conform to expectations of performance as SNR decreases, but it does
%   not seem to be a major issue.  Overall performance is strong and, given
%   strategic mic placement close to the kit and tuning of the thresholding
%   parameters (as well as other parameters), noise should not be a major
%   concern (for reference, 0 dB SNR means the mic would be detecting
%   external room noise to be on average as loud as the drum kit - this is
%   highly unlikely since the mic's diaphragm will be pointed appropriately
%   at the kit, the mic's polar response leads to sounds that are
%   off-center to be relatively attenuated, loud sounds shouldn't be as
%   close to the mic as the drum kit is anyways, and sound intensity
%   follows the Inverse Square Law.

%-----------Step #1: Read in drum training samples and test beat-----------

%[y_b, fs] = audioread("recordings/test_beat_lab.wav", 'native');
[y_h, fs] = audioread("recordings/hihat.wav",  'native');
[y_k, fs] = audioread("recordings/kick.wav",  'native');
[y_s, fs] = audioread("recordings/snare.wav",  'native');
[y_t, fs] = audioread("recordings/tom.wav", 'native');

%y_b = cast(y_b, 'double');
y_h = cast(y_h, 'double');
y_k = cast(y_k, 'double');
y_s = cast(y_s, 'double');
y_t = cast(y_t, 'double');

% BELOW: For SNR Evaluation
% [noise,fs] = audioread('recordings/ross_2.wav');
% h_missed_perc = zeros(1,27);
% k_missed_perc = zeros(1,27);
% s_missed_perc = zeros(1,27);
% h_false_num = zeros(1,27);
% k_false_num = zeros(1,27);
% s_false_num = zeros(1,27);
% snrs = zeros(1,27);
% step = sqrt(10^.2);
% noise = noise * (step^-16);
% Above: For SNR Evaluation

%--------------------------------End Step #1-------------------------------

%----------Step #2: Take the STFT of each sample and the test beat---------
FFT_length = 2048;
Hop_length = FFT_length / 4;
wlen = 1024;

win = hamming(wlen);
[totspec_h, totf_h, t] = stft(y_h, win, Hop_length, FFT_length, fs);

[totspec_k, totf_k, t] = stft(y_k, win, Hop_length, FFT_length, fs);

[totspec_s, totf_s, t] = stft(y_s, win, Hop_length, FFT_length, fs);

[totspec_t, totf_t, t] = stft(y_t, win, Hop_length, FFT_length, fs);

half_f = ((length(totf_h)-1)/2); % Middle of spectrogram (Hz = 0)

% we only want the positive half of the spectrogram
prespec_h = totspec_h(half_f+1:end, :);
prespec_k = totspec_k(half_f+1:end, :);
prespec_s = totspec_s(half_f+1:end, :);
prespec_t = totspec_t(half_f+1:end, :);

% we need power representation, not magnitude
spec_h = abs(prespec_h);
spec_k = abs(prespec_k);
spec_s = abs(prespec_s);
spec_t = abs(prespec_t);

f = totf_h(half_f:end); % Frequencies at which the STFTs are evaluated
%--------------------------------End Step #2-------------------------------

%--------Step #3: Get the band-wise sums of the pieces' spectrograms-------
% Define the frequency band partitioning
arrayOfFrequencies = [100; 200; 300; 400; 510; 630; 770; 920; 1080; 1270;
                        1480; 1720; 2000; 2320; 2700; 3150; 3700; 4400;
                        5300; 6400; 7700; 9500; 12000; 15500; 22050];
% arrayOfFrequencies = [100; 200; 300; 400; 500; 1000; 5000; 10000; 22050];

% Initialize frequency band matricies for each piece of the kit
specband_h = zeros(length(arrayOfFrequencies), size(spec_h, 2));
specband_k = zeros(length(arrayOfFrequencies), size(spec_k, 2));
specband_s = zeros(length(arrayOfFrequencies), size(spec_s, 2));
specband_t = zeros(length(arrayOfFrequencies), size(spec_t, 2));

% spec_h = ones(size(prespec_h));
% Perform frequency band summation
for samp = 1:half_f
    % Hi-hat
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_h(i, :) = specband_h(i, :) + spec_h(samp, :);
            break
        end
    end

    % Kick
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_k(i, :) = specband_k(i, :) + spec_k(samp, :);
            break
        end
    end

    % Snare
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_s(i, :) = specband_s(i, :) + spec_s(samp, :);
            break
        end
    end
    
    % Snare
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_t(i, :) = specband_t(i, :) + spec_t(samp, :);
            break
        end
    end
end
%--------------------------------End Step #3-------------------------------

%--------Step #4: Calculate basis matrix B for each piece of the kit-------
% Approach 4A: Take the average of each spectrograms' Hz bands over time
    B_h = mean(specband_h.').';
    B_k = mean(specband_k.').';
    B_s = mean(specband_s.').';
    B_t = mean(specband_t.').';
% End Approach 4A

% Approach 4B: Compute basis matrix using NNMF
%     k = 1; % Rank of NNMF, we only need to use k=1
%     % Initialize Bs and Gs to ones
%     B_h = ones(length(arrayOfFrequencies), k);
%     G_h = ones(k, size(specband_h, 2));
%     B_k = ones(length(arrayOfFrequencies), k);
%     G_k = ones(k, size(specband_k, 2));
%     B_s = ones(length(arrayOfFrequencies), k);
%     G_s = ones(k, size(specband_s, 2));
%     % Ones matricies of same size as training sample spectrograms
%     ones_h = ones(size(specband_h));
%     ones_k = ones(size(specband_k));
%     ones_s = ones(size(specband_s));
%
%     % Approach 4B-1: Do NNMF using a set number of iterations
%     num_iters = 100;
%     % Hi-Hat
%     for iter = 1:num_iters
%         B_h = B_h.*( ((specband_h./(B_h*G_h))*G_h.')./ (ones_h*G_h.') );
%         G_h = G_h.*((B_h.'*(specband_h./(B_h*G_h)))./(B_h.'*ones_h));
%     end
%
%     % Kick
%     for iter = 1:num_iters
%         B_k = B_k.*( ((specband_k./(B_k*G_k))*G_k')./(ones_k*G_k') );
%         G_k = G_k.*((B_k'*(specband_k./(B_k*G_k)))./(B_k'*ones_k));
%     end
%
%     % Snare
%     for iter = 1:num_iters
%         B_s = B_s.*( ((specband_s./(B_s*G_s))*G_s')./(ones_s*G_s') );
%         G_s = G_s.*((B_s'*(specband_s./(B_s*G_s)))./(B_s'*ones_s));
%     end
% End Approach 4B-1
%
% % Approach 4B-2: Do NNMF using cost threshold evaluation
% cost_threshold = 1;
% % Hi-Hat
% div_h = 1000; %The change in the value returned by the cost function
% cost_h_old = 0; % The previous value returned by the cost function
% while abs(div_h) > cost_threshold
%     B_h = B_h.*( ((specband_h./(B_h*G_h))*G_h.')./ (ones_h*G_h.') );
%     G_h = G_h.*((B_h.'*(specband_h./(B_h*G_h)))./(B_h.'*ones_h));
%     approx_h = B_h*G_h;
%     cost_h_new = sum(sum(specband_h.*log10(specband_h./approx_h) - ...
%                          specband_h + approx_h));
%     div_h = cost_h_new - cost_h_old;
%     cost_h_old = cost_h_new;
% end
%
% % Kick
% div_k = 1000; %The change in the value returned by the cost function
% cost_k_old = 0; % The previous value returned by the cost function;
% while abs(div_k) > cost_threshold
%     B_k = B_k.*( ((specband_k./(B_k*G_k))*G_k')./(ones_k*G_k') );
%     G_k = G_k.*((B_k'*(specband_k./(B_k*G_k)))./(B_k'*ones_k));
%     approx_k = B_k*G_k;
%     cost_k_new = sum(sum(specband_k.*log10(specband_k./approx_k) - ...
%                          specband_k + approx_k));
%     div_k = cost_k_new - cost_k_old;
%     cost_k_old = cost_k_new;
% end
%
% % Snare
% div_s = 1000; %The change in the value returned by the cost function
% cost_s_old = 0; % The previous value returned by the cost function
% while abs(div_s) > cost_threshold
%     B_s = B_s.*( ((specband_s./(B_s*G_s))*G_s')./(ones_s*G_s') );
%     G_s = G_s.*((B_s'*(specband_s./(B_s*G_s)))./(B_s'*ones_s));
%     approx_s = B_s*G_s;
%     cost_s_new = sum(sum(specband_s.*log10(specband_s./approx_s) - ...
%                          specband_s + approx_s));
%     div_s = cost_s_new - cost_s_old;
%     cost_s_old = cost_s_new;
% end
% % End Approach 4B-2

% concatenate basis matricies
B_fixed = [B_h B_k B_s B_t];

fileID = fopen('b_fixed.h','w');
    fprintf(fileID,'#include <stdint.h>\n\n');
    fprintf(fileID,'float32_t b_fixed[][4] = {\n');
    for j=1:size(B_fixed, 1)
        fprintf(fileID, '{');
        for i=1:size(B_fixed, 2)
            fprintf(fileID,'%f,',B_fixed(j,i));
        end
        if j == 25
            fprintf(fileID, '}\n');
        else
            fprintf(fileID, '},\n');
        end
    end
    fprintf(fileID,'};');
    fclose(fileID);

% [y, fs] = audioread("recordings/click2.wav", 'native');
% % [y_real, fs] = audioread("recordings/click2.wav");
% y = cast(y, 'int32');
% fileID = fopen('click.h','w');
% fprintf(fileID,'#include <stdint.h>\n\n');
% fprintf(fileID,'int32_t click[11025] = ');
% fprintf(fileID, '{');
% for i=1:size(y, 1)
%     fprintf(fileID,'%i', y(i,1));
%     if i == 11025
%         fprintf(fileID, '};');
%     else
%         fprintf(fileID, ', ');
%     end
% end
% fclose(fileID);
%
% sound(y(:,1), fs);
