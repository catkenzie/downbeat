
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

[y_b, fs] = audioread("recordings/beat.wav", 'native');
[y_h, fs] = audioread("recordings/hihat.wav",  'native');
[y_k, fs] = audioread("recordings/kick.wav",  'native');
[y_s, fs] = audioread("recordings/snare.wav",  'native');
[y_t, fs] = audioread("recordings/tom.wav", 'native');
[y_c, fs] = audioread("recordings/crash.wav", 'native');

y_b = cast(y_b, 'double');
y_h = cast(y_h, 'double');
y_k = cast(y_k, 'double');
y_s = cast(y_s, 'double');
y_t = cast(y_t, 'double');
y_c = cast(y_c, 'double');

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
FFT_length = 1024;
Hop_length = FFT_length / 2;
[totspec_h, totf_h] = stft(y_h, fs, 'Window', hamming(FFT_length), ...
                           'OverlapLength', Hop_length, ...
                           'FFTLength', FFT_length);
[totspec_k, totf_k] = stft(y_k, fs, 'Window', hamming(FFT_length), ...
                           'OverlapLength', Hop_length,  ...
                           'FFTLength', FFT_length);
[totspec_s, totf_s] = stft(y_s, fs, 'Window', hamming(FFT_length), ...
                           'OverlapLength', Hop_length,  ...
                           'FFTLength', FFT_length);
[totspec_t, totf_t] = stft(y_t, fs, 'Window', hamming(FFT_length), ...
                           'OverlapLength', Hop_length,  ...
                           'FFTLength', FFT_length);
[totspec_c, totf_c] = stft(y_c, fs, 'Window', hamming(FFT_length), ...
                           'OverlapLength', Hop_length,  ...
                           'FFTLength', FFT_length);

half_f = (length(totf_h)/2); % Middle of spectrogram (Hz = 0)

% we only want the positive half of the spectrogram
prespec_h = totspec_h(half_f+1:end, :);
prespec_k = totspec_k(half_f+1:end, :);
prespec_s = totspec_s(half_f+1:end, :);
prespec_t = totspec_t(half_f+1:end, :);
prespec_c = totspec_c(half_f+1:end, :);

% we need power representation, not magnitude
spec_h = abs(prespec_h);
spec_k = abs(prespec_k);
spec_s = abs(prespec_s);
spec_t = abs(prespec_t);
spec_c = abs(prespec_c);

f = totf_h(half_f:end); % Frequencies at which the STFTs are evaluated
%--------------------------------End Step #2-------------------------------

%--------Step #3: Get the band-wise sums of the pieces' spectrograms-------
% Define the frequency band partitioning
% arrayOfFrequencies = [100; 200; 300; 400; 510; 630; 770; 920; 1080; 1270;
%                         1480; 1720; 2000; 2320; 2700; 3150; 3700; 4400;
%                         5300; 6400; 7700; 9500; 12000; 15500; 22050];
arrayOfFrequencies = [44; 88; 132; 176; 220; 264; 308; 352; 396; 440;
                        510; 630; 770; 920; 1080; 1380; 1740; 2580;
                        4250; 6400; 7700; 9500; 12000; 15500; 22050];
% arrayOfFrequencies = [100; 200; 300; 400; 500; 1000; 5000; 10000; 22050];

% Initialize frequency band matricies for each piece of the kit
specband_h = zeros(length(arrayOfFrequencies), size(spec_h, 2));
specband_k = zeros(length(arrayOfFrequencies), size(spec_k, 2));
specband_s = zeros(length(arrayOfFrequencies), size(spec_s, 2));
specband_t = zeros(length(arrayOfFrequencies), size(spec_t, 2));
specband_c = zeros(length(arrayOfFrequencies), size(spec_c, 2));

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
    
    % Tom
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_t(i, :) = specband_t(i, :) + spec_t(samp, :);
            break
        end
    end
    
    % Crash
    for i = 1:length(arrayOfFrequencies)
        if f(samp) < arrayOfFrequencies(i)
            specband_c(i, :) = specband_c(i, :) + spec_c(samp, :);
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
    B_c = mean(specband_c.').';
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
B_fixed = [B_h B_k B_s B_t B_c];
%--------------------------------End Step #4-------------------------------

%-----Step #5: Onset detection of input signal using fixed basis matrix----

% BELOW: For SNR Evaluation
% nSim = 1000;
% for snr_step = -15:11
%     % Necessary code for SNR Plots (Step #7)
%     noise = noise * step;
%     noise_wrap = cat(1, noise, noise);
%     pSignal = rms(y_b)^2;
%     pNoise = rms(noise)^2;
%     ston = 10*log10(pSignal/pNoise);
%     h_missed = 0;
%     h_false = 0;
%     k_missed = 0;
%     k_false = 0;
%     s_missed = 0;
%     s_false = 0;
%     fprintf('Signal-to-noise ratio: %f\n',ston);
%     snrs(snr_step+16) = ston;
%     for iSim = 1:nSim
%         offset = randi(length(noise));
%         noise_seg = noise_wrap(offset:offset+size(y_b)-1);
% ABOVE: For SNR Evaluation

        noise_seg = zeros(length(y_b), 1); % For running w/out added noise

        [totspec_b, totf_b] = stft(y_b + noise_seg, fs, ...
                               'Window', hamming(FFT_length), ...
                               'OverlapLength', Hop_length, ...
                               'FFTLength', FFT_length);
        prespec_b = totspec_b(half_f+1:end, :);
        spec_b = abs(prespec_b);

        % Initialize the amplitude envelopes of each piece of the kit
        Amp_env_h = zeros(1, size(spec_b, 2));
        Amp_env_k = zeros(1, size(spec_b, 2));
        Amp_env_s = zeros(1, size(spec_b, 2));
        Amp_env_t = zeros(1, size(spec_b, 2));
        Amp_env_c = zeros(1, size(spec_b, 2));

        num_nnmf_iters = 100; % For approach 5A
        nnmf_cost_threshold = 0.001; % For approach 5B

        % Evaluate the input signal one frame at a time
        for frame = 1:size(spec_b, 2)
            % Initialize frequency band matrix for the input signal
            specband_b = zeros(length(arrayOfFrequencies), 1);
            % Perform frequency band summation
            for samp = 1:half_f
                for i = 1:length(arrayOfFrequencies)
                    if f(samp) < arrayOfFrequencies(i)
                        specband_b(i) = specband_b(i) + spec_b(samp, frame);
                        break
                    end
                end
            end

            % Initialize the gain matrix
            G_b = ones(5, 1);
            % Ones matrix of same size as input signal spectrogram
            ones_b = ones(size(specband_b));

            % Approach 5A: Do NNMF using a set number of iterations
            for i = 1:num_nnmf_iters
                G_b = G_b.*((B_fixed'*(specband_b./(B_fixed*G_b)))./(B_fixed'*ones_b));
            end
            % End Approach 5A

        %     % Approach 5B: Do NNMF using cost threshold evaluation
%             div_b = 1000;
%             cost_b_old = 0;
%             while abs(div_b) > .001
%                 G_b = G_b.*((B_fixed'*(specband_b./(B_fixed*G_b)))./(B_fixed'*ones_b));
%                 approx_b = B_fixed*G_b;
%                 cost_b_new = sum(sum(specband_b.*log10(specband_b./approx_b) - ...
%                                      specband_b + approx_b));
%                 div_b = cost_b_new - cost_b_old;
%                 cost_b_old = cost_b_new;
%             end
        %     % End Approach 5B

            % Save gain values of each peice of the kit at each frame
            Amp_env_h(frame) = G_b(1);
            Amp_env_k(frame) = G_b(2);
            Amp_env_s(frame) = G_b(3);
            Amp_env_t(frame) = G_b(4);
            Amp_env_c(frame) = G_b(5);
        end
%--------------------------------End Step #5-------------------------------

        % Plot the superimposed amplitude envelopes
%          for j = 1:2
%             for i = 2:length(Amp_env_k) - 1
%                 if Amp_env_k(i-1) > Amp_env_k(i) && ...
%                     Amp_env_k(i+1) > Amp_env_k(i)
%                     Amp_env_k(i) = (Amp_env_k(i-1) + Amp_env_k(i+1))/2;
%                 end
%             end
%          end
        figure;
        plot(Amp_env_k);
        axis tight;
        
        dog_k = diff(log(Amp_env_k));
        dog_k(dog_k < 0) = 0;
        figure;
        plot([0 dog_k]);
        axis tight;
        
%         figure;
%         plot(log(Amp_env_t));
%         axis tight;
%         hold on
%         plot(Amp_env_k);
%         axis tight;
%         plot(Amp_env_s);
%         axis tight;
%         figure
%         plot([0 diff(Amp_env_s)])
%         hold on
%         plot(zeros(1,length(Amp_env_s)))
%         axis tight;
%         skews = zeros(1, length(Amp_env_h));
%         for i = 3:length(derv) - 2
%             skews(i) = 3*(mean(Amp_env_h(i-2:i+2)) - ...
%                             median(Amp_env_h(i-2:i+2))) / ...
%                             std(Amp_env_h(i-2:i+2));
%         end

%         figure;
%         plot([0 diff(Amp_env_h)]);
%         axis tight;
%         hold on
%         plot(diff(Amp_env_k));
%         axis tight;
%         plot(diff(Amp_env_s));
%         axis tight;

% %--------------------Step #6: Onset Detection Thresholding-----------------
%
%         h_detected = false;
%         h_onset_frame = 0;
%         h_onset_level = 0;
%         k_detected = false;
%         k_onset_frame = 0;
%         k_onset_level = 0;
%         s_detected = false;
%         s_onset_frame = 0;
%         s_onset_level = 0;
%
%         % BELOW: For Producing CSV Output
%         fileID = fopen('transcription.csv','w');
%         fprintf(fileID, 'instrument,time(s)\n');
%         % Above: For Producing CSV Output
%
%         % BELOW: For SNR Evaluation
% %         h_onset_frames_noise = zeros(1, size(spec_b, 2));
% %         k_onset_frames_noise = zeros(1, size(spec_b, 2));
% %         s_onset_frames_noise = zeros(1, size(spec_b, 2));
%         % BELOW: For SNR Evaluation
%
%         for frame = 1:size(spec_b, 2)
%             % For Hi-Hat
%             if Amp_env_h(frame) > 2
%                 h_detected = true;
%                 if Amp_env_h(frame) > h_onset_level
%                     h_onset_frame = frame;
%                     h_onset_level = Amp_env_h(frame);
%                 end
%             end
%             if h_detected && frame - h_onset_frame == 4
%                 if Amp_env_h(frame) > 0.002
%                     % BELOW: For Producing CSV Output
%                     h_onset_time = ((h_onset_frame - 1)*512)/44100;
%                     fprintf(fileID, 'hi-hat,%u\n',h_onset_time);
%                     % ABOVE: For Producing CSV Output
%
%                     % BELOW: For SNR Evaluation
% %                     h_onset_frames_noise(h_onset_frame) = 1;
%                     % ABOVE: For SNR Evaluation
%                 end
%                 h_detected = false;
%                 h_onset_level = 0;
%             end
%
%             % For Kick
%             if Amp_env_k(frame) > 3.3
%                 k_detected = true;
%                 if Amp_env_k(frame) > k_onset_level
%                     k_onset_frame = frame;
%                     k_onset_level = Amp_env_k(frame);
%                 end
%             end
%             if k_detected && frame - k_onset_frame == 4
%                 if Amp_env_k(frame) > 0.5
%                     % BELOW: For Producing CSV Output
%                     k_onset_time = ((k_onset_frame - 1)*512)/44100;
%                     fprintf(fileID, 'kick,%u\n',k_onset_time);
%                     % ABOVE: For Producing CSV Output
%
%                     % BELOW: For SNR Evaluation
% %                     k_onset_frames_noise(k_onset_frame) = 1;
%                     % ABOVE: For SNR Evaluation
%                 end
%                 k_detected = false;
%                 k_onset_level = 0;
%             end
%
%             % For Snare
%             if Amp_env_s(frame) > 7.7
%                 s_detected = true;
%                 if Amp_env_s(frame) > s_onset_level
%                     s_onset_frame = frame;
%                     s_onset_level = Amp_env_s(frame);
%                 end
%             end
%             if s_detected && frame - s_onset_frame == 4
%                 if Amp_env_s(frame) > 4.3
%                     % BELOW: For Producing CSV Output
%                     s_onset_time = ((s_onset_frame - 1)*512)/44100;
%                     fprintf(fileID, 'snare,%u\n',s_onset_time);
%                     % ABOVE: For Producing CSV Output
%
%                     % BELOW: For SNR Evaluation
% %                     s_onset_frames_noise(s_onset_frame) = 1;
%                     % ABOVE: For SNR Evaluation
%                 end
%                 s_detected = false;
%                 s_onset_level = 0;
%             end
%         end
%--------------------------------End Step #6-------------------------------

% %------------------------Step #7: Develop SNR Plots------------------------
%
%       % For Hi-Hat
%         if sum(h_onset_frames_noise) > 7
%             h_false = h_false + (sum(h_onset_frames_noise) - 7);
%         elseif sum(h_onset_frames_noise) < 7
%             h_missed = h_missed + (7 - sum(h_onset_frames_noise));
%         end
%
%       % For Kick
%         if sum(k_onset_frames_noise) > 5
%             k_false = k_false + (sum(k_onset_frames_noise) - 5);
%         elseif sum(k_onset_frames_noise) < 5
%             k_missed = k_missed + (5 - sum(k_onset_frames_noise));
%         end
%
%       % For Snare
%         if sum(s_onset_frames_noise) > 7
%             s_false = s_false + (sum(s_onset_frames_noise) - 7);
%         elseif sum(s_onset_frames_noise) < 7
%             s_missed = s_missed + (7 - sum(s_onset_frames_noise));
%         end
%     end
%
%     h_missed_perc(snr_step + 16) = 1 - ((h_missed/nSim)/7);
%     k_missed_perc(snr_step + 16) = 1 - ((k_missed/nSim)/5);
%     s_missed_perc(snr_step + 16) = 1 - ((s_missed/nSim)/7);
%     h_false_num(snr_step + 16) = (h_false/nSim);
%     k_false_num(snr_step + 16) = (k_false/nSim);
%     s_false_num(snr_step + 16) = (s_false/nSim);
% end
%
% figure
% plot(snrs, h_missed_perc);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Percent Correct Onset Detection');
% title('Performance of Hi-Hat Classification in Noise');
%
% figure
% plot(snrs, k_missed_perc);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Percent Correct Onset Detection');
% title('Performance of Kick Classification in Noise');
%
% figure
% plot(snrs, s_missed_perc);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Percent Correct Onset Detection');
% title('Performance of Snare Classification in Noise');
%
% figure
% plot(snrs, h_false_num);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Average Number False Onsets');
% title('Performance of Hi-Hat Classification in Noise');
%
% figure
% plot(snrs, k_false_num);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Average Number False Onsets');
% title('Performance of Kick Classification in Noise');
%
% figure
% plot(snrs, s_false_num);
% axis tight
% xlabel('SNR (dB)');
% ylabel('Average Number False Onsets');
% title('Performance of Snare Classification in Noise');
% %--------------------------------End Step #7-------------------------------
%
 fileID = fopen('workspace_452\classification-master.zip_expanded\classification-master\display_fft\src\b_fixed.h','w');
    fprintf(fileID,'#include <stdint.h>\n\n');
    fprintf(fileID,'float32_t b_fixed[][5] = {\n');
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
