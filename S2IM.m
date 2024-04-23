clear all;
close all;


%% Algorithm Parameters;
itr = 10; % The number of iterations
NA = 0.1; % The numerical aperture of the collection point spread function
STEPSIZE = 0.1; % The gradient descent rate
SAVE_RESULTS = 1;

lambda = 0.605; k=2*pi/lambda;  % wavelength (micons) and wave number
F1=16;
F2=75;
F3=100;
F4=100;

mag = (F2/F1)*(F4/F3); % magnification
pscrop = 6.5/(mag); % Pixels size (microns)
NAs = 0.5

resolution = lambda /  NA
resolution_pixel = resolution /pscrop
psf_size = 2*floor(resolution_pixel) +1
psf = fspecial('gaussian',[psf_size,psf_size], resolution_pixel/2.355);


%%
loadfile = 'Input.mat';
load(loadfile);

%%
frames_alligned = (clean_frames +1)/max(clean_frames(:));
clear clean_frames;
average_frame = average_frame;
MIN = min(frames_alligned(:))

xlim = size(frames_alligned,1);
ylim = size(frames_alligned,2);
Nill = size(frames_alligned, 3);
xcut = max(xshift(:))-min(xshift(:))
ycut = max(yshift(:))-min(yshift(:))
xycat = round(max(xcut,ycut))+1
xmax =xycat;
ymax = xycat;

guess = mean(frames_alligned,3);
guess = sqrt(guess);
illumination = padarray( guess,[2*xmax 2*ymax],MIN);

%% Super-resolution algorithm.

STEPSIZE_1 = STEPSIZE;
STEPSIZE_2 = STEPSIZE;

decay_rate = 0.0;

diff_Iobj = zeros(size(average_frame));
diff_speckle = zeros(size(average_frame));
err = zeros(1,itr);
for i = 1:itr
    i

    for jjj = 1:size(frames_alligned,3)
        I_image_current = frames_alligned(:,:,jjj);

        illumination =  fraccircshift(illumination,[xshift(jjj) yshift(jjj)]); 


        I_speckle_current = illumination(1+2*xmax:xlim+2*xmax,1+2*ymax:ylim+2*ymax);
        I_m = I_speckle_current.*guess;
        I_est = conv2(I_m, psf,'same');
        I_est = max(I_est,MIN);
        

        ratrat = I_image_current ./ I_est;
        lograt = log(ratrat);
        I_diff = I_image_current - I_est;
 

        % kl =  1 -lograt -ratrat; 
        kl = -real(lograt);
        kl_conv =conv2(kl, psf,'same');
        

  

        grad_Iobj = I_speckle_current .* kl_conv;
        grad_speckle = guess .* kl_conv; 
        diff_Iobj = decay_rate * diff_Iobj - STEPSIZE_1 * grad_Iobj;
        diff_speckle = - STEPSIZE_2 * grad_speckle;

        guess = guess +diff_Iobj;
        guess = max(guess,0);
        illumination(1+2*xmax:xlim+2*xmax,1+2*ymax:ylim+2*ymax) = illumination(1+2*xmax:xlim+2*xmax,1+2*ymax:ylim+2*ymax)+ diff_speckle;

        illumination =  fraccircshift(illumination,[-xshift(jjj) -yshift(jjj)]); 

        illumination = max(illumination,0);


    end
    % KL = I_diff.* lograt;
    KL = I_est.* log(I_est ./ I_image_current) + I_image_current - I_est;
    err(i) =  (mean(mean(KL)));
    % err(i) =  mean(mean(abs(I_diff)));
end

myresult = guess;
myresult = myresult - min(myresult(:));
myresult = myresult / max(myresult(:));

avv = average_frame;
avv = avv - min(avv(:)); 
avv = avv / max(avv(:));

if SAVE_RESULTS

    one_frame = frames_alligned(:,:,1);
    Filename = strjoin(['result_NA_' string(NA) '_itr_' string(itr) '_STEPSIZE_' string(STEPSIZE) '_.mat'],'')
    save(Filename, 'one_frame', 'guess', 'illumination', 'xshift','yshift', "average_frame")
end




