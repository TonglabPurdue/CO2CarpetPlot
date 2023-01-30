



close all;
clear;
clc;
   
%% Stuff that won't appear in final code, just necessary for testing

load('test_fmri_data.mat');
imag_funcA = image_func_full(:,:,:,1:510);
imag_funcB = image_func_full;

%% Enter parameters
TR=1; %TR of fMRI scan
N0=1; %Starting time point (for fMRI time series cropping, if different length than end-tidal CO2 time series)
N=510; %Ending time point (for fMRI time series cropping)

%% Compute global mean (averaged over whole brain) time series

func_mat = zeros(size(imag_funcA,1)*size(imag_funcA,2)*size(imag_funcA,3), size(imag_funcA,4));
for i=1:size(imag_funcA,4)
    func_mat(:,i) = reshape(imag_funcA(:,:,:,i), [], 1);
end
GMean_brain = nanmean(func_mat);
    
%% Compute CO2 Delay Map

ratio=0.1; %Controls ratio for oversampling time series
lp=0.001;  %Lower frequency bound for time series bandpass filter
hp=0.02;   %Upper frequency bound for time series bandpass filter
ref_ts=GMean_brain';   %Time series to be used for cross-correlation with voxel time series
imag_data=imag_funcA;  %fMRI time series data
window=30; %Cross-correlation window width
num=2;     %Order of Butterworth filter

delay_map = Compute_Delay_Map(imag_data,ref_ts,lp, hp, TR,ratio,window,num, N0, N);
   
%% Create sorted carpet plot

%Vectorize delay map
dmap_vec = reshape(delay_map, [], 1);

%Sort fMRI time series according to voxel delays
func_mat = reshape(imag_funcB, [], size(imag_funcB,4));
func_mat_ex = [dmap_vec func_mat];
ind_valid = find(func_mat(:,100)>0); %Make sure there are no time series of all 0s
func_mat_ex = func_mat_ex(ind_valid, :);
sorted_mat_ex = sortrows(func_mat_ex, 1, 'descend');
sorted_mat = sorted_mat_ex(:, 2:end);

%Demean and bandpass filter the time series (0.001-0.02 Hz)
sorted_mat_filt = zeros(size(sorted_mat));
sorted_mat_dt = sorted_mat_filt;
parfor i=1:size(sorted_mat,1)
    sorted_mat_dt(i,:) = demean(sorted_mat(i,:));
    sorted_mat_filt(i,:) = filtf(double(sorted_mat_dt(i,:)), 0.001, 0.02, 1/TR, 2);
end

%Isolate middle carpet plot section
dm_vec = sorted_mat_ex(:,1);
middle_val = median(dm_vec);
start_point = nnz(dm_vec>(middle_val+10));
stop_point = length(dm_vec) - nnz(dm_vec<(middle_val-10));
perc_cropped = (stop_point-start_point)/length(dm_vec);
disp(strcat('Percent of voxels included in middle carpet plot section is: ', perc_cropped));
im1 = sorted_mat_filt(start_point:stop_point, N0:N);
    
    
%% Find CO2 Edge

%Choose number of lines to look for
num_lines = 1;

%Choose to calc on rising (black to white) or falling (white to black) edges
%-1 = rising edge, 1 = falling edge
edge = -1;

%Set contrast thereshold for limiting accepted edges
contrast_threshold = 0.2;

%% Smooth image

scale_factor = 1; %Used to scale computation of edge angles if desired
height = size(im1,1);

orig_avg = mean(im1);

sort_avg = sort(orig_avg, 'descend');

%Create filter h for image blurring - meant to smooth out image for
%clearer edges
h = 1/12*[1 1.5 1;
    1.5 2 1.5;
    1 1.5 1];
h = conv2(h,h);

%Filter data for several repititions
smoothed_data = filter2(h, im1);
for i=1:5
smoothed_data = filter2(h, smoothed_data);
end
old_im1 = im1;
im1 = smoothed_data;
smoothed_avg = mean(im1(:,:));

%% Derivative
%Apply derivative filter
h = 1/2*[1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_avg = edge * filter2(h, smoothed_avg);
%Adjust the ends of the derivatived data since the front and back will be
%skewed
for i=1:filtsize
    derivatived_avg(filtsize-(i-1)) = derivatived_avg(filtsize-(i-2));
    derivatived_avg(size(derivatived_avg,2)-filtsize+(i)) = derivatived_avg(size(derivatived_avg,2)-filtsize+(i-1));
end

%% Find Peaks
max_ind = zeros(1,num_lines);
temp = derivatived_avg(:,1:300)-0.0003;
% figure
% plot(temp);

%Here we find the n=num_lines highest peaks of the derivative data,
%representing locations where we want to draw a line
for i=1:num_lines
    is_line_on_back_edge = 1;
    is_line_on_front_edge = 1;
    while is_line_on_back_edge == 1 || is_line_on_front_edge == 1
        [~, max_ind(i)] = max(temp);
        a=0;
        while (max_ind(i) + a <= size(im1,2)) && (temp(max_ind(i) + a) > 0)
            temp(max_ind(i) + a) = 0;
            a= a+1;
        end
        if max_ind(i) + a <= size(im1,2)
            is_line_on_back_edge = 0;
        else
            is_line_on_back_edge = 1;
        end
        a=-1;
        while (max_ind(i) + a >= 1) && (temp(max_ind(i) + a) > 0)
            temp(max_ind(i) + a) = 0;
            a= a-1;
        end
        if max_ind(i) + a >= 1
            is_line_on_front_edge = 0;
        else
            is_line_on_front_edge = 1;
        end
    end
    if nnz(temp<=0) == length(temp)
        disp('Not enough locations for a line!');
        num_lines = i
        break;   
    end
end

neur_assignment = zeros(1, num_lines);

max_ind = sort(max_ind(max_ind>0), 'ascend');
%Plot derivative over data to illustrate where the algorithm is deciding to
%draw lines
% figure
% imagesc(im1);
% hold on;
% colormap(gray);
% caxis([-1 1]);
% plot(derivatived_avg*20000, 'r', 'LineWidth', 2);
% hold off;

%% Draw Lines
%Initialize space to store line data
linfit = zeros(num_lines, 2);
linfitline = zeros(num_lines, height-1+1);

%Apply derivative filter to all data
h = 1/8*[1 0 -1; 2 0 -2; 1 0 -1];
filtsize = floor(size(h,2) / 2);
derivatived_data = edge * filter2(h, im1);
for b=1:filtsize
    derivatived_data(:,filtsize-(b-1)) = derivatived_data(:,filtsize-(b-2));
    derivatived_data(:,size(derivatived_data,2)-filtsize+(b)) = derivatived_data(:,size(derivatived_data,2)-filtsize+(b-1));
end

temp = derivatived_avg - 0.0003;
avg_contrast = zeros(1, num_lines);
for i=1:num_lines
    location = max_ind(i);
    fin(i)=0;
    %For now we are considering the range around the peak where the
    %derivative is greater than 0, plus adding an extra 2 points on either
    %side
    while (max_ind(i) + fin(i) <= size(im1,2)) && (temp(max_ind(i) + fin(i)) > 0)
        temp(max_ind(i) + fin(i)) = 0;
        fin(i)= fin(i)+1;
    end
    st(i)=-1;
    while (max_ind(i) + st(i) >= 1) && (temp(max_ind(i) + st(i)) > 0)
        temp(max_ind(i) + st(i)) = 0;
        st(i)= st(i)-1;
    end
    
    %Compute averaged contrast
    cont_st = smoothed_avg(st(i)+max_ind(i));
    cont_end = smoothed_avg(fin(i)+max_ind(i));
    avg_contrast(i) = abs(cont_end - cont_st);
    
    
    fin(i) = fin(i)+ 2;
    st(i) = st(i) - 2;
    
    st(i)=st(i)+max_ind(i);
    if st(i) < 1
        st(i) = 1;
    end
    fin(i)=fin(i)+max_ind(i);
    if fin(i) > size(im1,2)
        fin(i) = size(im1,2);
    end
    data = derivatived_data(1:height, st(i):fin(i));
    %We can display the considered data if we want, but it's commented out
    %as seen below
%     figure
%     imshow(data);
%     colormap(gray);    
    
    %Find horizontal point where derivative is maximized for each row
    max_ind2 = zeros(1, size(data, 1));
    for b=1:size(data,1)
        [~, max_ind2(b)] = max(data(b, :));
    end
    %Account for extreme values that will probably throw off algorithm
    max_ind2(length(max_ind2)) = max_ind2(length(max_ind2)-1);
    max_ind2(1) = max_ind2(2);
    
    %Calculate the best fit line for our derivative peak locations
    datayax = 1:size(max_ind2,2);
%     figure
    %plot(flip(max_ind2), datayax);
%     scatter(flip(max_ind2), datayax);
%     hold on;
    linfit(i,:) = polyfit(datayax, max_ind2, 1);
    linfitline(i,:) = linfit(i,1)*datayax + linfit(i,2);
%     plot(flip(linfitline), datayax);
%     hold off;
%     title('Finding Best Fit Line')
    
end

%Calculate slopes and angles
slopes = -1./linfit(:,1);
angles = atan(slopes * scale_factor)*180/pi;
times = ones(size(slopes)) * size(im1,1) ./ slopes;

%% Main Fig 1: Create figure with first CO2 edge and transit time for this subject

transit_time = times*TR;
figure
imagesc(sorted_mat_dt(:, 2:end))
colormap(gray);
caxis([-1 1]);
name = strcat('Subject', {', '}, 'Transit Time =', {' '}, string(transit_time));
title(name);
ylabel('Voxels');
xlabel('Time (TRs)');
hold on;
yline(start_point, 'b', 'Linewidth', 3)
yline(stop_point, 'b', 'Linewidth', 3)
hold on;
plot_mat = [];
for i=1%:num_lines
    if avg_contrast(i) > contrast_threshold
        if slopes(i) < 0
            plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax-10, 'g', 'LineWidth', 2);
        else
            plot((linfitline(i,:)-mean(linfitline(i,:))+max_ind(i)), datayax+start_point, 'r', 'LineWidth', 2);
        end
        plot_mat = [plot_mat; [times(i)*TR max_ind(i)]];
    end
end

%%
function   delay_map = Compute_Delay_Map(imag_data,ref_ts,lp, hp, TR,ratio,window,num, N0, N)
% This function is to generate delay maps
% by Jinxia Yao on 21Feb2020, compiled for CO2 carpet plot paper on
% 29Dec2022 by Bradley Fitzgerald

    ratioNtoOld=floor(TR/ratio); %upsample ratio; TR/ratioNtoOld is the new temporal resolution
    xcorr_range=floor(window/(TR/ratioNtoOld));

    %Demean and filter reference (CO2) time series
    ref_ts = ref_ts(N0:N);
    ref_ts=demean(oversample_ts(double(ref_ts),ratioNtoOld));
    ref_ts=filtf(ref_ts,lp,hp,1/(TR/ratioNtoOld),num);

    %Set up filtered time series
    AA1=size(imag_data);
    Total=AA1(1)*AA1(2)*AA1(3); 
    ts_raw = double(reshape(imag_data,[],AA1(4)));  %each row is a timeseries
    ind_valid = find(ts_raw(:,100)>0);
    ts_detrend=demean(oversample_ts(double(ts_raw(ind_valid,:)'),ratioNtoOld));  % #timepoint * #voxel
    ts_filtered = filtf(ts_detrend,lp,hp,1/(TR/ratioNtoOld),num) ; % #voxel *#timepoint

    %Compute delays for each voxel time series
    coef=[];delay=[];
    bg_1d_2_3d = zeros(Total,1);
    delay_whole = bg_1d_2_3d;
    for w=1:length(ind_valid)   
            [r,p]=xcorr(ts_filtered(:,w),ref_ts,xcorr_range,'coeff');
            [cc,ind]=mmax(r);   % cc is the maximum cc 
             T=p(ind)*TR/ratioNtoOld;
            coef = [coef cc];
            delay = [delay T];
    end
    delay_whole(ind_valid) = delay;
    delay_map = reshape(delay_whole,AA1(1),AA1(2),AA1(3),1);
end

%%
function y=filtf(input_signal,lpass, hpass,freq,num)
    % modified by Jinxia Yao on 02/21/2020.
    % input is a vector / matrix (calculated by columns)
    L=length(input_signal(:,1));
    N=floor(L/10);
    signal_padding = padarray(input_signal,N,'symmetric','both'); %padding 
    [b,a] = butter(num, [lpass*2/freq, hpass*2/freq], 'bandpass');
    y_padding = filtfilt(b, a, signal_padding);
    y=y_padding(N+1:N+L,:);
end

%%
function [demean_ts]=demean(ts)
    % change to the following code in 20200220 by Jinxia Yao;
    % input can be a vector/matrix;calculated by columns;
    detrend_ts = detrend(ts);  % by column
    demean_ts = detrend_ts./std(detrend_ts);
end

%%
%This one was pulled from NIFTI folder in Functions, make sure it is
%properly attributed
function TS=oversample_ts(ts,ratioNtoOld)
%% oversample
    %ratioNtoOld=5 mean oversample 5 times
    [m,n]=size(ts);
    kk=round(m*ratioNtoOld);
    TS=zeros(kk,n);
    for i=1:n
        TS(:,i)=interpft (ts(:,i), kk);
    end 
end

%%
function [y,ind]=mmax(A)
    % Modified by Jinxia Yao on 02/21/2020
    % return the value with the max absolute value and its corresponding index.
    [y1,ind] = max(abs(A));
    y=A(ind);
end