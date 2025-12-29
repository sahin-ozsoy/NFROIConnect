function [EEG, source_voxel_data] = roi_activityOnlyActivityFeedback(EEG,fReset,roiIndices)
% roi_activityOnlyActivityFeedback - compute connectivity between ROIs
%  modified from EEGLAB roi_activity to remove the welch function.
%  optimized for speed for NFB


% Cortex mesh or volume
% ---------------------
persistent cortex;
persistent leadfield;
persistent P_eloreta;


if fReset
    s = load('nfbcortex.mat');
    cortex = s.cortex;   

    % leadfield matrix (fieldtrip)
    % ------------------------------------------
    leadfield = EEG.dipfit.sourcemodel;
    leadfield.gain = reshape( [ leadfield.leadfield{:} ], [length(leadfield.label) 3 length(leadfield.leadfield)]);
    leadfield.gain = permute(leadfield.gain, [1 3 2]);
    leadfield = leadfield.gain;
end

nvox = size(cortex.Vertices, 1);

% common average reference transform
nbchan = EEG.nbchan;
% H = eye(nbchan) - ones(nbchan) ./ nbchan;

% apply to data and leadfield
% tmpData = reshape(H*EEG.data(g.chansel, :), nbchan, EEG.pnts, EEG.trials);
% leadfield = reshape(H*leadfield(:, :), nbchan, nvox, 3);
tmpData = EEG.data;

%% source reconstruction
if fReset
    C = cov(tmpData(:, :)');
    alpha = 0.05*trace(C)/length(C);
    Cr = C + alpha*eye(nbchan);
    P_eloreta = lcmvFast(Cr, leadfield,fReset,struct('alpha', 0, 'onedim', 0));
end
% [~, P_eloreta] = lcmv(Cr, leadfield,struct('alpha', 0, 'onedim', 0));
source_voxel_data = 10^3 .* reshape(tmpData(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
% source_voxel_data = 10^3*source_voxel_data; % the units are nA*m
    
% number of ROIs in the Desikan-Killiany Atlas
nROI = length(cortex.Atlas.Scouts);

% Keep only the first strongest component for each ROI   
source_roi_data = zeros(size(source_voxel_data,1), nROI);    

% Calculate only required rois
for iROI = roiIndices        
    ind_roi = cortex.Atlas.Scouts(iROI).Vertices;
    [source_roi_data_tmp, ~] = roi_getactFeedback(source_voxel_data, ind_roi, 1);
    source_roi_data(:, iROI) = source_roi_data_tmp;        
end

source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);

% Output paramters
EEG.roi.source_roi_data = single(source_roi_data);
EEG.roi.srate = EEG.srate; % add srate to the roi

