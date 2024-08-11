% roi_activityOnlyActivity - call roi_connectivity_process to compute
%                                connectivity between ROIs
% Usage:
%  modified from EEGLAB roi_activity to remove the welch function
%
% Inputs:
%  EEG - EEGLAB dataset
%
% Required inputs:
%  'leadfield'   - [string or struct] leadfield structure (Fieldtrip)
%                  or file containing Fieldtrip or Brainstrom leadfield
%                  matrix.
%  'sourcemodel' - [string] source model file also containing Atlas info.
% 
% Optional inputs:
%  'model'            - ['eLoretaFieldtrip'|'lcmvFieldtrip'|'eLoreta'|'lcmv'] distributed
%                       source localization method. Default is 'lcvm'
%                       (Beamforming) from ROIconnect. 'eLoretaFieldtrip' and 
%                       'lcmvFieldtrip' are alternative implementations in
%                       Fieldtrip that should return similar results.
%  'sourcemodel2mni'  - [9x float] homogeneous transformation matrix to convert
%                       sourcemodel to MNI space.
%  'sourcemodelatlas' - [string] name of Atlas to use (must be contained
%                       in Atlas field of the sourcemodel file.
%  'nPCA'             - [interger] Number of PCA component for each ROI. Each ROI
%                       is made of many voxel. Instead of averaging their activity,
%                       this function takes the x first PCA components, then use
%                       these to compute connectivity (default is 3)
%  'naccu'            - [interger] For bootstrap, number of accumulation. Default is 
%                        none.
%  'eloretareg'       - [float] regularization term for eLoreta. Default is 0.05.
%  'roiactivity'      - ['on'|'off'] compute ROI activity. Default is on. If
%                       you just need voxel activity, you can set this option to
%                '      off'.
%  'chansel'          - [cell array of string] channel selection. Default is all.
%  'exportvoxact'     - ['on'|'off'] export voxel activity in EEG.roi.source_voxel_data
%                       These array are huge, so the default is 'off'.
%  'fooof'            - ['on'|'off'] enable FOOOF analysis (this method can be used to parameterize neural power spectra and is described here: https://fooof-tools.github.io/fooof/). Default is 'off'.
%  'fooof_frange'     - [ ] FOOOF fitting range. Default is [1 30] like in the MATLAB example: 
%                       https://github.com/fooof-tools/fooof_mat/blob/main/examples/fooof_example_one_spectrum.m.
%  'freqresolution'   - [integer] Desired frequency resolution (in number of frequencies). If
%                       specified, the signal is zero padded accordingly.
%                       Default is 0 (means no padding).
%  'lowmemory'        - ['on'|'off'] Option to run the code with low memory, though, it might take significantly longer to complete. When turned on, the estimation of voxel-wise spectral power 
%                       will require less memory.
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%  source_voxel_data - voxel data activity (voxels x times x trials).
%                      Usually several Gb in size.
%
% Author: Stefan Haufe and Arnaud Delorme
%
% Example: call pop_roi_activity instead because it will
% compute the leadfield matrix automatically using Dipfit information.

% Copyright (C) Arnaud Delorme, arnodelorme@gmail.com
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [EEG, source_voxel_data] = roi_activity(EEG, varargin)

if nargin < 2
    help roi_activity;
    return
end

% decode input parameters
% -----------------------
g = finputcheck(varargin, { ...
    'leadfield'         {'struct' 'string'}  {{} {}}         '';
    'headmodel'         'string'             { }             ''; % sometimes useful when loading volume to see which voxels are inside/outside
    'sourcemodel'       'string'             { }             '';
    'sourcemodel2mni'   'real'               { }             [];
    'sourcemodelatlas'  'string'             { }             '';
    'modelparams'       'cell'               { }         { 0.05 };
    'model'             'string'             { 'eLoretaFieldtrip' 'lcmvFieldtrip' 'eLoreta' 'lcmv' } 'lcmv';
    'nPCA'              'integer'            { }              3;
    'downsample'        'integer'            { }              1;
    'chansel'           {'integer' 'cell'}   { {} {} }       {};
    'roiactivity'       'string'             { 'on' 'off' }  'on';
    'channelpower'      'string'             { 'on' 'off' }  'off';
    'exportvoxact'      'string'             { 'on' 'off' }  'off';
    'fooof'             'string'             { 'on' 'off'}   'off';
    'fooof_frange'      ''                   {}              [1 30];
    'freqresolution'    'integer'            {}               0;
    'lowmemory'         'string'             { 'on' 'off'}    'off';              
    'outputdir'         'string'  { }              '' }, 'roi_activity');
if ischar(g), error(g); end
if isempty(g.leadfield), error('Leadfield is mandatory parameter'); end

% Creating result folder
if ~isempty(g.outputdir)
    mkdir(fullfile( g.outputdir, 'data'));
end

% Cortex mesh or volume
% ---------------------
[~,~,ext] = fileparts(g.sourcemodel);

if strcmpi(ext, '.head')
    [~, grid, labels, strlabels ] = load_afni_atlas(g.sourcemodel, g.headmodel, g.sourcemodel2mni, g.downsample);
    uniqueROIs = unique(labels);
    nROI = length(uniqueROIs);
    cortex.Atlas(1).Name = g.sourcemodelatlas;
    for iROI = 1:nROI
        indVertices = find(labels == uniqueROIs(iROI));
        cortex.Atlas(1).Scouts(iROI).Label    = strlabels{iROI};
        cortex.Atlas(1).Scouts(iROI).Vertices = indVertices;
    end
    cortex.Vertices = grid;
else
    cortex = load(g.sourcemodel);
    if isfield(cortex, 'Faces')
        % make brainstorm coordinate system consistent with MNI coordinates for
        % plotting (in terms of axis directions)
        disp('Brainstorm cortex mesh detected - transforming to MNI coordinates');
        tf = traditionaldipfit(g.sourcemodel2mni);
        pos      = tf*[cortex.Vertices ones(size(cortex.Vertices,1),1)]';
        pos      = pos';
        cortex.Vertices = pos(:,1:3);
    elseif isfield(cortex, 'cortex') && isfield(cortex, 'atlas')
        hm = cortex;
        clear cortex;
        % align with MNI coordinates
        if ~isempty(g.sourcemodel2mni)
            tf = traditionaldipfit(g.sourcemodel2mni);
            pos      = tf*[hm.cortex.vertices ones(size(hm.cortex.vertices,1),1)]';
            pos      = pos';
        else
            pos = hm.cortex.vertices;
        end
        cortex.Vertices = pos(:,1:3);
        cortex.Faces = hm.cortex.faces;

        % make Alejandro Atlas definition compatible with Brainstrom one
        nROI = length(hm.atlas.label);
        cortex.Atlas(1).Name = g.sourcemodelatlas;
        for iROI = 1:nROI
            indVertices = find(hm.atlas.colorTable == iROI);
            cortex.Atlas(1).Scouts(iROI).Label    = hm.atlas.label{iROI};
            cortex.Atlas(1).Scouts(iROI).Vertices = indVertices;
        end
    elseif isnumeric(cortex) && mod(size(cortex,1),3) == 0 && size(cortex,2) == 6
        % NFT matrix
        cortextmp = cortex;
        clear cortex;
        cortex.Vertices = cortextmp(:,1:3);
        cortex.Atlas(1).Name = g.sourcemodelatlas;
    elseif ~isfield(cortex, 'Vertices')
        % code below is functional to load a mesh
        % However, need to align with an Atlas
        % This can be achieve with Fieldtrip functions
        sourcemodelOriOld = ft_read_headshape(fullfile(ftPath, 'template', 'sourcemodel', 'cortex_20484.surf.gii'));

        error('Unknown mesh format')
    end
end

% Select Atlas
% ------------
found = false;
for iAtlas = 1:length(cortex.Atlas)
    if strcmpi(cortex.Atlas(iAtlas).Name, g.sourcemodelatlas)
        cortex.Atlas = cortex.Atlas(iAtlas);
        found = true;
        break
    end
end
if ~found
    error('Atlas not found');
end

% leadfield matrix (Brainstorm or Fieldtrip)
% ------------------------------------------
if ~isstruct(g.leadfield)
    leadfield = load(g.leadfield, '-mat');
else
    leadfield = g.leadfield;
end
if isstruct(leadfield) && isfield(leadfield, 'roiconnectleadfield') 
    leadfield = leadfield.roiconnectleadfield;
elseif isstruct(leadfield) && isfield(leadfield, 'Gain') 
    % brainstorm
    % make format compatible with Stefan's routines
    leadfield = permute(reshape(leadfield.Gain, [], 3, nvox), [1 3 2]);
elseif isstruct(leadfield) && isfield(leadfield, 'leadfield') 
    % fieldtrip
    oldLeadfield = leadfield;
    leadfield.gain = reshape( [ leadfield.leadfield{:} ], [length(leadfield.label) 3 length(leadfield.leadfield)]);
    leadfield.gain = permute(leadfield.gain, [1 3 2]);
    leadfield = leadfield.gain;
elseif isfield(leadfield, 'LFM')
    % NFT
    leadfield = leadfield.LFM;
else
    disp('Warning: unknown leadfield matrix format, assuming array of gain values');
end

nvox = size(cortex.Vertices, 1);
nvox2 = size(leadfield,2);
if ~isequal(nvox, nvox2)
    error('There must be the same number of vertices/voxels in the leadfield and source model');
end
if isempty(g.chansel)
    g.chansel = [1:EEG.nbchan];
%     if isfield(EEG.dipfit, 'chansel')
%         g.chansel = EEG.dipfit.chansel;
%     else
%         g.chansel = 1:EEG.nbchan;
%     end
elseif iscell(g.chansel)
    g.chansel = eeg_decodechan(EEG.chanlocs, g.chansel);
end
if ~isequal(size(leadfield,1), length(g.chansel))
    error('There must be the same number of channels in the leadfield and in the list of selected channels');
end

fres = EEG.pnts/2;

% from the MVGC toolbox, compute frequencies in Hz for a
frqs = sfreqs(fres, EEG.srate);

% common average reference transform
nbchan = length(g.chansel);
H = eye(nbchan) - ones(nbchan) ./ nbchan;

% apply to data and leadfield
tmpData = reshape(H*EEG.data(g.chansel, :), nbchan, EEG.pnts, EEG.trials);
leadfield = reshape(H*leadfield(:, :), nbchan, nvox, 3);

%% source reconstruction
if strcmpi(g.model, 'eLoreta')
    % eLORETA inverse projection kernel
    disp('Computing eLoreta...');
    P_eloreta = mkfilt_eloreta_v2(leadfield, g.modelparams{:});
    
    % project to source space
    source_voxel_data = reshape(tmpData(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
elseif strcmpi(g.model, 'LCMV')
    C = cov(tmpData(:, :)');
    if length(g.modelparams) == 1
        lcmv_reg = g.modelparams{1};
    end
    alpha = lcmv_reg*trace(C)/length(C);
    Cr = C + alpha*eye(nbchan);
    [~, P_eloreta] = lcmv(Cr, leadfield, struct('alpha', 0, 'onedim', 0));
    source_voxel_data = reshape(tmpData(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
    source_voxel_data = 10^3*source_voxel_data; % the units are nA*m
else
    % transform the data to continuous so we can get an estimate for each sample
    EEG2 = EEG;
    EEG2.data = EEG2.data(:,:);
    EEG2.pnts = size(EEG2.data,2);
    EEG2.trials = 1;
    EEG2 = eeg_checkset(EEG2);
    dataPre = eeglab2fieldtrip(EEG2, 'preprocessing', 'dipfit');  
    
    % prepare data
    cfg = [];
    cfg.channel = {'all', '-EOG1'};
    cfg.reref = 'yes';
    cfg.refchannel = {'all', '-EOG1'};
    dataPre = ft_preprocessing(cfg, dataPre);

    % load head model and prepare leadfield matrix
    vol = load('-mat', g.headmodel);

    % source reconstruction
    cfg             = [];
    if lower(g.model(1)) == 'e'
        cfg.method      = 'eLoreta';
    else
        cfg.method      = 'lcmv';
    end
    try
        cfg.(g.sourcemethod) = struct(g.modelparams{:});
    catch, end
    cfg.sourcemodel = oldLeadfield;
    cfg.headmodel   = vol.vol;
    cfg.keeptrials  = 'yes';
    source          = ft_sourceanalysis(cfg, dataPre);  % compute the source
    
    % reformat for ROI analysis below
    source_voxel_data = reshape([ source.avg.mom{:} ], 3, size(source.avg.mom{1},2), length(source.avg.mom));
    source_voxel_data = permute(source_voxel_data, [2 3 1]);
end
    
% number of ROIs in the Desikan-Killiany Atlas
nROI  = length(cortex.Atlas.Scouts);

% ROI labels
labels = {cortex.Atlas.Scouts.Label};

% keep only the first nPCA strongest components for each ROI
clear tmpData
if strcmpi(g.roiactivity, 'on')
    nfreq = fres + 1;
    tmpData = reshape(source_voxel_data, EEG.pnts, EEG.trials*size(source_voxel_data,2)*size(source_voxel_data,3));
    source_roi_data = [];
    data_pnts = EEG.pnts;

    % zero padding if necessary
    if g.freqresolution ~= 0
        required_zeros = g.freqresolution - fres;
        if required_zeros < 0
            error('Desired frequency resolution cannot be lower than the actual resolution of the signal.')
        end
        pad = zeros(required_zeros, size(tmpData, 2));
        tmpData = cat(1,pad,tmpData,pad);
        frqs = sfreqs(g.freqresolution, EEG.srate);
        data_pnts = size(tmpData, 1);
        nfreq = data_pnts/2 + 1;
    end
    source_roi_power = zeros(nfreq, nROI);
    
    % compute power using the Welch method
    if strcmpi(g.lowmemory, 'off')
        % default version
        % [tmpWelch, ftmp] = pwelch(tmpData, data_pnts, floor(data_pnts/2), data_pnts, EEG.srate); % ftmp should be equal frqs 
    else
        % version for less memory, can be very slow
        warning('The sequential spectral estimation can take a very long time. It is recommended to run the faster version on a computer/cluster with enough memory.')
        tmpWelch = zeros(nfreq, size(tmpData, 2));
        fprintf('Progress of %d:', size(tmpData, 2));
        for ivox = 1:size(tmpData, 2)
            if mod(ivox, 10000) == 0
                fprintf('%d', ivox);
            elseif mod(ivox, 1000) == 0
                fprintf('.');
            end
            [tmpWelch(:, ivox), ftmp] = pwelch(tmpData(:, ivox), data_pnts, floor(data_pnts/2), data_pnts, EEG.srate); % ftmp should be equal frqs 
        end
        fprintf('\n');
    end
    % tmpWelch = reshape(tmpWelch, size(tmpWelch,1), EEG.trials, size(source_voxel_data,2), size(source_voxel_data,3));
    % tmpWelch = squeeze(mean(tmpWelch,2)); % remove trials size freqs x voxels x 3
    % tmpWelch = squeeze(mean(tmpWelch,3)); % remove 3rd dim size freqs x voxels
    % source_voxel_power = tmpWelch;

    % fooof settings
    if strcmpi(g.fooof, 'on')
        f_range = g.fooof_frange; % freq range where 1/f should be fitted 
        settings = struct(); % use defaults
        slope = zeros(1, nROI);
        PS_corrected = zeros(size(frqs, 1), size(frqs, 2), nROI);
    end
    
    source_roi_data = zeros(size(source_voxel_data,1), g.nPCA*nROI);
    fprintf('Computing ROI activity:');
    source_roi_power_norm = [];
    
    for iROI = 1:nROI
        if mod(iROI, 5) == 0
            fprintf('.');
        end
        ind_roi = cortex.Atlas.Scouts(iROI).Vertices;
        % [~, source_roi_power_norm(iROI)] =  roi_getpower(source_voxel_data, ind_roi); 
        % source_roi_power(:,iROI) = mean(tmpWelch(:, ind_roi),2); % shape: (101, nROI)

        [source_roi_data_tmp, ~] = roi_getact(source_voxel_data, ind_roi, g.nPCA);
        source_roi_data(:, (iROI-1)*g.nPCA+1:iROI*g.nPCA) = source_roi_data_tmp;
        if strcmpi(g.fooof, 'on')
            ps1 = source_roi_power(:,iROI);
            fooof_result = fooof(frqs, ps1, f_range, settings, true);
            
            offset = fooof_result.aperiodic_params(1);
            slope(iROI) = fooof_result.aperiodic_params(2);
            y = (-slope(iROI) .* log10(frqs)) + offset;
            PS_corrected(:,:,iROI) = 10*log10(ps1)-10*y;
        end
    end
    fprintf('\n');

    % version with nPCA components
    source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);
%     source_roi_data = permute(reshape(source_roi_data, data_pnts, EEG.trials, []), [3 1 2]); % error when fres is chosen to be 400 Hz

else
    source_roi_data = [];
    source_roi_power = [];
    source_roi_power_norm = [];
end
disp('Done');

% Output paramters
EEG.roi.cortex    = cortex;
EEG.roi.atlas     = cortex.Atlas.Scouts;
if strcmpi(g.exportvoxact, 'on')
    EEG.roi.source_voxel_data     = source_voxel_data; % large (takes lots of RAM)
end
% EEG.roi.source_voxel_power    = single(source_voxel_power);
EEG.roi.source_roi_data       = single(source_roi_data);
EEG.roi.source_roi_power      = source_roi_power; % used for plotting
EEG.roi.source_roi_power_norm = source_roi_power_norm; % used for cross-sprectum
EEG.roi.freqs     = frqs;
EEG.roi.nPCA      = g.nPCA;
EEG.roi.nROI      = nROI;
EEG.roi.pnts      = EEG.pnts;
EEG.roi.srate     = EEG.srate;
EEG.roi.atlas     = cortex.Atlas;
EEG.roi.srate     = EEG.srate;
EEG.roi.leadfield = g.leadfield;
EEG.roi.headmodel = g.headmodel;
EEG.roi.parameters = varargin;
if exist('P_eloreta', 'var')
    EEG.roi.P_eloreta = single(P_eloreta);
end
if strcmpi(g.fooof, 'on')
    EEG.roi.fooof_results.slope = slope;
    EEG.roi.fooof_results.offset = offset;
    EEG.roi.fooof_results.PS = ps1;
    EEG.roi.fooof_results.PS_corrected = squeeze(PS_corrected);
end

% get channel power for comparison
if strcmpi(g.channelpower, 'on')
    tmpdata = permute(EEG.data(g.chansel, :, :), [2 1 3]); % pnts trials channels
    tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
    [tmpWelch,ftmp] = pwelch(tmpdata, data_pnts, data_pnts/2, data_pnts, data_pnts/2); % ftmp should be equal frqs 
    tmpWelch = reshape(tmpWelch, size(tmpWelch,1), EEG.nbchan, EEG.trials);
    tmpWelch = squeeze(mean(tmpWelch,3)); % remove trials size freqs x voxels x 3
    EEG.roi.channel_power = tmpWelch;
end


