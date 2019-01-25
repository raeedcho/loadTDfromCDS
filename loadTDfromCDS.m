function trial_data = loadTDfromCDS(filename,params)
    % loadtdfromcds Loads trial_data structure from CDS file
    %   trial_data = loadtdfromcds(filename,params)
    %   Inputs:
    %       filename - (string) location and name of CDS file to load
    %       params - (struct) indicating which signals to load into trial_data
    %           array_name - (string) name of array, e.g. 'S1'
    %           cont_signal_names - (cell array of strings) list of signal names to extract.
    %               Could be one of:
    %                   'pos'
    %                   'vel'
    %                   'acc'
    %                   'force'
    %                   'motor_control'
    %                   'markers'
    %                   'joint_ang'
    %                   'joint_vel'
    %                   'muscle_len'
    %                   'muscle_vel'
    %                   'opensim_hand_pos'
    %                   'opensim_hand_vel'
    %                   'opensim_hand_acc'
    %                   'opensim_elbow_pos'
    %                   'opensim_elbow_vel'
    %                   'opensim_elbow_acc'
    %               List needs to be in row vector form (default: {})
    %           extract_emg - (bool) whether or not to extract emg signals (default: false)
    %           event_names - (cell array of string) list of event names to extract.
    %               Only supports events that end in 'Time' (default: {'startTime','endTime'})
    %           bin_size - (numeric) bin size at which to load trial_data (default: 0.01)
    %           trial_meta - (cell array of strings) meta information about each trial (like target_direction)
    %               to be loaded in from CDS trial table (names should be desired column names of trial table)
    %           meta - (struct) meta information for given file to be put into loaded trial data, like epoch name
    %               Note: most typical meta information (e.g. monkey, date, task) are automatically extracted from
    %               the CDS, so you don't need to add them here.

    %% default variables
    cont_signal_names = {};
    extract_emg = false;
    event_names = {'startTime','endTime'};
    bin_size = 0.01;
    trial_meta = {};
    meta = [];
    array_name = 'S1';

    assignParams(who,params)

    % check filename
    assert(ischar(filename),'filename must be a string')

    %% parameter integrity checks
    assert(iscell(cont_signal_names),'cont_signal_names needs to be a cell array')
    assert(islogical(extract_emg),'extract_emg needs to be a bool')
    assert(iscell(event_names),'event_names needs to be a cell array')
    assert(isnumeric(bin_size),'bin_size needs to be a number')
    assert(iscell(trial_meta),'trial_meta needs to be a cell')
    assert(isempty(meta) || isstruct(meta),'meta needs to be a struct')

    meta.task={'TRT'}; % for the loading of cds
    
    %% Make TD
    cont_signal_labels = get_signal_labels(cont_signal_names);
    spike_routine = @processCDSspikes;
    cds_routine = @processCDS;
    
    if extract_emg
        emg_signal_names = get_emg_names();
    else
        emg_signal_names = {};
    end
    
    % trial_data loading parameters...
    if ~isempty(meta)
        td_params = struct('bin_size',bin_size,'meta',meta);
    else
        td_params = struct('bin_size',bin_size);
    end
    
    %% load it in
    % get signal info
    signal_info = { ...
        initSignalStruct( ...
            'filename',filename, ...
            'routine',spike_routine, ...
            'params',struct(), ... % can pass arbitrary parameters. Shouldn't need much with CDS
            'name',array_name, ... % it gets stored under this name... in case of spikes, this gives S1_spikes
            'type','spikes', ... % which type... see documentation of initSignalStruct
            'label',''), ... % label can be indices or strings
        initSignalStruct( ... % continuous data
            'filename',filename, ...
            'routine',cds_routine, ...
            'params',struct('trial_meta',{trial_meta}), ...
            'name',[...
                cont_signal_names,...
                emg_signal_names,...
                event_names,...
                ], ... % stored in this name, matched to the labels below which correspond to the output of the processing routine
            'type',[...
                repmat({'generic'},1,length(cont_signal_names)),...
                repmat({'emg'},1,length(emg_signal_names)),...
                repmat({'event'},1,length(event_names)),...
                ],...
            'label',[...
                cont_signal_labels,...
                strcat('EMG_',emg_signal_names),...
                event_names,...
                ], ... % can also pass [1 2],[3 4] etc if you know the arrangment of the signals in the data struct
            'operation',[]), ...
        };
    
    % load trial_data (will result in warning for no meta info, but we're
    % taking most meta info from the CDS anyway)
    trial_data = convertDataToTD(signal_info,td_params);
    
    % add some meta information
    if any(strcmpi(cont_signal_names,'markers'))
        trial_data.marker_names = sort(get_marker_labels());
    end
    if any(contains(cont_signal_names,'joint'))
        trial_data.joint_names = get_joint_labels();
    end
    if any(contains(cont_signal_names,'muscle'))
        trial_data.muscle_names = get_muscle_labels();
    end
    if any(strcmpi(cont_signal_names,'motor_control'))
        trial_data.motorcontrol_names = {'MotorControlSho','MotorControlElb'};
    end
    
    % make it pretty
    trial_data = reorderTDfields(trial_data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = get_signal_labels(signal_names)
    % returns labels given a list of signal names
    labels = cell(1,length(signal_names));
    for sig_idx = 1:length(signal_names)
        switch lower(signal_names{sig_idx})
            case 'pos'
                labels{sig_idx} = {'x','y'};
            case 'vel'
                labels{sig_idx} = {'vx','vy'};
            case 'acc'
                labels{sig_idx} = {'ax','ay'};
            case 'force'
                labels{sig_idx} = {'fx','fy','fz','mx','my','mz'};
            case 'motor_control'
                labels{sig_idx} = {'MotorControlSho','MotorControlElb'};
            case 'markers'
                labels{sig_idx} = sort(get_marker_labels());
            case 'joint_ang'
                labels{sig_idx} = strcat(get_joint_labels(),'_ang');
            case 'joint_vel'
                labels{sig_idx} = strcat(get_joint_lables(),'_vel');
            case 'muscle_len'
                labels{sig_idx} = strcat(get_muscle_labels(),'_len');
            case 'muscle_vel'
                labels{sig_idx} = strcat(get_muscle_labels(),'_muscVel');
            case 'opensim_hand_pos'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_handPos'});
            case 'opensim_hand_vel'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_handVel'});
            case 'opensim_hand_acc'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_handAcc'});
            case 'opensim_elbow_pos'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_elbowPos'});
            case 'opensim_elbow_vel'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_elbowVel'});
            case 'opensim_elbow_acc'
                labels{sig_idx} = strcat({'X','Y','Z'},{'_elbowAcc'});
        end
    end
end

function names = get_emg_names()
    % function to return the names of EMGs
    
    names = {...
        'BiMed',...
        'FCR',...
        'FCU',...
        'FDS',...
        'DeltAnt',...
        'DeltMid',...
        'DeltPos',...
        'Trap',...
        'Lat',...
        'TerMaj',...
        'InfSpin',...
        'TriMid',...
        'TriLat',...
        'TriMed',...
        'Brad',...
        'ECRb',...
        'ECU',...
        'EDC',...
        'PecSup',...
        'PecInf',...
        'Brach',...
        'BiLat',...
        };
end

function names = get_joint_labels()
    % function to return the names of joints
    
    names = {...
        'shoulder_adduction',...
        'shoulder_rotation',...
        'shoulder_flexion',...
        'elbow_flexion',...
        'radial_pronation',...
        'wrist_flexion',...
        'wrist_abduction',...
        };
end

function names = get_marker_labels()
    % function to return the names of markers, as they appear in raw cds format
    
    temp = {...
        'Marker_1',...
        'Marker_2',...
        'Marker_3',...
        'Marker_4',...
        'Marker_5',...
        'Marker_6',...
        'Marker_7',...
        'Marker_8',...
        'Shoulder_JC',...
        'Pronation_Pt1',...
        };
    names = cell(1,3*length(temp));
    
    for i = 1:length(temp)
        names(((i-1)*3+1):(i*3)) = strcat(temp(i),{'_y','_z','_x'});
    end
end

function names = get_muscle_labels()
    % function to return the names of muscles
    
    names = {...
        'abd_poll_longus',...
        'anconeus',...
        'bicep_lh',...
        'bicep_sh',...
        'brachialis',...
        'brachioradialis',...
        'coracobrachialis',...
        'deltoid_ant',...
        'deltoid_med',...
        'deltoid_pos',...
        'dorsoepitrochlearis',...
        'ext_carpi_rad_longus',...
        'ext_carp_rad_brevis',...
        'ext_carpi_ulnaris',...
        'ext_digitorum',...
        'ext_digiti',...
        'ext_indicis',...
        'flex_carpi_radialis',...
        'flex_carpi_ulnaris',...
        'flex_digit_profundus',...
        'flex_digit_superficialis',...
        'flex_poll_longus',...
        'infraspinatus',...
        'lat_dorsi_sup',...
        'lat_dorsi_cen',...
        'lat_dorsi_inf',...
        'palmaris_longus',...
        'pectoralis_sup',...
        'pectoralis_inf',...
        'pronator_quad',...
        'pronator_teres',...
        'subscapularis',...
        'supinator',...
        'supraspinatus',...
        'teres_major',...
        'teres_minor',...
        'tricep_lat',...
        'tricep_lon',...
        'tricep_sho',...
        };
end

function out = processCDSspikes(filename,params)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loads a CDS file and returns a field
    spiking_chans  = 1:96;
    exclude_units  = 255; % sort id of units to exclude
    assignParams(who,params); % overwrite parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_flag = false;
    
    % load the CDS
    if ~isempty(filename)
        load(filename);
    else
        error_flag = true;
        disp(['ERROR: ' mfilename ': no filename provided']);
    end
    
    [data, wf] = deal(cell(1,length(cds.units)));
    tmax = 0;
    for unit = 1:length(cds.units)
        data{unit} = cds.units(unit).spikes.ts;
        wf{unit} = cds.units(unit).spikes.wave; % waveforms aren't supported right now in convertDataToTD
        tmax = max(tmax,cds.units(unit).spikes.ts(end)); % find max timestamp for time vector
    end
    
    % assume right now that the blackrock sampling is 30kHz
    t = (0:1/double(30000):tmax)';
    
    labels = [vertcat(cds.units.chan) vertcat(cds.units.ID)];
    
    % remove bogus units
    bad_idx = ~ismember(labels(:,1),spiking_chans) | ismember(labels(:,2),exclude_units);
    labels = labels(~bad_idx,:);
    data = data(~bad_idx);
    wf = wf(~bad_idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.meta   = [];
    out.data   = data;
    out.wf     = wf;
    out.labels = labels;
    out.t      = t;
    out.error_flag = error_flag;
end

function out = processCDS(filename,signal_info)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loads a CDS file and returns all continuous signals it can extract, along
    % with all the event data (trial table entries ending in 'Time').
    
    % params
    trial_meta = {}; % list of fields from the trial table to import into TD
    assignParams(who,signal_info.params); % overwrite parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_flag = false;
    
    % load the CDS
    if ~isempty(filename)
        load(filename);
    else
        error_flag = true;
        disp(['ERROR: ' mfilename ': no filename provided']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first the events...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    event_names = cds.trials.Properties.VariableNames(endsWith(cds.trials.Properties.VariableNames,'Time'));
    
    % if events are time, it expects them in cells like spiking data
    %   you can also give it already-binned events and it just passes them
    %   along
    event_data = cell(1,length(event_names));
    for eventnum = 1:length(event_names)
        event_data{eventnum} = cds.trials.(event_names{eventnum});
    end
    
    % this part could be automated
    event_labels = event_names;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now the continuous signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure out which signals to extract
    signal_names = {};
    if cds.meta.hasKinematics
        signal_names = [signal_names {'kin'}];
    end
    if cds.meta.hasForce
        signal_names = [signal_names {'force'}];
    end
    if cds.meta.hasEmg
        signal_names = [signal_names {'emg'}];
    end
    if cds.meta.hasLfp
        % signal_names = [signal_names {'lfp'}];
    end
    if cds.meta.hasAnalog
        % analog could be a lot of stuff...
        for analog_idx = 1:length(cds.analog)
            header = cds.analog{analog_idx}.Properties.VariableNames;
            % Figure out if we have motor control data
            if any(contains(header,'MotorControl'))
                signal_names = [signal_names {'motorcontrol'}];
                motorcontrol_idx = analog_idx;
            end
            if any(endsWith(header,'_ang'))
                signal_names = [signal_names {'joint_ang'}];
                joint_ang_idx = analog_idx;
            end
            if any(endsWith(header,'_vel'))
                signal_names = [signal_names {'joint_vel'}];
                joint_vel_idx = analog_idx;
            end
            if any(endsWith(header,'_len'))
                signal_names = [signal_names {'muscle_len'}];
                muscle_len_idx = analog_idx;
            end
            if any(endsWith(header,'_muscVel'))
                signal_names = [signal_names {'muscle_vel'}];
                muscle_vel_idx = analog_idx;
            end
            if any(contains(header,'_handPos'))
                signal_names = [signal_names {'hand_pos'}];
                hand_pos_idx = analog_idx;
            end
            if any(contains(header,'_handVel'))
                signal_names = [signal_names {'hand_vel'}];
                hand_vel_idx = analog_idx;
            end
            if any(contains(header,'_handAcc'))
                signal_names = [signal_names {'hand_acc'}];
                hand_acc_idx = analog_idx;
            end
            if any(contains(header,'_elbowPos'))
                signal_names = [signal_names {'elbow_pos'}];
                elbow_pos_idx = analog_idx;
            end
            if any(contains(header,'_elbowVel'))
                signal_names = [signal_names {'elbow_vel'}];
                elbow_vel_idx = analog_idx;
            end
            if any(contains(header,'_elbowAcc'))
                signal_names = [signal_names {'elbow_acc'}];
                elbow_acc_idx = analog_idx;
            end
            if any(contains(header,'Frame')) || any(contains(header,'Marker'))
                signal_names = [signal_names {'markers'}];
                markers_idx = analog_idx;
            end
        end
    end
    
    % extract signals
    samp_rate = zeros(1,length(signal_names));
    [timevec_cell,cont_data_cell,signal_labels] = deal(cell(1,length(signal_names)));
    for signum = 1:length(signal_names)
        switch lower(signal_names{signum})
            case 'kin'
                % do kin stuff
                data_table = cds.kin;
                data_cols = startsWith(data_table.Properties.VariableNames,{'x','y','vx','vy','ax','ay'});
            case 'force'
                % do force stuff
                data_table = cds.force;
                data_cols = startsWith(data_table.Properties.VariableNames,{'fx','fy','fz','mx','my','mz'});
            case 'emg'
                % do emg stuff (mostly will be processed by special 'emg' tag in convertDataToTD
                data_table = cds.emg;
                data_cols = startsWith(data_table.Properties.VariableNames,'EMG');
            case 'motorcontrol'
                % do motor control stuff
                data_table = cds.analog{motorcontrol_idx};
                data_cols = startsWith(data_table.Properties.VariableNames,'MotorControl');
            case 'joint_ang'
                % do opensim stuff
                data_table = cds.analog{joint_ang_idx};
                data_cols = endsWith(data_table.Properties.VariableNames,'_ang');
            case 'joint_vel'
                % do opensim stuff
                data_table = cds.analog{joint_vel_idx};
                data_cols = endsWith(data_table.Properties.VariableNames,'_vel');
            case 'muscle_len'
                % do opensim stuff
                data_table = cds.analog{muscle_len_idx};
                data_cols = endsWith(data_table.Properties.VariableNames,'_len');
            case 'muscle_vel'
                % do opensim stuff
                data_table = cds.analog{muscle_vel_idx};
                data_cols = endsWith(data_table.Properties.VariableNames,'_muscVel');
            case 'hand_pos'
                % do opensim stuff
                data_table = cds.analog{hand_pos_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_handPos'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'hand_vel'
                % do opensim stuff
                data_table = cds.analog{hand_vel_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_handVel'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'hand_acc'
                % do opensim stuff
                data_table = cds.analog{hand_acc_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_handAcc'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'elbow_pos'
                % do opensim stuff
                data_table = cds.analog{elbow_pos_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_elbowPos'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'elbow_vel'
                % do opensim stuff
                data_table = cds.analog{elbow_vel_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_elbowVel'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'elbow_acc'
                % do opensim stuff
                data_table = cds.analog{elbow_acc_idx};
                data_cols = find(contains(data_table.Properties.VariableNames,'_elbowAcc'));
                % don't really need guide since it's X Y Z (as long as the user requests those signals by label, in order)
            case 'markers'
                % do marker stuff
                assert(cds.meta.hasKinematics,'CDS has no kinematics!')
                
                marker_table = cds.analog{markers_idx};
                marker_cols = marker_table.Properties.VariableNames(...
                    ~strcmpi(marker_table.Properties.VariableNames,'Frame') &...
                    ~strcmpi(marker_table.Properties.VariableNames,'t')...
                    );
                marker_names = cell(1,length(marker_cols)*3);
                % get actual labels (assuming that the raw marker format is still y,z,x in lab coordinates...)
                for i = 1:length(marker_cols)
                    marker_names(((i-1)*3+1):(i*3)) = strcat(marker_cols(i),{'_y','_z','_x'});
                end
                
                marker_data = marker_table{:,3:end};
                t = marker_table.t;
                
                % interpolate to uniform sampling rate (treat nan as missing)
                dt = (t(end)-t(1))/length(t);
                tGrid = (t(1):dt:t(end))';
                marker_data_interp = zeros(length(tGrid),size(marker_data,2));
                for i=1:size(marker_data,2)
                    real_idx = find(~isnan(marker_data(:,i)));
                    marker_data_interp(:,i) = interp1(t(real_idx),marker_data(real_idx,i),tGrid);
                end
                
                % set in data_table
                data_table = array2table([tGrid marker_data_interp],'VariableNames',[{'t'} marker_names]);
                data_cols = 2:width(data_table);
            otherwise
                error('No idea what this signal is (%s)',signal_names{signum})
        end
        signal_labels{signum} = data_table.Properties.VariableNames(data_cols);
        cont_data_cell{signum} = data_table{:,data_cols};
        timevec_cell{signum} = data_table.t;
        samp_rate(signum) = 1/mode(diff(data_table.t));
    end
    
    % resample everything to the highest sampling rate
    [final_rate,maxrate_idx] = max(samp_rate);
    t_end = 0;
    for signum = 1:length(signal_names)
        % [P,Q] = rat(final_rate/samp_rate(signum),1e-7);
        [P,Q] = rat(final_rate/samp_rate(signum)); % sometimes P and Q are very big, so use default tolerance...
        if P~=1 || Q~=1
            assert(signum~=maxrate_idx,'Something went wrong with the resample code...')
            
            % figure out where the NaNs are before resampling (mostly for markers, to see where they cut out)
            nan_spots = isnan(cont_data_cell{signum});
            if any(any(nan_spots))
                nanblock_thresh = 0.3/mode(diff(timevec_cell{signum})); % tolerate nan blocks up to 0.3 seconds long
                nan_transitions = diff([zeros(1,size(nan_spots,2));nan_spots;zeros(1,size(nan_spots,2))]);
                nanblock_endpoints = cell(1, size(nan_spots,2));
                for i = 1:size(nan_spots,2)
                    nan_starts = find(nan_transitions(:,i)==1);
                    nan_stops = find(nan_transitions(:,i)==-1);
                    nanblock_lengths = nan_stops-nan_starts;
                    
                    nan_starts(nanblock_lengths<nanblock_thresh) = [];
                    nan_stops(nanblock_lengths<nanblock_thresh) = [];
                    
                    % save times for nanblocks
                    nanblock_endpoints{i} = [timevec_cell{signum}(nan_starts) timevec_cell{signum}(nan_stops-1)];
                end
            else
                nanblock_endpoints = {};
            end
            
            % need to resample
            % need to detrend first...
            % detrend first because resample assumes endpoints are 0
            a = zeros(2,size(cont_data_cell{signum},2));
            dataDetrend = zeros(size(cont_data_cell{signum},1),size(cont_data_cell{signum},2));
            for i = 1:size(cont_data_cell{signum},2)
                % in case start or end are nans
                nanners = isnan(cont_data_cell{signum}(:,i));
                data_poly = cont_data_cell{signum}(~nanners,i);
                
                t_poly = timevec_cell{signum}(~nanners);
                a(1,i) = (data_poly(end)-data_poly(1))/(t_poly(end)-t_poly(1));
                a(2,i) = data_poly(1);
                
                dataDetrend(:,i) = cont_data_cell{signum}(:,i)-polyval(a(:,i),timevec_cell{signum});
            end
            temp=resample(dataDetrend,P,Q);
        
            % interpolate time vector
            % using upsample -> downsample to save memory (it's the same thing
            % as the reverse) but it adds extra points at the end that aren't
            % in the resampled data
            resamp_vec = ones(size(cont_data_cell{signum},1),1);
            resamp_vec = upsample(downsample(resamp_vec,Q),P);
            ty=upsample(downsample(timevec_cell{signum},Q),P);
            ty=interp1(find(resamp_vec>0),ty(resamp_vec>0),(1:length(ty))');
            
            % get rid of extrapolated points at the end
            extrap_idx = isnan(ty);
            ty(extrap_idx) = [];
            temp(extrap_idx(1:size(temp,1)),:) = [];
    
            % retrend...
            dataResampled = zeros(size(temp,1),size(temp,2));
            for i=1:size(dataDetrend,2)
                dataResampled(:,i) = temp(:,i)+polyval(a(:,i),ty(:,1));
            end
            
            % set nan blocks back into resampled data...
            if ~isempty(nanblock_endpoints)
                for i=1:size(dataResampled,2)
                    for j=1:size(nanblock_endpoints{i},1)
                        inan = ty>=nanblock_endpoints{i}(j,1) & ty<=nanblock_endpoints{i}(j,2);
                        dataResampled(inan) = NaN;
                    end
                end
            end
            
            % assign back into cell
            cont_data_cell{signum} = dataResampled;
            timevec_cell{signum} = ty;
        end
        
        % collect max time for time extension
        t_end = max(t_end,timevec_cell{signum}(end));
    end
    
    % extend time vectors to be the same length and interpolate to unified time
    dt = mode(diff(timevec_cell{maxrate_idx}));
    t = (0:dt:t_end)';
    for signum = 1:length(signal_names)
        % interpolate to new time vector (fill extrapolated points with NaNs)
        cont_data_cell{signum} = interp1(timevec_cell{signum},cont_data_cell{signum},t);
    end
    
    % try horizontally concatenating...If everything went well, things should
    % be the right length...
    cont_data = horzcat(cont_data_cell{:});
    cont_labels = horzcat(signal_labels{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now the meta...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meta_info = struct(...
        'monkey',cds.meta.monkey,...
        'task',cds.meta.task,...
        'date_time',cds.meta.dateTime,...
        'trialID',cds.trials.number',...
        'result',cds.trials.result'...
        );
    
    % import the trial table meta info
    for i=1:length(trial_meta)
        if ischar(trial_meta{i}) && ismember(trial_meta{i},cds.trials.Properties.VariableNames)
            meta_info.(trial_meta{i}) = cds.trials.(trial_meta{i})';
        else
            warning('Element %d of trial_meta was not imported',i)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.meta   = meta_info;
    out.cont_data   = cont_data;
    out.cont_labels = cont_labels;
    out.event_data   = event_data;
    out.event_labels = event_labels;
    out.t      = t;
    out.error_flag = error_flag;
end
