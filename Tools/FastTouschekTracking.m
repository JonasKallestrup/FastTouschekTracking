function TL = FastTouschekTracking(settingFilePath)
addpath('/psi/home/kallestrup_j/at/atmat/')
addpath('/psi/home/kallestrup_j/at_jk/');
atpath;

if isstring(settingFilePath) || ischar(settingFilePath)
    load(settingFilePath,'Settings')
else
    Settings = settingFilePath;
end
t1 = tic;% start time taking to compute total duration

ring = Settings.ring;
G = Settings.G;

%% Compute DAs for the various values of dp/p for the lattice
[u_da,up_da,x_da,xp_da,Settings] = computeDA(Settings);
dpoffsets_DA = Settings.dpoffsets_DA;

%% perform interface method
% localize indices in lattice with non-zero length
index = 1:length(ring);
length_elems = atgetfieldvalues(ring,index,'Length');
index(length_elems==0) = [];% don't track where there's no length

% Compute closed orbit for on-momentum particle
if ~isempty(findcells(ring,'Frequency'))
    orbit0 = findorbit6(ring,1:length(ring)+1);
else
    orbit0 = findorbit4(ring,0,1:length(ring)+1);
    orbit0(5:6,:)= 0 ;
end

MA_pos = InterfaceMethod_binary(+1);% run interface method for positive energies
MA_neg = InterfaceMethod_binary(-1);% run interface method for negative energies
duration_total = toc(t1);% finish time taking

[~,LinData] = atlinopt4(ring,index);
dpp= [MA_pos',MA_neg'];
params = Settings.params;

[TL,~]=TouschekPiwinskiLifeTime_inputoptics(ring,dpp,Settings.Ib,index,'lo',LinData,'pa',params);
TL = TL/60/60;

save(Settings.outputName,'Settings','MA_pos','MA_neg','duration_total','index','TL','LinData','params','x_da','xp_da','u_da','up_da')

    function MA_found = InterfaceMethod_binary(posNeg)
        % use binary search to find local MA
        MA_lost = posNeg*Settings.dpmax*ones(1,numel(index));
        MA_check = MA_lost;
        MA_survives = zeros(1,numel(index));
        MA_found = MA_survives;
        index_binary = index; % copy of lattices indexes to be checked

        while ~isempty(index_binary)
            Lout = trackAround(ring,MA_check,index_binary);
            toDelete = [];
            for a = 1:numel(index_binary)
                isIn = isWithin(Lout(1:2,a),MA_check(a));% check if particle is within polygon
                if isIn % particle survives
                    dMA(a) = (MA_lost(a)-MA_check(a))/2;
                    MA_survives(a) = MA_check(a);

                else % particle is lost
                    dMA(a) = (MA_survives(a)-MA_check(a))/2;
                    MA_lost(a) = MA_check(a);
                end

                if abs(dMA(a)) < Settings.dpResolution
                    if MA_survives(a) == 0
                        MA_found(a) = NaN;
                    else
                        MA_found(a) = MA_survives(a);
                    end
                    toDelete = [toDelete,a];
                else
                    MA_check(a) = MA_check(a)+dMA(a);
                end
            end

            index_binary(toDelete) = [];
            MA_lost(toDelete) = [];
            MA_survives(toDelete) = [];
            MA_check(toDelete) = [];
            dMA(toDelete) = [];
        end
    end

    function Lout = trackAround(ring,dp,indices)
        Lin = orbit0(:,indices(1))+[1e-9;0;1e-9;0;dp(1);0];
        for a = 1:numel(indices)
            if a < numel(indices)
                Lin = linepass(ring(indices(a):indices(a+1)),Lin,numel(ring(indices(a):indices(a+1))));
                Lin = [Lin,orbit0(:,indices(a+1))+[1e-9;0;1e-9;0;dp(a+1);0]];
            else
                Lout = linepass(ring(indices(a):end),Lin,numel(ring(indices(a):end))+1);
            end
        end
    end

    function isIn = isWithin(xxp,dpval)
        % Check if particle is within DA-polygon. Do linear interpolation 
        % between DA slices if necessary. 

        if any(dpval == dpoffsets_DA)
            % DA slice with requested dp exists; interpolation not needed.
            idx_correctdp = dpval==dpoffsets_DA;
            ui = u_da{idx_correctdp};
            upi = up_da{idx_correctdp};
        else 
            % interpolate between two DA slices
            idx_low = find(dpval>dpoffsets_DA,1,'last');
            idx_high = find(dpval<dpoffsets_DA,1,'first');
            if or(isempty(idx_low),isempty(idx_high))
                ui = 0;
            else
                % calculate interpolated polygon
                lambda = (dpval-dpoffsets_DA(idx_low))/(dpoffsets_DA(idx_high)-dpoffsets_DA(idx_low));
                ui = (1-lambda)*u_da{idx_low}+lambda*u_da{idx_high};
                upi = (1-lambda)*up_da{idx_low}+lambda*up_da{idx_high};
            end
        end
        if all(ui == 0)
            isIn = 0;
        else
            % transfer to physical phase-space
            P = inv(G)*[ui;upi];
            isIn = inpolygon(xxp(1),xxp(2),P(1,:),P(2,:));
        end
    end



end
