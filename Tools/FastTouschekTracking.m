function [TL,MA_pos,MA_neg] = FastTouschekTracking(settingFilePath)
%{
Main execution of FastTouschekTracking method.
First, the dynamic apertures for the specified range of dp-values
is computed.
Then, the momentum acceptance is found using the interface-method
Lastly, the Touschek lifetime is computed using Piwinski's formula 

A path to the file location of a .mat file containing the
"Settings" struct can be given as input instead of the struct itself


REQUIRED fields of "Settings" struct:
Settings.ring           - lattice
Settings.outputName     - location and filename of output of results
Settings.dpmax          - maximum value of dp/p for momentum acceptance search
                        - also maximum dp/p used for dynamic aperture search
                        (see computeDA())
Settings.dpResolution   - resolution of dp/p in momentum acceptance search
                        - also boundary-resolution of x,x',dp/p 3d volume
                        (see computeDA())
Settings.dpSliceStep    - dp/p separation of two adjecent x,x' DA slices
Settings.Ib             - single bunch current (Touschek Calc.)
Settings.params         - beam parameters (computed using, e.g., atx())
    Settings.params.modemittance(1:2)   - horiz. and vert. mode-emittance
    Settings.params.espread             - bunch energy spread
    Settings.params.blength             - bunch length
    

OPTIONAL field of "Settings" struct:
Settings.G              - matrix to normalize x,x' phase space 


NOTE: see also requirements for "Settings" in the chosen dynamic aperture
computation method!

%}


fprintf('\n\n ### BEGINNING FAST TOUSCHEK TRACKING ###\n\n')
if isstring(settingFilePath) || ischar(settingFilePath)
    load(settingFilePath,'Settings')
else
    Settings = settingFilePath;
end
t1 = tic;% start time taking to compute total duration

ring = Settings.ring;

% convert xmax and x-resolution into an maximum ampltidue and resolution
if isfield(Settings,'G')
    G = Settings.G;
else
    [~,LinData] = atlinopt2(ring,1);
    G = [1/sqrt(LinData.beta(1)),0;LinData.alpha(1)/sqrt(LinData.beta(1)),sqrt(LinData.beta(1))];
    Settings.G = G;
end


%% Compute DAs for the various values of dp/p for the lattice
fprintf(['   Launching DA computation. Method: ',Settings.DAmethod,'\n'])
[u_da,up_da,x_da,xp_da,Settings] = computeDA(Settings);
dpoffsets_DA = Settings.dpoffsets_DA;
fprintf('   DA computation finished\n')

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

fprintf('   Computing positive momentum acceptance\n')
MA_pos = InterfaceMethod_binary(+1);% run interface method for positive energies
% MA_pos = InterfaceMethod_line(+1);% run interface method for positive energies

fprintf('   Computing negative momentum acceptance\n')
MA_neg = InterfaceMethod_binary(-1);% run interface method for negative energies
% MA_neg = InterfaceMethod_line(-1);% run interface method for negative energies

duration_total = toc(t1);% finish time taking

[~,LinData] = atlinopt4(ring,index);
dpp= [MA_pos',MA_neg'];
params = Settings.params;

fprintf('   Computing Touschek lifetime\n')
[TL,~]=TouschekPiwinskiLifeTime_inputoptics(ring,dpp,Settings.Ib,index,'lo',LinData,'pa',params);
TL = TL/60/60;
fprintf(['       Result: \n',num2str(TL),' hours\n'])

save(Settings.outputName,'Settings','MA_pos','MA_neg','duration_total','index','TL','LinData','params','x_da','xp_da','u_da','up_da')


    function MA_found = InterfaceMethod_binary(posNeg)
        % Interface method based on a binary search at each element location
        values = 0:posNeg*Settings.dpResolution:posNeg*Settings.dpmax;
        low = ones(1,numel(index));
        high = numel(values).*ones(1,numel(index));
        index_binary = index;

        while ~isempty(index_binary)
            whichMembers = find(ismember(index,index_binary));
            k = floor((high-low)/2+low);
            Lout = trackAround(ring,values(k),index_binary);

            isIn = [];
            for a =  1:numel(index_binary)
                isIn(a) = isWithin(Lout(1:2,a),values(k(a)));% check if particle is within polygon
            end
            isIn = logical(isIn);
            high(~isIn) = k(~isIn);
            low(isIn) = k(isIn);

            isAtResolutionLimit = (high-low)<=1;
            MA_found(whichMembers(isIn)) = values(k(isIn));
            
            index_binary(isAtResolutionLimit)=[];
            high(isAtResolutionLimit)=[];
            low(isAtResolutionLimit)=[];
        end

    end


    function MA_found = InterfaceMethod_line(posNeg)
        % Interface method based on a line search at each element location
        values = 0:posNeg*Settings.dpResolution:posNeg*Settings.dpmax;
        index_line = index;
        for b =  1:numel(values)
            whichMembers = find(ismember(index,index_line));
            Lout = trackAround(ring,values(b)*ones(size(index_line)),index_line);

            isIn = [];
            for a =  1:numel(index_line)
                isIn(a) = isWithin(Lout(1:2,a),values(b));% check if particle is within polygon
            end
            isIn = logical(isIn);
            MA_found(whichMembers(isIn)) = values(b);
            index_line(~isIn) = [];
            if isempty(index_line)
                break;
            end
        end
    end


    function Lout = trackAround(ring,dp,indices)
        %Function to track particles from multiple, different indicies to the reference location
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
            % transfer from normalized to physical phase-space
            P = inv(G)*[ui;upi];
            isIn = inpolygon(xxp(1),xxp(2),P(1,:),P(2,:));
        end
    end



end
