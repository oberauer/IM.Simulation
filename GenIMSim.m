% Simulation implementation of the IM - Simulates several experiments
% This is the master script for selecting models, experiments, and parameter values
% Version using the KO/Manohar mediated-binding mechanism

clear all
close all

global E
global C
global P
%path(strcat(pwd,'\Compiling'),path);

model = 1;  % 1 = IMSim

% Experiments
% 1: Continuous reproduction, sequential presentation in separate locations, cue by location, SPC, n responses (Oberauer & Lin, 2023)
% 2: Continuous reproduction, sequential presentation in one location, cue by feature, SPC, 1 response (Gorgoraptis et al.)
% 3: Continuous reproduction, forward serial recall
% 4: Continuous reproduction, sequential vs. simultaneous presentation, set size (+ CDA and alpha power)
% 5: Continuous reproduction, setsize + RI (Pertzov)
% 6: Continuous reproduction, seq. vs. sim. encoding of bindings between 2 features (Shepherdson et al., 2022)
% 7: Continuous reproduction, sim. vs seq. with varying presentation rates: UZH students (Zepp & Oberauer, in prep)
% 8: Continuous reproduction, sim. vs seq. with varying presentation rates: Prolific (Zepp & Oberauer, in prep)
% 9: Continuous reproduction, 2 arrays with varying set sizes and inter-array time (Zepp & Oberauer, in prep)
% 10: Continuous reproduction, consolidation time and mask SOA (Ricker & Sandry, 2018)
% 11: Continuous reproduction, setsize + mask-SOA (Bays et al., 2011); 
% 12: Continuous reproduction, masking (Agaoglu et al., 2015)
% 13: Continuous reproduction, CDA, individual differences 
% 14: Continuous reproduction, channel tuning functions from IEM
% 15: Continuous reproduction, track CDA, alpha power, and IEM-CTF over time
% 16: Continuous reproduction, multi-feature objects and CDA
% 17: Continuous reproduction, simultaneous presentation, set size + pre-cue
% 18: Continuous reproduction, simultaneous presentation, set size + retro-cue
% 19: Continuous reproduction, cue-target interval
% 20: Continuous reproduction, wheel-attraction effect (Souza et al., 2016, JEP:HPP, Exp. 6)
% 21: Continuous reproduction, guided refreshing (Souza et al., 2015, ANYAS)
% 22: Continuous reproduction, removal (Williams & Woodman, 2012)
% 23: Continuous reproduction, removal (Gunseli et al., 2015)
% 24: Continuous reproduction, varying RI and absence of test interference (Hautekiet & Oberauer, 2026)
% 25: Continuous reproduction, double-cue with IEM-CFT of feature content
% 26: Continuous reproduction, interruption with CTF (van Moorselaar et al. 2017)
% 27: CD, setsize, retrocue
% 28: CD, set size and continuous degrees of change
% 29: CD with retro-cueing and re-loading (Souza et al., 2014, JEP:HPP, Exp. 2)
% 30: CD with retro-cue and response selection delay
% 31: CD with 2-cue intrusion (Rerko & Oberauer, 2013, Exp. 2)
% 32: CD with 3-cue ABA vs. CBA (Rerko & Oberauer, 2013, Exp. 3)
% 33: CD with SOA (array, cue) to test sensory memory (Pratte & Greene, 2023, Exp. 1)
% 34: Change localization with variation of response set size (He, Kellen & Singmann, 2026)

% 40 = Test retro-cue mechanisms individually
% 41 = Test retro-cue mechanisms fully crossed, for single and dual retro-cue
% 42 = Retro-cue & strengthening
% 43 = Setsize and Alpha oscillations of spatial attention
% 44 = Effect of Nb on accuracy
% 45 = Generic Parameter-Sensitivity simulation for continuous reproduction (sequential, set size 4)
% 46 = Generic Parameter-Sensitivity simulation for continuous reproduction (simultaneous, set size 6, fit MM and SDM)
% 47 = Generic Parameter-Sensitivity simulation for change detection (simultaneous, set-size 6)

saveResults = 0;
Exp = 34;

% Exp's     1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
Material = [1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]; % 1 = 360 degrees, 2 = discrete items, 3 = 180 degrees 

Setsize = 6;  % default value (can be overwritten later)
fitMM = 0;   % fit mixture model?
fitIMSim = 0; % fit IM?

%%% Experimental Constants/Defaults

E.ntrials = 200;     % number of trials to run per subject and condition
E.nsubj = 20;        % number of subjects
E.ngroups = 1;       % number of groups of subjects
E.material = Material(Exp);  % 1 = features on a continuous circular dimension (e.g., color wheel); 2 = highly distinct features
E.targetDim = 1;     % feature dimension of the target stimuli: 1 = color, 2 = orientation, 3 = spatial location
E.test = 1;          % 1 = continuous reproduction, 2 = change detection (here: set default value, can be overwritten later)
E.context = 1;       % 1 = spatial location, 2 = another feature
E.calibrateAmp = 1;  % level to calibrate amplification factor: 1 = population, 2 = individual
E.maxsetsize = 8;    % maximum set size
E.prestime = 0.5;    % presentation time before onset of another attended stimulus that interrupts consolidation (e.g., another item, or an attention-demanding distractor)
E.ISI = 0;           % inter-item interval (for seq. presentation): time between offset of previous and onset of next stimulus (-> SOA = E.prestime + E.ISI)
E.outsize = 1;       % output size = number of item tested (default: only 1 item tested)
E.forwardrecall = 0; % default: random recall order
E.presentation = 1;  % 1 = simultaneous, 2 = sequential
E.nfeat = 1;         % number of features of each object that need to be remembered
E.layout = 1;        % 1 = spatially distributed on a virtual circle, 2 = all centrally presented
E.PreRetro = 2;      % 1 = pre-cue, 2 = retro-cue
E.CTI = [0, 1, 1, 0.5, 0, 1]; % cue-target interval for no-cue, valid cue, invalid cue, refreshing, delay of response selection, multi-cueing
E.CQI = [0, 1, 1, 0,   1, 1]; % cue-question interval for no-cue, valid cue, invalid cue, refreshing, delay of response selection, multi-cueing
E.RI = 1;            % retention interval up to onset of the cue (default setting)
E.wheel = 1;         % 1 = color-wheel, 2 = grey wheel
E.MaskSOA = 100;     % default: no mask, so the SOA is huge
E.nDistr = [0, 0];   % number of to-be-attended distractors before/after the retro-cue
E.saveResults = saveResults;

%%% Constants of the memory system and the simulations

if E.material == 1, C.nstim = 360; end    % number of stimuli - here: 360 different colors (or orientations)
if E.material == 2, C.nstim = 12; end     % number of stimuli - here: 12 different highly distinct features (maximally spaced on the circular feature space)
if E.material == 3, C.nstim = 180; end % orientation of bars: 180 orientations
C.nfeatures = 1;  %number of feature dimensions (default)
C.nloc = 13;      %number of possible object locations (on a virtual circle around fixation)
C.nc = 360;       %number of units to represent color space (or orientation space)
C.nCat = 8;       % number of content categories
C.nLocCat = 8;    % number of context categories
C.x = pi*(1:C.nc)./180;  % x axis for the distributions of population codes in circular color space / orientation space
C.tstep = 0.05;    % time steps for simulation (in s)
C.nsamples = 100;  % number of samples in multinomial accumulator
C.ContentFun = @VonMises;
C.ContextFun = @VonMises;
C.accum = 1; 
C.CDA = 1;         % 1: use number of committed binding units; 2 = use summed weight in weight matrices
C.spatInhibResolution = 0; % 0 = model the dynamics of spatial attention in one step per consolidation event; 1 = model it iteratively by time steps
C.seqVariant = 2;  % 1 = consolidation and replenishment until strength is reached; 2 = consolid until strength is reached; replenishment until time is over; 3 = both until time is over
C.consolidAttempt = 100; % default: all items are attempted to be consolidated
C.retroCueConsolid = 0;  % use retro-cue to consolidate: 0 = not at all, 1 = only once per item to reach C.strength, 2 = unlimited additional codes with separate sets of binding units

%%% Parameters

P.kappaf_feat = 25;  % precision of original stimuli (in the sensory layer), which is also the feature precision in the focus of attention, for content features
P.kappa_feat = 25;   % mean precision of categories, for content features
P.kappaf_ctx = 25;   % precision of original stimuli, and the focus of attention, for context (needs to be fairly high, otherwise CW intrusion becomes too big)
P.kappa_ctx = 25;    % precision of categories for context
P.kappaCatSD = 3;    % SD of precision values of categories (variability across categories)
P.mCatSD = 5;        % SD of deviation of category center from equal spacing
P.delta = 0.80;      % proportion of committed binding units that remain committed, and weights that remain, upon encoding of each new item
P.pMax = 1.0;        % the initial proportion of binding units recruited 
P.pBase = 0.3;       % minimal (base) strength of bindings (lower asymptote)
P.keepFocus = 0.3;   % probability of keeping the last-presented item in the FoA until test
P.a = 0.1;           % strength of item memory - implemented as "C.locationnoise" in CreateStimuli: all location cues receive some baseline activation
P.nb = 100;          % number of units in the binding layer
P.nbNorm = sqrt(P.nb); % normalization constant depends on mean P.nb, not on individual P.nb (and not on manipulation of P.nb in simulation 33)
P.spatinhib = 0.0;   % global inhibition of spatial attention distribution
P.TopDownSpatAttn = 0; % strength of top-down modulation of spatial attention from the FoA in WM (AfocusLoc)
P.maskWindow = 0.05;  % mean of time window within which a mask or a cue is integrated with the current feature Map
P.maskWindowSD = 0.75; % SD (as porportion of mean) of time window of integration
P.SDstrengthFX = 0.1; % SD of encoding strength into FX
P.selfactFX = 1;     % self-activation of FX
P.inhibFX = 0.002;   % global inhibition in FX that causes decay
P.eraseFX = 0.2;     % degree to which FX is erased by onset of a new attended stimulus (1 = not at all, 0 = completely)
P.cRate = 10;        % rate of short-term consolidation (gain in strength of bindings)
P.rRate = 4;         % rate of release of BP units
P.cRateFactor = 1;   % proportional reduction of cRate for Ricker's dots on a ring
P.cRateSD = 0.5;     % 0.7 - SD of consolidation rates (as proportion of mean)
P.cStrength = 0.9;   % proportion of maximal strength that consolidation aims for - when that strength is reached, consolidation stops
P.cBallistic = 0.5;  % probability of consolidation being ballistic
P.filter = [0.1, 0.1, 0.1]; % strength of encoding of the test display (colorwheel or probe) when attended (with probability P.eraseFX) 
P.rad1 = 0.7;        % proportion of radius of memory array to radius of color wheel (for computation of color-wheel interference as a function of distance between wheel and target location)
P.outputinterference = 0; % proportion of reduction of W
P.wnoise = 0.03;     % noise added to W at each time step to implement decay
P.cueingStrength = 1;     % amount of strengthening by re-encoding in response to retro-cue (0 = none)
P.removalThreshold = 0.3; % binding units with abs(strength) below the threshold are removed - now expressed as proportion of the maximal absolute retrieved binding strength
P.removalTau = 0.8;  % threshold for the logistic translating cue validity into removal strength
P.removalGain = 10;  % gain for the logistic translating cue validity into removal strength
P.inhib = 0;         % global inhibition of activation during retrieval
P.cuerate = 5;       % rate of using the cue
%P.dnoise = 2.5;     % SD of noise added to each accumulator in recall/recollection of a feature
P.driftnoise = 1.0;  % SD of noise added to each accumulator for recognition decision
%P.boundary = [70, 20]; % boundary for response for [recall/recollection, recognition decision] 

P.dnoise = 1.5;      % SD of noise added to each accumulator in recall/recollection of a feature
P.boundary = [30, 10]; % boundary for response for [recall/recollection, recognition decision] 

P.sz = 0;             % starting point variability for yes/no accumulators for recognition
P.kappacrit = 1;      % proportion of kappa_feat: meta-cognitive estimate of average precision -> used in Bayesian optimal decision rule for same-change decision in Change Detection

C.indVar = {'nb', 'dnoise', 'delta', 'pBase', 'kappa_feat', 'kappaf_feat', 'kappa_ctx', 'kappaf_ctx',  ...
           'eraseFX', 'inhibFX', 'cRate', 'rRate', 'removalThreshold'};  % parameters to be varied in individual-differences simulation
       
%             nb    dnoise delta pBase  k_f   kf_f  k_c   kf_c  eraseFX  inhibFX  cRate rRate remThr        
C.maxIndVar = [500,   20,   1,   0.9,    30,   30,   30,   50,    1,     0.1,     20    20    10];   % max. value of parameters varied in individual-differences simulation
C.SDfactor =  [0.25,  0.25, 0.25, 0.25,  0.25  0.25, 0.25, 0.25,  0.25,  0.25,    0.25  0.25  0.25]; % SD as fraction of mean
C.logistVar = [0      0     1     1      0     0     0     0      1      0        0     0     0];    % whether or not the parameter's individual differences are normal or on a logit scale
C.nullVar =   [NaN,   0,    1     0      0     0     0     0      0      0        NaN   NaN   0];    % value of the parameter that neutralizes an effect -> don't create individual differences because then the values depart from the null value!
C.groupVar = ['none'];                                                                               % grouping of participants in individual-differences simulation - max. one parameter name as string variable!
C.maxGroupVar = repmat(inf, 1, length(C.indVar)); 

if model==1, Model = @IMSim; end          % depending on the model chosen, define Model as the function to be used for simulating data

% Choose a function to run the experiment selected
if Exp == 1, SetsizeSPC(Model, 1, Setsize, 0); end  % sequential presentation, continuous reproduction with serial-position effects
if Exp == 2, SetsizeSPC(Model, 2, Setsize, 0); end  % sequential presentation, continuous reproduction with serial-position effects
if Exp == 3, SetsizeSerialRecall(Model, Setsize); end
if Exp == 4, SimSeqAlphaCDA(Model, Setsize); end % sequential or simultaneous presentation, CDA and Alpha power suppression
if Exp == 5, SetsizeRI(Model, 1, fitMM, fitIMSim); end  % set size and retention-interval variation
if Exp == 6, SimSeq(Model, fitMM); end % sequential/simultaneous encoding of 2 objects with 2 features
if Exp == 7, SimSeqPresentationRate(Model, P.cRate, fitMM); end % simultaneous encoding vs. sequential encoding with varying presentation rates (2nd parameter: consolidation rate)
if Exp == 8, P.cStrength = 0.7; SimSeqPresentationRate(Model, P.cRate, fitMM); end
if Exp == 9, TwoArrayISI; end    % two successive arrays with varying set sizes, varying inter-array interval
if Exp == 10, Consolidation(Model, fitMM); end % consolidation time with SOA variation
if Exp == 11, SetsizeMaskSOA(Model, fitMM); end  % Bays et al. (2011)
if Exp == 12, Masking(Model, fitMM); end  % Agaoglu et al., (2015)
if Exp == 13, [Kestimate, MMSD, MMguessing, MMtranspos, Mwact, Decayrate] = SetsizeIndDiff(Model, 0); end %Continuous reproduction, individual differences and CDA
if Exp == 14, [MMSD, MMguessing, MMtranspos] = SetsizeCTF(Model, 0, 0); end % Continuous reproduction, channel tuning functions
if Exp == 15, SetsizeTrackNeuralSignals(Setsize); end
if Exp == 16, MultiFeatureCDA(Model, Setsize); end  % Variation of number of features and setsize
if Exp == 17, SetsizeCueing(Model, 1, fitMM, fitIMSim); end  % Continuous reproduction, Pre-cue
if Exp == 18, SetsizeCueing(Model, 2, fitMM, fitIMSim); end  % Continuous reproduction, Retro-cue
if Exp == 19, CueTargetInterval(Model, 6, 0, 0); end  % Souza et al., 2016
if Exp == 20, WheelAttraction(Model, Setsize, fitMM, fitIMSim); end % Continuous reproduction, retro-cue and wheel attraction
if Exp == 21, Refreshing(Model, Setsize, fitMM, fitIMSim); end  % Continuous reproduction, guided refreshing
if Exp == 22, Removal(Model, 2, 0, 0); end  % Woodman & Williams (2012)
if Exp == 23, Removal(Model, 4, 0, 0); end  % Gunseli et al (2015)
if Exp == 24, RetroCueDecayInterference(Model, 5, 0, 0); end % Hautekiet & Oberauer (2026)
if Exp == 25, DoubleCueCTF(Model); end
if Exp == 26, InterruptCTF; end  % van Moorselaar et al. (2017)
if Exp == 27, E.test = 2; E.wheel = 0; SetsizeCueingCD(Model); end  % CD for set-size and retro-cue manipulation
if Exp == 28, E.test = 2; E.wheel = 0; SetsizeDeltaCD(Model); end  % CD for set-size and degree-of-change manipulation
if Exp == 29, E.test = 2; E.wheel = 0; Reloading(Model); end % CD for retro-cue and re-loading experiment
if Exp == 30, E.test = 2; E.wheel = 0; DelayRS(Model); end  % CD for retro-cue and delay of response selection
if Exp == 31, E.test = 2; E.wheel = 0; MultiCueIntrusion(Model); end  % 2-cues (last always valid), with intrusion probes sometimes matching the first-cued item
if Exp == 32, E.test = 2; E.wheel = 0; MultiCueABA(Model); end  % 3-cues (last always valid), with CBA vs. ABA cueing sequence
if Exp == 33, E.test = 2; E.wheel = 0; SensoryMemoryCD(Model); end  % CD with varying SOA from array to probe
if Exp == 34, E.test = 3; E.wheel = 0; ROC(Model, 3); end  % reconstruction of ROC curves from change localization with variable response set size. Second parameter = probe type of change

if Exp == 40, RetroCueSeparateMechanisms(Model, [1,6], 1, 1:2, fitMM); end  % Retro-cue exploration. Arguments are Mechanisms, Tasks (1=CR, 2=CD), Cueing conditions (1=neutral, 2=valid, 3=invalid)
if Exp == 41, RetroCueFullDesign(Model, C.indVar, C.maxIndVar, fitMM); end
if Exp == 42, RetroCueStrength(Model, fitMM); end  % Retro-cue exploration
if Exp == 43, SetsizeAlpha(@IMSimAlpha, 8); end
if Exp == 44, NbindingCapacity(Model); end
if Exp == 45, ParameterSensitivity(Model, 'dnoise', [1:5]); end
if Exp == 46, ParameterSensitivity2(Model, 'a', [0, 0.05, 0.1, 0.2, 0.3, 0.4]); end
if Exp == 47, ParameterSensitivityCD(Model, 'a', [0, 0.05, 0.1, 0.2, 0.3, 0.4]); end

