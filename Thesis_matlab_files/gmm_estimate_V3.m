function gmm_estimate_V3(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 4;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1 : block.NumInputPorts
    block.InputPort(i).Dimensions        = 1 ;
    block.InputPort(i).DatatypeID  = 0 ;  % double
    block.InputPort(i).Complexity  = 'Real' ;
    block.InputPort(i).DirectFeedthrough = true ;
end
% Override output port properties
block.OutputPort(1).Dimensions       = 7 ;
block.OutputPort(1).DatatypeID  = 0 ; % double
block.OutputPort(1).Complexity  = 'Real' ;

% Register parameters
block.NumDialogPrms     = 0;
Sample_Time = 0.05 ; % Change in S-function as well !!!!
% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [Sample_Time 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
% block.RegBlockMethod('Derivatives', @Derivatives);
% block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)

  block.NumDworks = 1;  
  block.Dwork(1).Name            = 'x0';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;


%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%% Initialize Dwork
%   block.Dwork(1).Data = block.InputPort(2).Data;%block.DialogPrm(1).Data;
%   clear getTime3   % To reset persistent variable !
%   clear DefineHypotheses3

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)

block.Dwork(1).Data = 0; %block.InputPort(2).Data; % 0
 clear getTime3   % To reset persistent variable !
 clear DefineHypotheses3



%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

NewPoint = [block.InputPort(1).Data; block.InputPort(2).Data; block.InputPort(3).Data] ;
% NewPoint = [ block.InputPort(1).Data; block.InputPort(2).Data ] ;

time3 = getTime3();
BestHypothesis_V3 = DefineHypotheses3( NewPoint, time3 ) ; 
save BestHypothesis_V3 BestHypothesis_V3

load  BestHypothesis_V1

d_ex_flag = 0 ;
surface_threshold = 1 ;
if block.InputPort(3).Data <= surface_threshold && block.InputPort(4).Data <= surface_threshold
    d_ex_flag = 1 ;   % data exchange flag
end

persistent gmm_Dist_3_1
if isempty(gmm_Dist_3_1)
    gmm_Dist_3_1 = 0 ;
else
    if  d_ex_flag && ~isempty(BestHypothesis_V1) && ~isempty( BestHypothesis_V3 )  
    gmm_Dist_3_1 = abs(Variational_distance(BestHypothesis_V3, BestHypothesis_V1)) ;
    end
end

if  ~isempty(BestHypothesis_V1) && ~isempty( BestHypothesis_V3 )
%      Just to see how it works, calculate two types of distances!
     Variational_distance_13= abs(Variational_distance(BestHypothesis_V1, BestHypothesis_V3)) ;
     Bhattacharyya_distance_31 = abs(Bhattacharyya_distance(BestHypothesis_V3, BestHypothesis_V1)) ;
     Bhattacharyya_distance_13 = abs(Bhattacharyya_distance(BestHypothesis_V1, BestHypothesis_V3)) ;
else 
    Variational_distance_13 = 0 ;
    Bhattacharyya_distance_31 = 0 ;
    Bhattacharyya_distance_13 = 0 ;    
end

block.OutputPort(1).Data = [ gmm_Dist_3_1, gmm_Dist_3_1, d_ex_flag,...
   BestHypothesis_V3.K,Variational_distance_13,Bhattacharyya_distance_31,Bhattacharyya_distance_13] ;

%end Outputs


%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)

block.Dwork(1).Data = block.InputPort(1).Data;

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

