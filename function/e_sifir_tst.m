function [Err, T_estimate] = e_sifir_tst(x_plus, Eamp, T, Param, Trim, varargin)
%E_SIFIR_TST Test a non-linear FIR EMG-torque model.
%
% [Err T_estimate] = e_sifir_tst(x_plus, Eamp, T, [Q D Tol ii], Trim)
%
% Computes the estimated torque(s) from a previously fit non-linear FIR
% EMG-torque model [presumably fit using function e_sifir_trn()]
% and the error(s) between the estimated and actual
% torque(s).  The default error measure is the RMS.  For multiple-output
% systems, errors may be computed along each output component (default)
% or in a distance sense. Full definitions of the required function
% inputs are provided in e_sifir_trn().  x_plus holds
% the model coefficients generated by e_sifir_trn().  Eamp and T are
% the data from the test trial(s), and Param and Trim are the same
% parameters used in training.  The Tol value of Param(3) is not used.
%
% If Eamp and T are NOT cell arrays,
% this function computes the EMG-predicted torque, T_estimate, using
% the entire USABLE sample range of Eamp.  This range will omit the
% first Q torque values (due to the lags needed in the FIR model) and
% ii additional values due to estimation into the future.  However, it is
% useful for T_estimate to have the same length as T.  Thus, the first
% useable T value is pre-pended (Q+ii) times so that T_estimate is
% returned to the caller with the same length as T, and so that each
% sample in T_estimate is properly time-aligned with samples in Eamp
% and T.
%
% After T_estimate is computed, the error (default = RMS) between it and T
% is computed [ignoring the first Q+ii+Trim(1) values and also the last
% Trim(2) values].  This error value is returned as the scalar Err,
% and T_estimate is returned as a vector. When T is an array (i.e., when
% there are multiple outputs), the error can be computed in two alternative
% manners.  By default, error is computed along each output component.
% If variable argument 'Edist' is set to 'distance, then a distance error
% is first computed at
% each time instant (treating each output as a dimension in the distance
% computation) and then the RMS or MAV is computed across time.
%
% When Eamp and T ARE cell arrays, then output T_estimate will be
% a cell array of the same dimensions.  Each element of T_estimate
% holds the corresponding EMG-predicted torque associated with the
% corresponding Eamp{} and T{} elements.  Output Err will be a regular
% matrix of the same dimensions as Eamp and T.  Each element of vector Err
% holds the error associated with the corresponding Eamp{}
% and T{} elements.
% 
% EMG Amplitude Estimation Toolbox - Edward (Ted) A. Clancy - WPI

% Copyright (c) Edward A. Clancy, 2015.
% This work is licensed under the Aladdin free public license.
% For copying permissions see license.txt.
% email: ted@wpi.edu

% 22 June 2015.

%%%%%%%%%%%%%%%%%%%%%%%%%% Process Command Line %%%%%%%%%%%%%%%%%%%%%%%%%%
% A few overall checks.
if nargin<5, error('Need >= 5 input arguments.'); end
if iscell(Eamp) == 0, Eamp = {Eamp}; end  % Coerce to cell array.
if iscell(T)    == 0, T    = {T};    end  % Coerce to cell array.
if sum(size(T)==size(Eamp))~=2, error('T and Eamp dimensions differ.'); end

% Eamp and T elements.
for el=2:numel(Eamp)  % Are all elements the same size?
  if sum(size(Eamp{1})==size(Eamp{el}))~=2, error('{Eamp} sizes differ.'); end
  if sum(size(T{1})   ==size(T{el}))   ~=2, error('{T} sizes differ.');    end
end
if length(Eamp{1})~=length(T{1}), error('Eamp and T time durations differ.'); end
for el=1:numel(Eamp)  % Coerce each Eamp as Eamp(ci,m); T as T(m,co).
  if size(Eamp{el},1) > size(Eamp{el},2), Eamp{el} = Eamp{el}'; end
  if size(T{el},1)    < size(T{el},2),    T{el}    = T{el}';    end
end

% Parameters.
if size(Param,1)~=1 && size(Param,2~=1), error('Param not a vector.'); end
if length(Param)~=4, error('length(Param) ~= 4.'); end
Q = Param(1);  D = Param(2);  ii = Param(4);
if  Q<0, error('Q less than zero.');  end
if  D<1, error('D less than one.');   end
if ii<0, error('ii less than zero.'); end

% Trim.
if size(Trim,1)~=1 && size(Trim,2~=1), error('Trim not a vector.'); end
if length(Trim)~=2, error('length(Trim) ~= 2.'); end
if Trim(1)<0, error('Trim(1) must be non-negative.'); end
if Trim(2)<0, error('Trim(2) must be non-negative.'); end

% Process optional commands, if any.
if rem(length(varargin),2)==1, error('Need even number of optional args.'); end
Edist = 'Component'; % Default method for multiple-output data.
Emeth = 'RMS';  % Default error computation method.

for i = 1:2:length(varargin)
  switch varargin{i}
      
    case 'Edist'  % Specify the method for multiple outputs.
      switch varargin{i+1}
        case {'Component', 'Distance'}, Edist = varargin{i+1};
        otherwise, error(['Bogus "Edist": "' varargin{i+1} '".'])
      end
    
    case 'Emeth'  % Specify the method for error computation.
      switch varargin{i+1}
        case {'MAV', 'RMS'}, Emeth = varargin{i+1};
        otherwise, error(['Bogus "Emeth": "' varargin{i+1} '".'])
      end
    
    otherwise, error(['Bogus PropertyName "' varargin{i} '".'])
      
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_estimate = cell(size(Eamp));  % Pre-allocate.
Err = zeros(size(Eamp));  % Pre-allocate.

for m = 1:numel(Eamp)  % Loop over the test data sets.
  % Build the test data design matix A using ALL available data.
  [A, b1] = e_sifir_ab(Eamp(m), T(m), Param);  % A  is (Nr-Q-ii) by (Q+1)*ci*D.
                                               % b1 is (Nr-Q-ii)*Ns by Nco.
  
  % Estimate last (Nr-Q-ii) torques.
  T_estimate{m} = A*x_plus;  % T_estimate is (Nr-Q-ii)*Ns by co.

  % Pre-pend Q+ii torques ==> T, T_estimate are aligned + have same length.
  T_estimate{m} = [ones(Q+ii,1)*T_estimate{m}(1,:); T_estimate{m}];
  b = [ones(Q+ii,1)*b1(1,:); b1];

  % Compute the desired error.
  Range = 1+Q+ii+Trim(1) : length(T{1})-Trim(2); % Computation index range.
  switch Edist  % Component vs. distance error for multiple channels.
    case 'Component'
      switch Emeth
        case 'MAV'
          Err(m) = mean(mean( abs(T_estimate{m}(Range,:)-b(Range,:)) ));
        case 'RMS'
          Err(m) = sqrt(  mean(mean( (T_estimate{m}(Range,:)-b(Range,:)).^2 ))  );
      end
    case 'Distance'
      E2 = (T_estimate{m}(Range,:)-b(Range,:)) .^ 2; % Array/vector: squared errors.
      if size(E2,2)>1, E2 = sum(E2,2); end  % Sum if there are multiple outputs.
      switch Emeth  % E2 is now a row vector of squared distance errors.
        case 'MAV',  Err(m) = mean(sqrt(E2));
        case 'RMS',  Err(m) = sqrt(mean(E2));
      end
  end
end

if length(Eamp)==1, T_estimate = T_estimate{1}; end  % If one test set, no cell.

return

