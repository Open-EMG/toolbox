function [x_plus, A, b] = e_sifir_trn(Eamp, T, Param, Trim)
%E_SIFIR_TRN Identify (i.e., train) a non-linear FIR EMG-torque model.
%
% x_plus = e_sifir_trn(EMGamp, T, [Q D Tol ii], Trim)
%
% Uses regularized (via the pseudo-inverse technique) linear least squares
% [Press et al., 1994] to compute one or more sets of fit parameters of the
% non-linear FIR EMG-torque model specified by Koirala et al. [in review],
% including use of multiple training sets as well as extension to multiple
% inputs and outputs.  Can also be used for fitting of other linear dynamic
% models.  The model is of the form:
%  T(m+ii,co) = SUM_c_in(  SUM_d( SUM_q (
%           fit(q,d,ci,co) .* (Eamp(m-q,ci) .^ d)) )  ) + error,
% where:
%   m: (scalar integer) sample (time) index.  Not actually used within this
%      function, but useful in understanding the model.
%   ii: (scalar integer) number of samples into the future at which
%      the torque should be estimated.  Set ii=0 to estimate torque
%      at the current time.  Parameter ii cannot be negative.  Note that as
%      ii increases, the number of points available to the least squares
%      fit decreases.  See Koirala at al. [in review] for more details, but
%      EMG can be used to estimate torque into the future, with excellent
%      performance out to 60 ms and very good performance out to at least
%      100 ms.
%   ci: (scalar integer) input channel number.
%   Nci: (scalar integer) number of input channels.
%   co: (scalar integer) output channel number.
%   Nco: (scalar integer) number of output channels.
%   T(m) or T(m,co): If a vector, gives the output torque values
%      at times m for the one output.  A column vector is the
%      preferred convention.  If a matrix, each row (or column)
%      gives the torque values for each output (facilitating multiple-
%      output models).  In either case, the longer dimension is
%      taken as the time dimension; nonetheless, time indexed along
%      the column dimension is preferred.  Internally, always converted to
%      T(m,co) orientation and to a cell variable, which is used to
%      facilitate multiple recordings (see below).  If a cell array,
%      each element of T must have the same dimensions.  A cell array
%      format (see below) is used to facilitate multiple recordings.
%   D: (scalar integer) Max polynomial power.  D=1 for linear model.
%   Q: (scalar integer) Maximum time lag, in samples.  Q=0 ==> static
%      model.
%   Eamp(ci,m): Matrix of input EMGsigma estimates.  There is no distinction
%     needed to distinguish between extension- and flexion-oriented
%     Eamps.  Each row (or column) gives the estimate across all times.
%     The longer dimension is taken as the time dimension, the other
%     dimension is taken as the channel.  The length of the time
%     dimension must match the T time length.  Internally, always converted
%     to Eamp(ci,m) orientation and to a cell variable (see below).
%     Nonetheless, the preferred convention is for one channel per row
%     (time along the row dimension).  
%     Can also be a cell variable, which is used to facilitate multiple
%     recordings (see below).  If so, each cell element
%     must have the same dimensions.  T and Eamp should NOT have
%     already had start-up portions removed.
%
%  Tol: Tolerance of removal of singular values, based on the ratio of
%     the largest singular value of the design matix.
%  Param: A convenient vector to pack the values of [Q D Tol ii].
%  Trim: Two-element vector, giving the number of time samples to remove
%     from the beginning, Trim(1), and end, Trim(2), of Eamp and T
%     due to the various filter startup transients.  These samples are
%     removed prior to any fitting.  Note that the next Q values in T
%     are also unused, due to the inherent startup of the FIR filter being fit.
%     The first Q subsequent values in Eamp ARE used, since they are lag
%     times to the first-used sample time: m=(Q=1), assuming m=1 indexes the
%     first sample.  Trim values must be non-negative and cannot
%     leave the resulting Eamp and T lengths too small for proper fitting.
%
%  x_plus: Fit coefficients, listed in the normal format corresponding
%     to processing by the pseudo-inverse approach.  See e_sifir_x()
%     for a full description of the coefficient ordering.
%  A: Design matrix for the least squares fit.
%  b: Output matrix for the least squares fit.
%
% Multiple Training Recordings:
% Multiple training recordings can be combined for use in one parameter
% fit.  It is common to do so when recording trials are comprised
% of multiple components and/or multiple repetitions.
% A component refers to a trial that does not provide complete,
% stand-alone information for a parameter fit.  For example, if EMG is being
% related to two degrees of freedom, distinct recordings may have been
% made while exciting each degree of freedom, respectively.  Each of these two recordings
% has insufficient information alone to fit a complete model.  But, the recording
% pair forms a complete information set.  In this case, the pair of
% recordings must always be used together when training a model.
% A repetition refers to a repeated trial, generally with
% identical (or equivalent) recording conditions.  Any one of a set of repeated
% trials, or any combination of them, could be used to fit a model.  For this
% training routine, all supplied repetitions are used to fit the model.
%
% When utilizing multiple training recordings, inputs T and
% Eamp are supplied as cell arrays of the same dimensions.  The
% torque data from one cell is associated with the EMG amplitude data from
% the corresponding cell in Eamp.  Using cell indexing, these inputs
% are indexed as T{Cmp, Rep} and Eamp{Cmp, Rep}, where
% Cmp indexes the components and Rep indexes the
% repetitions.  Each element of T{} contains a vector/matrix as
% specified above, each of the same dimensions.  Each element must
% represent the same output channels.  Thus, if a recording only excites one
% degree of freedom, the user must still supply data for the unused
% degrees of freedom.  Presumably, null vectors would be supplied.
% Each element of Eamp{} contains a vector/matrix as specified
% above, each of the same dimensions.  Each element must represent
% the same input channels.
%
% Note that parameter fitting (as implemented in this function) does not
% distinguish between components and repetitions---all distinct trials
% are combined into one data set for fitting.  The testing function,
% e_sifir_tst, also does not distinguish between components and
% repetitions.  It simply tests all distinct trials.  However, cross-validation
% does distinguish between components and repetitions.
% Components are always grouped together, while repetitions are used
% to form the train-test combinations.
%
% EMG Amplitude Estimation Toolbox - Edward (Ted) A. Clancy - WPI

% Copyright (c) Edward A. Clancy, 2014.
% This work is licensed under the Aladdin free public license.
% For copying permissions see license.txt.
% email: ted@wpi.edu

% 30 December 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%% Process Command Line %%%%%%%%%%%%%%%%%%%%%%%%%%
% A few overall checks.
if nargin~=4, error('Need exactly 4 input arguments.'); end
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
Q = Param(1);  D = Param(2);  Tol = Param(3);  ii = Param(4);
if  Q<0, error('Q less than zero.');  end
if  D<1, error('D less than one.');   end
if ii<0, error('ii less than zero.'); end

% Trim.
if size(Trim,1)~=1 && size(Trim,2~=1), error('Trim not a vector.'); end
if length(Trim)~=2, error('length(Trim) ~= 2.'); end
if Trim(1)<0, error('Trim(1) must be non-negative.'); end
if Trim(2)<0, error('Trim(2) must be non-negative.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the fit. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim startup transients from the sequences.
for el=1:numel(Eamp), Eamp{el} = Eamp{el}(:,1+Trim(1):end-Trim(2)); end
for el=1:numel(T),    T{el}    = T{el}(1+Trim(1):end-Trim(2),:);    end

% Build the design matix A and output vector/matrix b.
[A, b] = e_sifir_ab(Eamp, T, Param);  % Build design matrix A, output b.

% Perform the least squares and compute the coefficients.
Aplus = pinv(A, Tol*norm(A));  % norm(A) ==> largest singular value of A.
x_plus = Aplus * b;  % Matrix multiply.  x_plus is (Q+1)*ci*D by co.

return
