function [A, b] = e_sifir_ab(Eamp, T, Param)
%E_SIFIR_AB Assemble design (A) and output (b) matrices for non-linear FIR EMG-torque model.
%
% [A, b] = e_sifir_ab(EMGamp, T, [Q D Tol ii])
%
% Assembles design matrix "A" and output vector/matrix "b" for the
% non-linear FIR EMG-torque model.  Both the training and testing
% functions for this model must assemble these matrices, thus the
% assembly code is modularized to this one function.  Inputs Eamp(ci,m),
% T(m,co) and Param are as defined in e_sifir_trn() and e_sifir_tst()  --
% as cell arrays.  It is assumed that these inputs have already been
% error checked.  Having one function to assemble these matrices helps
% to ensure uniformity in the ordering of information within these
% matrices.  The ordering of training and testing must match.
%
% For a single input channel, D=1, ii=0 and only one set, the design matrix
% is of the form (Q=3 is used for this example):
%
%         | Eamp(4)      Eamp(3)     Eamp(2)     Eamp(1)   |
%         | Eamp(5)      Eamp(4)     Eamp(3)     Eamp(2)   |
%   A11 = | Eamp(6)      Eamp(5)     Eamp(4)     Eamp(3)   |
%         |                  ...                           |
%         | Eamp(Nr)    Eamp(Nr-1)  Eamp(Nr-2)  Eamp(Nr-3) |
%
% If two input channels are available, then the second channel is
% concatenated to the right of the first.  Using partitioned matrix
% notation, for two channels, construct A11 and A12 then combine as:
%
%   A_{Two Channels} = | A11 A12 |
% 
% If a second data set/recording (or more, but two will be used in this
% example) is used, they are combined as additional rows.  Let the two
% channels of the second data set have their individual "A" matrices be
% denoted A21 and A22.  Then,
%
%   A_{Two Sets, Two Channels} = | A11 A12 |
%                                | A21 A22 |
% Note that additional recordings may be either "components" (arranged
% along columns of cell arrays Eamp and T) or "repetitions" (arranged
% along rows of Eamp and T).  The additional recordings are added in the
% order specified by unitary indexing through Eamp and T.  By MATLAB
% convention, therefore, components are added first, then repetitions.
%
% Next, some models will include terms that raise the EMGamp to a power.
% Those terms are appended to the right.  For a model that includes
% second-degree polynomical terms, we would have:
%
%   A = | A11  A12  A11.^2  A12.^2|
%       | A21  A22  A21.^2  A22.^2|
%
% The output "b" will be a vector when there is only one channel, but a
% matrix when there are additional channels (one additional column per
% additional channel).  Note that the first Q output values are not used,
% due to the start-up  transient of the dynamic model.
%      If T1 is a column vector of the first output channel and T2 is a
% column vector of the second output channel, then a two-channel "b"
% matrix would be assembled as (still assuming ii=0):
%
%  b = | T1(Q+1)  T2(Q+1) |
%      | T1(Q+2)  T2(Q+2) |
%      | T1(Q+3)  T2(Q+3) |
%               ...
%      | T1(Nr)   T2(Nr)  |
%
% Finally, if ii>0, then EMG is being used to estimate torque at
% future times.  Thus, the first ii times (i.e., the first ii rows)
% within the torque vector/matrix must be removed.  In that way,
% all times are advanced by ii samples.  Accordingly, the last
% ii rows of the design matrix "A" must also be removed.
%
% EMG Amplitude Estimation Toolbox - Edward (Ted) A. Clancy - WPI

% Copyright (c) Edward A. Clancy, 2014.
% This work is licensed under the Aladdin free public license.
% For copying permissions see license.txt.
% email: ted@wpi.edu

% 30 December 2014.

%%%%%%%%%%%%%%%%%%%%%%% Build the design matrix A. %%%%%%%%%%%%%%%%%%%%%%%
Q = Param(1);  D = Param(2);  ii = Param(4);
Nci = size(Eamp{1},1);  % Number of input channels.
Nco = size(T{1}, 2);    % Number of output channels.
Ns  = numel(T);         % Number of data sets (a.k.a., recordings).
Nr  = length(T{1});     % Length of one input or output set.

A = zeros( (Nr-Q-ii)*Ns, (Q+1)*Nci*D );  % Full design matrix.

% Insert design matrix values for polynomial power equal to one.
for el=1:Ns                    % Loop over data sets.
  Row = 1 + (el-1)*(Nr-Q-ii);  % Reset row index.
  Col = 1;                     % Reset col index.
  for ci=1:Nci                 % Loop over input channels.
    for qq=0:Q                 % Loop over the lags.
      % Assemble one column for D=1, this input channel, this set.
      A(Row:Row+Nr-Q-ii-1, Col+qq) = Eamp{el}(ci, Q+1-qq:Nr-qq-ii)';
    end
    Col = Col + Q+1;           % Increment for next input channel.
  end
end

% Insert design matrix polynomial powers above 1.
for dd=2:D
  A(:,(Q+1)*Nci*(dd-1)+1:(Q+1)*Nci*dd) = A(:,1:(Q+1)*Nci) .^ dd;
end

% Assemble the output matrix.
b = zeros((Nr-Q-ii)*Ns,Nco);  % Pre-allocate.
for el=1:Ns                   % Assemble, removing start-up.
  b( 1+(el-1)*(Nr-Q-ii):el*(Nr-Q-ii), : ) = T{el}(Q+1+ii:end,:);
end

return
