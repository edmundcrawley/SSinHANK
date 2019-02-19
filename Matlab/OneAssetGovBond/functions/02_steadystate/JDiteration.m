function [joint_distr]=JDiteration(m_star,P_H,mpar,grid)
%% JDITERATION:
%
% Description: Iterates the joint distribution over m,k,h using a transition matrix
% obtained from the house distributing the households optimal choices. It
% distributes off-grid policies to the nearest on grid values.

% Copyright (c) 2014-02-28
% Christian Bayer, Ralph L�tticke, Lien Pham-Dao, and Volker Tjaden
% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Precautionary Savings, Illiquid Assets, and the Aggregate Consequences of
%  Shocks to Household Income Risk', Bonn mimeo
% http://wiwi.uni-bonn.de/hump/wp.html
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================


%% find next smallest on-grid value for money and capital choices

weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

% Adjustment case
[Dist_m,idm] = genweight(m_star,grid.m);

idm=repmat(idm(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));

for hh=1:mpar.nh
    
    %Corresponding weights
    weight11_aux = (1-Dist_m(:,hh));
    weight12_aux =  (Dist_m(:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P_H(hh,:);
end

weight11=permute(weight11,[1 3 2]);
weight12=permute(weight12,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);

H=sparse(rowindex,[index11(:); index12(:)],...
    [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), a without interest


%% Joint transition matrix and transitions


distJD=9999;
countJD=1;
[joint_distr,~]=eigs(H',1);
joint_distr=joint_distr'./sum(joint_distr);

while (distJD>1e-14 || countJD<50) && countJD<10000
    
    joint_distr_next=joint_distr*H;
    joint_distr_next=full(joint_distr_next);
    joint_distr_next=joint_distr_next./sum(joint_distr_next);
    
    distJD=max((abs(joint_distr_next(:)-joint_distr(:))));
    
    countJD=countJD+1;
    joint_distr=joint_distr_next;
    
end

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function
