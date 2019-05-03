function [ v_z_cc, v_zw_cc ] = fn_ww__ext__cc_weights( N1, a, b )
%fn_ww__ext__cc_weights: EXTERNAL Calc CC qudrature weights
%
%   [ v_z_cc, v_zw_cc  ] = fn_ww__ext__cc_weights( N1, a, b )
%
% Calculate Clenshaw-Curtis points and weights. Code is from Greg von
% Winckel, see license and original headers below.
%
% TAGS: EXT
%
% ORIGIN URL : https://se.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature
%
% LICENSE
%
% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
% ORIGINAL HEADER
%
% fclencurt.m - Fast Clenshaw Curtis Quadrature
%
% Compute the N nodes and weights for Clenshaw-Curtis
% Quadrature on the interval [a,b]. Unlike Gauss 
% quadratures, Clenshaw-Curtis is only exact for 
% polynomials up to order N, however, using the FFT
% algorithm, the weights and nodes are computed in linear
% time. This script will calculate for N=2^20+1 (1048577
% points) in about 5 seconds on a normal laptop computer.
%
% Written by: Greg von Winckel - 02/12/2005
% Contact: gregvw(at)chtm(dot)unm(dot)edu
%

N=N1-1; bma=b-a;
c=zeros(N1,2);
c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
v_zw_cc=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
v_z_cc=0.5*((b+a)+N*bma*f(1:N1,2));


end