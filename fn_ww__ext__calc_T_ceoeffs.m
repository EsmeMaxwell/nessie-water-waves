function [ B ] = fn_ww__ext__calc_T_ceoeffs( a_w )
%fn_ww__ext__calc_T_ceoeffs: EXTERNAL Calc Chebyshev coefficients
%
% [ B ] = fn_ww__ext__calc_T_ceoeffs( a_w )
%
% Takes array of eigenvectors and calculates the spectral coefficients.
% External library code by Greg von Winckel. See license and original
% headers below.
%
%
% ORIGIN URL : https://se.mathworks.com/matlabcentral/fileexchange/4591-fast-chebyshev-transform-1d?s_tid=prof_contriblnk
%
% TAGS : EXT
%
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
% 
% 
% ORIGINAL HEADERS
% 
% Fast Chebyshev Transform 
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
% direction - set equal to 1 for nodal to spectral
%             anything else for spectral to nodal
%
% Written by Greg von Winckel 03/08/04  
% Contact: gregvw@chtm.unm.edu


direction = 1;


[N,M]=size(a_w);

if direction==1 % Nodal-to-spectral
    F=ifft([a_w(1:N,:);a_w(N-1:-1:2,:)]);
    B=([F(1,:); 2*F(2:(N-1),:); F(N,:)]);
else            % Spectral-to-nodal
    F=fft([a_w(1,:); [a_w(2:N,:);a_w(N-1:-1:2,:)]/2]);
    B=(F(1:N,:));
end


end