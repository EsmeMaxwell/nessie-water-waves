function I = fn_ww__ext__dim__simpsons(f,dz)
%fn_ww__ext__dim__simpsons: EXTERNAL DIM calc Simpsons quadrature (from Juan Medina)
% 
%   [tc_kc,tc_EL] = fn_ww__ext__dim__calc_KC_EL(I,Ik,z,b_no_EL)
% 
% Calculate Simpson's rule integration. Used in Yan Li's DIM code but is
% from Juan Carlo Medina, origina license and headers included below.
%
%
% TAGS: EXT
% 
% ORIGIN URL : https://se.mathworks.com/matlabcentral/fileexchange/28726-simpson-s-rule-integration
%
% LICENSE
%
% Copyright (c) 2010, Juan Medina
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
% ORIGINAL HEADERS
%
% function I = simpsons(f,a,b,n)
% This function computes the integral "I" via Simpson's rule in the interval [a,b] with n+1 equally spaced points
% 
% Syntax: I = simpsons(f,a,b,n)
% 
% Where,
%  f= can be either an anonymous function (e.g. f=@(x) sin(x)) or a vector
%  containing equally spaced values of the function to be integrated
%  a= Initial point of interval
%  b= Last point of interval
%  n= # of sub-intervals (panels), must be integer
% 
%  Written by Juan Camilo Medina  - The University of Notre Dame
%  09/2010 (copyright Dr. Simpson)
% 
% 
% Example 1:
% 
% Suppose you want to integrate a function f(x) in the interval [-1,1].
% You also want 3 integration points (2 panels) evenly distributed through the
% domain (you can select more point for better accuracy).
% Thus:
%
% f=@(x) ((x-1).*x./2).*((x-1).*x./2);
% I=simpsons(f,-1,1,2)
% 
% 
% Example 2:
% 
% Suppose you want to integrate a function f(x) in the interval [-1,1].
% You know some values of the function f(x) between the given interval,
% those are fi= {1,0.518,0.230,0.078,0.014,0,0.006,0.014,0.014,0.006,0}
% Thus:
%
% fi= [1 0.518 0.230 0.078 0.014 0 0.006 0.014 0.014 0.006 0];
% I=simpsons(fi,-1,1,[])
%
% note that there is no need to provide the number of intervals (panels) "n",
% since they are implicitly specified by the number of elements in the
% vector fi
% 
% 
% 
% 
% See also
%   fn_ww__ext__dim__c()


% if numel(f)>1 % If the input provided is a vector
%     n=numel(f)-1; h=(b-a)/n;
%     I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end))+f(end));
% else % If the input provided is an anonymous function
%     h=(b-a)/n; xi=a:h:b;
%     I= h/3*(f(xi(1))+2*sum(f(xi(3:2:end-2)))+4*sum(f(xi(2:2:end)))+f(xi(end)));
% end

I= dz/3.*(f(1,:)+2*sum(f(3:2:end-2,:))+4*sum(f(2:2:end,:))+f(end,:)); 

%% test here
% n_iterval = 16;
% x = (0:1/n_iterval:1)';
% m_x = repelem(x,1,7);
% v_itg = m_x.^2;
% a = v_itg;
% sum1 = 1/n_iterval/3*(a(1,:)+2*sum(a(3:2:end-2,:))+4*sum(a(2:2:end,:))+a(end,:))
end