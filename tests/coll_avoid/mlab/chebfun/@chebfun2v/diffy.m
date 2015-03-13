function G = diffy( F, n )
%DIFFY   Differentiate a CHEBFUN2V with respect to its first argument
%   DIFFY(F) returns a CHEBFUN2V representing the derivative of F in its first
%   argument. This is the same as DIFF(F,1,2).
%
%   DIFFY(F,N) returns a CHEBFUN2V representing the Nth derivative of F in its
%   first argument. This is the same as DIFF(F,N,2).
%
%   This command is for convenience as the syntax for DIFF, inherited from the
%   DIFF command for matrices, can be confusing.
% 
% See also DIFFX, DIFF. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to first derivative. 
    n = 1; 
end

G = diff( F, n, 1 );

end
