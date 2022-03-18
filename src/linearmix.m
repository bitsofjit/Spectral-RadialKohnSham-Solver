function [cdin] =  linearmix(cdin, cdout, linearmixing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% 
% Inputs
% ====== 
% cdin(:)         : input cd to current iteration
%                   overwritten by next input
% cdout(:)        : output cd from current iteration
% linearmixing    : mixing parameter - 0.1 means
%                   mix 10% or output with input
%                   Check that linear mixing parameter is between 0 and 1 
%                   and that input and output charge densities have the 
%                   same dimension

if(linearmixing <= 0.0 || linearmixing >= 1.0)
   error(' Error in linearMix: linearmixing out of range');
end

if (size(cdin) ~= size(cdout) )
   error(' error in linearMix: mixmatch in cd array dimensions');
end

% Now just mix them and we're done
cdin = cdin + linearmixing*(cdout - cdin);

end