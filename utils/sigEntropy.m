function ent = sig_entropy(inSig,params)
%
% accepts a one-dimensional signal and computes the entropy
%
% params.normalize: 0 or 1. Whether or not to normalize by the
%                   length of the signal
%
% 12/11/2008 First Version, Stefan Tomic
%


nSig = length(inSig);

%remove values == 0
inSig(inSig == 0) = [];

inSig = inSig./sum(inSig);

ent = -sum(inSig.*log2(inSig))./log2(nSig);



