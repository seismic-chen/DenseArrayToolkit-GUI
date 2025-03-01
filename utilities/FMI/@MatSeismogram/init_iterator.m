function mySeis = init_iterator(mySeis, nPerRetrieve);

% INIT_ITERATOR Initialize iterator
%
% mySeis = init_iterator(mySeis, nPerRetrieve)
% Input arguments:
%    nPerRetrieve: number of seimograms per retrieve
%  
% Note: You have to initialize iterator before a multiple retrieving job.
%
% See alos enable_iterator, disable_iterator


mySeis.iterator.nTotal = 0;
mySeis.iterator.nPerRetrieve = nPerRetrieve;
mySeis.iterator.nRetrieved = 0;
mySeis.iterator.nRemained = 0;
mySeis.iterator.nIter = 0;
