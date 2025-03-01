function mySeis = disable_iterator(mySeis);

% DISABLE_ITERATOR Disable iterator 
%
% Note: The iterator is disabled initially when you new a MatSeismogram object. 
%
% See also enable_iterator, init_iterator

mySeis.iterator.enabled = 0;

