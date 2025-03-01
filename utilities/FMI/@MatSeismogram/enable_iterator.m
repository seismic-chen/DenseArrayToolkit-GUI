function mySeis = enable_iterator(mySeis);

% ENABLE_ITERATOR Enable iterator 
%
% Note: The iterator is disabled initially when you new a MatSeismogram object. 
%
% See also disable_iterator, init_iterator

mySeis.iterator.enabled = 1;
