function axh_return = get_axhandle(axh_array,tag)
% Returns an axes handle with a specified tag from a list of axes handles
%
% axh_return = get_axhandle(axh_array,tag)
%
% searches for an axes with a given tag (see axes tag property)
% from a list of axes specified by an axes handle array.
% This function returns -1 if an axes with the given tag was not found.
%
% INPUT
%  axh_array: an array of axes handles 
%  tag: the tag to search for
%
% OUTPUT
%  an axes handle for the axes with the given tag
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved 
%
% Author:
% Stefan Tomic 10/07

tagArray = get(axh_array,'tag');
axIdx = find(strcmp(tag,tagArray));

if(isempty(axIdx))
  axh_return = -1;
end

axh_return = axh_array(axIdx);
