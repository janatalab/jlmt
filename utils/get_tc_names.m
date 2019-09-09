function names = get_tc_names(times)
% Returns the names (strings) for time constants specified in times.  These
% are used to index into the .signals field of the li output structure.
%
% This function duplicates the functionality of the HalfDecayNameGen()
% subfunction in calc_li

% 07Sep2019 PJ

ntime = length(times);
names = cell(1,ntime);
for j=1:ntime
  t = times(j);
  if mod(t,fix(t))
    names{j} = regexprep(sprintf('tc_%1.4f',t),'\.','p');
  else
    names{j} = sprintf('tc_%d',t);
  end
end

return
 