function display_som(trained_som, hitdata, params)
% Displays SOM projection results using labels stored in trained_som.
%
% See also: label_som

% 08Sep2019 Petr Janata - adapted relevant secdtions of train_som_w_melody.m
labels = hitdata.labels;
[outdir,mapstub] = fileparts(params.somfname);

% Generate the output PostScript file name
figname = fullfile(outdir,sprintf('%s.ps',mapstub));
fprintf('Saving figures to %s\n', figname)

if params.add_hit_strength
  for ilabel = 1:length(labels)
    labels{ilabel} = sprintf('%s\n%3.0f', hitdata.labels{ilabel}, hitdata.hit_strength(ilabel));
  end
end

figure
som_show(trained_som,'empty','Best unit for each key')
som_show_add('label',labels)
if params.savefig, savefig(figname,''), end


figure
som_show(trained_som,'empty','Best unit for each key')
som_show_add('hit',hitdata.hits');
if params.savefig, savefig(figname,'-append'), end

figure
imagesc(hitdata.hits')
set(gca,'xtick',1:length(hitdata.categories),'xticklabel',hitdata.categories)
if params.savefig, savefig(figname,'-append'), end


return
end

function savefig(figname,append_str)
  print(figname,'-dpsc',append_str)
end