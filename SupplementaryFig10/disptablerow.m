function disptablerow(varlist, strlist)

for i=1:length(varlist)
    fprintf('%-12s\t',varlist{i});
end
for i=1:length(strlist)
    fprintf('%-20s\t',strlist{i});
end
fprintf('\n');
for i=1:length(varlist)
    fprintf('%12g\t',evalin('base',varlist{i}));
end
for i=1:length(strlist)
    fprintf('%-20s\t',evalin('base',strlist{i}));
end
fprintf('\n');

end