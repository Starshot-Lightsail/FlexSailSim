function mystr = getdisptablerow(varlist, strlist)

% for i=1:length(varlist)
%     fprintf('%-12s\t',varlist{i});
% end
% for i=1:length(strlist)
%     fprintf('%-12s\t',strlist{i});
% end
% fprintf('\n');
mystr = '';
for i=1:length(varlist)
    mystr = [mystr sprintf('%-15g\t',evalin('base',varlist{i}))];
end
for i=1:length(strlist)
    mystr = [mystr sprintf('%-20s\t',evalin('base',strlist{i}))];
end

end