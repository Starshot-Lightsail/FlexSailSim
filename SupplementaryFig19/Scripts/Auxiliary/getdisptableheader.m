function mystr = getdisptableheader(varlist, strlist)

mystr = '';
for i=1:length(varlist)
    mystr = [mystr sprintf('%-15s\t',varlist{i})];
end
for i=1:length(strlist)
    mystr = [mystr sprintf('%-20s\t',strlist{i})];
end