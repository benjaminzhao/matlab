%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2013-05-20
% find the index of string in name_array
% return: index if successful
%         error if failed
% ___________________________________________________
%%
function X = find_i(X_name, X_str) 

X = 0;
for i = 1:length(X_name);
    if isequal(X_str, X_name{i})
        X = i;
        break
    end
end
if X == 0
    error('name unmatched');
end







