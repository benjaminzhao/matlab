function test()
%team_vol = 4;
timer = tic();
a = [0 0 0 2 5 6 7 0 0 0];
b = [1 8 5 6 4 7 0 0 0 0];
max = -1;
max_id = 0;
for i = 1:1:10
    found = 0;
    for j = 1:1:10
        if i == b(j)
            found = 1;
            break;
        end
    end
    if found == 0
        if a(i) > max
            max = a(i);
            max_id = i;
        end
    end
    toc(timer);
end
disp(max_id);
%%re-order big->small
% for i = 1:team_vol;
%             for j = i+1:(team_vol+1);
%                 if a(i, 1) < a(j, 1)
%                     temp = a(i, 1);
%                     a(i, 1) = a(j, 1);
%                     a(j, 1) = temp;
%                 end
%             end
% end
%disp(a);
% out = fopen('out.txt','w');
% for i = 1:5
%     fprintf(out,'%d\r\n',a(i));
% end
% fclose(out);
