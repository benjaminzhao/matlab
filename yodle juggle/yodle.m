%% yodle programming test
%% juggler and circuits
function result = yodle()
timer = tic;

c_num = 0; % cline number = team number
j_num = 0; % jline number = player number in total
element_num = 3;
j_assign_num = 10; % each j assgin to c number


C_H = zeros();
C_E = zeros();
C_P = zeros();

J_H = zeros();
J_E = zeros();
J_P = zeros();
J_C = zeros();

%% input
%Open file for reading
testin = fopen('jugglefest.txt','r');
if testin < 0
    error('input file open error!');
end

while ~feof(testin)
    tline = fgetl(testin);
    
    if isempty(tline) || tline(1) == ' '
        continue;%skip empty line
    elseif tline(1) == 'C' || tline(1) == 'c'
        % this is a C line
        c_num = c_num + 1;
        temp = sscanf(tline, 'C C%d H:%d E:%d P:%d');
        temp(1) = temp(1) + 1;% C_ID
        if element_num == 3
            C_H(temp(1)) = temp(2); %col increase
            C_E(temp(1)) = temp(3);
            C_P(temp(1)) = temp(4);
        end
    elseif tline(1) == 'J' || tline(1) == 'j';
        % this is a J line
        j_num = j_num + 1;
        if j_assign_num == 3
            temp = sscanf(tline, 'J J%d H:%d E:%d P:%d C%d,C%d,C%d');
        elseif j_assign_num == 10
            temp = sscanf(tline, 'J J%d H:%d E:%d P:%d C%d,C%d,C%d,C%d,C%d,C%d,C%d,C%d,C%d,C%d');
        end
        temp(1) = temp(1) + 1; % J_ID
        if element_num == 3
            J_H(temp(1)) = temp(2); % col increase
            J_E(temp(1)) = temp(3);
            J_P(temp(1)) = temp(4);
        end
        J_C(temp(1),1) = temp(5)+1;
        J_C(temp(1),2) = temp(6)+1;
        J_C(temp(1),3) = temp(7)+1;
        if j_assign_num == 10
            J_C(temp(1),4) = temp(8)+1;
            J_C(temp(1),5) = temp(9)+1;
            J_C(temp(1),6) = temp(10)+1;
            J_C(temp(1),7) = temp(11)+1;
            J_C(temp(1),8) = temp(12)+1;
            J_C(temp(1),9) = temp(13)+1;
            J_C(temp(1),10) = temp(14)+1;
        end
    end
end
fclose(testin);
% disp(C_H);
% disp(C_E);
% disp(C_P);
% disp(J_H);
% disp(J_E);
% disp(J_P);

%% calc
team_num = c_num;
team_vol = j_num / c_num; % team volume = player number in one team

C = [C_H;C_E;C_P];
%disp(C);
J = [J_H',J_E',J_P'];
%disp(J);
%disp(J_C);

%J result = J dot product C
temp = J * C;
J_R = [temp,(1:j_num)'];
J_C_ASS_2 = zeros(j_num, c_num-j_assign_num);
%disp(J_R);

for j = 1:j_num
    j_temp = zeros(2, j_num);
    j_temp(1) = temp(j);
    j_temp(2,:) = 1:j_num;
    for i = 1:j_assign_num
        j_temp(1, J_C(i)) = NaN;
        j_temp(2, J_C(i)) = NaN;
    end
    m = 1;
    jtemp = zeros(2, c_num-j_assign_num);
    for i = 1:j_num
        if j_temp(1, i) ~= NaN
          jtemp(1,m) = j_temp(1,i);
          jtemp(2,m) = j_temp(2,i);
          m = m + 1;
        end
    end
    
    for m = 1:c_num-j_assign_num-1
        for n = m:c_num-j_assign_num
            if jtemp(1, m) < jtemp(1,n)
                t = j_temp(1, m);
                jtemp(1, m) = jtemp(1, n);
                j_temp(1, n) = t;
                t = jtemp(2, m);
                jtemp(2, m) = jtemp(2, n);
                jtemp(2, n) = t;
            end
        end
    end
    J_C_ASS_2(j, 1:c_num-j_assign_num) = jtemp(2,1:c_num-j_assign_num);
    disp(j);
    toc(timer);
end

J_C_ASS = zeros(j_num, c_num+1); % save j_id poped out c#
zr = ones(1,j_num);
J_C_ASS = [J_C, J_C_ASS_2, zr'];


%% proc
C_R = zeros(team_num, team_vol+1);

%row = j_num;
[row,col] = size(J_R);

for index_1 = 1:1:row
    J_ID = index_1;
    j_order_id = J_C_ASS(J_ID, c_num+1);

    % get the 1st J's C_ID
    J_C_ID = J_C_ASS(J_ID, j_order_id);
    
    j_order_id = j_order_id + 1;
    J_C_ASS(J_ID, c_num+1) = j_order_id;
    
%     if J_C_ID > c_num || J_C_ID < 1
%         error('J_C ID error!');
%     end
	C_R(J_C_ID, team_vol+1) = J_R(J_ID, c_num+1);% save ID to C_R buffer pos    
    
    % do re-order, big-->small
    for i = 1:1:team_vol
        for j = (i+1):1:(team_vol+1)
            if C_R(J_C_ID, i) == 0 && C_R(J_C_ID, j) ~= 0
                temp = C_R(J_C_ID, i);
                C_R(J_C_ID, i) = C_R(J_C_ID, j);
                C_R(J_C_ID, j) = temp;  
            elseif C_R(J_C_ID, i) ~= 0 && C_R(J_C_ID, j) ~= 0
                if J_R(C_R(J_C_ID, i), J_C_ID) < J_R(C_R(J_C_ID, j), J_C_ID)
                    temp = C_R(J_C_ID, i);
                    C_R(J_C_ID, i) = C_R(J_C_ID, j);
                    C_R(J_C_ID, j) = temp;
                end
            end
        end
    end 

    %if pop out a J
    while C_R(J_C_ID, team_vol+1) ~= 0
        
        % a J pop out from team, get its J_ID
        J_ID = C_R(J_C_ID, team_vol+1);
    
        j_order_id = J_C_ASS(J_ID, c_num+1);
        % get the J's C_ID
        J_C_ID = J_C_ASS(J_ID, j_order_id);
        
        j_order_id = j_order_id + 1;
        J_C_ASS(J_ID, c_num+1) = j_order_id;

        % save ID to C_R buffer pos
        C_R(J_C_ID, team_vol+1) = J_R(J_ID, c_num+1);

        % do re-order, big-->small
        for i = 1:team_vol
            for j = i+1:(team_vol+1)
                % if empty, fill it
                if C_R(J_C_ID, i) == 0 && C_R(J_C_ID, j) ~= 0
                    temp = C_R(J_C_ID, i);
                    C_R(J_C_ID, i) = C_R(J_C_ID, j);
                    C_R(J_C_ID, j) = temp;
                %if not empty, 
                elseif C_R(J_C_ID, i) ~= 0 && C_R(J_C_ID, j) ~= 0
                    if J_R(C_R(J_C_ID, i), J_C_ID) < J_R(C_R(J_C_ID, j), J_C_ID)
                        temp = C_R(J_C_ID, i);
                        C_R(J_C_ID, i) = C_R(J_C_ID, j);
                        C_R(J_C_ID, j) = temp;
                    end
                end
            end
        end

    end %% while
    disp(index_1);
    toc(timer);
end




%% out
% Open or create new file for reading and writing.
% Append data to the end of the file. 
testout = fopen('testout.txt','w');
if testout < 0
    error('output file open error!');
end
for i = 1:1:c_num
    fprintf(testout,'C%d ',i-1);
    for j = 1:1:team_vol-1
        fprintf(testout,'J%d ',C_R(i,j)-1);
        for k = 1:1:j_assign_num-1
            fprintf(testout,'C%d:%d ',J_C_ASS(C_R(i,j),k)-1,    J_R( C_R(i,j), J_C_ASS(C_R(i,j),k)) );
        end
        k = j_assign_num;
        fprintf(testout,'C%d:%d,',J_C_ASS(C_R(i,j),k)-1, J_R(C_R(i,j), J_C_ASS(C_R(i,j),k)) );
    end
    j = team_vol;
    fprintf(testout,'J%d ',C_R(i,j)-1);
    for k = 1:1:j_assign_num-1
        fprintf(testout,'C%d:%d ',J_C_ASS(C_R(i,j),k)-1, J_R(C_R(i,j), J_C_ASS(C_R(i,j),k)) );
    end
    k = j_assign_num;
    fprintf(testout,'C%d:%d\r\n',J_C_ASS(C_R(i,j),k)-1, J_R(C_R(i,j), J_C_ASS(C_R(i,j),k)) );
end
fclose(testout);

%%result = sum( C_R(1970+1,1:team_vol) ) - 6;

