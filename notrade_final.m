%%%计算本城市群的贸易
clear all
clc

City_type=xlsread('城市群分类.xlsx','Sheet1','D2:D314');

%%%2012
load('2012_Chinese_City_MRIO.mat');
CO2_2012=xlsread('data_input','CO2_2012');
CO2_12=reshape(CO2_2012',[313*42,1]);
F_12=MRIO.F;
Z_12=MRIO.Z;
VA_12=MRIO.VA;
IM_12=MRIO.Import;
EX_12=MRIO.Export;
X_12=sum(F_12,2)+sum(Z_12,2)+EX_12;

intensity12=CO2_12./X_12;

intensity12(isnan(intensity12))=0;
intensity12(isinf(intensity12))=0;


for i=1:313

    F2_12(:,i)=sum(F_12(:,(i-1)*5+1:(i-1)*5+5),2);

end

%%城市出口计算

for i= 1:313

    for s= 1:42
      
     E12(i,s)=sum(Z_12((i-1)*42+s,:),2)+sum(F2_12((i-1)*42+s,:),2)-sum(Z_12((i-1)*42+s,(i-1)*42+1:i*42),2)-F2_12((i-1)*42+s,i);

    end
    
end


% 预处理City_type的唯一类型
unique_city_types = unique(City_type, 'rows'); 

for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_12((i-1)*42+s,(j-1)*42+1:j*42),2) + F2_12((i-1)*42+s,j);
                q = q + t;
            else
                t=0;
                q = q + t;
            end
        end
        E12_ua(i, s) = q;
    end
end
    

%%城市进口计算
for i= 1:313

    for s= 1:42
     
       M12(i,s)=sum(sum(Z_12(s:42:312*42+s,(i-1)*42+1:i*42)))-sum(sum(Z_12((i-1)*42+s,(i-1)*42+1:i*42)))+sum(F2_12(s:42:312*42+s,i),1)-F2_12((i-1)*42+s,i);

    end

end

for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_12((j-1)*42+s,(i-1)*42+1:i*42)) + F2_12((j-1)*42+s,i);
                q = q + t;
                            else
                t=0;
                q = q + t;
            end
        end
        M12_ua(i, s) = q;
    end
end



%%本地供给本地
for i= 1:313

    for j= 1:42

       EtoL12(i,j)=sum(Z_12((i-1)*42+j,(i-1)*42+1:i*42),2)+F2_12((i-1)*42+j,i);

    end
    
end

EtoL12_2=sum(EtoL12,2);

%%反事实情景2012城市间的

NE12=E12-M12;

for i= 1:313

  A=Z_12((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_12((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity12((i-1)*42+1:i*42))*L;

  Notrade_EX_sector_12(:,i)=Multiplier*NE12(i,:)';

end

Notrade_NE_12=sum(Notrade_EX_sector_12,1)';


%%反事实情景2012城市群之间的

NE_ua12=E12_ua-M12_ua;

for i= 1:313

  A=Z_12((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_12((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity12((i-1)*42+1:i*42))*L;

  Notrade_sector_12(:,i)=Multiplier*NE_ua12(i,:)';

end

Notrade_NE_ua_12=sum(Notrade_sector_12,1)';


%%%2015
load('313 MRIO_2015.mat');
CO2_2015=xlsread('data_input','CO2_2015');
CO2_15=reshape(CO2_2015',[313*42,1]);
F_15=MRIO.F;
Z_15=MRIO.Z;
VA_15=MRIO.VA;
IM_15=MRIO.IM;
EX_15=MRIO.EX;
X_15=sum(F_15,2)+sum(Z_15,2)+EX_15;

intensity15=CO2_15./X_15;

intensity15(isnan(intensity15))=0;
intensity15(isinf(intensity15))=0;

for i=1:313

    F2_15(:,i)=sum(F_15(:,(i-1)*5+1:(i-1)*5+5),2);

end

P = zeros(42, 313); %初始化P矩阵

%%城市群城市出口计算

for i= 1:313

    for j= 1:42

       E15(i,j)=sum(Z_15((i-1)*42+j,:),2)+sum(F2_15((i-1)*42+j,:),2)-sum(Z_15((i-1)*42+j,(i-1)*42+1:i*42),2)-F2_15((i-1)*42+j,i);

    end
end

%%城市群城市进口计算

for i= 1:313

    for s= 1:42
     
    M15(i,s)=sum(sum(Z_15(s:42:312*42+s,(i-1)*42+1:i*42)))-sum(sum(Z_15((i-1)*42+s,(i-1)*42+1:i*42)))+sum(F2_15(s:42:312*42+s,i),1)-F2_15((i-1)*42+s,i);

    end

end

%%2015城市群内出口
for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_15((i-1)*42+s,(j-1)*42+1:j*42),2) + F2_15((i-1)*42+s,j);
                q = q + t;
                            else
                t=0;
                q = q + t;
            end
        end
        E15_ua(i, s) = q;
    end
end
    
%%2015城市群内进口
for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_15((j-1)*42+s,(i-1)*42+1:i*42)) + F2_15((j-1)*42+s,i);
                q = q + t;
                            else
                t=0;
                q = q + t;
            end
        end
        M15_ua(i, s) = q;
    end
end


%%本地供给本地

for i= 1:313

    for j= 1:42

       EtoL15(i,j)=sum(Z_15((i-1)*42+j,(i-1)*42+1:i*42),2)+F2_15((i-1)*42+j,i);

    end
    
end

EtoL15_2=sum(EtoL15,2);

%%反事实情景2015城市间的

NE15=E15-M15;

for i= 1:313

  A=Z_15((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_15((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity15((i-1)*42+1:i*42))*L;

  Notrade_EX_sector_15(:,i)=Multiplier*NE15(i,:)';

end

Notrade_NE_15=sum(Notrade_EX_sector_15,1)';

%%反事实情景2015城市群之间的

NE_ua15=E15_ua-M15_ua;

for i= 1:313

   A=Z_15((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_15((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity15((i-1)*42+1:i*42))*L;

  Notrade_sector_15(:,i)=Multiplier*NE_ua15(i,:)';

end

Notrade_NE_ua_15=sum(Notrade_sector_15,1)';


%2017构建反事实情景，城市群城市不进口、不出口，

load('China city MRIO_2017.mat');
CO2_2017=xlsread('data_input','CO2_2017');
CO2_17=reshape(CO2_2017',[42*313,1]);

F_17=MRIO.F;
Z_17=MRIO.Z;
VA_17=MRIO.VA;
IM_17=MRIO.IM;
EX_17=MRIO.Export;

X_17=sum(F_17,2)+sum(Z_17,2)+EX_17;

F2_17=zeros(313*42,313);

intensity17=CO2_17./X_17;
intensity17(isnan(intensity17))=0;
intensity17(isinf(intensity17))=0;

for i=1:313

    F2_17(:,i)=sum(F_17(:,(i-1)*5+1:(i-1)*5+5),2);

end

P = zeros(42, 313); % 初始化 P 矩阵

%%城市群城市出口计算

for i= 1:313

    for j= 1:42

       E17(i,j)=sum(Z_17((i-1)*42+j,:),2)+sum(F2_17((i-1)*42+j,:),2)-sum(Z_17((i-1)*42+j,(i-1)*42+1:i*42),2)-F2_17((i-1)*42+j,i);

    end
end

%%城市群城市进口计算

for i= 1:313

    for s= 1:42
     
    M17(i,s)=sum(sum(Z_17(s:42:312*42+s,(i-1)*42+1:i*42)))-sum(sum(Z_17((i-1)*42+s,(i-1)*42+1:i*42)))+sum(F2_17(s:42:312*42+s,i),1)-F2_17((i-1)*42+s,i);

    end

end

%2017年城市群内出口
for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_17((i-1)*42+s,(j-1)*42+1:j*42)) + F2_17((i-1)*42+s,j);
                q = q + t;
                            else
                t=0;
                q = q + t;
            end
        end
        E17_ua(i, s) = q;
    end
end
    
%2017年城市群内进口
for i = 1:313
    m = City_type(i,:); % 当前城市类型
    % 找到所有与当前城市类型相同的城市索引
    same_type_cities = find(all(City_type == m, 2)); 

    for s = 1:42
        q = 0;
        for j = same_type_cities'
            if j ~= i
                t = sum(Z_17((j-1)*42+s,(i-1)*42+1:i*42)) + F2_17((j-1)*42+s,i);
                q = q + t;
                            else
                t=0;
                q = q + t;
            end
        end
        M17_ua(i, s) = q;
    end
end



%%本地供给本地

for i= 1:313

    for j= 1:42

       EtoL17(i,j)=sum(Z_17((i-1)*42+j,(i-1)*42+1:i*42),2)+F2_17((i-1)*42+j,i);

    end
    
end

EtoL17_2=sum(EtoL17,2);

%%反事实情景

NE17=E17-M17;

for i= 1:313

  A=Z_17((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_17((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity17((i-1)*42+1:i*42))*L;

  Notrade_EX_sector_17(:,i)=Multiplier*NE17(i,:)';

end

Notrade_NE_17=sum(Notrade_EX_sector_17,1)';


%%反事实情景2017城市群内的

NE_ua17=E17_ua-M17_ua;

for i= 1:313

   A=Z_17((i-1)*42+1:i*42,(i-1)*42+1:i*42)./X_17((i-1)*42+1:i*42)';

  A(isnan(A))=0;
  A(isinf(A))=0;

  I=eye(size(A));

  L=(I-A)^-1;

  Multiplier=diag(intensity17((i-1)*42+1:i*42))*L;

  Notrade_sector_17(:,i)=Multiplier*NE_ua17(i,:)';

end

Notrade_NE_ua_17=sum(Notrade_sector_17,1)';


%% ==== 聚合城市层面结果为城市群 ====

for i = 1:17
    idx = find(City_type == i);
    Notrade1_sum_au_12(i) = sum(Notrade_NE_12(idx));
    Notrade1_sum_au_15(i) = sum(Notrade_NE_15(idx));
    Notrade1_sum_au_17(i) = sum(Notrade_NE_17(idx));

    Notrade2_sum_au_12(i) = sum(Notrade_NE_ua_12(idx));
    Notrade2_sum_au_15(i) = sum(Notrade_NE_ua_15(idx));
    Notrade2_sum_au_17(i) = sum(Notrade_NE_ua_17(idx));

    t_NE12_ua=sum(NE_ua12,2);
    t_NE15_ua=sum(NE_ua15,2);
    t_NE17_ua=sum(NE_ua17,2);
    NE2_ua12(i)=sum(t_NE12_ua(idx));
    NE2_ua15(i)=sum(t_NE15_ua(idx));
    NE2_ua17(i)=sum(t_NE17_ua(idx));

    t_E12=sum(E12,2);
    t_E15=sum(E15,2);
    t_E17=sum(E17,2);
    T2_E12(i)=sum(t_E12(idx));
    T2_E15(i)=sum(t_E15(idx));
    T2_E17(i)=sum(t_E17(idx));
    
    t_M12=sum(M12,2);
    t_M15=sum(M15,2);
    t_M17=sum(M17,2);
    T2_M12(i)=sum(t_M12(idx));
    T2_M15(i)=sum(t_M15(idx));
    T2_M17(i)=sum(t_M17(idx));

    t_E12_ua=sum(E12_ua,2);
    t_E15_ua=sum(E15_ua,2);
    t_E17_ua=sum(E17_ua,2);
    E2_ua12(i)=sum(t_E12_ua(idx));
    E2_ua15(i)=sum(t_E15_ua(idx));
    E2_ua17(i)=sum(t_E17_ua(idx));
    
    t_M12_ua=sum(M12_ua,2);
    t_M15_ua=sum(M15_ua,2);
    t_M17_ua=sum(M17_ua,2);
    M2_ua12(i)=sum(t_M12_ua(idx));
    M2_ua15(i)=sum(t_M15_ua(idx));
    M2_ua17(i)=sum(t_M17_ua(idx));

end

Total_notrade=[Notrade1_sum_au_12',Notrade1_sum_au_15',Notrade1_sum_au_17'];
intra_ua_natrade=[Notrade2_sum_au_12',Notrade2_sum_au_15',Notrade2_sum_au_17'];
inter_ua_natrade=aa-bb;
dd=[NE2_ua12',NE2_ua15',NE2_ua17'];

ee=[E2_ua12',E2_ua15',E2_ua17'];
mm=[M2_ua12',M2_ua15',M2_ua17'];

trade_ua_import=[T2_M12',T2_M15',T2_M17'];
trade_ua_export=[T2_E12',T2_E15',T2_E17'];
