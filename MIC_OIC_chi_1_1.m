function [MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_OIC_chi_1_1(data,B,c,n,threshold_a)
% OIC used the chi2test restrict segment by compare with two chi2test
% table (Local compare), the x_fix also use chi2test restrict (global compare)
% each table must campare with the 2*2 table; But all chi2test used adjust,
% include the df > 1
% n is the numbers of sample
% mutual_I_2 represent fixed the first column of data(x1)
% mutual_I_1 represent fixed the second column of data(x2)
% B=n^0.6
% c is max clumps (c*x)
%  see also MIC_2variable_chi.m : using Global chi2test compare
%           MIC_2variable_optimization2: using Local compare and Backtracking
%           MIC_2variable_optimization2_2: used the chi2test restrict segment by compare with two chi2test
%                table (Local compare),and replace log2(min(x,y)) with log2(x_fix)
%           MIC_OIC_chi_2: used the chi2test restrict segment by compare with two chi2test
%                        table (Local compare), Backtracking with Local
%                        compare
if nargin<5
    threshold_a=0.01;
end
mutual_I_2=zeros(round(B/2)-1,round(B/2)-1);
mutual_I_1=zeros(round(B/2)-1,round(B/2)-1);
mutual_chi_2=zeros(round(B/2)-1,round(B/2)-1);
mutual_chi_1=zeros(round(B/2)-1,round(B/2)-1);
last_c2=zeros(round(B/2)+1,round(B/2)+1);
last_c1=zeros(round(B/2)+1,round(B/2)+1);
% randnum=randperm(n)';
% data=data(randnum,:);
n=int32(n);
[value1,pos1]=sortrows(data,1);
[value2,pos2]=sortrows(data,2);
D1=nan(n,2);
D1=int32(D1);
for i=2:round(B/2)
    Q_x1=equipartitionYaxis2c(value1(:,1),int32(i),n);
    Q_x2=equipartitionYaxis2c(value2(:,2),int32(i),n);
    D1(pos1,1)=Q_x1;
    vector1=D1(pos2,1);
    avg=int32(n/(round(B/i)*c));
    c_x2=getsuper2var(vector1,avg,n);
    D1(pos2,2)=Q_x2;
    vector2=D1(pos1,2);
    c_x1=getsuper2var(vector2,avg,n);
    len1=int32(length(c_x1)-2);
    len2=int32(length(c_x2)-2);
    sub_max_seg=int32(round(B/i));
    chi2value=chi2inv(1-threshold_a,i-1);
    [temI2,tem_c2,chi2]=getmutualI2var_fix4(vector1,c_x2,sub_max_seg,int32(i),n,len2,chi2value);
    [temI1,tem_c1,chi1]=getmutualI2var_fix4(vector2,c_x1,sub_max_seg,int32(i),n,len1,chi2value);
    mutual_I_2(i-1,1:length(temI2))=temI2';
    mutual_chi_2(i-1,1:length(temI2))=chi2';
    last_c2(i-1,1:length(tem_c2))=tem_c2;
    mutual_I_1(1:length(temI1),i-1)=temI1;
    mutual_chi_1(1:length(temI1),i-1)=chi1;
    last_c1(1:length(tem_c1),i-1)=tem_c1';
end
% [ind_row,ind_col] = ind2sub(size(mutual_chi_2),find(mutual_chi_2~=0));
% pdist([ind_row,ind_col]);
for i=2:(round(B/2)-1)^2
    if mutual_chi_2(i)~=0
        [ind_row,ind_col] = ind2sub(size(mutual_chi_2),i);
%         if ind_row<=ind_col
            if (mutual_chi_2(i)-mutual_chi_2(1))<chi2inv(1-threshold_a,(ind_row)*(ind_col)-1)
                mutual_I_2(i)=0;
            end
%         else
%             mutual_I_2(i)=0;
%         end
    end
    if mutual_chi_1(i)~=0
        [ind_row,ind_col] = ind2sub(size(mutual_chi_1),i);
%         if ind_col<=ind_row
            if (mutual_chi_1(i)-mutual_chi_1(1))<chi2inv(1-threshold_a,(ind_row)*(ind_col)-1)
                mutual_I_1(i)=0;
            end
%         else
%             mutual_I_1(i)=0;
%         end
    end
end
MIC=max(max([mutual_I_1,mutual_I_2]));
position=find(mutual_I_2==MIC, 1);
if ~isempty(position)
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    mycertain=1;
    certain_seg=pos_row+1;
    bestc=last_c2(pos_row,2:pos_col+1);
    position2=find(mutual_I_1==MIC, 1);
    if ~isempty(position2)
        pos_row2=mod(position2(1)-1,round(B/2)-1)+1;
        pos_col2=floor((position2(1)-1)/(round(B/2)-1))+1;
        mycertain2=2;
        certain_seg2=pos_col2+1;
        bestc2=last_c1(2:pos_row2+1,pos_col2)';
        if certain_seg2*(length(bestc2)+1)<certain_seg*(length(bestc)+1)
            mycertain=mycertain2;
            certain_seg=certain_seg2;
            bestc=bestc2;
        end
    end
else
    position=find(mutual_I_1==MIC, 1);
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    mycertain=2;
    certain_seg=pos_col+1;
    bestc=last_c1(2:pos_row+1,pos_col)';
end
end
