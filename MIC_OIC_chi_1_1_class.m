function [MIC,mybestc]=MIC_OIC_chi_1_1_class(data,B,c,threshold_a)
% OIC used the chi2test restrict segment by compare with two chi2test
% table (Local compare), the x_fix also use chi2test restrict (global compare)
% each table must campare with the 2*2 table; But all chi2test used adjust,
% include the df > 1
% n is the numbers of sample
% mutual_I_2 represent fixed the first column of data(x1)
% mutual_I_1 represent fixed the second column of data(x2)
% B=0.6
% c is max clumps (c*x)
% y=data(:,1;)
if nargin<4
    threshold_a=0.01;
end
n=length(data(:,1));
class_num=length(unique(data(:,1)));
if size(data,2)<3
    max_seg=round(n^B/class_num);
    [~,I2]=sortrows(data,2);
    vector1=data(I2,1);
    vector1=int32(vector1);
    n=int32(n);
    avg=int32(n/(max_seg*c));
    best_c=getsuper2var(vector1,avg,n);
    len2=int32(length(best_c)-2);
    chi2value=chi2inv(1-threshold_a,class_num-1);
    [mutual_I_2,tem_c,~]=getmutualI2var_fix4(vector1,best_c,int32(max_seg),int32(class_num),n,len2,chi2value);
    [MIC,pos]=max(mutual_I_2);
    mybestc=tem_c(2:pos+1)';
end
end