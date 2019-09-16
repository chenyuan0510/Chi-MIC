# Chi-MIC
A new method, ChiMIC, to calculate the MIC values

1） MATLAB is the tool of NDC;

2)  Make sure that c++ has installed in your computer;

3)  Running program “make.m” to compile the equipartitionYaxis2.c, getsuper2var.c and getmutualI2var_fix4.c to mex files;

4)  num=randperm(size(data,1)); 

    data=data(num',:);# scramble the samples
    
    the first column of data is Y (dependent variable), the rest of the columns (X) independent variable;
    
    # while Y is numerical data
    
    [MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_OIC_chi_1_1(data,B,c,sample_num);
    
    # sample_num=size(data,1);
    
    #[MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_OIC_chi_1_1(data,sample_num^0.55,5,sample_num);
    
    # while Y is numerical data discrete data
    
    [MIC,~]=MIC_OIC_chi_1_1_class(data,0.55,5)
    
