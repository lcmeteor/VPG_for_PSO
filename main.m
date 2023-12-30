 clear
 clc
% mex cec13_func.cpp -DWINDOWS

%% Benchmark settings
func_num=28; %1-5umi,6-20mul,21-28com
Xmin=-100;
Xmax=100; 
foptimal=[-1400:100:-100 100:100:1400];
fhd=str2func('cec13_func');

runs=50; 

%% run algorithm
for D=30
    switch D
        case 10, pop_size=100; iter_max=1000;
        case 30, pop_size=100; iter_max=3000;
        case 50, pop_size=100; iter_max=5000;
        case 100, pop_size=100; iter_max=10000;
    end
    
    fbest=zeros(func_num,runs); %records the gbestval in each run on each function
    fb_mean=zeros(1,func_num); %records the mean value of all runs on each function
    fb_best=zeros(1,func_num); %records the best value of all runs on each function

    for i=1:20
        func_no=i;
        for j=1:runs
            T=cputime;
            [gbest,gbestval,gbestlog,xdlog,pdlog]=VPG_PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,foptimal(func_no),func_no);
            fbest(i,j)=gbestval;
            fprintf('func==%d, run==%d,gbestval==%g\n',i,j,gbestval);
        end
        fb_mean(i)=mean(fbest(i,:));
        fb_best(i)=min(fbest(i,:));
    end
end
  
