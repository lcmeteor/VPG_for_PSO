function [gbest,gbestval,gbestlog,xdlog,pdlog]= VPG_QPSO_func(fhd,Dimension,Particle_Number,Max_Gen,Rmin,Rmax,foptimal,varargin)
% QPSO enabled by the VPG strategy
% Reference: J. Sun, W. Fang, X. Wu, V. Palade, W. Xu, Quantum-behaved particle swarm optimization: analysis of individual particle behavior and parameter selection, Evolutionary Computation 20.3 (2012) 349-393.

%% initilization
Me=Max_Gen; 
Ps=Particle_Number;
D=Dimension; 
gbestlog=zeros(1,Me); 
xdlog=zeros(1,Me-1);
pdlog=zeros(1,Me-1);
if length(Rmin)==1
    RminM=repmat(Rmin,Ps,D); 
    RmaxM=repmat(Rmax,Ps,D);
end
pos=RminM+(RmaxM-RminM).*rand(Ps,D);
fit=feval(fhd,pos',varargin{:})-foptimal; 
fitcount=Ps; 

pbest=pos; 
pbestval=fit; 
[gbestval,gbestid]=min(pbestval); 
gbestlog(1)=gbestval; 
gbest=pbest(gbestid,:); 
gbestrep=repmat(gbest,Ps,1); 

%% parameter settings  
alpha_start=1.0;
alpha_end=0.5;

%% iterative process
for t=2:Me

    %% calculate mbest
    pavg=mean(pbest);
    pavgrep=repmat(pavg,Ps,1);
    
    %% Using the VPG strategy to set virtual positions
    [Vpbest,Vpavg,xdlog(t-1),pdlog(t-1)] = VPG(pos,pbest,gbest,pavgrep,t,Me,RmaxM,RminM);
     
    %% update particles
    alpha=(alpha_start-alpha_end)*(1-t/Me)+alpha_end;
    fi=rand(Ps,D);
    p=fi.*Vpbest+(1-fi).*gbestrep; %Vpbest replace pbest
    b=alpha*abs(Vpavg-pos); %Vpavg replace pavg
    u=rand(Ps,D);
    pos=p+((-1).^ceil(0.5+rand(Ps,D))).*b.*(-log(u)); 
    pos=((pos>=RminM)&(pos<=RmaxM)).*pos...
        +((pos<RminM)|(pos>RmaxM)).*(RminM+(RmaxM-RminM).*rand(Ps,D)); 

    %% update pbests and gbest
    fit=feval(fhd,pos',varargin{:})-foptimal; 
    fitcount=fitcount+Ps; 
    tmp=(pbestval<fit); 
    temp=repmat(tmp',1,D); 
    pbest=temp.*pbest+(1-temp).*pos; 
    pbestval=tmp.*pbestval+(1-tmp).*fit; 
    [gbestval,gid]=min(pbestval); 
    gbestlog(t)=gbestval;
    gbest=pbest(gid,:);
    gbestrep=repmat(gbest,Ps,1);
end