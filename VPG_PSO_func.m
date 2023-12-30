function [gbest,gbestval,gbestlog,xdlog,pdlog] = VPG_PSO_func(fhd,Dimension,Particle_Number,Max_Gen,Rmin,Rmax,foptimal,varargin)
% Canonical PSO enabled by the VPG strategy

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

Vmax=0.5*(RmaxM-RminM); 
Vmin=-Vmax;
vel=Vmin+2.*Vmax.*rand(Ps,D); 

pbest=pos; 
pbestval=fit; 
[gbestval,gbestid]=min(pbestval); 
gbestlog(1)=gbestval; 
gbest=pbest(gbestid,:); 
gbestrep=repmat(gbest,Ps,1); 

%% parameter setting
w_start=0.729;
w_end=0.729;
cc1=2.0412;
cc2=0.9477;

%% iterative process
for t=2:Me
    %% using the VPG strategy to set virtual positions
    [Vpbest,Vvel,xdlog(t-1),pdlog(t-1)] = VPG(pos,pbest,gbest,vel,t,Me,RmaxM,RminM);
    
    %% update particles
    w=(w_start-w_end)*(1-t/Me)+w_end;
    aa=cc1.*rand(Ps,D).*(Vpbest-pos)+cc2.*rand(Ps,D).*(gbestrep-pos); %Vpbest replace pbest
    vel=w.*Vvel+aa;  %Vvel replace vel
    pos=pos+vel;                 
    pos=((pos>=RminM)&(pos<=RmaxM)).*pos...
        +((pos<RminM)|(pos>RmaxM)).*(RminM+(RmaxM-RminM).*rand(Ps,D));
    fit=feval(fhd,pos',varargin{:})-foptimal; 
    fitcount=fitcount+Ps;
    
    %% update pbests and gbest
    tmp=(pbestval<fit); 
    temp=repmat(tmp',1,D); 
    pbest=temp.*pbest+(1-temp).*pos; 
    pbestval=tmp.*pbestval+(1-tmp).*fit; 
    [gbestval,gid]=min(pbestval); 
    gbestlog(t)=gbestval; 
    gbest=pbest(gid,:); 
    gbestrep=repmat(gbest,Ps,1);
end