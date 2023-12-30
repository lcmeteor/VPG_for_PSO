function [gbest,gbestval,gbestlog,xdlog,pdlog]= VPG_GLPSO_func(fhd,Dimension,Particle_Number,Max_Gen,Rmin,Rmax,foptimal,varargin)
% GLPSO enabled by the VPG strategy
% Reference: Y.J. Gong, J.J. Li, Y. Zhou, Y. Li, H.S.H. Chung, Y.H. Shi, J. Zhang, Genetic learning particle swarm optimization, IEEE Transactions on Cybernetics 46.10 (2015) 2277-2290.

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

%% parameter settings 
offp=pos;
exep=pos;
fe=fit;
w=0.7298;
cc=1.49618;
pm=0.01;
sg=7;
ci=0;

%% iterative process
for t=2:Me
    %% Using the VPG strategy to set virtual positions
    [Vpbest,Vvel,xdlog(t-1),pdlog(t-1)] = VPG(pos,pbest,gbest,vel,t,Me,RmaxM,RminM);

    %% update particles
    for n=1:Ps 
        for d=1:D
            k=randi(Ps);
            if fit(n)<fit(k)
                fi=rand();
                offp(n,d)=fi*Vpbest(n,d)+(1-fi)*gbest(d); %Vpbest replace pbest
            else
                offp(n,d)=Vpbest(k,d); %Vpbest replace pbest
            end
        end
        for d=1:D
            if rand<pm
                offp(n,d)=Rmin+rand*(Rmax-Rmin);
            end
        end
    end
    fo=feval(fhd,offp',varargin{:})-foptimal;
    for n=1:Ps
        if fo(n)<fe(n)
            exep(n,:)=offp(n,:);
            fe(n)=fo(n);
            ci=0;
        else
            ci=ci+1;
        end
        if ci>sg
            cop=randperm(Ps);
            for m=1:floor(0.2*Ps)
                if exep(cop(m),:)<exep(n,:)
                    exep(n,:)=exep(cop(m),:);
                    fe(n)=fe(cop(m));
                end
            end
        end
    end
    vel=w.*Vvel+cc.*rand(Ps,D).*(exep-pos); %Vvel replace vel
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
    [gbestval,gbestid]=min(pbestval);
    gbestlog(t)=gbestval;
    gbest=pbest(gbestid,:);
end
