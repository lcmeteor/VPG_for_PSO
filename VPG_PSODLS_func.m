function [gbest,gbestval,gbestlog,xdlog,pdlog]= VPG_PSODLS_func(fhd,Dimension,Particle_Number,Max_Gen,Rmin,Rmax,foptimal,varargin)
% PSODLS enabled by the VPG strategy
% Reference: W. Ye, W. Feng, S. Fan, A novel multi-swarm particle swarm optimization with dynamic learning strategy, Applied Soft Computing 61 (2017) 832-843.

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
w_start=0.9;
w_end=0.4;
c1=1.499445;
c2=1.499445;
M=10;
subpopsize=Ps/M;

%% initilization of subswarms
subxdlog=zeros(M,Me-1);
subpdlog=zeros(M,Me-1);
for m=1:M
    subswarm(m,:,:)=pos(((m-1)*subpopsize+1):(m*subpopsize),:);
    subswarmval(m,:)=fit(((m-1)*subpopsize+1):(m*subpopsize));
    [~,index]=min(subswarmval(m,:));
    lbest(m,:)= subswarm(m,index,:);
    lbestval(m)= subswarmval(m,index);
end

%% iterative process
for t=2:Me
    %% Using the VPG strategy to set virtual positions
    [Vpbest,Vvel,xdlog(t-1),pdlog(t-1)] = VPG(pos,pbest,gbest,vel,t,Me,RmaxM,RminM);
    for m=1:M
        [subVpbest(m,:,:),subVvel(m,:,:),subxdlog(m,t-1),subpdlog(m,t-1)] =...
            VPG(pos(((m-1)*subpopsize+1):(m*subpopsize),:),...
            pbest(((m-1)*subpopsize+1):(m*subpopsize),:),...
            lbest(m,:),vel(((m-1)*subpopsize+1):(m*subpopsize),:),t,Me,...
            RmaxM(((m-1)*subpopsize+1):(m*subpopsize),:),...
            RminM(((m-1)*subpopsize+1):(m*subpopsize),:));
    end

    %% update positions and pbests
    w=(w_start-w_end)*(1-t/Me)+w_end;
    mlbest=mean(lbest);   
    for m=1:M
        for n=1:subpopsize
            p=m*(M-1)+n;
            if rand()<(t/Me)
                subvirpbest=squeeze(subVpbest(m,n,:))';
                subvirvel=squeeze(subVvel(m,n,:))';
                vel(p,:)=w*subvirvel+c1.*rand(1,D).*(subvirpbest-pos(p,:))+c2.*rand(1,D).*(lbest(m,:)-pos(p,:)); %using virtual positions of subswarm
            else
                vel(p,:)=w*Vvel(p,:)+c1.*rand(1,D).*(Vpbest(p,:)-pos(p,:))+c2.*rand(1,D).*(mlbest-pos(p,:)); %using virtual positions of whole swarm
            end
            pos(p,:)=pos(p,:)+vel(p,:);
            pos(p,:)=((pos(p,:)>=RminM(p,:))&(pos(p,:)<=RmaxM(p,:))).*pos(p,:)...
                +((pos(p,:)<RminM(p,:))|(pos(p,:)>RmaxM(p,:))).*(RminM(p,:)+(RmaxM(p,:)-RminM(p,:)).*rand(1,D));
            fit(p)=feval(fhd,pos(p,:)',varargin{:})-foptimal;
            fitcount=fitcount+1;
            if fit(p)<pbestval(p)
                pbestval(p)=fit(p);
                pbest(p,:)=pos(p,:);
            end
            if fit(p)<lbestval(m)
                lbestval(m)=fit(p);
                lbest(m,:)=pos(p,:);
            end
        end
    end

    %% update gbest
    [gbestval,gbestid]=min(pbestval);
    gbestlog(t)=gbestval;
    gbest=pbest(gbestid,:);
end