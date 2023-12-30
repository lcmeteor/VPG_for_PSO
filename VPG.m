function [Vpbest,Vvel,xdiv,pdiv] = VPG(pos,pbest,gbest,vel,t,Me,RmaxM,RminM)

[Ps,D]=size(pos);

%% parameter setting of VPG
pdiv_start=0.5; %initial value of the baseline 
pdiv_end=0; %end value of the baseline 
cube=10; %decline rate of the baseline 

%% calculate normalized swarm diversity
npos=(pos-RminM)./(RmaxM-RminM);
npbest=(pbest-RminM)./(RmaxM-RminM);
nxavg=mean(npos);
xdis=zeros(1,Ps);
for n=1:Ps
    xdis(n)=sqrt(sum((npos(n,:)-nxavg).^2));
end
xdiv=mean(xdis)/sqrt(D); %D(X)
npavg=mean(npbest);
pdis=zeros(1,Ps);
for n=1:Ps
    pdis(n)=sqrt(sum((npbest(n,:)-npavg).^2));
end
pdiv=mean(pdis)/sqrt(D); %D(P)

%% switch phase in each iteration and set virtual positions
Vpbest=pbest;
Vvel=vel;
pdiv_base=(pdiv_start-pdiv_end)*(1-t/Me)^cube+pdiv_end;
if pdiv<pdiv_base
    if xdiv<pdiv_base
        %% divergence phase
        for n=1:Ps
            for d=1:D
                if pbest(n,d)>=gbest(d)
                    edgeR=RmaxM(n,d);
                else
                    edgeR=RminM(n,d);
                end
                ind=(1-((pbest(n,d)-gbest(d))/(edgeR-gbest(d)+eps))^(1/cube))*Me;
                if ind>t
                    pbasevp=(edgeR-gbest(d))*(1-t/Me)^cube+gbest(d);
                    Vpbest(n,d)=rand()^((Me-t)/(ind-t))*(edgeR-pbasevp)+pbasevp;
                else
                end
            end
        end
        Vvel=vel;
    else 
        %% normal phase
        Vpbest=pbest;
        Vvel=vel;
    end
else 
    %% acceleration phase
    for n=1:Ps
        for d=1:D
            if pbest(n,d)>=gbest(d)
                edgeR=RmaxM(n,d);
            else
                edgeR=RminM(n,d);
            end
            ind=(1-((pbest(n,d)-gbest(d))/(edgeR-gbest(d)+eps))^(1/cube))*Me;
            if ind<t
                pbasevp=(edgeR-gbest(d))*(1-t/Me)^cube+gbest(d);
                Vpbest(n,d)=pbasevp-rand()^(t/(t-ind))*(pbasevp-gbest(d));
            else
            end
            if (pbest(n,d)-gbest(d))==0
                Vvel(n,d)=(Vpbest(n,d)-gbest(d))/(pbest(n,d)-gbest(d)+eps)*vel(n,d);
            else
                Vvel(n,d)=(Vpbest(n,d)-gbest(d))/(pbest(n,d)-gbest(d))*vel(n,d);
            end
        end
    end
end
