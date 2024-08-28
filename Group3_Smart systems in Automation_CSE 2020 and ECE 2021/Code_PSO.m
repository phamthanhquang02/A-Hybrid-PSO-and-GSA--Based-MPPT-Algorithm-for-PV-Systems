function D = PSO(Vpv,Ipv)
%%%CHECKEDDDDDD%%%%%%
persistent u;%u-th particle
persistent dcurrent;%store current duty cycle
persistent pbest;%store local best dc for power
persistent p;%  power for each particle
persistent dc;%store duty cycle ~ position 
persistent v;%velocity
persistent counter;%delay iteration
persistent gbest;%store global best dc for power
persistent convergence;%check dc convergence status
persistent P_prev;%store previous power
%initialization
num_agency = 3;
if(isempty(counter))
    counter = 0;
    dcurrent=0.5;
    gbest=0.9;
    p = zeros(num_agency,1);
    v=zeros(num_agency,1);
    pbest = zeros(num_agency,1);
    u=0;
    dc=zeros(num_agency,1);
    dc(1)=0.05;
    dc(2)=0.4;
    dc(3)=0.8;
    P_prev = Ipv*Vpv;
    convergence = 0;
end
if(convergence == 1)
        if(abs(Vpv*Ipv-P_prev)/(P_prev + 0.0001) >= 0.3)
            counter = 0;
            dcurrent=0.5;
            gbest=0.9;
            p = zeros(num_agency,1);
            v=zeros(num_agency,1);
            pbest = zeros(num_agency,1);
            u=0;
            dc=zeros(num_agency,1);
            convergence = 0;
            dc(1)=0.05;
            dc(2)=0.65;
            dc(3)=0.8;
        end
end
%%3800
if(counter >=1 && counter < 4000)
    D=dcurrent;
    counter= counter+1;
    return;
end

if(u>=1 && u<= num_agency)
    if((Vpv*Ipv)>p(u))
        p(u) = Vpv*Ipv;
        pbest(u)=dcurrent;
    end
end
u=u+1;
if(u==num_agency+2)
    u=1;
end
if(u>=1 && u<= num_agency)
    D=dc(u);    
    dcurrent=D;
    counter=1;
    return;
elseif(u==num_agency+1)
    [~,i]=max(p);
    gbest=pbest(i);
    D=gbest;
    dcurrent=D;
    counter=1;
    for i=1:num_agency
     v(i)=updatevelocity(v(i),pbest(i),dc(i),gbest);
      dc(i)=updateduty(dc(i),v(i));
    end 
     P_prev = Ipv*Vpv;
    convergence = checkconvergence(dc(1),dc(2),dc(3));
    return;
    
else
    D=0.5
end
end

function vfinal=updatevelocity(velocity,pobest,d,gwbest)

w=0.1;
c1=1.2;
c2=1.2;

vfinal = (w*velocity)+(c1*rand(1)*(pobest-d))+(c2*rand(1)*(gwbest-d));
end
function status = checkconvergence(d1,d2,d3)
    status = 0;
    rate = 0.05;
    if(abs(d1-d2) < rate)
        if(abs(d1-d3) < rate)
            if(abs(d2-d3) < rate)
                status = 1;
            end
        end
    end
end
function dfinal=updateduty(d,velocity)
dup=d+velocity;
if(dup>1)
    dfinal=1;
elseif(dup<0.1)
    dfinal=0,1;
else
    dfinal=dup;
end
end