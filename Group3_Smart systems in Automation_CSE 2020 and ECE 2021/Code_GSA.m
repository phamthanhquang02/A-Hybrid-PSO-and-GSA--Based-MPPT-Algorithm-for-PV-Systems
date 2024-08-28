function D = GSA(Vpv,Ipv)
%%CHECKED%%
persistent u;%u-th particle
persistent dcurrent;%store current duty cycle
persistent pbest;%store local best dc for power
persistent force; %store force
persistent acceleration; %store acceleration
persistent mass; % mass
persistent q; %  strength of mass
persistent p; %  power for each particle
persistent p_current; %  power current for each particle
persistent p_min; %  power min for each particle
persistent worse;   %store best worse of each particle 
persistent dc; %store duty cycle ~ position 
persistent v;  %velocity
persistent counter; %delay iteration
persistent iteration;
persistent gbest;%store global best dc for power
%initialization
num_agency = 3;
max_iter = 500;
if(isempty(counter))
    counter = 0;
    dcurrent = 0.5; 
    gbest = 0.5;
    pbest = zeros(num_agency,1);
    worse = zeros(num_agency,1);
    v = zeros (num_agency,1);
    force = zeros(num_agency,1);
    mass = zeros(num_agency,1);
    q = zeros(num_agency,1);
    p = zeros(num_agency,1);
    p_current = zeros(num_agency,1);
    p_min=zeros(num_agency,1);
    p_min(1)= Vpv*Ipv;
    p_min(2)= Vpv*Ipv;
    p_min(3)= Vpv*Ipv;
    acceleration=zeros(num_agency,1);
    u = 0;
    dc = zeros (num_agency,1); 
    iteration = 1;
    
    %initialize position for each particle
    dc(1) = 0.2;
    dc(2) = 0.5;
    dc(3) = 0.8;
end

if(counter >=1 && counter < 100)
    D=dcurrent;
    counter= counter+1;
    return;
end

if(u>=1 && u<=num_agency)
    p_current(u) = Vpv*Ipv;
    if((Vpv*Ipv)>=p(u))
        p(u) = Vpv*Ipv;
        pbest(u)=dcurrent;
    end
    if(Vpv*Ipv < p_min(u))
        p_min(u) = Vpv*Ipv;
        worse(u) = dcurrent;
    end
end
u=u+1;
if(u== num_agency + 2)
    u=1;
end
if(u >= 1 && u <= num_agency)
    %Avoid over shooting
    if(iteration < max_iter)
        D=dc(u);
        dcurrent=D;
        counter=1;
    return;
    else
        D = dcurrent;
        return
    end
elseif(u==num_agency+1)
    iteration = iteration +1;
    [~,i]=max(p);
    gbest=pbest(i);
    D=gbest;
    dcurrent=D;
    counter=1;
     %Calculate strength of mass
    for i = 1:num_agency
     q(i) = (p_current(i) - worse(i))/(pbest(i)-worse(i));
    end
     %Calculate sum of strength of mass

     sum_strength_of_mass = q(1) + q(2) + q(3);
    
    %Calculate mass   
    for i = 1:num_agency
     mass(i) = q(i)/sum_strength_of_mass;
    
    end
    %Calculate force
    alpha = 20;
    G0 = 1;
    G = G0 * exp(-alpha*iteration/max_iter);
    %G = 6.67430 * 10^-13; %gravitational constant
    e = 2.2204*10^-16;
    force(1) = rand()*G*(mass(3)*mass(1)*(dc(3)-dc(1))/(Euclidian_distance(dc(3),dc(1))+e) + mass(2)*mass(1)*(dc(3)-dc(1))/(Euclidian_distance(dc(3),dc(1))+e));
    force(2) = rand()*G*(mass(3)*mass(2)*(dc(3)-dc(2))/(Euclidian_distance(dc(3),dc(2))+e) + mass(1)*mass(2)*(dc(1)-dc(2))/(Euclidian_distance(dc(1),dc(2))+e));
    force(3) = rand()*G*(mass(2)*mass(3)*(dc(2)-dc(3))/(Euclidian_distance(dc(2),dc(3))+e) + mass(1)*mass(3)*(dc(1)-dc(3))/(Euclidian_distance(dc(1),dc(3))+e));
    %Avoid over shooting
    if(iteration > max_iter)
        disp("should done")
        D=dcurrent;
        return;
    end
    %Calculate acceleration 
    for i = 1:num_agency
        acceleration(i) = force(i)/mass(i);
         disp(force(i))
    end
    
    for i=1:num_agency  
     v(i)=updatevelocity(v(i),acceleration(i));
     dc(i)=updateduty(dc(i),v(i));
     
    end  
    return;
    
else
    D=0.5
end
end
function d = Euclidian_distance(d1,d2)
    d = sqrt(d1^2+d2^2);
end
function vfinal=updatevelocity(velocity,acceleration)
    vfinal = rand()*velocity + acceleration;
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