function D = PSOGSA2(Vpv,Ipv)
    %%CHECKED TESTING%%
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
    persistent gbest;%store global best dc for power
    persistent convergence;%check dc convergence status
    persistent P_prev;%store previous power
    %initialization
    num_agency = 3;
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
        P_prev = Ipv*Vpv;
        convergence = 0;

        %initialize position for each particle
        dc(1) = 0.2;
        dc(2) = 0.6;
        dc(3) = 0.9;
    end

    if(convergence == 1)
            if(abs(Vpv*Ipv-P_prev) >= 0.3*P_prev)
                disp("help")
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
                P_prev = Ipv*Vpv;
                convergence = 0;

                %initialize position for each particle
                dc(1) = 0.2;
                dc(2) = 0.6;
                dc(3) = 0.9;
            end
    end
    %Delay for more visualization
    if(counter >=1 && counter < 500)
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
        %Calculate strength of mass
        disp('chia 1');
        for i = 1:num_agency
            q(i) = (p_current(i) - 0.99*worse(i))/(pbest(i)-worse(i)+ 0.0001);
        end
         %Calculate sum of strength of mass

        sum_strength_of_mass = q(1) + q(2) + q(3);

        %Calculate mass   
        disp('chia 2');
        for i = 1:num_agency
            mass(i) = q(i)*5/sum_strength_of_mass;
        end
        %Calculate force
        G = 6.67430 * 10^-13; %gravitational constant
        e = 2.2204*10^-16;
        disp('chia 3');
        force(1) = rand()*G*(mass(3)*mass(1)*(dc(3)-dc(1))/(Euclidian_distance(dc(3),dc(1))+e) + mass(2)*mass(1)*(dc(3)-dc(1))/(Euclidian_distance(dc(3),dc(1))+e));
        force(2) = rand()*G*(mass(3)*mass(2)*(dc(3)-dc(2))/(Euclidian_distance(dc(3),dc(2))+e) + mass(1)*mass(2)*(dc(1)-dc(2))/(Euclidian_distance(dc(1),dc(2))+e));
        force(3) = rand()*G*(mass(2)*mass(3)*(dc(2)-dc(3))/(Euclidian_distance(dc(2),dc(3))+e) + mass(1)*mass(3)*(dc(1)-dc(3))/(Euclidian_distance(dc(1),dc(3))+e));
        %Calculate acceleration 
        disp('chia 4');
        for i = 1:num_agency
            acceleration(i) = force(i)/mass(i);
        end

        disp('chia 5');
        for i=1:num_agency
            v(i)=updatevelocity(v(i),acceleration(i),dc(i),gbest);
            dc(i)=updateduty(dc(i),v(i));
        end  
        %disp(P_prev)
        %disp(convergence)

        P_prev = Ipv*Vpv;
        convergence = checkconvergence(dc(1),dc(2),dc(3));
        return;

    else
        D=0.5;
    end
end

function status = checkconvergence(d1,d2,d3)
    status = 0;
    rate = 0.002;
    if(abs(d1-d2) < rate)
        if(abs(d1-d3) < rate)
            if(abs(d2-d3) < rate)
                status = 1;
            end
        end
    end
end
function vfinal=updatevelocity(velocity,acceleration,d,gwbest)
    w=0.1;
    c1=1.2;
    c2=1.5;

    vfinal = w*velocity + (rand(1)*c1*(acceleration))+(c2*rand(1)*(gwbest-d));
end
function d = Euclidian_distance(d1,d2)
    d = sqrt(d1^2+d2^2);
end
function dfinal=updateduty(d,velocity)
    dup=d+velocity;
    if(dup>1)
        dfinal=1;
    elseif(dup<0.1)
        dfinal=0.1;
    else
        dfinal=dup;
    end
end