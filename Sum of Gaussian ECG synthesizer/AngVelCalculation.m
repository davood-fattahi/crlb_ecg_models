function omega=AngVelCalculation(phase,Ts)


y3=phase(2:end)-phase(1:end-1);
b=find(y3<0);
y3(b)=y3(b+1);
omega=[y3(1)  y3]./Ts; % angular velocity