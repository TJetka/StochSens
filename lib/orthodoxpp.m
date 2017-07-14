function orthodoxpp(iteration)
figure;
hold off;
count=0;
for i=1:iteration,
    params=[rand(1,4)*4,rand(1,4)*10,rand(1,4)*4,rand(1,4)*10];
    [t,x]=ode45(@(t,x) fpp(t,x,params), [0,100],[0;0;0]);
    V1=params(1);
    V2=params(2);
    if V1<V2;
        plot(t,x,'-')
        hold on;
        t1=1;
        t2=2;
        t3=3;
        D1=feval(@fpp,t1,x(t1,:),params);
        D2=feval(@fpp,t1,x(t2,:),params);
        D3=feval(@fpp,t1,x(t3,:),params);
        D0=[D1(3,1),D2(3,1),D3(3,1)];
        deri=fpp(D0);
        for j=1:2;
            if deri(j)-deri(j+1)>0; 
                count=count+1;
            end
        end
    end
end
b=count;
title('HKP, AP and RRP')
xlabel ('t'), ylabel ('HKP, AP and RRP')
legend('HKP','AP','RRP'), grid