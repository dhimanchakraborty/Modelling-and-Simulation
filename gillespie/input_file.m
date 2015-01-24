function input_file()
import methods.*
clear all;
p.k1 = 2;
p.k2 = 3;
p.k3 = 1;
p.k4 = 2;
p.k5 = 3;
p.k6 = 3;
tspan = [0, 10000];
x0 = [1000, 0, 0, 0, 0, 0];
%% Specify reaction network
prop_func = @propensity;
react_matrix = [-1  1  0  0  0  0
		 	    -1  0  1  0  0  0
			     0 -1  0  1  0  0
		 	     0  0 -1  0  1  0
		 	     0  0  0 -1  0  1
		 	     0  0  0  0 -1  1];
[t,x] = m_GillespieAlgo(react_matrix, prop_func, tspan, x0, p);
%[t,x] = m_FirstReaction(react_matrix, prop_func, tspan, x0, p);

%% Plot time course
figure();
plot(t,x);
%plot(t,x(:,6));
xlabel('time (s)');
ylabel('Copies');
legend({'A','B1','B2','C1','C2','D'});

end

function a = propensity(x, p)
A = x(1);
B1 = x(2);
B2 = x(3);
C1 = x(4);
C2 = x(5);
D = x(6);

a = [p.k1*A;
     p.k2*A;
     p.k3*B1;
     p.k4*B2;
     p.k5*C1;
     p.k6*C2];
end