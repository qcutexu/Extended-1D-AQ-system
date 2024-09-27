% This subroutine aims for plotting

Aval = zeros(locdim, variables.N);
Qval = zeros(locdim, variables.N);
R0val = zeros(locdim, variables.N);
partialR0val = zeros(locdim, variables.N);
X = zeros(locdim, variables.N);

for iel = 1 : variables.N
    [basephi, basephiD, ~, xx] = DGelemcalc(variables.pdeg, [-1, 1], variables.cardiogrid(iel), variables.cardiogrid(iel+1));
    [A, ~] = basiseval(iel, locdim, basephi, basephiD, Asol);
    [Q, ~] = basiseval(iel, locdim, basephi, basephiD, Qsol);
    [R0, partialR0, ~, ~] = Stenosis(xx, variables.type, variables);
    Aval(:, iel) = A';
    Qval(:, iel) = Q';
    R0val(:, iel) = R0';
    partialR0val(:, iel) = partialR0';
    X(:, iel) = xx';
end

subplot(2, 2, 1);
plot(X, Aval, 'b', LineWidth=1.5);
xlabel('x','FontSize',16);ylabel('A','FontSize',16);

subplot(2, 2, 2);
plot(X, Qval, 'r', LineWidth=1.5);
xlabel('x','FontSize',16);
ylabel('Q','FontSize',16);
drawnow;

p = variables.G0_./(0.18^2).*(sqrt(Aval)-R0val)...
    + variables.nu*variables.rho*(variables.alpha/(variables.alpha-1))*Qval./Aval./R0val.*partialR0val;

subplot(2,2,3);
plot(X, p, 'r-', LineWidth=1.5);   
xlabel('x','FontSize',16);
ylabel('pressure','FontSize',16);
drawnow;

subplot(2,2,4);
plot(X, Qval./Aval, 'b', LineWidth=1.5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%