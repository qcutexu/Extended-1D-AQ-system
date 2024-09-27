% Shihan Guo, 03/07/2024
% East China Normal University
% Driver script for solving the 1D AQ equation

clear;
clc;
tic
Variables;
time = 0;

% initialize Asol, Qsol and Atmp, Qtmp
% Atemp and Qtemp are intermediate stage in RK scheme
locdim = variables.pdeg + 1;
Asol = zeros(locdim, variables.N);
Qsol = zeros(locdim, variables.N);
Atmp1 = zeros(locdim, variables.N);
Qtmp1 = zeros(locdim, variables.N);
Atmp2 = zeros(locdim, variables.N);
Qtmp2 = zeros(locdim, variables.N);
Atmp3 = zeros(locdim, variables.N);
Qtmp3 = zeros(locdim, variables.N);


for i = 1 : variables.N
    Asol(:, i) = DGL2proj(variables.pdeg, i, 1, variables);
    Qsol(:, i) = DGL2proj(variables.pdeg, i, 2, variables);
end

% obtian Jacobian
[~, ~, De, ~] = DGelemcalc(variables.pdeg, 1, variables.cardiogrid(1), variables.cardiogrid(2));

% loop over time
while (time < FinalTime)

    if (time + variables.dt > FinalTime)
        variables.dt = FinalTime - time;
    end
    
    % RK3
    % specify boundary data on boundary nodes    
    Bext = ones(1, 4);
    [Bext(1), Bext(2), Bext(3), Bext(4)] = cardiobdryext_V(Asol, Qsol, variables);

    parfor iel = 1 : variables.N
        [Atmp1(:, iel), Qtmp1(:, iel)] = RKrhs(iel, Asol, Qsol, Bext, variables);
        Atmp1(:, iel) = Asol(:, iel) + variables.dt/De*Atmp1(:, iel);
        Qtmp1(:, iel) = Qsol(:, iel) + variables.dt/De*Qtmp1(:, iel);
    end

    [Atmp1, Qtmp1] = limiter(Atmp1, Qtmp1, variables);

    [Bext(1), Bext(2), Bext(3), Bext(4)] = cardiobdryext_V(Atmp1, Qtmp1, variables);
    parfor iel = 1 : variables.N
        [Atmp2(:, iel), Qtmp2(:, iel)] = RKrhs(iel, Atmp1, Qtmp1, Bext, variables);
        Atmp2(:, iel) = (3*Asol(:, iel) + Atmp1(:, iel) + variables.dt/De*Atmp2(:, iel))/4;
        Qtmp2(:, iel) = (3*Qsol(:, iel) + Qtmp1(:, iel) + variables.dt/De*Qtmp2(:, iel))/4;
    end
    
    [Atmp2, Qtmp2] = limiter(Atmp2, Qtmp2, variables);

    [Bext(1), Bext(2), Bext(3), Bext(4)] = cardiobdryext_V(Atmp2, Qtmp2, variables);
    parfor iel = 1 : variables.N
        [Atmp3(:, iel), Qtmp3(:, iel)] = RKrhs(iel, Atmp2, Qtmp2, Bext, variables);
        Asol(:, iel) = (Asol(:, iel) + 2*Atmp2(:, iel) + 2*variables.dt/De*Atmp3(:, iel))/3;
        Qsol(:, iel) = (Qsol(:, iel) + 2*Qtmp2(:, iel) + 2*variables.dt/De*Qtmp3(:, iel))/3;
    end

    [Asol, Qsol] = limiter(Asol, Qsol, variables);

    time = time + variables.dt;
    TIME = sprintf('Time=====%d',time); disp(TIME);

    % Plotting results
    if mod(floor(time/variables.dt),2000) == 1
        plot_subroutine;
    end

end
toc
save("modified.mat")