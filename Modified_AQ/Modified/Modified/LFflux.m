% input: iel = element number
%        pdeg = polynomial order
%        Asol
%        Qsol = vector of coefficients for the basis function
%        Bext = specified boundary data
% 
% output: flux is a (pdeg + 1, 2) matrix
%         flux(j, 1) = flux evaluated at xe times values of jth basis
%         flux(j, 2) = flux evaluated at xe1 times values of jth basis 
%         
%         fluxA corresponds the area equation
%         fluxQ corresponds the mass flux equation
% 
%         xe---(element iel)---xe1
% 
% description: returns the Lax-Friedrichs flux values for element iel

function [fluxA, fluxQ] = LFflux(iel, pdeg, Asol, Qsol, Bext, variables)
    
    N = variables.N;
    cardiogrid = variables.cardiogrid;
    type = variables.type;
    G0_ = variables.G0_;
    rho = variables.rho;
    alpha = variables.alpha;
    locdim = pdeg + 1;
    fluxA = zeros(locdim, 2);
    fluxQ = zeros(locdim, 2);

    %%%%%% INLET %%%%%%
    if iel == 1 
        %%%%%%%% FOR RIGHT NODE xe1 %%%%%%%%
        xlocL = 1; % limit from the left --> xe1
        [basephiL, basephiLD, ~, xe] = DGelemcalc(pdeg, xlocL, cardiogrid(iel), cardiogrid(iel+1));

        xlocR = -1;  % limit from the right --> xe1
        [basephiR, basephiRD, ~, ~] = DGelemcalc(pdeg, xlocR, cardiogrid(iel+1), cardiogrid(iel+2));

        % compute A0 at xe
        [~, partialR0, ~, ~] = Stenosis(xe, type, variables);

        % compute AL, QL, AR, QR at xe
        [AL, ~] = basiseval(iel, locdim, basephiL, basephiLD, Asol);
        [QL, ~] = basiseval(iel, locdim, basephiL, basephiLD, Qsol);

        [AR, ~] = basiseval(iel+1, locdim, basephiR, basephiRD, Asol);
        [QR, ~] = basiseval(iel+1, locdim, basephiR, basephiRD, Qsol);

        % compute eigenvalues at xe
        alpha_c = -2/35*(partialR0^2);
        lamL1 = (alpha+alpha_c)*QL/AL + sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamL2 = (alpha+alpha_c)*QL/AL - sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamR1 = (alpha+alpha_c)*QL/AR + sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));
        lamR2 = (alpha+alpha_c)*QL/AR - sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));

        % compute F(U) at xe
        FAL = QL;
        FAR = QR;
        FQL = (alpha+alpha_c)*(QL^2)/AL + G0_/3/rho*(AL^1.5)/(0.18^2);
        FQR = (alpha+alpha_c)*(QR^2)/AR + G0_/3/rho*(AR^1.5)/(0.18^2);

        % obtain Lax-Friedrich flux
        fluxA(:, 2) = basephiL .* (0.5*(FAL + FAR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(AR - AL));
        fluxQ(:, 2) = basephiL .* (0.5*(FQL + FQR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(QR - QL));
        
        %%%%%%%% FOR LEFT NODE xe %%%%%%%%
        %specified boundary data
        Ab = Bext(1); 
        Qb = Bext(2);

        xlocR = -1; % limit from the right --> xe
        [basephibdry, ~, ~, xxB] = DGelemcalc(pdeg, xlocR, cardiogrid(iel), cardiogrid(iel+1));

        % compute A0 at xxB
        [~, partialR0, ~, ~] = Stenosis(xxB, type, variables);
        alpha_c = -2/35*(partialR0^2);

        fluxA(:, 1) = Qb * basephibdry;
        fluxQ(:, 1) = ((alpha+alpha_c)*(Qb^2)/Ab + G0_/3/rho*(Ab^1.5)/(0.18^2)) * basephibdry;

    %%%%%% OUTLET %%%%%%
    elseif iel == N
        %%%%%%%% FOR LEFT NODE xe %%%%%%%%
        xlocL = 1; % limit from the left --> xe
        [basephiL, basephiLD, ~, xe] = DGelemcalc(pdeg, xlocL, cardiogrid(iel-1), cardiogrid(iel));

        xlocR = -1;  % limit from the right --> xe
        [basephiR, basephiRD, ~, ~] = DGelemcalc(pdeg, xlocR, cardiogrid(iel), cardiogrid(iel+1));

        % compute A0 at xe
        [~, partialR0, ~, ~] = Stenosis(xe, type, variables);

        % compute AL, QL, AR, QR at xe
        [AL, ~] = basiseval(iel-1, locdim, basephiL, basephiLD, Asol);
        [QL, ~] = basiseval(iel-1, locdim, basephiL, basephiLD, Qsol);

        [AR, ~] = basiseval(iel, locdim, basephiR, basephiRD, Asol);
        [QR, ~] = basiseval(iel, locdim, basephiR, basephiRD, Qsol);

        % compute eigenvalues at xe
        alpha_c = -2/35*(partialR0^2);
        lamL1 = (alpha+alpha_c)*QL/AL + sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamL2 = (alpha+alpha_c)*QL/AL - sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamR1 = (alpha+alpha_c)*QL/AR + sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));
        lamR2 = (alpha+alpha_c)*QL/AR - sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));
        
        % compute F(U) at xe
        FAL = QL;
        FAR = QR;
        FQL = (alpha+alpha_c)*(QL^2)/AL + G0_/3/rho*(AL^1.5)/(0.18^2);
        FQR = (alpha+alpha_c)*(QR^2)/AR + G0_/3/rho*(AR^1.5)/(0.18^2);

        % obtain Lax-Friedrich flux
        fluxA(:, 1) = basephiR .* (0.5*(FAL + FAR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(AR - AL));
        fluxQ(:, 1) = basephiR .* (0.5*(FQL + FQR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(QR - QL));
        
        %%%%%%%% FOR RIGHT NODE xe1 %%%%%%%%
        %specified boundary data
        Ab = Bext(3); 
        Qb = Bext(4);

        xlocR = 1; % limit from the left --> xe1
        [basephibdry, ~, ~, xxB] = DGelemcalc(pdeg, xlocR, cardiogrid(iel), cardiogrid(iel+1));

        % compute A0 at xxB
        [~, partialR0, ~, ~] = Stenosis(xxB, type, variables);
        alpha_c = -2/35*(partialR0^2);

        fluxA(:, 2) = Qb * basephibdry;
        fluxQ(:, 2) = ((alpha+alpha_c)*(Qb^2)/Ab + G0_/3/rho*(Ab^1.5)/(0.18^2)) * basephibdry;
    
    %%%%% INTERIOR POINTS %%%%%% 
    else
        %%%%%%% FOR LEFT NODE xe %%%%%%%%
        xlocL = 1; % limit from the left --> xe
        [basephiL, basephiLD, ~, xe] = DGelemcalc(pdeg, xlocL, cardiogrid(iel-1), cardiogrid(iel));

        xlocR = -1;  % limit from the right --> xe
        [basephiR, basephiRD, ~, ~] = DGelemcalc(pdeg, xlocR, cardiogrid(iel), cardiogrid(iel+1));

        % compute A0 at xe
        [~, partialR0, ~, ~] = Stenosis(xe, type, variables);

        % compute AL, QL, AR, QR at xe
        [AL, ~] = basiseval(iel-1, locdim, basephiL, basephiLD, Asol);
        [QL, ~] = basiseval(iel-1, locdim, basephiL, basephiLD, Qsol);

        [AR, ~] = basiseval(iel, locdim, basephiR, basephiRD, Asol);
        [QR, ~] = basiseval(iel, locdim, basephiR, basephiRD, Qsol);

        % compute eigenvalues at xe
        alpha_c = -2/35*(partialR0^2);
        lamL1 = (alpha+alpha_c)*QL/AL + sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamL2 = (alpha+alpha_c)*QL/AL - sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamR1 = (alpha+alpha_c)*QL/AR + sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));
        lamR2 = (alpha+alpha_c)*QL/AR - sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));

        % compute F(U) at xe
        FAL = QL;
        FAR = QR;
        FQL = (alpha+alpha_c)*(QL^2)/AL + G0_/3/rho*(AL^1.5)/(0.18^2);
        FQR = (alpha+alpha_c)*(QR^2)/AR + G0_/3/rho*(AR^1.5)/(0.18^2);

        % obtain Lax-Friedrich flux
        fluxA(:, 1) = basephiR .* (0.5*(FAL + FAR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(AR - AL));
        fluxQ(:, 1) = basephiR .* (0.5*(FQL + FQR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(QR - QL));

        %%%%%%% FOR RIGHT NODE xe1 %%%%%%%%
        xlocL = 1; % limit from the left --> xe1
        [basephiL, basephiLD, ~, xe] = DGelemcalc(pdeg, xlocL, cardiogrid(iel), cardiogrid(iel+1));

        xlocR = -1;  % limit from the right --> xe1
        [basephiR, basephiRD, ~, ~] = DGelemcalc(pdeg, xlocR, cardiogrid(iel+1), cardiogrid(iel+2));

        % compute A0 at xe
        [~, partialR0, ~, ~] = Stenosis(xe, type, variables);

        % compute AL, QL, AR, QR at xe
        [AL, ~] = basiseval(iel, locdim, basephiL, basephiLD, Asol);
        [QL, ~] = basiseval(iel, locdim, basephiL, basephiLD, Qsol);

        [AR, ~] = basiseval(iel+1, locdim, basephiR, basephiRD, Asol);
        [QR, ~] = basiseval(iel+1, locdim, basephiR, basephiRD, Qsol);

        % compute eigenvalues at xe
        alpha_c = -2/35*(partialR0^2);
        lamL1 = (alpha+alpha_c)*QL/AL + sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamL2 = (alpha+alpha_c)*QL/AL - sqrt(((alpha+alpha_c)*QL/AL)^2 - (alpha+alpha_c)*(QL^2)/(AL^2) + G0_/rho/2*sqrt(AL)/(0.18^2));
        lamR1 = (alpha+alpha_c)*QL/AR + sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));
        lamR2 = (alpha+alpha_c)*QL/AR - sqrt(((alpha+alpha_c)*QR/AR)^2 - (alpha+alpha_c)*(QR^2)/(AR^2) + G0_/rho/2*sqrt(AR)/(0.18^2));

        % compute F(U) at xe
        FAL = QL;
        FAR = QR;
        FQL = (alpha+alpha_c)*(QL^2)/AL + G0_/3/rho*(AL^1.5)/(0.18^2);
        FQR = (alpha+alpha_c)*(QR^2)/AR + G0_/3/rho*(AR^1.5)/(0.18^2);

        % obtain Lax-Friedrich flux
        fluxA(:, 2) = basephiL .* (0.5*(FAL + FAR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(AR - AL));
        fluxQ(:, 2) = basephiL .* (0.5*(FQL + FQR) - 0.5*max(abs([lamL1, lamL2, lamR1, lamR2]))*(QR - QL));

    end    

    return