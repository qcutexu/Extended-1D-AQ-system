%clc; clear;
%load('../tex/Case123/Workspace_case3_steadystatemembrane.mat')
clear  t x
global  c_R0 c_partialR0 LL c_Re c_uzdata c_Urdata dx N c_dAdzdata dt 

Uz=Q./A;
RR=sqrt(A);

Re=zeros(size(A));
Urdata=zeros(size(A));
uzdata=zeros(size(A));
u_rSol_All=zeros(101,2,21);

Re = rho* Q ./ A / nu .*R0 *2;
Urdata = Q ./ A *3/2 .* R0 / L;
uzdata = Q./A *3/2; %parabolic profile. maximum is 3/2 times.

loop=1;
%case1: 35, 50, 80
%case2:

for iii= [36,47]*1%case3[36,42,48]%case2 [50,63,70]%case1[46,51,54]
    
    LL=RR(1,iii);
    c_R0=R0(1,iii);
    c_partialR0=partialR0(1,iii);
    c_Re=0.1;%Re(1,iii);
    c_Urdata=Urdata(1,iii);
    c_uzdata=uzdata(1,iii);
    c_dAdzdata =0;%urR(1,iii);
    [r_coord,u_rSol]=solve_bvp();
    u_rSol_All(:,1,loop)=r_coord(:);
    u_rSol_All(:,2,loop)=u_rSol(1,:);

loop=loop+1;
end

saved(:,1)=u_rSol_All(:,1,1);
saved(:,2)=u_rSol_All(:,2,1);
saved(:,3)=u_rSol_All(:,1,2);
saved(:,4)=u_rSol_All(:,2,2);
saved(:,5)=u_rSol_All(:,1,3);
saved(:,6)=u_rSol_All(:,2,3);


function [rcoord,Dsol]=solve_bvp()
global c_R0 c_partialR0 LL c_Re c_uzdata c_Urdata
    xmesh=linspace(1e-4,LL,101);
    options=bvpset('stats','off','RelTol',1e-6);
    solinit=bvpinit(xmesh,[0 0]);
    sol = bvp4c(@OdeBVP, @OdeBC, solinit, options);
    Dsol=deval(sol,xmesh).*c_Urdata;
    rcoord=xmesh;
    plot(sol.x,sol.y(1,:).*c_Urdata,'+');
   
    hold on;

end


function f=OdeBVP(x,y)
global c_R0 c_partialR0 LL c_Re c_uzdata c_Urdata c_dAdzdata 

f =[y(2) 
    (y(1)*c_Re/c_Urdata - c_R0/x)*y(2) + (2*y(1)*c_Re*c_uzdata/(c_R0 * c_Urdata) * (1 - x^2 / LL^2)- 4 * c_uzdata / LL) * c_partialR0 - c_R0 * y(1) / x/x];

end

function res=OdeBC(ya,yb)
global A dx N c_dAdzdata

res=[ya(1)-0
     -yb(1)+c_dAdzdata    
];
end
    


