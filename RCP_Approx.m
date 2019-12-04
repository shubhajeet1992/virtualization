function [RCP] = RCP_Approx(kappa,X,Delta,BSLocationX,BSLocationY)
%compute RCP for given decision variables x and delta for given kappa.

%   parameters
% BSLocationX = [1.6,0.6,1.1,1.4,1.8,2,1.1,0.3,0.3,0.6]; 
% BSLocationY = [1.7,0.6,1.7,0.5,1.9,0.7,0.4,0.6,1.3,1]; 
Vertc1=[0,0,0];
Vertc2=[0,2,0];
Vertc3=[2,0,0];
Vertc4=[2,2,0];
mu=10^(46/100)/ 1000 ;
muNdB=46*ones(1,length(BSLocationX));
muN=10.^(muNdB/10) ./1000;
alpha=4;
W=20*10^6;
sigma=10^(-174/ 10) *W / 1000 ;
lambda=20;
A=4;
%kappa=1*10^6;
%q=0.7;

% compute distance between BSs and radius of BSs
D=[];
for i = 1: length(BSLocationX)
  Dist=[];
  Dist=sqrt((BSLocationX(i)-BSLocationX).^2 + (BSLocationY(i)-BSLocationY).^2); 
  pt=[BSLocationX(i),BSLocationY(i),0];
  D_edge=zeros(1,4);
  D_edge(1)= norm(cross(Vertc1-Vertc2,pt-Vertc2))/norm(Vertc1-Vertc2);
  D_edge(2)= norm(cross(Vertc2-Vertc3,pt-Vertc3))/norm(Vertc3-Vertc3); 
  D_edge(3)= norm(cross(Vertc3-Vertc4,pt-Vertc4))/norm(Vertc3-Vertc4); 
  D_edge(4)= norm(cross(Vertc4-Vertc1,pt-Vertc1))/norm(Vertc4-Vertc1); 
  D(i,:)= nonzeros(Dist);
  Q(i)=min([nonzeros(X.*Dist);nonzeros(D_edge)]);
  AB(i)=pi*Q(i)^2;
end
%display(Q)
%display(D)



% using integral ()
% R1= @(rho,C,U)integral2(@(w1,theta) 2.^(rho./(delta*W)).*(mu.*(U.^alpha).*(sigma+C)) ...
%      .*(exp(-mu.*(U.^alpha).*(2.^(rho./(delta*W)) -1)))...
%      .*exp(-1i.*w1.*C).*(1- x(1)*(1- (mu.*(U.^2 + d.^2 - 2.*U.*d.*cos(theta)).^(alpha/2)) ...
%            ./((2*pi.*((mu.*(U.^2 + d.^2 - 2.*U.*d.*cos(theta)).^(alpha/2))) - (1i.*w1)))))...
%            .*(2.*U./(q^2)),-Inf,Inf,0,2*pi);
% R= integral3(@(rho,C,U)arrayfun(R1,rho,C,U), 0,kappa ,0,Inf,0,q);
% rcp=1-abs(R);

% % using trapz()
Max=10^50;
Min=-Max;
Gap=2*Max;
y=10^20;
theta=0:2*pi/y:2*pi;
w=Min:Gap/y:Max;
C=0:Max/y:Max;
Rho=0:kappa/y:kappa;

for xx=1:length(X)
    if X(xx)==0|| Delta(xx)==0
        R_pdf(xx)=0;
    else
            x=X([1:xx-1 xx+2:end]);
            delta=Delta(xx);
            d=D(xx,:);
            q=Q(xx);
            Ab=AB(xx);
            U=0:q/y:q;
            for m=1:length(Rho)
                rho=Rho(m);
                for mm=1:length(C)
                    c=C(mm);
                    for mmm=1:length(U)
                        u=U(mmm);
                        for n=1:length(w)
                            w1=w(n);
                            for nn=1:length(x) % or, mu
                                dn=d(nn);
                                mun=muN(nn);
                                for nnn=1:length(theta)
                                    Theta=theta(nnn);
                                    r= (u^2 + dn^2 - 2*u*dn*cos(Theta))^(alpha/2);
                                    if r==0
                                        T1(nnn)=0;
                                    else
                                        T1(nnn)= mun*r/(abs(mun*r - 1i*w1)); %mun*r*exp(-c*mun*r + 1i*w1*c);
                                        if T1(nnn)<10^-10
                                            T1(nnn)=0;
                                        end
                                    end
                                end
                                phi_Ij(nn)=1- x(nn)*(1-((trapz(theta,T1))/2*pi));
                            end
                            T2(n)=exp(-1i*w1*c)*(prod(phi_Ij));
                        end
                        Temp1= abs(trapz(T2));
                        R(mmm)=(2^(rho*lambda*Ab/(delta*W)))*(mu*(u^(alpha+1)))*...
                                    (exp(-mu*(u^alpha)*(2^(rho*lambda*Ab/(delta*W)) -1)*(sigma+c)))* Temp1;        
                   end
                   T3(mm)= trapz(R);
                end
                T4(m)=trapz(T3);
            end
            T5=trapz(T4);
            % since we approximate the true integral...
            r_pdf=T5*Ab*lambda*log(2)/A/delta/W;
            if r_pdf >1
                R_pdf(xx)=1;
            else
                R_pdf(xx)=r_pdf;
            end
        end
        
end
  % display(R_pdf)    
        RCP= 1- sum(R_pdf);
        if RCP>1
            RCP=1;
        end
        if RCP<0
            RCP=0;
        end
end

