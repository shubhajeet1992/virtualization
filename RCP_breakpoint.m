function [breakpoints] = RCP_breakpoint(kappa, BS_No, X,BSLocationX,BSLocationY)
%RCP_BREAKPOINT returns breakpoints of RCP curve of BS mentioned in BS_No. X and Y are the
%coordinates of the coverage region of BS. Lx and Ly are the corrdinates of
%the BS location. d is the distance of neighbor BSs.
Delta=[0.2,0.4,0.6,0.8,1];
% BS_No=BS_No+1; % C++ to Matlab conversion
if (X(BS_No)==0)
    breakpoints=zeros(1,5);
else
%         BSLocationX = [1.6,0.6,1.1,1.4,1.8,2,1.1,0.3,0.3,0.6]; 
%         BSLocationY = [1.7,0.6,1.7,0.5,1.9,0.7,0.4,0.6,1.3,1]; 
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
        x=X([1:BS_No-1 BS_No+2:end]);
        Dist= (sqrt((BSLocationX(BS_No)-BSLocationX).^2 + (BSLocationY(BS_No)-BSLocationY).^2)); 
        pt=[BSLocationX(BS_No),BSLocationY(BS_No),0];
        D_edge=zeros(1,4);
        D_edge(1)= norm(cross(Vertc1-Vertc2,pt-Vertc2))/norm(Vertc1-Vertc2);
        D_edge(2)= norm(cross(Vertc2-Vertc3,pt-Vertc3))/norm(Vertc3-Vertc3); 
        D_edge(3)= norm(cross(Vertc3-Vertc4,pt-Vertc4))/norm(Vertc3-Vertc4); 
        D_edge(4)= norm(cross(Vertc4-Vertc1,pt-Vertc1))/norm(Vertc4-Vertc1); 
        d= nonzeros(Dist);
        q=min([nonzeros(X.*Dist);nonzeros(D_edge)]);
        % Ab=pi*q^2;
        if sum(X)==1
            Ab=A;
        elseif sum(X)==2
            Ab=A/2;
        else
            BSLocationX_new=nonzeros(X.*BSLocationX);
            BSLocationY_new=nonzeros(X.*BSLocationY);
            X_nz= find(X);
            ind=  find(X_nz== BS_No);
            AB=zeros(1,length(X));
            v=[];
            [v,cl] =(voronoin([BSLocationX_new, BSLocationY_new]));
            v1=[];
            v2=[];
            v1 = v(cl{ind},1)';
            v2 = v(cl{ind},2)';
            for ii=1:length(v1)
               if v1(ii)<0
                   v1(ii)=0;
               end
               if v1(ii)>2
                   v1(ii)=2;
               end
            end
            for ii=1:length(v2)
               if v2(ii)<0
                   v2(ii)=0;
               end
               if  v2(ii)>2
                   v2(ii)=2;
               end
            end
            if length(v1)<3
               v1(3)=2;
               v2(3)=2/2;
            end
            Ab = polyarea(v1,v2);
        end

        % 
        Max=10^50;
        Min=-Max;
        Gap=2*Max;
        y=10^20;
        theta=0:2*pi/y:2*pi;
        w=Min:Gap/y:Max;
        U=0:q/y:q;
        C=0:Max/y:Max;
        Rho=0:kappa/y:kappa;

        RR=[];
        phi_Ij=[];
        T1=[];
        T2=[];
        T3=[];
        T4=[];
        breakpoints=[];

        for xx=1:length(Delta)
            delta=Delta(xx);
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
                                    end
                                end
                                phi_Ij(nn)=1- x(nn)*(1-((trapz(theta,T1))/2*pi));
                            end
                            T2(n)=exp(-1i*w1*c)*(prod(phi_Ij));
                        end
                        Temp1= abs(trapz(T2));
                        %syms N;
        %                R(mmm)=double(symsum((N.*(lambda.^N)*(Ab.^N)./factorial(N)).*...
        %                                     (2.^(rho.*N./(delta*W))).*(mu*(u^(alpha+1)))*...
        %                                     (exp(-mu*(u^alpha).*(2.^(rho.*N/(delta*W)) -1)*(sigma+c)))* Temp1,N,0,Inf));
                        for N=1:10:200
                            RR(N)=(N.*(lambda.^N)*(Ab.^N)./factorial(N)).*...
                                            (2.^(rho.*N./(delta*W))).*(mu*(u^(alpha+1)))*...
                                            (exp(-mu*(u^alpha).*(2.^(rho.*N/(delta*W)) -1)*(sigma+c)))* Temp1;
                        end
                        R(mmm)=sum(RR);
                        end
                   T3(mm)= trapz(R);
                end
                T4(m)=trapz(T3);
            end
        T5=trapz(T4);
        point=1- (T5*exp(-Ab*lambda)*log(2)/A/delta/W);
        if point>1
          breakpoints(xx)=1;
        else
          breakpoints(xx)=point;
        end
        end
 end
end

