%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    FINITE DIFFERENCE SOLUTION OF POISSON'S EQUATION:
%                   Vxx + Vyy = G		
%    USING THE METHOD OF SUCCESSIVE OVER-RELAXATION		
% 		
%    NX     :  NO. OF INTERVALS ALONG X-AXIS		
%    NY     :  NO. OF INTERVALS ALONG Y-AXIS		
%    A X B  :  DIMENSION OF THE SOLUTION REGION		
%    V(I,J) :  POTENTIAL AT GRID POINT (X,Y) = H*(I,J)		
%              WHERE I = 0,1,...,NX, J = 0,1,....,NY		
%    H      :  MESH SIZE		
%  *******************************************************		
  		
    A=1;B=1;
    V1=0;V2=10;V3=20;V4=-10;
    		
    % SPECIFY BOUNDARY VALUES AND NECESSARY PARAMETERS		
    NX= 20; %4 12 20	
    NY= NX;	
    H = A/NX;
    % SET INITIAL GUESS EQUAL TO ZEROS OR TO AVERAGE OF FIXED VALUES

    for I=1:NX-1
        for J=1:NY-1
            V(I+1,J+1)=(V1 + V2 + V3 + V4)/4.0;
        end
    end
    % SET POTENTIALS AT FIXED NODES		
    for I = 1:NX-1
        V(I+1,1)=V1;
        V(I+1,NY+1)=V3;
    end
    for J=1:NY-1
        V(1,J+1)=V4;
        V(NX+1,J+1)=V2;
    end
    V(1,1)=(V1 + V4)/2.0;
    V(NX+1,1)=(V1 + V2)/2.0;
    V(1,NY+1)=(V3 + V4)/2.0;
    V(NX+1,NY+1)=(V2 + V3)/2.0;
    % FIND THE OPTIMUM OVER-RELAXATION FACTOR		
    T = cos(pi/NX) + cos(pi/NY);
    W = ( 8 - sqrt(64 - 16*T^2))/(T^2);
    disp(['SOR Factor Omega = ',num2str(W)])    
    W4 = W/4;
    % ITERATION BEGINS		
    NCOUNT = 0;
       
    loop = 1;
    while loop == 1;
        RMIN = 0;
        for I =1:NX-1
            X = H*I;
            for J = 1:NY-1
                Y = H*J;
                G = -36.0*pi*X*(Y - 1.0);
                R = W4*( V(I+2,J+1) + V(I,J+1) + V(I+1,J+2) + V(I+1,J)-4.0*V(I+1,J+1) - G*H*H  );	
                RMIN = RMIN + abs(R);
                V(I+1,J+1) =  V(I+1,J+1) + R;
            end
        end
        RMIN = RMIN/(NX*NY);
        if(RMIN>=0.0001)            
            NCOUNT = NCOUNT + 1;
            if(NCOUNT>100)                
                loop = 0;
                disp('SOLUTION DOES NOT CONVERGE IN 100 ITERATIONS')
            end
        else
            %Then RMIN is less than .0001 and then solution has converged
            loop = 0; 
            disp(['Solution Converges in ',num2str(NCOUNT),' iterations'])
            disp(['h = ', num2str(H)])
        end 
    end

Vnum = V;

%Grab original points a through i
        abc = zeros(1,9);
        a_tic = 1;
        vec = [0:H:1];        
        for ii = .25:.25:.75            
            for jj = .25:.25:.75
                xind = find(vec==ii);
                yind = find(vec==jj);
                %disp([xind,yind])
                abc(a_tic) = Vnum(xind,yind);
                a_tic = a_tic + 1;
            end
        end


%     OUTPUT THE FINITIE DIFFERENCE APPROX. RESULTS		

%  ---------------------------------------------------------		
%      CALCULATE THE EXACT SOLUTION		
% 		
%      POISSON'S EQUATION WITH HOMOGENEOUS BOUNDARY CONDITIONS		
%      SOLVED BY SERIES EXPANSION		
% 		
    for I =1:NX-1
        X = H*I;
        for J = 1:NY-1
            Y = H*J;
            SUM = 0;
            for M = 1:10   % TAKE ONLY 10 TERMS OF THE SERIES	
                FM = M;
                for N = 1:10
                    FN = N;
                    FACTOR1 = (FM*pi/A)^2  +  (FN*pi/B)^2;
                    FACTOR2 = ( (-1)^(M+N) )*144*A*B/(pi*FM*FN);
                    FACTOR3 = 1 - (1 - (-1)^N)/B;
                    FACTOR = FACTOR2*FACTOR3/FACTOR1;
                    SUM = SUM + FACTOR*sin(FM*pi*X/A)*sin(FN*pi*Y/B);
                end
            end
            VH = SUM;
		
%      LAPLACE'S EQUATION WITH INHOMOGENEOUS BOUNDARY CONDITIONS		
%      SOLVED USING THE METHOD OF SEPARATION OF VARIABLES		
		
            C1=4*V1/pi;
            C2=4*V2/pi;
            C3=4*V3/pi;
            C4=4*V4/pi;
            SUM=0;
            for K =1:10  % TAKE ONLY 10 TERMS OF THE SERIES	
                N=2*K-1;
                AN=N;	
                A1=sin(AN*pi*X/B);
                A2=sinh(AN*pi*(A-Y)/B);
                A3=AN*sinh(AN*pi*A/B);
                TERM1=C1*A1*A2/A3;	
                B1=sinh(AN*pi*X/A);	
                B2=sin(AN*pi*Y/A);	
                B3=AN*sinh(AN*pi*B/A);
                TERM2=C2*B1*B2/B3;	
                D1=sin(AN*pi*X/B);	
                D2=sinh(AN*pi*Y/B);	
                D3=AN*sinh(AN*pi*A/B);
                TERM3=C3*D1*D2/D3;	
                E1=sinh(AN*pi*(B-X)/A);
                E2=sin(AN*pi*Y/A);	
                E3=AN*sinh(AN*pi*B/A);
                TERM4=C4*E1*E2/E3;	
                TERM = TERM1 + TERM2 + TERM3 + TERM4;
                SUM=SUM + TERM;
            end
            VI = SUM;	
            Vexact(I+1,J+1) = VH + VI;	
        end 
    end

%Grab original points a through i
        abc2 = zeros(1,9);
        a_tic = 1;
        vec = [0:H:1];        
        for ii = .25:.25:.75            
            for jj = .25:.25:.75
                xind = find(vec==ii);
                yind = find(vec==jj);
                %disp([xind,yind])
                abc2(a_tic) = Vexact(xind,yind);
                a_tic = a_tic + 1;
            end
        end
    
figure(1),
        imagesc(flipud(Vnum')),
        colorbar
        ylabel('y'),        xlabel('x')
        title('Example 3.4: Poisson PDE')
        
        
format short g
disp('     numerical     exact')
disp([abc' abc2'])    
