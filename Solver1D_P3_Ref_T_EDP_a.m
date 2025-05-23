classdef Solver1D_P3_Ref_T_EDP_a
    properties
        a;b; %les limites de domaine
        h; % le pas de l'espace
        ht; % le pas du temps
        hadj; % le pas pour le derienrs deux elements au cas de l'adjustement.
        n; % le nombre des noeuds de l'espace
        nt; %le nombre des noeuds du temps
        Ue; % la solution exacte
        fe; % la fonction seconde membre
        T; % le matrice des elements (matrice de connectivit�)
        X; % matrice column des cordonn�es
        U;% la fonction approch� 
        K;% matrice global 
        F;% matrice seconde membre
        er; % matrice des erreurs
        ti; % l'index d'un element
        t0;% l'instant initiale
        tf; % l'instant finale
        cl1;% column contient le type de condition au limite gauche et sa valeur 
        cl2;% column contient le type de condition au limite droite et sa valeur 
        timeindex; % l'indice de l'instant
        tInt; % matrice columns des instants
        mrf; % la condition de l'�lement de r�ference
        fref; % la matrice de l'element ref�rence de type (phi*phi)
        dfref; % la matrice de l'element ref�rence de type (phi'*phi')
        mfref; % la matrice de l'element ref�rence de type (phi*phi')
        dimEle; % les dimensions de l'element
        gorder; % l'order du l'approximation du gauss
        tstepSave; % parametre pour sauvegarder le pas de temps
         %% ******************************** equation paramters*******************
        v;
        D;
        sg;
        u0;
        alpha;
        %% end of block
        
        %% helpful function
        compare;
    end
    
    methods
        function obj = Solver1D_P3_Ref_T_EDP_a(a,b,h,t0,tf,ht,Ue,mrf,order,alpha)
            disp('this work is a result of work of three memebres : AIT ALI EL HOSAYN, MOHAMED ABIDAR, HAITAM KHLIFI TAGHZOUTI');
            disp('the problem here can be solved using shape functions of type 3 with or without reference elements');
            obj.compare = @(n1,n2,n_dcp) round(n1*10^n_dcp) < round(n2*10^n_dcp); %% compare two doubles to avoid normal error comparing
            obj.a = a;
            obj.b = b; 
            obj.h = h; 
            obj.ht = ht;
            obj.hadj=0;
            obj.t0 = t0;
            obj.tf = tf;
            obj.Ue = Ue;
            obj.mrf = mrf;
            obj.gorder = order;
            obj = obj.descritize(); % calcule le nombre des noeuds pour P2 (nombre des noeuds impair)
            obj = obj.Meshing(); % faire le maillage
            obj = obj.adjust();
            obj.U = zeros(obj.n,obj.nt);
            obj.fref =(3/2)*obj.h.* obj.ref();
            obj.dfref = (2/3)*(1/obj.h).*obj.dref();
            obj.mfref = obj.mref();
            obj.dimEle = size(obj.T);
           %% scope for equation parameters
            obj.v = 1;
            obj.D= 1;
            obj.sg = 1;
            obj.alpha = alpha;
            obj.tstepSave = obj.ht;
            if(alpha == 0)
                obj.ht =1;
            end
            obj.fe = obj.getf();
            %% helfpul function
            
        end
        
       %% block for generating the mesh of space and time
        function obj = descritize(obj)
            obj.n = floor((obj.b - obj.a)/obj.h) +1 ;
            obj.nt =  floor((obj.tf - obj.t0)/obj.ht) +1;
        end
        function obj = descHelp(obj)
             sz = obj.n-1;
             if(rem(sz,3)~=0)
                obj.n = obj.n+1;  
                obj.h = (obj.b-obj.a)/obj.n;
                obj = obj.Meshing(); % comment me;
            end
        end
        
         function obj = Meshing(obj)
            obj.X = obj.a + obj.h*(0:obj.n-1)';
            obj.T =  [(1:3:obj.n-3)', (4:3:obj.n)'];
            obj.tInt = obj.t0 + obj.ht*(0:obj.nt-1)';
            obj = obj.descHelp(); %comment me
         end
        
         function obj = adjust(obj)
            %sometimes the final boundry is not reached so add 2 points one
            %in the middle of last element and one as the last element
            if(obj.compare(obj.X(obj.n,1),obj.b,2))
                obj.n = obj.n+3;
                obj.hadj =  (obj.b-obj.X(obj.n-3))/3;
                obj.X(obj.n-2,1) = obj.X(obj.n-3,1) + obj.hadj;
                obj.X(obj.n-1,1) = obj.X(obj.n-3,1) + 2*obj.hadj;
                obj.X(obj.n,1) = obj.b;
                obj.T =  [(1:3:obj.n-3)', (4:3:obj.n)']; % add one more element
            end
        end
      %% end of block  
      
       %% Helpful Function to get second membre 
            function fe = getf(obj)
                syms x t;
                syms f1(x,t) f2(x,t);
                ue = sym(obj.Ue);
                f1(x,t) =obj.alpha*diff(ue,t) + obj.v*diff(ue,x) - obj.D*diff(ue,x,2) +obj.sg*ue;
                f2(x,t) =int(ue,t,0,t);
                f1(x,t) = f1(x,t) + f2(x,t);
                fe = matlabFunction(f1) ;
            end
            
       %% block where building of reference matrices for only one time 
       %calculate refrence element matrice
         function mi = ref(obj)
            [xi,wi]  = obj.gauss(obj.gorder);
            mi = zeros(4,4);
            for i=1:4
                    for j=1:4
                        f = polyval(obj.Cphi(i),xi).*polyval(obj.Cphi(j),xi);
                        mi(i,j) = sum(wi.*f);              
                    end
            end
         end
         %calculate reference element for derevitive
         function ki = dref(obj)

             [xi,wi]  = obj.gauss(obj.gorder);
             ki = zeros(4,4);
             for i=1:4
                     for j=1:4
                         f = obj.derv_ref(obj.Cphi(i),xi).*obj.derv_ref(obj.Cphi(j),xi);
                         ki(i,j) = sum(wi.*f);              
                     end
             end
         end
         function mki = mref(obj)

             [xi,wi]  = obj.gauss(obj.gorder);
             mki = zeros(4,4);
             for i=1:4
                     for j=1:4
                         f = polyval(obj.Cphi(i),xi).*obj.derv_ref(obj.Cphi(j),xi);
                         mki(i,j) = sum(wi.*f);              
                     end
            end
         end
         %% end of block
      
      %% loop through time
        function obj = LoopTime(obj)
            for index=2:obj.nt
                obj.timeindex = index;
                obj.K = zeros(obj.n);    %matrice of elements
                obj.F = zeros(obj.n,1); %matrice of second membre
                obj = obj.assembly();
                obj = obj.CL(obj.cl1,obj.cl2);
                obj = obj.Solve();
            end           
        end
        %% end of looping through time
        
        %% building of the first part of the equation
         function obj = globalMatrice(obj)
            %loop through each element
            t = size(obj.T); %number of elements
             for k=1:t
                  elementKi = obj.localMatrice(k);
                  
                for i=1:4  %boucle sur les num�ros locaux
                     for j=1:4  %boucle sur les num�ros locaux
                          I=3*k+i-3;                            % num�ros globaux dans K
                          J=3*k+j-3;   
                          obj.K(I,J) = obj.K(I,J) + elementKi(i,j);
                      end
               end
             end
     
         end
        
         function  ki =  localMatrice(obj,k)
                mref = obj.fref;
                dmref = obj.dfref;
                mxref = obj.mfref;
                if(obj.mrf == 1)
                    if(obj.hadj>0) %if the space step is modified at the last element then calcualte the reference matrices based on the new step
                        if(k == obj.dimEle(1,1))
                             mref = (1/obj.h).*obj.hadj.*obj.fref;
                             dmref =(1/obj.hadj).*obj.h.*obj.dfref;
                        end
                    end
                    A = (obj.alpha+obj.sg.*obj.ht + obj.tInt(obj.timeindex).*0.5.*obj.ht ).*mref;
                    B = obj.ht.*obj.D.*dmref;
                    C = obj.v.*obj.ht.*mxref;
                    elseif(obj.mrf ==0)
                    obj.ti = k;
                    [x1, ~ ,~,x4] = obj.nodes(k);
                    A = zeros(4);
                    B = zeros(4);
                    C = zeros(4);
                    [tildex,tildec]  = obj.Agauss(obj.gorder,x1,x4);  
                    for i=1:4
                        for j=1:4
                             f = (obj.alpha+obj.sg.*obj.ht + obj.tInt(obj.timeindex).*0.5.*obj.ht).*polyval(obj.phi_(i),tildex).*polyval(obj.phi_(j),tildex);
                             A(i,j) = sum(tildec.*f);
                             f = (obj.v.*obj.ht).*polyval(obj.phi_(i),tildex).*obj.derv_phi(obj.phi_(j),tildex);
                             B(i,j) = sum(tildec.*f);
                             f = (obj.ht.*obj.D).*obj.derv_phi(obj.phi_(i),tildex).*obj.derv_phi(obj.phi_(j),tildex);
                             C(i,j) = sum(tildec.*f);
                        end     
                    end
                end
                ki = A + B + C;
                
         end
       %% end of block of first part of equation
       
        %% buliding the second part of the equation
         function obj =  SecondMatrice(obj)
             t = size(obj.T);
             for k=1:t        
                 
                  elemFi= obj.LocalSecondMatrice(k);
                  
              for i=1:4  %boucle sur les num�ros locaux
                    I=3*k+i-3;           
                   obj.F(I) = obj.F(I) + elemFi(i);
             end
             end
         end
         
         function elemFi =  LocalSecondMatrice(obj,k)      
             obj.ti = k;
             [x1, ~ ,~, x4] = obj.nodes(k) ;
             [u1,u2,u3,u4] = obj.getLastU(k,obj.timeindex-1);
             [tildex,tildec]  = obj.Agauss(obj.gorder ,x1,x4);
             A = zeros(4,1);
             B = zeros(4,1);
             C = zeros(4,1);
             for i=1:4                
                    f =obj.alpha.*(u1.*polyval(obj.phi_(1),tildex) + u2.*polyval(obj.phi_(2),tildex) + u3.*polyval(obj.phi_(3),tildex) +  u4.*polyval(obj.phi_(4),tildex)).*polyval(obj.phi_(i),tildex);
                    A(i) =  sum(tildec.*f);
                    f =obj.ht.*obj.fe(tildex,obj.tInt(obj.timeindex)).*polyval(obj.phi_(i),tildex);
                    B(i) = sum(tildec.*f);
                    f =  obj.ht.*0.5*obj.tInt(obj.timeindex).*obj.u0.*(polyval(obj.phi_(1),tildex) + polyval(obj.phi_(2),tildex) + polyval(obj.phi_(3),tildex) + polyval(obj.phi_(4),tildex)).*polyval(obj.phi_(i),tildex);
                    C(i) = sum(tildec.*f);
             end

             elemFi = A + B - C;         
         end
         %% end of building    

         %% construction of shape functions // derevative of these functions
         function abc =phi_(obj,i)
             [x1, x2 , x3, x4] = obj.nodes(obj.ti);
             if(i==1)
                 abc=[x1^3 x1^2 x1 1;x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1]\[1; 0 ;0 ; 0];
             elseif(i==2)
                 abc=[x1^3 x1^2 x1 1;x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1]\[0; 1 ;0 ; 0];
             elseif(i==3)
                 abc=[x1^3 x1^2 x1 1;x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1]\[0; 0 ;1 ; 0];   
             elseif(i==4)
                 abc=[x1^3 x1^2 x1 1;x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1]\[0; 0 ;0 ; 1];				 
             end   
         end
         
         function y = derv_phi(obj,abc,x)
            %Retourne la d�riv�e d'un polynome �valu�e en un point
            y =polyval(polyder(abc),x);
         end
         %% end of block of construction
         
         %% construction of hat functions with their derevitives 
          %Calcul des phi chapeaux
         function abc = Cphi(obj,i)
              if(i==1)
                 abc=[-9/16 9/16 1/16 -1/16];
             elseif(i==2)
                 abc=[27/16 -9/16 -27/16 9/16];
             elseif(i==3)
                  abc=[-27/16 -9/16 27/16 9/16];     
             elseif(i==4)
                 abc=[9/16 9/16 -1/16 -1/16];     
             end        
         end
         
         function y = derv_ref(obj,abc,x)
            %Retourne la d�riv�e d'un polynome �valu�e en un point
            y =polyval(polyder(abc),x);
         end
         
         %% loock up functions to get (U or nodes     
         function [u1, u2, u3, u4] = getLastU(obj,k,timeindex)
             u1 = obj.U(obj.T(k,1),  timeindex);
             u2 = obj.U(obj.T(k,1)+1,timeindex);
             u3 = obj.U(obj.T(k,1)+2,timeindex); 
			 u4 = obj.U(obj.T(k,1)+3,timeindex);
         end
         
        function [x1, x2, x3, x4] = nodes (obj,k)
             x1 = obj.X(obj.T(k,1));
             x2 = obj.X(obj.T(k,1)+1);
             x3 = obj.X(obj.T(k,1)+2); 
			 x4 = obj.X(obj.T(k,1)+3);
        end      
         %% end of block
         
         %% gauss block to get gauss point with weights
         
          function [tildex,tildec] = Agauss(obj,order,a,b)
            [xi,wi] = obj.gauss(order);
            tildec = ((b-a)/2)*wi;
            tildex = ((b-a)*xi)/2 + (b+a)/2;
         end
          function [xi,wi] = gauss(obj,order)
              %not used only the step, because it might be a non homogene
              %step at the last element
              if(order ==2)
                 xi = [1/sqrt(3) -1/sqrt(3) ];
                 wi = [1 1];
              elseif(order ==3)
                 xi = [0.7745966692414834, 0, -0.7745966692414834];
                 wi = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
              elseif(order ==4)
                  xi= [-0.3399810435848563, 0.3399810435848563,-0.8611363115940526, 0.8611363115940526] ; 
                  wi=[ 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538] ;
              elseif(order ==5)
                 xi = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
                 wi = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
              elseif(order ==6)
                 xi = [0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969,-0.9324695142031521, 0.9324695142031521 ];
                 wi = [0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 	0.1713244923791704, 0.1713244923791704];
              end
          end
          
           %% block to handle the boundry conditions
          function obj = CLS(obj,cl1,cl2)
              obj.cl1 = cl1;
              obj.cl2 = cl2;
          end
             %condition limite section **************
         function obj = CL(obj,cl1,cl2)
             %test the CL in the first point
             if(cl1(1)==0)
             %element 1 de type direchlet
                 obj = obj.drchlet(1,cl1);
             elseif(cl1(1)==1)
             %element 1 de type neuman
                 obj = obj.neuman(1,cl1); %the last element must be substracted from F (worked with +du2/dt not -du2/dt
             end
             
             if(cl2(1)==0)
             %element n de type direchlet
                 obj = obj.drchlet(obj.n,cl2);
             elseif(cl2(1)==1)
             %element n de type neuman
                obj = obj.neuman(obj.n,cl2.*-1); %the last element must be added to F (worked with +du2/dt not -du2/dt
 
             end
         end
         
         function obj = drchlet(obj,line,cl)
               obj.K(line,1:obj.n) = 0;
               obj.K(line,line) = 1;
               obj.F(line) = cl(2);
         end
         
         function obj = neuman(obj,line,cl)
               obj.F(line) =obj.F(line) - cl(2).*obj.ht.*obj.D;
         end
        %% end of block of CLS
        
         %% block to set the initial condtion
        %intilze the frist time solution   
        function obj = initialize(obj,u0)
            obj.U(:,1) = ones(obj.n,1).*u0;
            obj.u0 = u0;

        end
        %% end of block
        
        %% block of function of assembling the matrices
         function obj = assembly(obj)
            obj = obj.globalMatrice();
            obj = obj.SecondMatrice();
         end
         %% end of block
        
         %% solve the linear system
        function obj = Solve(obj)
            obj.U(:,obj.timeindex) =obj.K\obj.F;
        end
         %% end of block
         
         %% block where to calculate the error for each instance
        function obj =  erreur(obj)
           obj.er = zeros(obj.n,obj.nt)
            for tt = 1:obj.nt
            obj.er(:,tt) =obj.U(:,tt) -  obj.Ue(obj.X(:),obj.tInt(tt));
            end
        end
        %% end of block
                %% Calculate the mean root square for each instant 
        function mrsr =  mrs(obj)
           mrsr = zeros(1,obj.nt);
            for tt = 1:obj.nt
            mrsr(1,tt) =immse(obj.U(:,tt) , obj.Ue(obj.X(:),obj.tInt(tt)));
            end
        end
        %% end of block
    %% display the graph of excate solution with the approximated one
        function dispf(obj)
            figure('name', 'Comparaision: solutions exacte et approch�e ');

            i = 1;
            for t = obj.t0:obj.tstepSave:obj.tf
                title(['t =  ' num2str(t)])
                clf;
                plot(obj.X,obj.Ue(obj.X(:),t),'b');
                hold on; 
                plot(obj.X,obj.U(:,i),'r.-')
                i = i+1;
                axis( [min(obj.X), max(obj.X) min(min(obj.U)) max(max(obj.U))]);
                
                title(['t =  ' num2str(t)])
                ylabel('U/Ue');
                xlabel('noeuds');
                pause( 0.1);    
            end

        end
        
        function dispe(obj)
            figure('name','variation du l erreur avec le temps et l espace');

            i = 1;
            for t = obj.t0:obj.tstepSave:obj.tf
                clf;
                plot(obj.X,obj.er(:,i),'r.-')
                i = i+1;
                axis( [min(obj.X), max(obj.X) min(min(obj.er)) max(max(obj.er))]);
                title(['t =  ' num2str(t)])
                ylabel('erreur');
                xlabel('noeuds');   
                pause(.01);
            end
        end
        %% end of block
        %% get columns of U exacte
       function uev = UeV(obj)
            uev = zeros(obj.n,obj.nt);
            for tt = 1:obj.nt
            uev(:,tt) =obj.Ue(obj.X(:),obj.tInt(tt));
            end   
       end
    end
end
%%


