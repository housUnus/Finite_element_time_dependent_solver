classdef Solver1D_P2_Ref_T_EDP_a
    properties
        a;b; %the domain boundaries
        h; %  space step
        ht; % time step
        hadj; % step for the last two elements in case of adjustement
        n; % number of noeuds in space
        nt; % number of noeuds in time
        Ue; % exacte solution function
        fe; % seconde membre function-expression
        T; % the matrix of elements (connectivity matrix)
        X; % Column matrix of coordinates
        U;% The approximate function
        K;% Global matrix
        F;% Second member matrix
        er; % Error matrix
        ti; % The index of an element
        t0;% The initial moment
        tf; % Final moment
        cl1;% Column holding the type of left boundary condition and its value
        cl2;% Column holding the type of right boundary condition and its value
        timeindex; % The index of the current moment
        tInt; % Matrix columns of instants
        mrf; % The condition of the reference element
        fref; % The matrix of the type reference element (phi*phi)
        dfref; % The matrix of the type reference element (phi'*phi')
        mfref; % The matrix of the type reference element (phi*phi')
        dimEle; % The dimensions of the element
        gorder; % The order of the approximation of the gauss
        tstepSave; % Parameter to save the time step
         %% ******************************** equation paramters*******************
        v;
        D;
        sg;
        u0;
        alpha;
        %% end of block
 
    end
    
    methods
        function obj = Solver1D_P2_Ref_T_EDP_a(a,b,h,t0,tf,ht,Ue,mrf,order,alpha)
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
            obj = obj.descritize(); % calcule the number of nodes for P2 (odd number of nodes)r)
            obj = obj.Meshing(); % do the meshing
            obj = obj.adjust(); %
            obj.U = zeros(obj.n,obj.nt);
            obj.fref =obj.h.*obj.ref();
            obj.dfref = (1/obj.h).*obj.dref();
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
        end
        
       %% block for generating meshing of space and time
        function obj = descritize(obj)
            obj.n = floor((obj.b - obj.a)/obj.h) +1 ;
            if(rem(obj.n,2)==0)
                obj.n = obj.n+1; %if n is pair convert into impair and then recalcuate the step h 
                obj.h = (obj.b-obj.a)/obj.n;
            end
            obj.nt =  floor((obj.tf - obj.t0)/obj.ht) +1; 
        end
        
         function obj = Meshing(obj)
            obj.X = obj.a + obj.h*(0:obj.n-1)';
            obj.T =  [(1:2:obj.n-2)', (3:2:obj.n)'];
            obj.tInt = obj.t0 + obj.ht*(0:obj.nt-1)';
         end
        
        function obj = adjust(obj)
            %sometimes the final boundry is not reached so add 2 points one
            %in the middle of last element and one as the last element
            if(obj.X(obj.n)<obj.b)
                obj.n = obj.n+2;
                obj.hadj =  (obj.b-obj.X(obj.n-2))/2;
                obj.X(obj.n-1,1) = obj.X(obj.n-2,1) + obj.hadj;
                obj.X(obj.n,1) = obj.b;
                obj.T =  [(1:2:obj.n-2)', (3:2:obj.n)']; % add a one more element
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
            mi = zeros(3,3);
            for i=1:3
                    for j=1:3
                        f = polyval(obj.Cphi(i),xi).*polyval(obj.Cphi(j),xi);
                        mi(i,j) = sum(wi.*f);              
                    end
            end
         end
         %calculate reference element for derevitive
         function ki = dref(obj)

             [xi,wi]  = obj.gauss(obj.gorder);
             ki = zeros(3,3);
             for i=1:3
                     for j=1:3
                         f = obj.derv_ref(obj.Cphi(i),xi).*obj.derv_ref(obj.Cphi(j),xi);
                         ki(i,j) = sum(wi.*f);              
                     end
             end
         end
         function mki = mref(obj)

             [xi,wi]  = obj.gauss(obj.gorder);
             mki = zeros(3,3);
             for i=1:3
                     for j=1:3
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
       
                for i=1:3  %boucle sur les numéros locaux
                     for j=1:3  %boucle sur les numéros locaux
                           I=2*k+i-2;                            % numéros globaux dans K
                          J=2*k+j-2;   
                          obj.K(I,J) = obj.K(I,J) + elementKi(i,j);
                      end
               end
             end
     
         end
        
         function  ki =  localMatrice(obj,k)
                 
                %calcule la matrice élémentaire dans l'élément Ti
                %               Ti
                %       |-------|--------|
                %       x1      x2        x3     
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
                    [x1, ~ ,x3] = obj.nodes(k);
                    A = zeros(3);
                    B = zeros(3);
                    C = zeros(3);
                    [tildex,tildec]  = obj.Agauss(obj.gorder,x1,x3);  
                    for i=1:3
                        for j=1:3
                             f = (obj.alpha+obj.sg.*obj.ht + obj.tInt(obj.timeindex).*0.5.*obj.ht).*polyval(obj.phi_(i),tildex).*polyval(obj.phi_(j),tildex);
                             A(i,j) = sum(tildec.*f);
                             f = (obj.v.*obj.ht).*polyval(obj.phi_(i),tildex).*obj.derv_phi(obj.phi_(j),tildex);
                             B(i,j) = sum(tildec.*f);
                             f = (obj.ht.*obj.D).*obj.derv_phi(obj.phi_(i),tildex).*obj.derv_phi(obj.phi_(j),tildex);
                             C(i,j) = sum(tildec.*f);
                        end
                    end       
                end
                ki = A + B + C ;              
         end
       %% end of block of first part of equation
       
        %% buliding the second part of the equation
         function obj =  SecondMatrice(obj)
             t = size(obj.T);
             for k=1:t        
                 
                  elemFi= obj.LocalSecondMatrice(k);
                  
              for i=1:3  %boucle sur les numéros locaux
                    I=2*k+i-2;           
                   obj.F(I) = obj.F(I) + elemFi(i);
             end
             end
         end
         
         function elemFi =  LocalSecondMatrice(obj,k)      
             obj.ti = k;
             [x1, ~ , x3] = obj.nodes(k) ;
             [u1,u2,u3] = obj.getLastU(k,obj.timeindex-1);
             [tildex,tildec]  = obj.Agauss(obj.gorder ,x1,x3);
             A = zeros(3,1);
             B = zeros(3,1);
             C= zeros(3,1);
             
             for i=1:3
                    f =obj.alpha.*(u1.*polyval(obj.phi_(1),tildex) + u2.*polyval(obj.phi_(2),tildex) + u3.*polyval(obj.phi_(3),tildex)).*polyval(obj.phi_(i),tildex);
                    A(i) =  sum(tildec.*f);
                    f =obj.ht.*obj.fe(tildex,obj.tInt(obj.timeindex)).*polyval(obj.phi_(i),tildex);
                    B(i) = sum(tildec.*f);
                    f =  obj.ht.*0.5*obj.tInt(obj.timeindex).*obj.u0.*(polyval(obj.phi_(1),tildex) + polyval(obj.phi_(2),tildex) + polyval(obj.phi_(3),tildex)).*polyval(obj.phi_(i),tildex);
                    C(i) = sum(tildec.*f);             
             end
             elemFi = A + B - C;              
         end
         %% end of building    

         %% construction of shape functions // derevative of these functions
         function abc =phi_(obj,i)
             
             [x1, x2 , x3] = obj.nodes(obj.ti) ;
             if(i==1)
                 abc=[x1^2 x1 1;x2^2 x2 1; x3^2 x3 1]\[1; 0 ;0];
             elseif(i==2)
                 abc=[x1^2 x1 1;x2^2 x2 1; x3^2 x3 1]\[0; 1; 0];
             elseif(i==3)
                 abc=[x1^2 x1 1;x2^2 x2 1; x3^2 x3 1]\[0; 0; 1];       
             end   
         end
         
         function y = derv_phi(obj,abc,x)
            %Retourne la dérivée d'un polynome évaluée en un point
            y =polyval(polyder(abc),x);
         end
         
         %% construction of hat functions with their derevitives 
          %Calcul des phi chapeaux
         function abc = Cphi(obj,i)
              if(i==1)
                 abc=[1/2 -1/2 0];
             elseif(i==2)
                 abc=[-1 0 1];
             elseif(i==3)
                 abc=[1/2 +1/2 0];     
             end  
         end
         
         function y = derv_ref(obj,abc,x)
            %Retourne la dérivée d'un polynome évaluée en un point
            y =polyval(polyder(abc),x);
         end
         
         %% end of block of construction

         
         %% loock up functions to get (U or nodes     
         function [u1, u2, u3] = getLastU(obj,k,timeindex)
             u1 = obj.U(obj.T(k,1),timeindex);
             u2 = obj.U(obj.T(k,1)+1,timeindex);
             u3 = obj.U(obj.T(k,2),timeindex); 
         end
         
        function [x1, x2, x3] = nodes (obj,k)
             x1 = obj.X(obj.T(k,1));
             x2 = obj.X(obj.T(k,1)+1);
             x3 = obj.X(obj.T(k,2)); 
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
           obj.er = zeros(obj.n,obj.nt);
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
            figure('name', 'Comparaision: solutions exacte et approchée ');
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


