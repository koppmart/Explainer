classdef lightTree
    
    properties
        vars %list of used variables
        rules% list of used rules 1 = <   0 = >=       
    end
    
    methods
        function LT=lightTree(data, anomaly)            
            C=bsxfun(@minus,data,anomaly); 
            A=C<0;
            B=C>0;
            
            i=1;
         while size(A,1)>0             
            a=mean(A,1);
            b=mean(B,1);
            [ia,Ia]=max(a);%choose the maximal split feature
            [ib,Ib]=max(b);%choose the maximal split feature
            
            if ia==0 && ib==0
                break
            end
            
            if ia>ib
                LT.vars(i)=Ia; %store feature
                LT.rules(i)=0; %store rules
                
                B=B(A(:,Ia)==0,:);
                A=A(A(:,Ia)==0,:);                
                A(:,Ia)=0; % disable already used feature
                B(:,Ia)=0;
            else
                LT.vars(i)=Ib;  %store feature
                LT.rules(i)=1; %store rules
                
                A=A(B(:,Ib)==0,:); % chose samples deviating the rule              
                B=B(B(:,Ib)==0,:);
                A(:,Ib)=0; % disable already used feature
                B(:,Ib)=0;
            end
            i=i+1; 
            
            
         end
        
    end
    end
end

