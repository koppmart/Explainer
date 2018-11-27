classdef lightTreeRange
    
    properties
        vars %list of used variables
        rules% list of used rules 1 = <   0 = >=   
        thresholds % list of threshold per used variable, stored on vars index
     end
    
    methods
        function LT=lightTreeRange(data, anomaly)            
            X=[data;anomaly];
            w=range(X);
            mu=mean(X);
            C=bsxfun(@minus,data,anomaly); 
            X=bsxfun(@minus,X,mu);
            X=bsxfun(@rdivide,X,w);
            A=C<0;
            B=C>0;
            LT.thresholds=nan(size(data,2),1);
             
            i=1;
            try
         while size(A,1)>0             
            a=mean(A,1);
            b=mean(B,1);
            ia=max(a);%choose the maximal split feature
            ib=max(b);%choose the maximal split feature
            
            if ia==0 && ib==0 % should not happen but it happend.
                break
            end
            if ia>=ib
                Ia=find(a==ia);
                if length(Ia)>1 % if there is more than one with the same quality choose random one
                   Y=bsxfun(@minus,X(1:end-1,Ia),X(end,Ia));
                   Y=abs(Y);
                   Y=min(Y);
                   [~,idx]=max(Y);
                   Ia=Ia(idx);
                end
                LT.vars(i)=Ia; %store feature
                LT.rules(i)=0; %store rules
                Y=anomaly(Ia)-data(:,Ia);
                Y=sort(Y,'ascend');
                LT.thresholds(Ia)=anomaly(Ia)-Y(floor(size(data,1)*0.05));              
                
                
                B=B(A(:,Ia)==0,:);
                A=A(A(:,Ia)==0,:);                
                A(:,Ia)=0; % disable already used feature
                B(:,Ia)=0;
            else
                Ib=find(b==ib);
                if length(Ib)>1 % if there is more than one with the same quality choose random one                    
                   Y=bsxfun(@minus,X(1:end-1,Ib),X(end,Ib));
                   Y=abs(Y);
                   Y=min(Y);
                   [~,idx]=max(Y);
                   Ib=Ib(idx);
                end
                LT.vars(i)=Ib;  %store feature
                LT.rules(i)=1; %store rules
                Y=data(:,Ib)-anomaly(Ib); %count threshold on real data
                Y=sort(Y,'descend');
                LT.thresholds(Ib)=Y(floor(size(data,1)*0.05))+anomaly(Ib); % and set threshold at 90
                
                
                A=A(B(:,Ib)==0,:); % chose samples deviating the rule              
                B=B(B(:,Ib)==0,:);
                A(:,Ib)=0; % disable already used feature
                B(:,Ib)=0;
            end
            i=i+1; 
                       
         end
            catch Me
                display(Me)
            end
        
        end
    
        function y=eval(T,data)
            
            rule=T.constructRule;
            y=eval(rule);
            
        end
    
        function rule=constructRule(T)
            try
                rule='';
                
            if length(T.vars)>0
                for si=1:length(T.vars)
                            if (T.rules(si)==0)
                                rules(si)='>';
                            else
                                rules(si)='<';
                            end
                        end
                
                for i=1:length(T.vars)
                    if i==1
                        rule=sprintf('data(:,%d)%s%f',T.vars(i),rules(i),T.thresholds(T.vars(i)));
                    else
                        rule=sprintf('%s & data(:,%d)%s%f',rule,T.vars(i),rules(i),T.thresholds(T.vars(i)));
                    end
                end
            else
                rule='';
                
            end
            catch e
                disp(e);
            end
        end
    end
end

