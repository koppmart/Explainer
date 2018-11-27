classdef lightTreeRangeRobustMin
    
    properties
        vars %list of used variables
        rules% list of used rules 1 = <   0 = >=
        thresholds % list of threshold per used variable, stored on vars index
    end
    
    methods
        function LT=lightTreeRangeRobustMin(data, anomaly)
            X=[data;anomaly];
            data(isnan(data))=0;
            anomaly(isnan(anomaly))=0;
            w=range(X);
            mu=nanmean(X);
            C=bsxfun(@minus,data,anomaly);
            X=bsxfun(@minus,X,mu);
            X=bsxfun(@rdivide,X,w);
            A=C<0;
            B=C>0;
            LT.thresholds=nan(size(data,2),1);
            
            i=1;
            try
                while size(A,1)>0
                    a=nanmean(A,1);
                    b=nanmean(B,1);
                    ia=max(a);%choose the maximal split feature
                    ib=max(b);%choose the maximal split feature
                    
                    if ia==0 && ib==0 % should not happen but it happend.
                        break
                    end
                    if ia>=ib
                        Ia=find(a==ia);
                        if length(Ia)>1 % if there is more than one with the same quality choose the one with maximal robust minimum
                            Y=bsxfun(@minus,X(1:end-1,Ia),X(end,Ia));
                            Y=sort(abs(Y),'ascend');
                            y=Y(ceil(size(Y,1)*0.05),:);
                            [~,idx]=max(y);
                            Ia=Ia(idx);
                        end
                        LT.vars(i)=Ia; %store feature
                        LT.rules(i)='<'; %store rules x1>anomaly
                        Y=anomaly(Ia)-data(A(:,Ia),Ia);
                        Y=sort(Y,'ascend');
                        LT.thresholds(Ia)=anomaly(Ia)-0.25*Y(floor(length(Y)*0.05)+1);
                        
                        
                        B=B(A(:,Ia)==0,:);
                        A=A(A(:,Ia)==0,:);
                        A(:,Ia)=0; % disable already used feature
                        B(:,Ia)=0;
                    else
                        Ib=find(b==ib);
                        if length(Ib)>1 % if there is more than one with the same quality choose random one
                            Y=bsxfun(@minus,X(1:end-1,Ib),X(end,Ib));
                            Y=sort(abs(Y),'ascend');
                            y=Y(ceil(size(Y,1)*0.05),:);
                            [~,idx]=max(y);
                            Ib=Ib(idx);
                        end
                        LT.vars(i)=Ib;  %store feature
                        LT.rules(i)='>'; %store rules
                        Y=data(B(:,Ib),Ib)-anomaly(Ib); %count threshold on real data
                        Y=sort(Y,'ascend');
                        LT.thresholds(Ib)=0.25*Y(floor(length(Y)*0.05)+1)+anomaly(Ib); % and set threshold at
                        
                        
                        A=A(B(:,Ib)==0,:); % chose samples deviating the rule
                        B=B(B(:,Ib)==0,:);
                        A(:,Ib)=0; % disable already used feature
                        B(:,Ib)=0;
                    end
                    i=i+1;
                    
                    %          sprintf('%d - %d%d%d%d%d%d%d%d',LT.eq,LT.eqcnt)
                    
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
            rule='';
            if ~isempty(T.vars)
                for i=1:length(T.vars)
                    if i==1
                        rule=sprintf('data(:,%d)%s%f',T.vars(i),T.rules(i),T.thresholds(T.vars(i)));
                    else
                        rule=sprintf('%s & data(:,%d)%s%f',rule,T.vars(i),T.rules(i),T.thresholds(T.vars(i)));
                    end
                end
            else
                rule='';
            end
        end
        
    end %methods
end

