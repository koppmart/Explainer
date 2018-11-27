classdef explainerMax % with taking away succesfull dimensions.
    
    properties
        Trees  % T-1 cell array of tree stumps
        nvars  % #of features
        nTrees % num of trees
        stops %stop index for each anomaly
        ht % heights of tree
        q % anomal point quantile
        eventList
    end
    
    methods
        function AF=explainerMax(data, labels, subN, alpha,minTrees)
            if nargin<5
                minTrees=3;
            end
            
            %if k is bigger than data length, k=2/3 data size.
            if subN>size(data,1)
                subN=size(data,1)-floor(size(data,1)/3);
            end
            
            anomals=data(labels==1,:);
            normal=data(labels==0,:);
            la=size(anomals,1);
            AF.nvars=size(data,2);
            AF.nTrees=0;
            AF.stops=zeros(1,la);
            
            for i=1:la
                T={}; %empty tree
                q=1;
                fv=true(1,size(data,2)); % feature vector - default use all features
                out=randperm(size(normal,1));
                IDX=out(1:subN);
                try
                    gstk=0;
                    while (q>alpha && sum(fv)>1 || gstk<minTrees) % train trees on reduces dimension until anomaly point is in (1-alpha) quntile, distance wise
                        gstk=gstk+1;
                        %           if gstk>size(data,2)
                        %               T.vars
                        %           end
                        transfdata=zeros(subN+1,size(data,2));%Transform data from reduced dimensionality to full dimension
                        transfdata(:,fv)=[normal(IDX,fv);anomals(i,fv)];
                        
                        T=lightTreeRangeRobustMin(transfdata(1:subN,:),transfdata(subN+1,:));
                        
                        %check tree quality
                        mask=rand(1,subN)>0.5; %split the neighbourhood to halves
                        N1=transfdata(mask,:);
                        N2=transfdata(~mask,:);
                        da=knn_dist(transfdata(subN+1,:),N1); %count avg distance from anomaly to nearest points in N1
                        dn=knn_dist(N2,N1); %count distance from half of normal point to their nearest points in N1
                        
                        [~,V]=sort([da dn'],'ascend'); % sort distance
                        real_possition=find(V==1);
                        
                        q=real_possition/(length(dn)+1); %count quantile
                        
                        if q>alpha || gstk<minTrees %if good enough
                            AF.nTrees=AF.nTrees+1; %include tree into forest
                            AF.q(AF.nTrees)=q;
                            AF.Trees{AF.nTrees}=T;
                            AF.ht(AF.nTrees)=length(T.vars);
                            fv(T.vars(1))=false;     %remove root feature
                        end
                        
                    end
                    AF.stops(i)=AF.nTrees;%On this index ends trees for anomaly i
                catch Me
                    display(Me);
                end
            end
        end
        
        function [fi,fipa]=featureImportance(AF,o) %featureImportance and feature importance per Anomaly
            if nargin<2
                o='b'; %binary output
            end
            
            fi=zeros(1,AF.nvars); %basic feature importance
            
            if nargout >1
                fih=fi; % feature importance history
                fipa=zeros(length(AF.stops),AF.nvars); % feature importance per anomaly
                stopvar=1; % for iterating throught anomalies
            end
            
            
            if strcmp(o,'w') %weighted extraction
                for i=1:AF.nTrees
                    %Extract used variables from i.th tree
                    vars=AF.Trees{i}.vars;
                    for g=1:length(vars) %weighted
                        fi(vars(g))=fi(vars(g))+1/g;
                    end
                    if nargout>1 && i==AF.stops(stopvar)
                        fipa(stopvar,:)=fi-fih;%current weights - weights at last stop
                        fih=fi;%update stop
                        stopvar=stopvar+1;
                    end
                end
            else %binar extraction
                for i=1:AF.nTrees
                    %Extract used variables from i.th tree
                    vars=AF.Trees{i}.vars;
                    fi(vars)=fi(vars)+1;
                    if nargout>1 && i==AF.stops(stopvar)
                        fipa(stopvar,:)=fi-fih;%current weights - weights at last stop
                        fih=fi;%update stop
                        stopvar=stopvar+1;
                    end
                end
            end
        end
   
        % This function need tp be changed to use AF.stops instead of
        % AF.tpa
  %      function ruleSet=extractPerAnomalyRules(AF,minConf)
            
%             if nargin<2
%                 minConf=0.1;
%             end
%             
%             thr=nan(AF.tpa,AF.nvars);
%             usedVars=zeros(AF.tpa,AF.nvars);
%             rulesGreater=usedVars;
%             rulesLesser=usedVars;
%             
%             for i=1:AF.nTrees
%                 T=AF.Trees{i};
%                 r='';
%                 nvar=length(T.vars);
%                 if (nvar>0)
%                     for j=1:nvar
%                         if (T.rules(1)==1)
%                             s='<';
%                             rulesLesser(mod(i,AF.tpa)+1,T.vars(1))=1;
%                         else s='>';
%                             rulesGreater(mod(i,AF.tpa)+1,T.vars(1))=1;
%                         end
%                         usedVars(mod(i,AF.tpa)+1,T.vars)=1;
%                         thr(mod(i,AF.tpa)+1,T.vars)=T.thresholds(T.vars);
%                     end
%                 end
%                 if(mod(i,AF.tpa)==0)
%                     try
%                         %Aggregate rules into one or more
%                         
%                         clearvars C ca cb cbi cbCount
%                         [C,ca,cb]=unique(usedVars,'rows');
%                         for cbi=1:max(cb)
%                             cbCount(cbi)=sum(cb==cbi);
%                         end
%                         [v,V]=sort(cbCount,'descend');
%                         fi=1;
%                         while(fi<=length(v) && v(fi)/AF.tpa>minConf)
%                             line=C(V(fi),:);
%                             idcs=find(line>0);
%                             trs=nanmedian(thr(:,idcs));
%                             for si=1:length(idcs)
%                                 if (sum(rulesGreater(cb==(V(fi)),idcs))>sum(rulesLesser(cb==(V(fi)),idcs)))
%                                     rules(si)='>';
%                                 else
%                                     rules(si)='<';
%                                 end
%                             end
%                             r{fi}=explainerMax.constructRule(idcs,trs,rules,v(fi)/AF.tpa);
%                             fi=fi+1;
%                         end
%                         
%                         ruleSet{i/AF.tpa}=r;
%                         %clear vars and start for another anomaly
%                         thr=nan(AF.tpa,AF.nvars);
%                         usedVars=zeros(AF.tpa,AF.nvars);
%                         rulesGreater=usedVars;
%                         rulesLesser=usedVars;
%                         
%                     catch e
%                         disp(e);
%                     end
%                 end
%             end
%         end% End of extractPerAnomaly function
%         
       
        function [ruleSet,cbSets]=extractRulesCompact(AF,minConf,readability)
            % readability - h== human m==machine  c==confidenceBoost (feature subsets are returned)
            % minConf is minimal rule occurence in explanation
            
            if nargin<3
                readability='h';
            end
            if nargin<2
                minConf=0.1;
            end
            
            nvar=AF.nvars;
            thr=nan(AF.nTrees,nvar);
            usedVars=false(AF.nTrees,nvar);
            rulesGreater=usedVars;
            ruleSet={};
            
            %% Filling the matrcies I would love to kill this part
            for i=1:AF.nTrees
                T=AF.Trees{i};
                rulesGreater(i,T.vars(T.rules=='>'))=true;
                usedVars(i,T.vars)=true;
                thr(i,T.vars)=T.thresholds(T.vars);
            end
            
            %% Make unique rule sets out of them.
            cbSets=zeros(100,2*nvar);
            ruleNum=0;
            [C,~,cb]=unique([usedVars rulesGreater],'rows');
            uniqueRuleCount=max(cb);
            cbCount=histc(cb,1:uniqueRuleCount);
            [v,V]=sort(cbCount,'descend');
            fi=1;
            
            %% Finally extract rules
            while(fi<=length(v) && v(fi)/AF.nTrees>minConf)
                line=C(V(fi),1:nvar); %% line contains line with bin vector of used features
                idcs=find(line);      %% idcs are indeces of used features
                rules=zeros(1,length(idcs));
                mask=cb==(V(fi)); %% mask of all lines with the rule n. fi
                trs=median(thr(mask,idcs),1);
                
                ruleNum=ruleNum+1;
                if size(cbSets,1)<ruleNum
                    cbSets=[cbSets;zeros(ruleNum,2*nvar)];
                end
                
                for si=1:length(idcs)
                    if sum(rulesGreater(mask,idcs(si)))>0
                        rules(si)='>';
                        cbSets(ruleNum,idcs(si))=1;
                    else
                        rules(si)='<';
                        cbSets(ruleNum,idcs(si)+AF.nvars)=1;
                    end
                end
                
                ruleSet{end+1,1}=explainerMax.constructRule(idcs,trs,rules,v(fi)/AF.nTrees,readability);
                
                %% include inverse rules, should be used however it looks weird.
                if length(idcs)>1
                    ruleSet{end+1,1}=explainerMax.constructRule(idcs,trs,explainerMax.invers(rules),v(fi)/AF.nTrees,readability);
                end
                fi=fi+1;
            end
            cbSets=cbSets(ruleNum,:);
            
        end% End of extractcompactrules function
        
    end
    
    methods (Static= true)
        
        function rule=constructRule(idcs,thr,rules,conf,readability)
            if readability=='m'
                rule='';
                f1='%sdata(:,%d)%s%f';
                f2='%s & data(:,%d)%s%f';
            else
                rule=sprintf('(%0.1f%%)',conf*100);
                f1='%s x%d%s%0.2g';
                f2='%s AND x%d%s%0.2g';
            end
            
            for i=1:length(idcs)
                if i==1
                    rule=sprintf(f1,rule,idcs(i),rules(i),thr(i));
                else
                    rule=sprintf(f2,rule,idcs(i),rules(i),thr(i));
                end
            end
        end
        
         function invRules=invers(rules)
            mask=rules=='<';
            invRules(mask)='>';
            invRules(~mask)='<';
        end
        
    end
end

