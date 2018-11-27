classdef explainerMin % number of trees per anomaly as parameter, using lightTrees with maxMargin
    
    properties
        Trees  % T-1 cell array of tree stumps
        nvars  % #of features
        nTrees % num of trees
        stops %stop index for each anomaly
        tpa % number of trees trained per anomaly
        ht %
    end
    
    methods
        function AF=explainerMin(data, labels, k, tpa)
            
            data(isnan(data))=0;
            anomals=data(labels==1,:);
            normal=data(labels==0,:);
            la=size(anomals,1);
            AF.nvars=size(data,2);
            AF.nTrees=0;
            AF.stops=zeros(1,la);
            AF.tpa=tpa;
            
            
            %if k is bigger than data length, k=3/4 data size.
            if k>size(normal,1)
                k=size(data,1)-floor(size(data,1)/3);
            end
            
            for i=1:la
                T={}; %empty tree
                
                for j=1:tpa
                    
                    out=randperm(size(normal,1));
                    IDX=out(1:k);
                    
                    %T=lightTree(normal(IDX,:),anomals(i,:));
                    %T=lightTreeRange(normal(IDX,:),anomals(i,:));
                    T=lightTreeRangeRobustMin(normal(IDX,:),anomals(i,:));
                    
                    AF.nTrees=AF.nTrees+1; %include tree into forest
                    AF.Trees{AF.nTrees}=T;
                    AF.ht(AF.nTrees)=length(T.vars);
                end
                AF.stops(i)=AF.nTrees;%On this index ends trees for anomaly i
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
        
        function fi=featureImportanceUsingCustomTreeNumber(AF,nt)
            fi=zeros(1,AF.nvars); %basic feature importance
            idx=[];
            tpa=AF.stops(1);
            for s=1:length(AF.stops)
                start=AF.stops(s)-tpa+1;
                idx=[idx start:start+nt-1];
            end
            
            for i=1:length(idx)
                %Extract used variables from i.th tree
                vars=AF.Trees{idx(i)}.vars;
                fi(vars)=fi(vars)+1;
            end
        end
        
        function [ruleSet,cbSets]=extractPerAnomalyRules(AF,minConf,readability)
            if nargin<3
                readability='h';
            end
            if nargin<2
                minConf=0.1;
            end
            
            cbSets=zeros(1,2*AF.nvars);
            thr=nan(AF.tpa,AF.nvars);
            usedVars=zeros(AF.tpa,AF.nvars);
            rulesGreater=usedVars;
            rulesLesser=usedVars;
            
            for i=1:AF.nTrees
                T=AF.Trees{i};
                r='';
                nvar=length(T.vars);
                if (nvar>0)
                    for j=1:nvar
                        if (T.rules(j)=='<')
                            rulesLesser(mod(i,AF.tpa)+1,T.vars(1))=1;
                        else
                            rulesGreater(mod(i,AF.tpa)+1,T.vars(1))=1;
                        end
                    end
                    usedVars(mod(i,AF.tpa)+1,T.vars)=1;
                    thr(mod(i,AF.tpa)+1,T.vars)=T.thresholds(T.vars);
                end
                if(mod(i,AF.tpa)==0)
                    try
                        %Aggregate rules into one or more
                        clearvars C ca cb cbi cbCount
                        [C,~,cb]=unique(usedVars,'rows');
                        for cbi=1:max(cb)
                            cbCount(cbi)=sum(cb==cbi);
                        end
                        [v,V]=sort(cbCount,'descend');
                        fi=1;
                        while(fi<=length(v) && v(fi)/AF.tpa>minConf)
                            cbSets(end+1,:)=zeros(1,2*AF.nvars);
                            line=C(V(fi),:);
                            idcs=find(line>0);
                            trs=nanmedian(thr(:,idcs));
                            for si=1:length(idcs)
                                if (sum(rulesGreater(cb==(V(fi)),idcs))>sum(rulesLesser(cb==(V(fi)),idcs)))
                                    rules(si)='>';
                                    cbSets(end,idcs(si))=1;
                                else
                                    rules(si)='<';
                                    cbSets(end,idcs(si)+AF.nvars)=1;
                                end
                            end
                            
                            if readability=='m'
                                r{fi}=explainerMin.constructRuleForEval(idcs,trs,rules,readabilityA);
                            else
                                r{fi}=explainerMin.constructRule(idcs,trs,rules,v(fi)/AF.tpa,readability);
                            end
                            fi=fi+1;
                        end
                        
                        ruleSet{i/AF.tpa}=r;
                        %clear vars and start for another anomaly
                        thr=nan(AF.tpa,AF.nvars);
                        usedVars=zeros(AF.tpa,AF.nvars);
                        rulesGreater=usedVars;
                        rulesLesser=usedVars;
                        
                    catch e
                        disp(e);
                    end
                end
            end
            cbSets=cbSets(2:end,:);
        end% End of extractPerAnomaly function
        
       
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
                
                ruleSet{end+1,1}=explainerMin.constructRule(idcs,trs,rules,v(fi)/AF.nTrees,readability);
                
                %% include inverse rules, should be used however it looks weird.
                if length(idcs)>1
                    ruleSet{end+1,1}=explainerMin.constructRule(idcs,trs,explainerMin.invers(rules),v(fi)/AF.nTrees,readability);
                end
                fi=fi+1;
            end
            cbSets=cbSets(ruleNum,:);
            
        end% End of extractcompactrules function
        
        
        function [y,Y]=eval(AF,data)
            Y=zeros(size(data,1),AF.nTrees);
            for i=1:AF.nTrees;
                a=AF.Trees{i}.eval(data);
                Y(:,i)=a;
            end
            y=mean(Y,2);%>0.5;
        end
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

