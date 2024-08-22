function [weights, criteria, altweights, feedbackweights, combination] = ahp(M,C,varargin)
%
% Analytic Hierarchy Process (AHP) is a simple technique for organizing and analyzing complex multi-objective decisions.
% It combines both quantitative and qualitative analysis elements and it has particular application in group decision making.
% The philosophy of the technique is to decompose problem into a hierarchy of more easily understood sub-problems, each of which can be analyzed independently.
% Once the hierarchy is built, the decision makers systematically evaluate its various elements by comparing them to one another two at a time,
% with respect to their impact on an element above them in the hierarchy.
% The AHP converts these evaluations to numerical values that can be processed and compared over the entire range of the problem.
% A numerical weight is derived for each element of the hierarchy, allowing diverse and often incommensurable elements to be compared to one another in a rational and consistent way.
% In the final step of the process, numerical weights are calculated for each of the decision alternatives.
% These weights represent the alternatives' relative ability to achieve the goal.
%
% The function facilitates the following:
% - Simple AHP implementation
% - Multiple decision makers
% - Analytic Network Process (ANP): The generalization of the AHP, which incorporates dependences and feedbacks between decision criteria and options.
% - Fuzzy AHP and ANP: This are special versions of the simple AHP and ANP, which find application in fuzzy environments, where the relative importance of the decision criteria and the alternatives is uncertain.
% - Simulation: A Monte-Carlo simulation-based approach of AHP and ANP,
% which allows to compare distributions of weights and performs sensitivity
% analysis.
% - Cost-Benefit analysis: The benefit (AHP weights) in relationship with the cost of the respective option.
% - Optimization: In case of a resource allocation problem, the function estimates the optimal feasible combination of alternatives subject to the resources constraints.
% - Prediction combination: In case this is a forecasting combination problem, the function generates a weighted average forecast, using the combination weights and the individual forecasts as inputs.
%
% Author: Apostolos Panagiotopoulos
% Published: August 2020
%
%
% Inputs:
%
%   M:               An m-by-c matrix, where m is the pairwise comparison of all alternatives
%                    subject to an underline criterion, and c is the number of criteria. If there are more than one decision makers
%                    c is the number of criteria multiplied by the number of decision makers.
%                    e.g. If there are 3 criteria and 2 decisions makers:
%                    The first three columns are the pairwise comparison of
%                    the alternatives according to each criterion subject
%                    to the opinion of the first decision maker. The last
%                    three columns are the counterpart opinion of the second
%                    decision maker. If there is only one decision maker, c
%                    is the number of decision criteria.
%                    The decision makers have to pick from scale 1 to 9 how much more important is a criterion in comparison with its peers (1: equal importance to 9: extremely more important.
%                    If the underline criterion is less important than a peer, we use the reverse scale (1: equal importance to 1/9: extremely less important).
%                    d is the number of decision makers. If the decision
%                    maker is willing to use 'beneficial' or
%                    'nonbeneficial' inputs only (see optional inputs), M
%                    should be left empty [].
%
%   C:               An n-by-d matrix, where n is the pairwise comparison
%                    of the hierarchy criteria (in the same way with the alternatives
%                    comparison), and d is the number of decision makers.
%                    If there is only one decision maker, C is a vector
%                    n-by-1. If there is only one decision criterion, the input should be left
%                    empty [].
%
%
% Optional Inputs:
%
%   'beneficial':    If a criterion is expressed quantitatively (e.g. salary),
%                    decision makers can include the actual observed values,
%                    instead of doing a pairwise comparison
%                    for all the alternatives. 'beneficial' is used for criteria with
%                    beneficial value for the alternative options (the
%                    higher the better), e.g. storage space or a battery life expectancy.
%                    It should be a k-by-g matrix, where k is the
%                    value of a criterion for every option
%                    and g is the number of criteria.
%
%   'nonbeneficial': Like 'beneficial', 'nonbeneficial' should be
%                    used for quantitative criteria with non-beneficial value
%                    for the alternative options (the lower the better),
%                    e.g. maintenance cost or energy consumption. It should
%                    be a k-by-g matrix, where k is the value a criterion for every option
%                    and g is the number of criteria.
%
%   'feedback':      In case we want to use the Analytic Network Method.
%                    An n-by-k matrix, where n is the a pairwise comparison
%                    of the hierarchy criteria per alternative option,
%                    and k is the number of alternative options multiplied
%                    by the number of decision makers. If there is only one
%                    decision maker, k is the number of alternatives.
%
%   'fuzzy':         A true/false Boolean variable. If true, the function runs
%                    Fuzzy AHP or ANP instead of the simple counterpart. The default is false.
%                    Fuzzy options do not apply on quantitative ('beneficial' or
%                    'nonbeneficial') criteria, since the performance of
%                    the alternatives are neither subjective nor uncertain.
%                    However, it applies to their pairwise comparison of importance
%                    their pairwise comparison against an alternative in case of ANP.  
%
%   'names':         An o-by-1 cell array with the names of the alternative
%                    options/decisions. If empty the alternatives are named 'Option 1',
%                    'Option 2', etc.
%
%   'criterianames': An n-by-1 cell array with the names of the decision criteria. If empty the alternatives are named 'Criterion 1',
%                    'Criterion 2', etc.
%
%   'normalization': In case quantitative criteria ('beneficial' and 'nonbeneficial') are used,
%                    their normalization method should be specified. The
%                    alternative specifications are 'linear', 'relativelinear', 'vector', 'minmax', 'enhanced', and
%                    'logarithmic'. Default option is 'enhanced'. 'linear' is
%                    preferred when we are interested in the absolute values
%                    of the performance of an alternative option under a criterion, while
%                    'enhanced' is preferred when we are interested in the
%                    performance of an option, relatively to the
%                    counterpart performances of the rest of the options.
%                    All other specifications stand somewhere in between.
%
%   'simulation':    A true/false Boolean variable. If true, the model runs
%                    a simulation analysis around the original scores set by
%                    decision makers. The default is false.
%
%   'percentiles':   In case of 'simulation' this is a vertical vector of the
%                    calculation percentiles the decision makers want to include in the output.
%                    If empty, the function provides the average only.
%
%   'cost':          A vector o-by-a is the cost related with the
%                    alternative options. a is the number of costs for one option
%                    e.g. for a car we may have price and maintenance cost.
%                    If there is no 'cost' input, the
%                    function assumes that the cost of each option is 1.
%
%   'costbenefit':   A true/false Boolean variable. If true, the function shows
%                    a cost-benefit of the alternative options.
%                    i.e. the AHP weight estimation of each option divided
%                    by the counterpart cost. The cost-benefit analysis is
%                    performed only if the 'cost' input is used. The default is false.
%
%   'costweights':   In case of 'costbenefit' and there are more than one costs
%                    associated with the options, this is an a-by-d matrix with the
%                    importance of each cost. This can be expressed
%                    either as percentages (they should some up to 1) or
%                    pairwise comparisons. a is the costs comparison and d
%                    is the number of decisions makers. 'fuzzy' applies to
%                    the pairwise comparison of 'costweights' too.
%
%   'optimization':  In case of a resource allocation problem, this is
%                    an one-by-one matrix with the total resource availability. if it is
%                    not empty, the function solves an optimization problem (mixed-integer
%                    linear programming), where the objective function is the
%                    sum of the estimated AHP weights (reward), which is optimized subject to the cost of each option
%                    and the total resource availability.
%                    The output is the optimal combination of alternatives the decision
%                    makers should choose to maximize the total reward
%                    subject to the cost constraints. Default value is empty [].
%
%   'prediction':    In case this is a weighted average problem (e.g. forecasting combination).
%                    An t-by-o vector of the individual elements that are needed to be
%                    combined. o is the number of alternative predictions making unites,
%                    and t is the number of predictions (either periods or cross-sectional predictions).
%
%   'order':         The function puts weights table in hierarchical order according to one column.
%                    the alternatives are (i) 'weights': AHP/ANP results, (ii) 'use': optimization result,
%                    (iii) 'costbenefit': the output of cost-benefit analysis, and (iv) 'alternatives':
%                    the alternative options in alphabetical order. The default is empty.
%
%   'trials':        Number of trials for the sensitivity simulation. The default is 1000.
%
%   'sensitivity':   The sensitivity barrier of the simulation analysis.
%                    It should be greater than 1 The default is 1.
%
%   'senscriteria':  Instead of an overall level of sensitivity, the
%                    decision makers assign individual levels of sensitivities to the comparison
%                    each of the hierarchy criteria separately. This applies
%                    to qualitative criteria only.
%
%   'sensalts':      Similar to 'senscriteria', the decision makers can
%                    assign individual levels of sensitivity to the
%                    comparison of each alternative according to each of
%                    the selection criteria.
%
%   'sensfeedback':  Similar to 'sencriteria' and 'sensalts', the decision
%                    makers can assign individual levels of sensitivity to
%                    the comparison of the prevalence of each criterion to
%                    each alternative option.
%
%   'charts':        A true/false Boolean variable. If true, the model generates
%                    of sensitivity analysis charts. The default is false.
%
%   'optimadjust':   A true/false Boolean variable. If true and the cost
%                    is used, the model adjusts the weights subject to the
%                    optimization output. i.e. the weights of unused options
%                    is turned to 0. The default is false.
%
%   'randomcolor':   A true/false Boolean variable. If true, the color of the
%                    area in the charts is randomly selected. Else, the model assigns
%                    predetermined colors in each order. The default is false.
%
%   'tables':        A true/false Boolean variable. If true, 'tcriteria', 'altweights',
%                    and 'feedbackweights' outputs are presented in a table. Else
%                    they are presented in matrices. The default is true.
%
%
% Outputs:
%
%   weights:         A table, matrix or structured array, which is the combination and hierarchy weights
%                    estimated by the model. The output may also includes the results of the
%                    optimization and cost-benefit options if they have been selected. If
%                    option 'simulation' is selected, this output will
%                    include both the average weights and the weights of
%                    the selected percentiles.
%
%
% Optional outputs:
%
%   criteria:        A table or matrix of the estimated importance
%                    weights of the decision criteria.
%
%   altweights:      A table or matrix of the estimated importance
%                    weights of the alternatives per decision criterion
%
%   feedbackweights: A table or matrix of the estimated dominance
%                    weights of each of the alternatives per alternative option
%
%   cp:              In case we consider it a combination problem, the combination of
%                    the individual numerical inputs. If there are not
%                    individual projection inputs, cp will be empty.
%
%
   %% Input

    if (nargin == 1)
        C = [];
    end
    
    parseObj = inputParser;
    parseObj.addRequired('M');
    parseObj.addRequired('C');
    parseObj.addParameter('beneficial',[]);
    parseObj.addParameter('nonbeneficial',[]);
    parseObj.addParameter('feedback',[]);
    parseObj.addParameter('fuzzy',false);
    parseObj.addParameter('names',{});
    parseObj.addParameter('criterianames',{});
    parseObj.addParameter('normalization','enhanced');
    parseObj.addParameter('cost',[]);
    parseObj.addParameter('optimization',[]);
    parseObj.addParameter('costbenefit',false);
    parseObj.addParameter('costweights',[]);
    parseObj.addParameter('order',[]);
    parseObj.addParameter('simulation',false);
    parseObj.addParameter('prediction',[]);
    parseObj.addParameter('trials',1000);
    parseObj.addParameter('percentiles',[]);
    parseObj.addParameter('sensitivity',1);
    parseObj.addParameter('senscriteria',[]);
    parseObj.addParameter('sensalts',[]);
    parseObj.addParameter('sensfeedback',[]);
    parseObj.addParameter('charts',false);
    parseObj.addParameter('optimadjust',false);
    parseObj.addParameter('randomcolor',false);
    parseObj.addParameter('tables',true);
    
    parseObj.parse(M,C,varargin{:});
    
    M = parseObj.Results.M;
    C = parseObj.Results.C;
    B = parseObj.Results.beneficial;
    NB = parseObj.Results.nonbeneficial;
    feedback = parseObj.Results.feedback;
    fuzzy = parseObj.Results.fuzzy;
    names = parseObj.Results.names;
    criterianames = parseObj.Results.criterianames;
    normalization = parseObj.Results.normalization;
    cost = parseObj.Results.cost;
    optimization = parseObj.Results.optimization;
    costbenefit = parseObj.Results.costbenefit;
    costweights = parseObj.Results.costweights;
    order = parseObj.Results.order;
    simulation = parseObj.Results.simulation;
    prediction = parseObj.Results.prediction;
    trials = parseObj.Results.trials;
    percentiles = parseObj.Results.percentiles;
    sensitivity = parseObj.Results.sensitivity;
    senscriteria = parseObj.Results.senscriteria;
    sensalts = parseObj.Results.sensalts;
    sensfeedback = parseObj.Results.sensfeedback;
    charts = parseObj.Results.charts;
    optimadjust = parseObj.Results.optimadjust;
    randomcolor = parseObj.Results.randomcolor;
    tables = parseObj.Results.tables;
    
   %% Calculation

        M(M==0)=NaN;
        M(M<0)=-1./(M(M<0));
        
        if isempty(B) == false
            if isequal(normalization,'linear')
            elseif isequal(normalization,'relativelinear')
                B = B./max(B,[],1);
            elseif isequal(normalization,'minmax')
                B = (B-min(B,[],1))./(max(B,[],1)-min(B,[],1));
            elseif isequal(normalization,'vector')
                B = B./sqrt(sum(B.^2,1));
            elseif isequal(normalization,'enhanced')
                B = 1 - (max(B,[],1)-B)./sum((max(B,[],1)-B),1);
            elseif isequal(normalization,'logarithmic')
                B = log(B)./log(prod(B,1));
            end
            B = B./sum(B,1);
        end % isempty(B) == false
        
        if isempty(NB) == false
            if isequal(normalization,'linear')
                NB = 1./NB;
            elseif isequal(normalization,'relativelinear')
                NB = min(NB,[],1)./NB;
            elseif isequal(normalization,'minmax')
                NB = (max(NB,[],1)-NB)./(max(NB,[],1)-min(NB,[],1));
            elseif isequal(normalization,'vector')
                NB = 1 - NB./sqrt(sum(NB.^2,1));
            elseif isequal(normalization,'enhanced')
                NB = 1 - (NB-min(NB,[],1))./sum((NB-min(NB,[],1)),1);
            elseif isequal(normalization,'logarithmic')
                NB = 1-(1-log(NB)./log(prod(NB,1)))./(size(NB,1)-1);
            end
            NB = NB./sum(NB,1);
        end % isempty(NB) == false

        if isempty(C)
            Cmean = [];
            Mmean = geomean(M,2);
        else
            C(C==0)=NaN;
            C(C<0)=-1./(C(C<0));
            Cmean = geomean(C,2);
            Mmean = [];
            for j = 1:size(C,2):size(M,2)
                Mmean = [Mmean,geomean(M(:,j:j+size(C,2)-1),2)];
            end % 1:size(C,2):size(M,2)
        end % (nargin == 1)
        
        
        if isempty(feedback) == false

            feedback(feedback==0)=NaN;
            feedback(feedback<0)=-1./(feedback(feedback<0));
            Fmean = [];

            for j = 1:size(C,2):size(feedback,2)
                    Fmean = [Fmean,geomean(feedback(:,j:j+size(C,2)-1),2)];
            end % 1:size(C,2):size(M,2)
            
            if fuzzy % fuzzification
                Fmean = fuzz(Fmean);
            end % fuzzy
            
        end % isempty(feedback) == false
        
        
%% Names

        if isempty(names)
            if isempty(Mmean)
                J = 1:size([B,NB]);
            else
                J = max(roots([.5 -.5 -size(Mmean,1)]));
            end % isempty(Mmean)
            for j = 1:J
                names(j,1) = {['Option ',num2str(j)]};
            end % j = 1:size(w,1)
        end % (nargin <= 3) || isempty(names)

        if isempty(criterianames)
            if isempty(Cmean)
                criterianames(1,1) = 'Criterion';
            else
                for j = 1:max(roots([.5 -.5 -size(Cmean,1)]))
                    criterianames(j,1) = {['Criterion ',num2str(j)]};
                end % j = 1:size(w,1)
            end % isempty(Cmean)
        end % (nargin <= 3) || isempty(names)
        
        
 %% Regular AHP and ANP
 if simulation == false
        
        if fuzzy % fuzzification
           Cmean = fuzz(Cmean);
           Mmean = fuzz(Mmean);  
        end % fuzzy
        
        for l = 1:max(size(Cmean,2),size(Mmean,3))
        
            Cweights(:,l) = vectors(Cmean(:,l));
        
            if isempty(M) == false
                for j = 1:min(size(Mmean,2)) % a loop for every criterion
                    Mweights(:,j,l) = vectors(Mmean(:,j,min(l,size(Mmean,3)))); % Estimation of the decision weights vector per criterion using the "vectors" fuction below.
                end %j = 1:min(size(Cweights,1))
                if l == 1
                    mwsize = size(Mweights,2)+1:size(Mweights,2)+size([B NB],2);
                end
                    Mweights(:,mwsize,l) = [B NB];
            else
                Mweights(:,:,l) = [B NB];
            end % isempty(M) == false

        end % l = max(size(Cmean,3),size(Mmean,3))
        
        if isempty(feedback) == false

                for l = 1:size(Fmean,3)
                    
                    for j = 1:min(size(Mweights,1))
                        Fweights(:,j,l) = vectors(Fmean(:,j,l));
                    end % j = 1:min(size(Mweights,1))

                    weights(:,l) = networkweights(Cweights(:,l),Mweights(:,:,l),Fweights(:,:,l)); % calibration of weights through the analytic network process
                    
                end % l = size(Fmean,3)
                 FweightsP = mean(Fweights,3)./sum(mean(Fweights,3));   
        else
             for l = 1:size(Mweights,3)
                weights(:,l) = Mweights(:,:,l)*Cweights(:,l); % calibration of the hierarchy weights for all alternative decisions
             end % l = size(Mweights,3)
        end % isempty(feedback) == false
        
        weights = mean(weights,2)./sum(mean(weights,2));
        MweightsP = mean(Mweights,3)./sum(mean(Mweights,3),1);
        CweightsP = mean(Cweights,2)./sum(mean(Cweights,2));
        
    
    %% Simulation analysis

    else

        for j = 1:trials

            if (nargin == 1) || isempty(C)
            else
                if isempty(senscriteria)
                    senscriteria = ones(size(Cmean))*sensitivity;
                end % isempty(senscriteria)
                Csen = sen(Cmean,senscriteria);
                Cweights(:,j) = vectors(Csen);
            end % Cmean == []
            
            if isempty(M) == false
            
                if isempty(sensalts)
                    sensalts = ones(size(Mmean))*sensitivity;
                end % isempty(sensalts)
                Msen = sen(Mmean,sensalts);

                for l = 1:min(size(Mmean,2)) % a loop for every criterion
                    MweightsD(:,l) = vectors(Msen(:,l)); % Estimation of the decision weights vector per criterion using the "vectors" fuction below.
                end %l = 1:min(size(Cweights,1))
                
            else
                
                MweightsD = [];
                
            end % isempty(M) == false
                
            Mweights(:,:,j) = [MweightsD, B, NB];
            
            if isempty(feedback) == false
                
                if isempty(sensfeedback)
                    sensfeedback = ones(size(Fmean))*sensitivity;
                end % isempty(sensfeedback)
                Fsen = sen(Fmean,sensfeedback);

                for l = 1:min(size(Mweights,1))
                    Fweights(:,l,j) = vectors(Fsen(:,l));
                end % l = 1:min(size(Mweights,1))

                wsen(:,j) = networkweights(Cweights(:,j),Mweights(:,:,j),Fweights(:,:,j)); % calibration of weights through the analytic network process

            else
                
                wsen(:,j) = Mweights(:,:,j)*Cweights(:,j); % calibration of the hierarchy weights for all alternative decisions

            end % isempty(feedback) == false

        end % j = 1:trials

        weights = mean(wsen,2)/sum(mean(wsen,2));
        MweightsP = mean(Mweights,3)./sum(mean(Mweights,3),1);
        CweightsP = mean(Cweights,2)./sum(mean(Cweights,2));
        if isempty(percentiles) == false
            for j = 1:size(percentiles,1)
                weights(:,1+j) = prctile(wsen,percentiles(j),2)/sum(prctile(wsen,percentiles(j),2));
                for l = 1:size(Mweights,2)
                    MweightsP(:,l,1+j) = prctile(Mweights(:,l,:),percentiles(j),3)/sum(prctile(Mweights(:,l,:),percentiles(j),3));
                end % l = 1:size(Mweights,2)
                CweightsP(:,1+j) = prctile(Cweights,percentiles(j),2)/sum(prctile(Cweights,percentiles(j),2));
            end % j = 1:size(percentiles,1)
        end % isempty(percentiles) == false
                
        if isempty(feedback) == false            
            FweightsP = mean(Fweights,3)./sum(mean(Fweights,3),1);
            if isempty(percentiles) == false
                for j = 1:size(percentiles,1)
                    for l = 1:size(Fweights,2)
                        FweightsP(:,l,1+j) = prctile(Fweights(:,l,:),percentiles(j),3)/sum(prctile(Fweights(:,l,:),percentiles(j),3));
                    end % l = 1:size(Fweights,2)
                end % j = 1:size(percentiles,1)
            end % isempty(percentiles) == false
        end % isempty(feedback) == false
        

%% Charts
        
        if charts

            for j = 1:size(wsen,1)

               if randomcolor == true
                    xt = rand(1);
                    xt(1,2) = rand(1);
                    xt(1,3) = rand(1);
               else
                    if j == 1
                        xt = [0.3010 0.7450 0.9330];
                    elseif j == 2
                        xt = [0.75, 0.75, 0];
                    elseif j == 3
                        xt = [0.9290 0.6940 0.1250];
                    elseif j == 4
                        xt = [0.4660 0.6740 0.1880];
                    elseif j == 5
                        xt = [0 0.4470 0.7410];
                    elseif j == 6
                        xt = [0, 0.75, 0.75];
                    elseif j == 7
                        xt = [0.8500 0.3250 0.0980];
                    elseif j == 8
                        xt = [0.75, 0, 0.75];
                    else
                        xt = rand(1);
                        xt(1,2) = rand(1);
                        xt(1,3) = rand(1);
                    end
               end

                figure('Name',char(names(j)),'NumberTitle','off')
                [f,xi] = ksdensity(100*wsen(j,:)');
                area(xi,f,'FaceColor',xt)
                y = get(gca,'ylim');
                hold on
                plot(100*[weights(j,1) weights(j,1)],y,'black','LineWidth',1)
                txtels = ['\leftarrow','Mean: ',num2str(round(mean(wsen(j,:)')*100,2)),'%'];
                text(100*mean(wsen(j,:)'),max(y)*4/5,txtels,'Fontsize',15)            
                hold on
                plot(100*[median(wsen(j,:)') median(wsen(j,:)')],y,'--black','LineWidth',1)
                txtels = ['\leftarrow','Q2: ',num2str(round(median(wsen(j,:)')*100,2)),'% (median)'];
                text(100*median(wsen(j,:)'),max(y)/2,txtels,'Fontsize',15)
                hold on
                plot(100*[quantile(wsen(j,:)',0) quantile(wsen(j,:)',0)],y,'--black','LineWidth',1)
                txtels = ['\leftarrow','Q0: ',num2str(round(quantile(wsen(j,:)',0)*100,2)),'%'];
                text(100*quantile(wsen(j,:)',0),max(y)*0.5,txtels,'Fontsize',15)
                hold on
                plot(100*[quantile(wsen(j,:)',0.25) quantile(wsen(j,:)',0.25)],y,'--black','LineWidth',1)
                txtels = ['\leftarrow','Q1: ',num2str(round(quantile(wsen(j,:)',0.25)*100,2)),'%'];
                text(100*quantile(wsen(j,:)',0.25),max(y)*0.25,txtels,'Fontsize',15)
                hold on
                plot(100*[quantile(wsen(j,:)',0.75) quantile(wsen(j,:)',0.75)],y,'--black','LineWidth',1)
                txtels = ['\leftarrow','Q3: ',num2str(round(quantile(wsen(j,:)',0.75)*100,2)),'%'];
                text(100*quantile(wsen(j,:)',0.75),max(y)*0.75,txtels,'Fontsize',15)
                hold on
                plot(100*[quantile(wsen(j,:)',1) quantile(wsen(j,:)',1)],y,'--black','LineWidth',1)
                txtels = ['\leftarrow','Q4: ',num2str(round(quantile(wsen(j,:)',1)*100,2)),'%'];
                text(100*quantile(wsen(j,:)',1),max(y)/2,txtels,'Fontsize',15)
                set(gca,'Fontsize',15,'YTick',[])
                xlabel('%','Fontsize',15)
                title(['Weight distribution for ', char(names(j))],'Fontsize',15)
            end % j = 1:size(wsen,1)

        end % charts

    end % simulation

    %% Cost-Benefit analysis
    
    if costbenefit
        
        if isempty(cost)
            cost = ones(size(weights,1),1);
        end % isempty(cost)
        
        costper = cost./sum(cost,1);
        
        if isempty(costweights)
           costperMean = mean(costper,2);
        else
            Mcostweights = mean(costweights,2);
            if isequal(size(Mcostweights,1),size(cost,2))
                costweightsCalc = Mcostweights./sum(Mcostweights);
            else
                 if fuzzy
                     Mcostweights = fuzz(Mcostweights);
                 end % fuzzy
                
                 for l = 1:size(Mcostweights,3)
                     costweightsCalc(:,l) = vectors(Mcostweights(:,:,l));
                 end
                 
                 costweightsCalc(:,l) = mean(costweightsCalc,2)./sum(mean(costweightsCalc,2));
                
            end
            
            costperMean = costper*costweightsCalc;
            
        end
            
        cost_benefit = weights./costperMean(1:end,1);
            
    end

    %% Optimisation

    if isempty(optimization) == false
            
         for j = 1:size(weights,2)
             
             use(:,j) = optim(weights(:,j),[cost; optimization]);
             
             if optimadjust
                weights(:,j) = weights(:,j).*use(:,j)/sum(weights(:,j).*use(:,j),1);
             end % optimadjust
             
         end % j = 1:size(weights,2)
            
     end % isempty(optimization) == false            


    %% Prediction

     if isempty(prediction) == false

         for j = 1:size(weights,2)
            combination(:,j) = prediction*weights(:,j); % calibration of the combined projection if individual projections are available
         end
         
         combination = [];
         
     end % isempty(prediction) == false

%% Tables

     if tables
            
         if simulation
             
             Ws = weights;
             
             if isempty(optimization) == false && costbenefit
                 weights = table(weights(:,1),use(:,1),cost_benefit(:,1),...
                 'RowNames',names,'VariableNames',{'mean weights','mean use','mean cost-benefit'});
                 if isempty(percentiles) == false
                     for j = 2:size(Ws,2)
                         Var = table(Ws(:,j),use(:,j),cost_benefit(:,j),...
                         'RowNames',names,'VariableNames',{join(['weights at ',num2str(percentiles(j-1)),'%']),...
                         join(['use at ',num2str(percentiles(j-1)),'%']),join(['cost-benefit at ',num2str(percentiles(j-1)),'%'])});
                         weights = [weights Var];
                     end % j = 2:size(Ws,2)
                 end % isempty(percentiles) == false
             elseif isempty(optimization) == false
                 weights = table(weights(:,1),use(:,1),'RowNames',names,'VariableNames',{'mean weights','mean use'});
                 if isempty(percentiles) == false
                     for j = 2:size(Ws,2)
                         Var = table(Ws(:,j),use(:,j),'RowNames',names,'VariableNames',...
                         {join(['weights at ',num2str(percentiles(j-1)),'%']),join(['use at ',num2str(percentiles(j-1)),'%'])});
                         weights = [weights Var];
                     end % j = 2:size(Ws,2)
                 end % isempty(percentiles) == false
             elseif costbenefit
                 weights = table(weights(:,1),cost_benefit(:,1),'RowNames',names,'VariableNames',...
                 {'mean weights','mean cost-benefit'});
                 if isempty(percentiles) == false
                     for j = 2:size(Ws,2)
                         Var = table(Ws(:,j),cost_benefit(:,j),...
                         'RowNames',names,'VariableNames',{join(['weights at ',num2str(percentiles(j-1)),'%']),...
                         join(['cost-benefit at ',num2str(percentiles(j-1)),'%'])});
                         weights = [weights Var];
                     end % j = 2:size(Ws,2)
                 end % isempty(percentiles) == false
             else
                 weights = table(weights(:,1),'RowNames',names,'VariableNames',{'mean weights'});
                 if isempty(percentiles) == false
                     for j = 2:size(Ws,2)
                         Var = table(Ws(:,j),'RowNames',names,'VariableNames',{join(['weights at ',num2str(percentiles(j-1)),'%'])});
                         weights = [weights Var];
                     end % j = 2:size(Ws,2)
                 end % isempty(percentiles) == false
             end % isempty(optimization) == false && costbenefit

             criteria = table(CweightsP(:,1),'RowNames',criterianames,'VariableNames',{'mean criteria weights'});
             if isempty(percentiles) == false
                for j = 2:size(Ws,2)
                    Var = table(CweightsP(:,j),'RowNames',criterianames,'VariableNames',...
                    {join(['criteria weights at ',num2str(percentiles(j-1)),'%'])});
                    criteria = [criteria Var];
                end % j = 2:size(Ws,2)
             end % isempty(percentiles) == false
             
             altweights = table(MweightsP(:,1,1),'RowNames',names,'VariableNames',join([criterianames(1,1),'mean']));
             for j = 2:size(MweightsP,2)
                 Var = table(MweightsP(:,j,1),'VariableNames',join([criterianames(j,1),'mean']));
                 altweights = [altweights Var];
             end % j = 2:size(Mweights,2)
             if isempty(percentiles) == false
                 for j = 2:size(MweightsP,3)
                     for l = 1:size(MweightsP,2)
                        Var = table(MweightsP(:,l,j),'VariableNames',join([criterianames(l,1),'at',num2str(percentiles(j-1)),'%']));
                        altweights = [altweights Var];
                     end % l = 1:size(MweightsP,2)
                 end % j = 2:size(MweightsP,3)
             end % isempty(percentiles) == false
             
             if isempty(feedback) == false
                 feedbackweights = table(FweightsP(:,1,1),'RowNames',...
                     criterianames,'VariableNames',join([names(1,1),'mean']));
                 for j = 2:size(FweightsP,2)
                     Var = table(FweightsP(:,j,1),'VariableNames',join([names(j,1),'mean']));
                     feedbackweights = [feedbackweights Var];
                 end
                 if isempty(percentiles) == false
                    for j = 2:size(FweightsP,3)
                        for l = 1:size(FweightsP,2)
                            Var = table(FweightsP(:,l,j),'VariableNames',join([names(l,1),'at',num2str(percentiles(j-1)),'%']));
                            feedbackweights = [feedbackweights Var];
                        end % l = 1:size(FweightsP,2)
                    end % j = 2:size(FweightsP,3)
                end % isempty(percentiles) == false
             else
                 feedbackweights = [];
             end
                
         else
                
             if isempty(optimization) == false && costbenefit
                 weights = table(weights,use,cost_benefit,'RowNames',names,'VariableNames',{'weights','use','cost-benefit'});
             elseif isempty(optimization) == false
                 weights = table(weights,use,'RowNames',names,'VariableNames',{'weights','use'});
             elseif costbenefit
                 weights = table(weights,cost_benefit,'RowNames',names,'VariableNames',{'weights','cost-benefit'});
             else
                 weights = table(weights,'RowNames',names);
             end
                
             criteria = table(CweightsP,'RowNames',criterianames,'VariableNames',{'criteria weights'});

             altweights = table(MweightsP(:,1),'RowNames',names,'VariableNames',criterianames(1,1));

             for j = 2:size(CweightsP,1)
                 Var = table(MweightsP(:,j),'VariableNames',criterianames(j,1));
                 altweights = [altweights Var];
             end % j = 2:size(Mweights,2)

             if isempty(feedback) == false
                 feedbackweights = table(FweightsP(:,1),'RowNames',...
                     criterianames,'VariableNames',names(1,1));
                 for j = 2:size(FweightsP,2)
                     Var = table(FweightsP(:,j),'VariableNames',names(j,1));
                     feedbackweights = [feedbackweights Var];
                 end
             else
                 feedbackweights = [];
             end % isempty(feedback) == false
             
         end %simulation
         
        % the results can be in order
        if isempty(order) == false
            if isequal(order,'weights')
                weights = sortrows(weights,1,'descend');
            elseif isequal(order,'use')
                if simulation && isempty(optimization) == false
                    weights = sortrows(weights,'mean use','descend');
                elseif isempty(optimization) == false
                    weights = sortrows(weights,'use','descend');
                else
                    weights = sortrows(weights,1,'descend');
                end % simulation && isempty(optimization) == false
            elseif isequal(order,'costbenefit')
                if simulation && costbenefit
                    weights = sortrows(weights,'mean cost-benefit','descend');
                elseif costbenefit
                    weights = sortrows(weights,'cost-benefit','descend');
                else
                    weights = sortrows(weights,1,'descend');
                end % simulation && costbenefit
            elseif isequal(order,'alternatives')
                weights = sortrows(weights,'RowNames');
            end % isequal(order,'weights')
        end % isempty(order) == false
      
     else
            
          if isempty(optimization) == false && costbenefit
               weights2 = weights;
               clear weights
               weights.weights = weights2;
               weights.use = use;
               weights.cost_benefit = cost_benefit;
          elseif isempty(optimization) == false
               weights2 = weights;
               clear weights
               weights.weights = weights2;
               weights.use = use;
          elseif costbenefit
               weights2 = weights;
               clear weights
               weights.weights = weights2;
               weights.cost_benefit = cost_benefit;
          end
          criteria = CweightsP;
          altweights = MweightsP;
          if isempty(feedback) == false
              feedbackweights = FweightsP;
          else
              feedbackweights = [];
          end
            
     end % tables

        
    %% Functions

    function v = vectors(C) % The function that estimates the criteria and decision weights

        nc = max(roots([.5 -.5 -size(C,1)])); % calibration of the number of criteria or alternatives using the numbers of pairwise comparisons
        CM = eye(nc); % Creation of the comparison matrix from the comparison vectors
        
        for i = 1:nc-1 % Fill the zeroes in the area above the diagonal with the comparison scores and the area below the diagonal with their reverse value
            
            if i == 1
                sp = i;
            else
                sp=fp+1;
            end
            fp = sp+nc-1-i;
            c = C(sp:fp,1)';

            CM(i,i+1:end) = c;
            clear c
            CM(i+1:end,i) = 1./CM(i,i+1:end)';

        end % i = 1:nc-1
        
        v = geomean(CM^100,2); % estimation of the weights vector
        v = v/sum(v); % normalisation of the weights, dividing by the sum of weights

    end % v = vectors(C)

    function ff = fuzz(f) % fuzzification
        
           fL = f;
           fL(f>1) = f(f>1)-1;
           fL(f<1) = 1./((1./f(f<1))+1);
           fL(f==1) = 1/2;
           fH = f;
           fH(f>1) = f(f>1)+1;
           fH(f<1) = 1./((1./f(f<1))-1);
           fH(f==1) = 2;
           
           ff = cat(3,fL,f,fH);
           
    end % ff = fuzz(f)

    function w = networkweights(Cweights,Mweights,Fweights)

        wsm = zeros(size(Cweights,1)+size(Mweights,1)+1);
        wsm(2:size(Cweights,1)+1,1) = Cweights;
        wsm(size(Cweights,1)+2:end,2:size(Mweights,2)+1) = Mweights;
        wsm(2:size(Cweights,1)+1,size(Cweights,1)+2:end) = Fweights;
        wsm = wsm./sum(wsm);
        wsm = wsm^100;
        
        w = wsm(size(Cweights,1)+2:end,1);

    end % w = networkweights(Cweights,Mweights,Fweights)

    function Msen = sen(Mean,d)

        for i = 1:size(Mean,1)
            for k = 1:size(Mean,2)
                if Mean(i,k) > 1
                    Msen(i,k) = unifrnd(max(Mean(i,k)-d(i,k),d(i,k)+2-Mean(i,k)),Mean(i,k)+d(i,k));
                elseif Mean(i,k) < 1
                    Msen(i,k) = -1/unifrnd((-1/Mean(i,k))-d(i,k),min((-1/Mean(i,k)),1/(d(i,k)-(-1/Mean(i,k))+2)));
                elseif Mean(i,k) == 1
                    Msen(i,k) = unifrnd(1/(d(i,k)+1),(d(i,k)+1));
                end
            end % k = 1:size(Mean,2)
        end % i = 1:size(Mean,1)

    end % Msen = sen(Mean,d)

    function use = optim(w,c)

            w = -w;
            intcon = 1:size(w,1);
            A1 = eye(size(w,1));
            b1 = ones(size(w,1),1);
            A = [A1;c(1:end-1,:)'];
            b = [b1;c(end,:)'];
            lb = zeros(size(w,1),1);
            use = intlinprog(w,intcon,A,b,[],[],lb);

    end % use = optim(w,c)

end % End of function