% Main routine for running the MA-ES variants on the 2017 CEC Benchmarks
clear all, clc

global  initial_flag
initial_flag = 0;


% choose the problem dimensionality D in [10, 30, 50, 100]
D=10;

% bound constraint definitions for all 28 test functions
Xmin1=-100*ones(1,D);
Xmax1=+100*ones(1,D);
Xmin2=-100*ones(1,D);
Xmax2=+100*ones(1,D);
Xmin3=-100*ones(1,D);
Xmax3=+100*ones(1,D);
Xmin4=-10*ones(1,D);
Xmax4=+10*ones(1,D);
Xmin5=-10*ones(1,D);
Xmax5=+10*ones(1,D);
Xmin6=-20*ones(1,D);
Xmax6=+20*ones(1,D);
Xmin7=-50*ones(1,D);
Xmax7=+50*ones(1,D);
Xmin8=-100*ones(1,D);
Xmax8=+100*ones(1,D);
Xmin9=-10*ones(1,D);
Xmax9=+10*ones(1,D);
Xmin10=-100*ones(1,D);
Xmax10=+100*ones(1,D);
Xmin11=-100*ones(1,D);
Xmax11=+100*ones(1,D);
Xmin12=-100*ones(1,D);
Xmax12=+100*ones(1,D);
Xmin13=-100*ones(1,D);
Xmax13=+100*ones(1,D);
Xmin14=-100*ones(1,D);
Xmax14=+100*ones(1,D);
Xmin15=-100*ones(1,D);
Xmax15=+100*ones(1,D);
Xmin16=-100*ones(1,D);
Xmax16=+100*ones(1,D);
Xmin17=-100*ones(1,D);
Xmax17=+100*ones(1,D);
Xmin18=-100*ones(1,D);
Xmax18=+100*ones(1,D);
Xmin19=-50*ones(1,D);
Xmax19=+50*ones(1,D);
Xmin20=-100*ones(1,D);
Xmax20=+100*ones(1,D);
Xmin21=-100*ones(1,D);
Xmax21=+100*ones(1,D);
Xmin22=-100*ones(1,D);
Xmax22=+100*ones(1,D);
Xmin23=-100*ones(1,D);
Xmax23=+100*ones(1,D);
Xmin24=-100*ones(1,D);
Xmax24=+100*ones(1,D);
Xmin25=-100*ones(1,D);
Xmax25=+100*ones(1,D);
Xmin26=-100*ones(1,D);
Xmax26=+100*ones(1,D);
Xmin27=-100*ones(1,D);
Xmax27=+100*ones(1,D);
Xmin28=-50*ones(1,D);
Xmax28=+50*ones(1,D);

% number of constraint functions per problem
 
problem.gn=[1 1 1 2 2 0 0 0 1 0 1 2 3 1 1 1 1 2 2 2 2 3 1 1 1 1 2 2];
problem.hn=[0 0 1 0 0 6 2 2 1 2 1 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0];

% budget of function evaluations and generations depending on dimension D
MaxFES=D*20000;
MaxIter=D*2000;

%% Input values
% initial (fixed) strategy parameter setting
%
input.dim               = D;
input.budget            = MaxFES;
input.maxIter           = MaxIter;

input.delta             = 10^-4;                   % error margin for equality constraints
input.runs              = 25;                       % number of independent algorithm runs (standard -- 25 runs)

input.lambda            = 4*D;                       % population size
input.mu                = floor(input.lambda/3);    % parental ppopulation size

input.sigma             = 1;                        % initial mutation strength
      
% MA-ES specific strategy paprameters (standard parameter choices)
%
input.weights = log(input.mu+1/2)-log(1:input.mu)';     % muXone array for weighted recombination
input.weights = input.weights./sum(input.weights);      % normalize recombination weights array
input.mueff=1/sum(input.weights.^2);                    % variance-effectiveness of sum w_i x_i

input.cs = (input.mueff+2) / (D+input.mueff+5);  % t-const for cumulation for sigma control
input.c1 = 2 / ((D+1.3)^2+input.mueff);          % learning rate for rank-one update of M
input.cmu = min(1-input.c1, 2 * (input.mueff-2+1/input.mueff) / ((D+2)^2+input.mueff));  % and for rank-mu update of M
input.damps = 1 + 2*max(0, sqrt((input.mueff-1)/(D+1))-1) + input.cs; % damping for sigma usually close to 1

input.reps = 3; % number of repair repetitions within epsMAg-ES

problem.constr_fun_name = 'CEC2017';  % objective function class

% Selection of MA-ES variants for constrained optimization
ESvariant{1}='epsMAgES';	% original epsilonMAg-ES variant under investigation

    for jj=1:1    % for each ESvariant do

        strategy=ESvariant{jj};
        foldername = ['CEC17_RS_' strategy];
        efn = exist(foldername);
        if efn ~= 7
            mkdir(foldername);
        end

        disp(['ES variant ' strategy ' --- ' num2str(input.runs) ' independent runs!'])
        
        % selection of testfunctions (C01 to C28)
        funs=[1:28];  %choose subset if desired
        
        % on each of the selected CEC2017 test function do
        for kk=1:length(funs)  
            func_num=funs(kk);          % test function number
            
            initial_flag = 0;
            
            eval(['problem.lower_bounds=transpose(Xmin' int2str(func_num) ');']);
            eval(['problem.upper_bounds=transpose(Xmax' int2str(func_num) ');' ]);

            input.dim               = D;
            input.population_size   = 4*D;          
            
            for j=1:input.runs  % perform multiple runs on each test function

                eval(['[tab, global_best, dyn]=' strategy '(problem,input,func_num);']); % run strategy
                disp([func_num,j])
                             
                FitT(j,:)       = [tab(1,:) tab(2,:) tab(3,:)] % save fitness / constraint violation after 10%, 50% and 100% of the evaluation budget as recommended by the CEC2017 competition
                Dyn{j}          = dyn;               % store returned dynamics of run j
                GB{j}           = global_best;       % store globally best observation of run j
            end
            
            tab         = build_stats(FitT,input);    % convert matrix of fitnesses/constraintviolations FitT to the performance indicators recommended by the CEC2017 competition
            Stats(kk,:) = [tab(1,:) tab(2,:) tab(3,:)];

            save([foldername '/' strategy '_RS_on_fun' num2str(func_num) 'D' num2str(D) '2.mat'],'GB','Dyn','FitT','problem','input','-v7')

        end
    end
