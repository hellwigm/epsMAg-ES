function [out,global_best, dyn]=epsMAgES(problem,input,CEC_fun_no)
    %% Implementation of the epsilonMAg-ES for constrained optimiazation
    % Implementation according to
    % M. Hellwig and H.-G. Beyer, "A Matrix Adaptation Evolution Strategy for Constrained
    % Real-Parameter Optimization", 2018 IEEE Congress on Evolutionary
    % Computation (CEC), IEEE, 2018, https://dx.doi.org/10.1109/CEC.2018.8477950
    %
    
    dim         = length(problem.lower_bounds);
    sigma       = input.sigma;
    mu          = input.mu;
    lambda      = input.lambda;
    newpop.y    = zeros(dim,lambda);					% initialize new population matrix (n times NP)
	newpop.f    = 0;
	newpop.conv = 0;
	evals.fun   = 0;
	evals.con   = 0;
    
    g           = 0;
    termination = 0;
    
    % flags for saving the realizations after having consumed 10% and 50% of the
    % constrained function evaluation budget
    ll=0;
    lll=0;
    
    % Initialize dynamic (internal) strategy parameters and constants
    ps      = zeros(dim,1);                                 % evolution paths for sigma
    M       = eye(dim);                                     % transformation matrix
    sqrt_s  = sqrt(input.cs*(2-input.cs)*input.mueff);      % factor in path update
    
    % upper limit of admissilbe mutation strength values
    sigmax = max(problem.upper_bounds-problem.lower_bounds)/2; % 100;      NEW: set sigmax to half the maximum boundary width, formerly this parameter was simply set to 100.
    
    % Initialize random population of lambda individuals as done in the DE
    % variants
    for k=1:lambda
        % create initial population of uniformly distributed vectors within box-constraints 
        newpop.y(:,k)   = problem.lower_bounds...
                            +(problem.upper_bounds-problem.lower_bounds).*rand(dim,1); 
        % evaluate initial population
        [fval, gv, hv]  = feval(problem.constr_fun_name,newpop.y(:,k)',CEC_fun_no);
        newpop.f(k)     = fval;	  
        % calculate constraint violation according to CEC2017
        % recommendations
        newpop.conv(k)  = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./(problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no));
        % count constraint function evaluations
		evals.fun       = evals.fun + 1;
    end
       
    % initial parameters of the epsilon constraint handling approach
    TC=1000;
    n=ceil(0.9*size(newpop.conv,2));   
    index=eps_sort(newpop.f,newpop.conv,0); % initial sorting w.r.t. epsilon-level zero <-- lexicographical ordering 
    EPSILON=mean(newpop.conv(index(1:n)));
    Epsilon= EPSILON;    
    CP=max(3,(-5-log(EPSILON))/log(0.05));
    
    [ranking]           = eps_sort(newpop.f,newpop.conv,Epsilon);                          % epsilon Ranking
    ParentPop           = newpop.y(:,ranking(1:mu));
    yParent             = sum(ParentPop,2)./mu;
    
    best_ind            = ranking(1);
    best_val            = newpop.f(best_ind);				% best fitness value of current population
    best_y              = newpop.y(:,best_ind);
    best_conv           = newpop.conv(:,best_ind);
        
    evals.sum           = evals.fun;
    
    % best solution found so far
    global_best.y       = best_y; 				
	global_best.val     = best_val;
    global_best.conv    = best_conv;
    global_best.evals   = evals.fun;
    
    % Initialize algorithm dynamics
    dyn.gen(g+1)        = g;
    dyn.fev(g+1)        = evals.sum;
    dyn.fit(g+1)        = global_best.val;
    dyn.conv(g+1)       = global_best.conv;
    dyn.sigma(g+1)      = sigma;
    dyn.ynorm(g+1)      = norm(global_best.y);
     
    while ~termination
        % compute Moore-Penrose inverse of M for the back-calculation step
        % which is applied to correct the mutation vectors corresponding to repaired solutions
        Minv      = pinv(M,1e-12);
                
        % create new generation of offspring candidate solutions
        newpop.z    = randn(dim,lambda);
        newpop.d    = M*newpop.z;
        newpop.y    = repmat(yParent,1,lambda)+sigma.*newpop.d;
               
        for k=1:lambda          
            % initialization of repair count
            repi = 0;
            % check for bound constraint satisfaction and repair if
            % necessary
            newy =keep_range(newpop.y(:,k),problem.lower_bounds,problem.upper_bounds);
            % check whether repair has been carried out in order to apply
            % back-calculation,
            if ~isequal(newy,newpop.y(:,k))
                repi = 1;
            end
            % evaluation of offspring candiadate solution (in bounds) 
            [fval, gv, hv]  = feval(problem.constr_fun_name,newy',CEC_fun_no);
            fitval          = fval;                                                                
            % compute constraint violation
            convio          = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./(problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no));            
            evals.fun       = evals.fun + 1;
            % initialize individual repair step count
            h=1;            
            % apply gradient-based mutation step if conditions are satisfied
            if mod(g,dim)==0 && rand(1) <=0.2 
                while convio > 0 && h <= input.reps
                     new_mutant = gradientMutation(problem,newy,gv,hv,CEC_fun_no);
                     new_mutant = keep_range(new_mutant,problem.lower_bounds,problem.upper_bounds);
            
                     [fval, gv, hv]  = feval(problem.constr_fun_name,new_mutant',CEC_fun_no);
                     fitval = fval;                                                               
                     convio = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./(problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no));
                     evals.fun       = evals.fun + dim +1;
                     h=h+1;
                     newy = new_mutant;
                    if ~isequal(newy,newpop.y(:,k))
                        repi = 1;
                    end
                end 
                
            end
            newpop.y(:,k) = newy;
 
            % apply back-calculation if necessary
            if repi > 0
                newpop.d(:,k) = (newpop.y(:,k)-yParent)./sigma;
                newpop.z(:,k) = Minv*newpop.d(:,k);
            end
            newpop.f(k)     = fitval;                                                          
            newpop.conv(k)  = convio;
        end
           
          
        % Implementation of Epsilon Constraint Ordering 
        % feasible (constraint violation below epsilon value!!!) solutions dominate infeasible 
        % ones AND feasible solutions are sorted according to their fitness values
                
        [ranking]   = eps_sort(newpop.f,newpop.conv,Epsilon);
                
        best_ind    = ranking(1);
        best_val    = newpop.f(best_ind);             % best feasible fitness value of current population
        best_y      = newpop.y(:,best_ind);           % best feasible individual of current population
        best_conv   = newpop.conv(:,best_ind);
        
        % Sort by fitness and compute weighted mean into xmean
        parent_z = newpop.z(:,ranking(1:mu)) * input.weights;  % recombination
        parent_d = newpop.d(:,ranking(1:mu)) * input.weights;  % recombination
        yParent  = yParent + sigma *parent_d; % update population certroid
               
        % Cumulation: Update evolution paths
        ps = (1-input.cs) * ps + sqrt_s * parent_z; 

        % Update transformation matrix
        M = (1 - 0.5*input.c1 - 0.5*input.cmu) * M + (0.5*input.c1)*(M*ps)*ps';
        for m = 1:mu
            M = M + ((0.5*input.cmu*input.weights(m))*newpop.d(:,ranking(m)))*newpop.z(:,ranking(m))';
        end    
        
        % Reset step to prevent unstable pseudo inverse calculations
        liMM = M > 1e+12 | isnan(M) | isinf(M);
        siMM = M < -1e+12 | isnan(M) | isinf(M);
        if sum(sum(liMM))>1 || sum(sum(siMM))>1
            M = eye(input.dim);
            ps = ones(input.dim,1);
        end
                
        % Adapt mutation strength sigma
        sigma = min(sigma * exp((input.cs/2)*(norm(ps)^2/input.dim - 1)),sigmax); 
       
        % update the best solution found so far                        
        if (best_conv==0 && global_best.conv==0 && best_val < global_best.val) ||...
                (best_conv==global_best.conv && best_val < global_best.val) || best_conv<global_best.conv
            global_best.y     = best_y; 				
            global_best.val   = best_val;
            global_best.conv  = best_conv;
            global_best.evals = evals.fun; % number of function evaluatioan at which best so far was observed
        end
        evals.sum = evals.fun;%+evals.con;
                
        % termination criteria
        if g>=input.maxIter 
            termination = 1;
        elseif evals.sum>=input.budget
            termination = 1;
        end
        
        % update generation counter
        g=g+1; 
        
        % update epsilon-level threshold
        if(g>1 && g<TC)
          Epsilon=EPSILON*((1-g/TC)^CP);
        elseif(g+1>=TC)
          Epsilon=0;
        end   
        
        % logging the algorithm dynamics with respect to the global best
        % observation
        dyn.gen(g+1)    = g;
        dyn.fev(g+1)    = evals.sum;
        dyn.fit(g+1)    = global_best.val;
        dyn.conv(g+1)   = global_best.conv;
        dyn.sigma(g+1)  = sigma;
        dyn.ynorm(g+1)  = norm(global_best.y);
    
        if evals.sum>=input.budget*10/100 && ll==0
            fit10       = global_best.val;
            con10       = global_best.conv;
            [ff,gg,hh]  = feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
            c10_1       = sum(gg>1) + sum(abs(hh)>1);
            c10_2       = sum((gg>0.01) & (gg<1)) + sum(abs(hh)>0.01 & abs(hh)<1);
            c10_3       = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01);          
            ll          = 1;   % set flag to 1
        elseif evals.sum>=input.budget*50/100 && lll==0
            fit50       = global_best.val;
            con50       = global_best.conv;
            [ff,gg,hh]  = feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
            c50_1       = sum(gg>1) + sum(abs(hh)>1);
            c50_2       = sum((gg>0.01)&(gg<1)) + sum(abs(hh)>0.01 & abs(hh)<1);
            c50_3       = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01)  ;
            lll         = 1; % set flag to 1
        end
        
    end
    fit100=global_best.val;
    con100=global_best.conv;
    
    [ff,gg,hh]=feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
        c100_1    = sum(gg>1) + sum(abs(hh)>1);
        c100_2    = sum((gg>0.01)&(gg<1)) + sum(abs(hh)>0.01 & abs(hh)<1);
        c100_3    = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 &abs(hh)<0.01);

        
    out= [fit10 con10 c10_1 c10_2 c10_3; fit50 con50 c50_1 c50_2 c50_3; fit100 con100 c100_1 c100_2 c100_3];
    
end
