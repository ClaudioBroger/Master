function [cost,paraopt] = optimize_funobj(paramSpecs,para0,funobjstr,model,ExpDataPath, FlPlot, FlMode, PIdx,maxeval, threshold, noise)

% minimize ojective function funobjstr with MEIGO
    npara = length(para0);
    problem = [];
    opts = [];
    problem.f= funobjstr;
    problem.x_L= paramSpecs.bmin(PIdx);
    problem.x_U= paramSpecs.bmax(PIdx);
    problem.x_0 = para0;
    problem.vtr= 0; %1e-5; % or use threshold
    opts.tolc = 1e-7; %Claude 1e-1
    opts.iterprint = 1;
    opts.ndiverse = 100*npara; %Claude: 10
    opts.maxeval = maxeval;
    opts.maxtime = 30000; %20000
    opts.inter_save = 1;
    opts.local.solver  = 'fminsearch'; 
    opts.local.iterprint = 0;
    opts.local.n1  = 10;
    opts.local.n2  = 10;  %Claude:5
    opts.local.finish  = 0;
    opts.local.bestx  = 1;
    opts.local.tol    = 3; % Claude 2
    opts.log_var = find(paramSpecs.islog(PIdx));     
   

    Results = MEIGO(problem,opts,'ESS', model,ExpDataPath, FlPlot, FlMode, PIdx, noise);
    paraopt = Results.xbest;   

    cost = feval(funobjstr,paraopt,model,ExpDataPath, FlPlot, FlMode, PIdx, noise);
end
