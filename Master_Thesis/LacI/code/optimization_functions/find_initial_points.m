function [paraoptlist, costlist] = find_initial_points(funobjstr,model, ExpDataPath, FlPlot, FlMode, paramSpecs, PIdx, npoints, maxeval, threshold, noise)

    % create matrix of optimal parameter points with starting points sampled by latin hypercube
    nparaest = length(PIdx);
    bmin = paramSpecs.bmin(PIdx);
    bmax = paramSpecs.bmax(PIdx);
    islog = paramSpecs.islog(PIdx);
    lhsmatrix = lhsdesign(npoints,nparaest);
    paraestinit = zeros(npoints,nparaest);
    paraoptlist = zeros(npoints,nparaest);
    costlist = zeros(npoints,1);
    for ipoint = 1:npoints
        for iparaest = 1:nparaest
            if islog(iparaest)
                paraestinit(ipoint,iparaest) = 10^(lhsmatrix(ipoint,iparaest)*(log10(bmax(iparaest)) - log10(bmin(iparaest))) + log10(bmin(iparaest)));
            else
                paraestinit(ipoint,iparaest) = lhsmatrix(ipoint,iparaest)*(bmax(iparaest) - bmin(iparaest)) + bmin(iparaest);
            end
        end
    end
    
    parfor ipoint = 1:npoints %parfor on cluster  
        [cost, alphaopt] = optimize_funobj(paramSpecs,paraestinit(ipoint,:),funobjstr, model, ExpDataPath, FlPlot, FlMode, PIdx, maxeval, threshold, noise);
        costlist(ipoint) = cost;
        paraoptlist(ipoint,:) = alphaopt;    
    end

end
