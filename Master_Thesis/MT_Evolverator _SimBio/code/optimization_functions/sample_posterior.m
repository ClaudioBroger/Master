function OutV = sample_posterior(paraoptlist, funobjstr, paramSpecs, PIdx, nsamples, threshold, model,  data_path, FlPlot, FlMode, noise)

% sample the parameter posterior with hyperspace

funobj = str2func(funobjstr);
bmin = paramSpecs.bmin(PIdx)';
bmax = paramSpecs.bmax(PIdx)';
islog = logical(paramSpecs.islog(PIdx))';
bmin(islog) = log10(bmin(islog)); % put bounds in log space
bmax(islog) = log10(bmax(islog));

funlog = @(alpha)funobj(inlinspace(alpha,islog),model,  data_path, FlPlot, FlMode, PIdx, noise);
alphaoptlistlog = paraoptlist;
alphaoptlistlog(:,islog) = log10(paraoptlist(:,islog));


MCsamples = [];

for ipoint = 1:size(alphaoptlistlog,1)

    OutM = MCexp(funlog,threshold,alphaoptlistlog(ipoint,:),bmax,bmin,nsamples/size(alphaoptlistlog,1));
    MCsamples = [MCsamples; OutM.V];

end

OutE = ELexp(funlog,threshold,MCsamples,bmax,bmin,2*nsamples);

OutV = Volint(funlog,threshold,[OutM.V; OutE.V],bmax,bmin,nsamples);

end