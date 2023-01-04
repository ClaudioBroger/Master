function paralin = inlinspace(paraestlog, islog)
paralin = paraestlog;
paralin(:,islog) = 10.^paraestlog(:,islog);
end