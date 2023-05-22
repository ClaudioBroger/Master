special2 = table2array(rand_parameter(2,:));
special18 = table2array(rand_parameter(18,:));

average = rand_parameter;
average(2,:) = [];
average(18,:) = [];

average = table2array(average);

average = mean(average);


for k = 1:height(average)
    
    [h2(k), p2(k)] = ttest2(special2, average(k,:));

    [h18(k), p18(k)] = ttest2(special18, average(k,:));

end

[h2, p2] = ttest2(special2, average);

[h18, p18] = ttest2(special18, average);