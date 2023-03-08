

for k= 1:num_draws
        [hill(k,1) ec50(k,1)] = doseResponse(SimFluoValues1combined.Dose(SimFluoValues1combined.Draw == k), SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k));

end

for k= 1:num_draws
        [hill2(k,1) ec502(k,1)] = doseResponse(SimFluoValues2combined.Dose(SimFluoValues2combined.Draw == k), SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k));

end

for k= 1:num_draws
        [hill3(k,1) ec503(k,1)] = doseResponse(SimFluoValues3combined.Dose(SimFluoValues3combined.Draw == k), SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k));

end


plot(ec50)
plot(hill)
plot(ec502)
plot(hill2)
plot(ec503)
plot(hill3)



SimFluoValuescombined([SimFluoValuescombined.Repression_Coefficient == 'LacIrep1' SimFluoValuescombined.Draw == 3])


reps = ["LacIrep1" "LacIrep2" "LacIrep3"];
for k = 1:num_draws
    for j = 1:length(reps)
        Values = SimFluoValuescombined(SimFluoValuescombined.Repression_Coefficient == reps(j),:);
        [hillcombined(k+j,1) ec50combined(k+j,1)] = doseResponse(Values.Dose(Values.Draw == k), Values.SimFluoValues(Values.Draw==k));
    end
end
