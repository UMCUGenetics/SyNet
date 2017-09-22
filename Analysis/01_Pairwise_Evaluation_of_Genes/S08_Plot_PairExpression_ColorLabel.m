i = 7688;
j = 7832;
close all
plot(zData(Patient_Label==0,i), zData(Patient_Label==0, j), 'ob');
hold on
plot(zData(Patient_Label==1,i), zData(Patient_Label==1, j), '+r');
xlabel(Gene_Name{i});
ylabel(Gene_Name{j});