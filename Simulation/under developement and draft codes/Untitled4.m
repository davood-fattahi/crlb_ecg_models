clc
lineLength=fprintf('%s will be %d this year.\n',name,age);

for i=1:20
name = 'Alice';   
age = i;
fprintf(repmat('\b', 1, lineLength));
lineLength=fprintf('%s will be %d this year.\n',name,age);
end