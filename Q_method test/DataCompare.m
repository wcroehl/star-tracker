clear all, clear clc 
%script for comapring any similar two data

%Input data
FILENAME1=input("What is the first file name?: ","s");
FILENAME2=input("What is the second file name?: ","s");
level=0;
k=0;
data1=load(FILENAME1);
data2=load(FILENAME2);
[row,column]=size(data1);

%level adjustment
while k == 0
level=input("What's level? 1-easy, 2-normal, 3-hard: ");
   if level == 1 || level ==2 || level == 3
       k=1;
   end
end

switch level
    case 1
    data1 = round(data1,2);
    data2 = round(data2,2);
    
    case 2
    data1 = round(data1,4);
    data2 = round(data2,4);
end

%Data comparison 
for i=1:row
   for k=1:column
       if data1(i,k) ~= data2(i,k) 
           fprintf("(%.4f,%.4f)",i,k)
           error("data does not match")
       end
   end
end