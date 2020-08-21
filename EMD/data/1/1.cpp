A= xlsread('1.xls');
[~,B]= xlsread('1_name.xlsx');
[row,column]=size(A);
fid = fopen('histogram.txt','wt');
fid2 = fopen('histogram_name.txt','wt');

for i=1:row
count=0;
 for j=1:column 
  if(isnan(A(i,j))==0&&A(i,j)~=0)
   count=count+1;
   fprintf(fid,'%d,0,0,',j);  
   fprintf(fid,'%f ',A(i,j));  
   end
 end
  
 if(count>=1)
 fprintf(fid,'\n'); 
fprintf(fid2,'%s\n',string(B(i)));
end
end
fclose(fid);
fclose(fid2);

