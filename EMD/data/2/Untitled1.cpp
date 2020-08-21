A= xlsread('2.xlsx');

[row,column]=size(A);
fid = fopen('histogram.txt','wt');


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

end
end
fclose(fid);


