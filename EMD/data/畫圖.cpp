clear all;
cd C:\Users\88696\eclipse-workspace\test\data

C= dlmread("vector.txt", ' ')	;
C=transpose(C);
[numRows,numCols]=size(C);

for j=0:numCols-1
figure(j+1);
filename=sprintf("p=%d,1.txt",j);
A = dlmread(filename, ' ')	;
A=transpose(A);
[numRows2,numCols2]=size(A);
N=power(numCols2,0.5);
filename=sprintf("p=%d,2.txt",j);
B = dlmread(filename, ' ')	;
B=transpose(B);
[numRows3,numCols3]=size(B);
D=numCols3;
x=zeros(1,2*N);
y=zeros(1,2*N);

for i=1:N
hold on
plot(A(1,1+N*(i-1):i*N),A(2,1+N*(i-1):i*N))
end



for i=2:2:2*N

x(i-1)=A(1,i/2);
x(i)=A(1,i/2+N*(N-1));

y(i-1)=A(2,i/2);
y(i)=A(2,i/2+N*(N-1));

end

hold on
for i=2:2:2*N
hold on
plot(x(i-1:i),y(i-1:i))
end

hold on
plot(B(1,1:D),B(2,1:D),'o')
vector_title=sprintf("vector=(%f,%f,%f)",C(1,j+1),C(2,j+1),C(3,j+1));

title(vector_title)


end


