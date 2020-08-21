
package emd_package;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
public class PictureTotxt {
	 public static void main(String[] argv) throws IOException  {
		 FileWriter fw = new FileWriter("data\\histogram.txt");
		 FileWriter fw2 = new FileWriter("data\\histogram_name.txt");
		 int i,j,k,n=0,y=0;
		 data[] database;
		 feature_t[] f;
		 picture wa;
		 int[][][] pixel;
		 double[]  w;
			database=new data[5];
			for(y=0;y<database.length;y++) {
				n=0;
				f=new feature_t[100000]; 
				w= new double[100000];
				wa= new picture("data\\"+Integer.toString(y+1)+".jpg");
				 pixel=wa.value();
			 for(i=0;i<256;i++)
				 for(j=0;j<256;j++)
					 for(k=0;k<256;k++) {
						 if(pixel[i][j][k]!=0) {
						 f[n]=new feature_t(i,j,k);
						 	w[n]=pixel[i][j][k];
						 	
						 	n++;
					 }
						 }
			
			 for(i=0;i<256;i++)
				 for(j=0;j<256;j++)
					 for(k=0;k<256;k++) {
						 if(pixel[i][j][k]!=0) {
							 fw.write(String.valueOf(i));
							 fw.write(",");
							 fw.write(String.valueOf(j));
							 fw.write(",");
							 fw.write(String.valueOf(k));
							 fw.write(",");
							 fw.write(String.valueOf(pixel[i][j][k]));
							 fw.write(" ");
						 }
						 
					
					 }
			 fw.write("\r\n");
			 fw2.write(String.valueOf(Integer.toString(y+1)+".jpg"));
			 fw2.write("\r\n");
			 
			}
			System.out.print("complete");
			fw.close();
			fw2.close();
	 
	 }
	
	 }




