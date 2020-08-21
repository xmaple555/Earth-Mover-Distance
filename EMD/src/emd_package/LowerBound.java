package emd_package;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.special.Erf;

public class LowerBound {
	int bound_region_number = 20;
	int sub_number = 255*2;
	double[][] vector = { { 1, 0, 0 } ,{ 0, 1,0 },{ 0, 0, 1}  };
	double threshold = 10000;
	double threshold2 = 1000;
	int DEBUG_LEVEL = 1;
	int data_number = 5;
	int array_length_x = 256;
	int array_length_y = 256;
	int array_length_z = 256;
	final int MAX_SIG_SIZE = 3000;
	/////////////////
	int number_of_performance_of_EMD = 0;

	double tmin, tmax;
	double[] sub;
	double b_min, b_max, m_min, m_max;
	double[] SS = new double[2];
	double[] SN = new double[2];
	double[] SW = new double[2];
	double[] SE = new double[2];

	public static void main(String[] argv) throws IOException, CloneNotSupportedException {
		LowerBound a = new LowerBound();

		double temp = 0;
		double totoa_weight_Q = 0;
		double totoa_weight_P = 0;
		int flag = 1;
		data[] database;
		database = new data[a.data_number];
		float[][][] pixel;
		data_Q[] Q = new data_Q[a.data_number];
		int i, j, k, n = 0, n_Q = 0, x = 0, y = 0, z = 0;
		feature_t[] f, f_Q;
		double[] w, w_Q;
		signature_t s1, s2;
		double[] Q_mean;
		double[] P_mean;
		String[] histogram_name = new String[a.data_number];
		Multiplie_region[] G = new Multiplie_region[a.data_number];
		for (i = 0; i < a.data_number; i++) {
			G[i] = new Multiplie_region();
			Q[i] = new data_Q();
			Q[i].Q = new data[a.vector.length];
			Q[i].sub_Q = new double[a.vector.length][a.sub_number];
			Q[i].no_intersection_error = new double[a.vector.length];
		}

		for (x = 0; x < a.vector.length; x++) {
			for (y = 0; y < database.length; y++) {
				n = 0;
				f = new feature_t[a.MAX_SIG_SIZE];
				w = new double[a.MAX_SIG_SIZE];

				pixel = a.histogram(y);

				for (i = 0; i < a.array_length_x; i++)
					for (j = 0; j < a.array_length_y; j++)
						for (k = 0; k < a.array_length_z; k++) {
							if (pixel[i][j][k] != 0) {
								f[n] = new feature_t(i, j, k);
								w[n] = pixel[i][j][k];

								n++;
							}
						}

				pixel = null;
				database[y] = new data(f, w, a.vector[x], n);
				database[y].name = y + 1;

			}
			a.init(database);
			a.Bounding_Region(G, x, database);
			for (i = 0; i < a.data_number; i++) {

				Q[i].Q[x] = database[i];
				if (Q[i].Q[x].standard_deviation_is_zero == 0) {
					Q[i].Q[x].get_feature_and_weight();
					Q[i].sub_Q[x] = a.process_Q_error(database[i], x);
					Q[i].no_intersection_error[x] = a.Err_integration(Q[i].Q[x], a.tmax, a.tmin);
					Q[i].Q[x].f = null;
					Q[i].Q[x].w = null;
				}
			}
		}
		histogram_name = a.histogram_name();
		for (x = 0; x < a.data_number; x++) {

			n_Q = 0;
			f_Q = new feature_t[a.MAX_SIG_SIZE];
			w_Q = new double[a.MAX_SIG_SIZE];
			totoa_weight_Q = 0;
			pixel = a.histogram(x);
			for (i = 0; i < a.array_length_x; i++)
				for (j = 0; j < a.array_length_y; j++)
					for (k = 0; k < a.array_length_z; k++) {
						if (pixel[i][j][k] != 0) {
							f_Q[n_Q] = new feature_t(i, j, k);
							w_Q[n_Q] = pixel[i][j][k];
							totoa_weight_Q = totoa_weight_Q + w_Q[n_Q];
							n_Q++;

						}
					}

			Q_mean = a.mean(w_Q, n_Q);
			for (i = 0; i < n_Q; i++) {
				w_Q[i] = w_Q[i] / totoa_weight_Q;

			}

			pixel = null;
			for (y = 0; y < a.data_number; y++) {
				if (x == y) {

				} else {
					flag = 1;
					f = new feature_t[a.MAX_SIG_SIZE];
					n = 0;
					f = new feature_t[a.MAX_SIG_SIZE];
					w = new double[a.MAX_SIG_SIZE];
					pixel = a.histogram(y);
					totoa_weight_P = 0;
					for (i = 0; i < a.array_length_x; i++)
						for (j = 0; j < a.array_length_y; j++)
							for (k = 0; k < a.array_length_z; k++) {
								if (pixel[i][j][k] != 0) {
									f[n] = new feature_t(i, j, k);
									w[n] = pixel[i][j][k];
									totoa_weight_P = totoa_weight_P + w[n];
									n++;
								}
							}
					P_mean = a.mean(w, n);

					for (i = 0; i < n; i++) {
						w[i] = w[i] / totoa_weight_P;

					}

					pixel = null;
					for (j = 0; j < a.vector.length; j++) {

						if (a.DEBUG_LEVEL == 1)
							System.out.println("Q=" + histogram_name[x] + " P=" + histogram_name[y]);

						if (G[y].g[j] != null && Q[x].Q[j].standard_deviation_is_zero == 0) {
							temp = a.EMDLowerBound(G[y].g[j], Q[x], j);
							if (a.DEBUG_LEVEL == 1)
								System.out.println("vector=(" + a.vector[j][0] + "," + a.vector[j][1] + ","
										+ a.vector[j][2] + ")" + ",EMDLB=" + temp);

							if (temp >= a.threshold) {
								flag = 0;
								break;
							}
						}
					}
					if (flag == 1) {
						s1 = new signature_t(n, f, w);
						s2 = new signature_t(n_Q, f_Q, w_Q);

						if (a.EMDJoin(s1, s2) == 1) {
							for (z = 0; z < P_mean.length; z++) {
								if (Math.abs(P_mean[z] - Q_mean[z]) > a.threshold2) {
									flag = 0;
									break;
								}
							}

							if (flag == 1)
								System.out
										.println("Q=" + histogram_name[x] + " P=" + histogram_name[y] + " are similar");

						}
					}

				}
			}
		}
		System.out.println("after pruned,number of performance of EMD:  " + a.number_of_performance_of_EMD);
		System.out.println("before pruned,number of performance of EMD: " + (a.data_number * (a.data_number - 1)));
	}

	double[] mean(double[] w, int n) {
		double[] mean = new double[vector.length]; // histogram 的維度
		int i, j;
		for (j = 0; j < mean.length; j++) {
			for (i = 0; i < n; i++) {
				mean[j] = mean[j] + w[i];
			}
			mean[j] = mean[j] / n;
		}
		return mean;
	}

	double CDF_normal(double x, double mean, double standard_deviation) {

		return (0.5 * (1 + Erf.erf((x - mean) / (standard_deviation * Math.pow(2, 0.5)))));
	}

	void init(data[] database) throws IOException {
		int i, j;

		sub = new double[sub_number + 1];

		double temp2;

		tmin = Double.MAX_VALUE;
		tmax = -Double.MAX_VALUE;
		for (j = 0; j < database.length; j++) {

			database[j].get_feature_and_weight();
			for (i = 0; i < database[j].n; i++) {

				if (database[j].f[i].x > tmax)
					tmax = database[j].f[i].x;
				if (database[j].f[i].x < tmin)
					tmin = database[j].f[i].x;
			}
		}

		temp2 = (tmin + tmax) / 2;

		for (j = 0; j < database.length; j++) {
			database[j].t_change(temp2);
			database[j].process_data();

		}
		tmax = tmax - (int) temp2;
		tmin = tmin - (int) temp2;
		m_min = database[0].m;
		m_max = database[0].m;
		b_min = database[0].b;
		b_max = database[0].b;
		for (j = 0; j < database.length; j++) {
			if (database[j].b > b_max)
				b_max = database[j].b;
			if (database[j].b < b_min)
				b_min = database[j].b;

			if (database[j].m > m_max)
				m_max = database[j].m;
			if (database[j].m < m_min)
				m_min = database[j].m;

		}

		for (i = 0; i <= sub_number; i++)
			sub[i] = tmin + i * (tmax - tmin) / sub_number;

	}

	double Err_integration(data P, double upper, double lower) throws IOException {
		double value;

		value = CDF_Integration(P, upper, lower);

		value = value - Normal_CDF_Integration(P.mean, P.standard_deviation, upper, lower);

		return value;
	}

	double CDF_Integration(data P, double upper, double lower) throws IOException {

		double value = 0;

		int i = 0;
		value = value + CDF(lower, P) * ((Math.ceil(lower) - lower));
		for (i = (int) Math.ceil(lower); i < (int) Math.floor(upper); i++) {
			value = value + CDF(i, P);
		}

		value = value + CDF(Math.floor(upper), P) * (((upper - Math.floor(upper))));

		return value;
	}

	float[][][] histogram(int x) throws IOException {
		int i;
		String[] keyword;
		String[] keyword2;
		String[] keyword3;
		String strFile = "data\\histogram.txt";
		BufferedReader input = new BufferedReader(new FileReader(strFile));
		StringBuffer sb = new StringBuffer();

		String strInput = input.readLine();
		while (strInput != null) {
			sb.append(strInput).append("\n");

			strInput = input.readLine();
		}

		keyword = sb.toString().split("\n");

		float[][][] pixel2 = new float[array_length_x][array_length_y][array_length_z];
		keyword2 = keyword[x].split(" ");

		for (i = 0; i < keyword2.length; i++) {
			keyword3 = keyword2[i].split(",");
			pixel2[Integer.valueOf(keyword3[0])][Integer.valueOf(keyword3[1])][Integer.valueOf(keyword3[2])] = Float
					.valueOf(keyword3[3]);

		}

		input.close();

		return pixel2;
	}

	String[] histogram_name() throws IOException {
		String[] keyword;
		String strFile = "data\\histogram_name.txt";
		BufferedReader input = new BufferedReader(new FileReader(strFile));
		StringBuffer sb = new StringBuffer();

		String strInput = input.readLine();
		while (strInput != null) {
			sb.append(strInput).append("\n");

			strInput = input.readLine();
		}

		keyword = sb.toString().split("\n");
		return keyword;
	}

	double PDF(double x, double mean, double standard_deviation) {

		return (Math.exp(-(Math.pow(x - mean, 2) / (2 * Math.pow(standard_deviation, 2))))
				/ (standard_deviation * Math.pow(2 * Math.PI, 0.5)));

	}

	double Normal_CDF_Integration(double mean, double standard_deviation, double upper, double lower) {
		double value = 0;

		value = value + CDF_normal(upper, mean, standard_deviation) * (upper - mean);
		value = value + PDF(upper, mean, standard_deviation);

		value = value - CDF_normal(lower, mean, standard_deviation) * (lower - mean);
		value = value - PDF(lower, mean, standard_deviation);

		return value;
	}

	double EMD_Normal(double mean_P, double standard_deviation_P, double mean_Q, double standard_deviation_Q,
			double tis) {

		if (Boolean.valueOf(tis >= tmin && tis <= tmax))
			return (Math
					.abs(Normal_CDF_Integration(mean_P, standard_deviation_P, tis, tmin)
							- Normal_CDF_Integration(mean_Q, standard_deviation_Q, tis, tmin))
					+ Math.abs(Normal_CDF_Integration(mean_P, standard_deviation_P, tmax, tis)
							- Normal_CDF_Integration(mean_Q, standard_deviation_Q, tmax, tis)));
		else {
			return Math.abs(Normal_CDF_Integration(mean_P, standard_deviation_P, tmax, tmin)
					- Normal_CDF_Integration(mean_Q, standard_deviation_Q, tmax, tmin));
		}
	}

	double get_tis(double mean_P, double standard_deviation_P, double mean_Q, double standard_deviation_Q) {
		return ((mean_P * standard_deviation_Q - mean_Q * standard_deviation_P)
				/ (standard_deviation_Q - standard_deviation_P));
	}

	double CDF(double n, data P) {
		int i;
		double value = 0;

		feature_t[] f = P.f;
		double[] w = P.w;

		if (n < P.t_max) {
			for (i = 0; i < P.n; i++) {

				if (f[i].x <= n)
					value = value + w[i];

			}

		} else
			return 1;
		return (value);
	}

	double Errmin(data P, double tj, double[] sub_min) throws IOException {
		if (tj < tmin || tj >= tmax) {
			P.get_feature_and_weight();
			double temp = Err_integration(P, tmax, tmin);
			P.f = null;
			P.w = null;
			return temp;

		} else {

			int i;
			for (i = 0; i < sub_number; i++) {
				if (tj >= sub[i] && tj < sub[i + 1]) {
					return sub_min[i];

				}
			}
			System.out.println("err error");
			return 0;

		}
	}

	double Errmax(data P, double tj, double[] sub_max) throws IOException {
		if (tj < tmin || tj >= tmax) {
			P.get_feature_and_weight();
			double temp = Err_integration(P, tmax, tmin);
			P.f = null;
			P.w = null;
			return temp;

		} else {

			int i;
			for (i = 0; i < sub_number; i++) {
				if (tj >= sub[i] && tj < sub[i + 1]) {
					return sub_max[i];

				}
			}
			System.out.println("err error");
			return 0;

		}
	}

	region[] Bounding_Region(Multiplie_region[] G, int p, data[] P) throws IOException {
		int n = bound_region_number;
		int i, j;
		double temp3;

		double sin1, cos1, tan1;
		double sin2, cos2, tan2;
		double[] m_plus, b_plus;
		m_plus = new double[n];
		b_plus = new double[n];
		double[][] position = new double[(n + 1) * (n + 1)][2];
		double[] temp2 = new double[P.length];
		region[] a = new region[(n) * (n)];

		SS = get_intersection(-tmin, m_max, b_min, -tmax, m_min, b_min);
		SN = get_intersection(-tmin, m_min, b_max, -tmax, m_max, b_max);
		SE = get_intersection(-tmin, m_max, b_min, -tmax, m_max, b_max);
		SW = get_intersection(-tmin, m_min, b_max, -tmax, m_min, b_min);

		temp3 = Euclid(SW, SS) / n;

		for (i = 0; i < n; i++) {

			m_plus[i] = Euclid(SW, SS) * (i + 1) / n;
			b_plus[i] = Euclid(SE, SS) * (i + 1) / n;
		}

		tan1 = Math.abs((SW[1] - SS[1]) / (SW[0] - SS[0]));
		cos1 = 1 / Math.pow((1 + tan1 * tan1), 0.5);
		sin1 = tan1 * cos1;
		tan2 = Math.abs((SE[1] - SS[1]) / (SE[0] - SS[0]));
		cos2 = 1 / Math.pow((1 + tan2 * tan2), 0.5);
		sin2 = tan2 * cos2;
		for (j = 0; j < (n + 1); j++) {
			if (j == 0) {
				temp2[0] = SS[0];
				temp2[1] = SS[1];

			}

			else {
				temp2[0] = SS[0] - m_plus[j - 1] * cos1;
				temp2[1] = SS[1] + m_plus[j - 1] * sin1;

			}

			for (i = 0; i < (n + 1); i++) {

				if (i == 0) {

					position[i + j * (n + 1)][0] = temp2[0];
					position[i + j * (n + 1)][1] = temp2[1];
					continue;
				}

				else {

					position[i + j * (n + 1)][0] = temp2[0] + b_plus[i - 1] * cos2;
					position[i + j * (n + 1)][1] = temp2[1] + b_plus[i - 1] * sin2;
				}

			}

		}

		for (i = 0; i < (n); i++) {
			for (j = 0; j < (n); j++) {

				a[j + i * (n)] = new region(position[j + i * (n + 1)], position[j + i * (n + 1) + (n + 1 + 1)],
						position[j + i * (n + 1) + 1], position[j + i * (n + 1) + (n + 1 + 1 - 1)]);

			}

		}

		for (i = 0; i < P.length; i++) {
			if (P[i].standard_deviation_is_zero == 0) {
				A: for (j = 0; j < ((n) * (n)); j++) {
					if (j < n) {

						for (; j < n; j++) {
							if (P_dominates_Q(P[i].m, P[i].b, a[j].Mu.m, a[j].Mu.b) == 1) {
								a[j].put_data(P[i]);
								a[j].tmax = tmax;
								a[j].tmin = tmin;
								a[j].sub = sub;
								a[j].SE = SE;
								a[j].SW = SW;
								a[j].SN = SN;
								a[j].SS = SS;

								G[i].g[p] = a[j];

								break A;
							}
						}

					}

					if (P_dominates_Q(P[i].m, P[i].b, a[j].Mu.m, a[j].Mu.b) == 1
							&& P_dominates_Q(a[j].Ml.m, a[j].Ml.b, P[i].m, P[i].b) == 1) {
						a[j].put_data(P[i]);
						a[j].tmax = tmax;
						a[j].tmin = tmin;
						a[j].sub = sub;
						a[j].SE = SE;
						a[j].SW = SW;
						a[j].SN = SN;
						a[j].SS = SS;

						G[i].g[p] = a[j];

						break A;

					}
					if (Boolean.valueOf(j < n * n && j >= (n * n - n))) {
						for (j = (n * n - 1); j >= (n * n - n); j--) {
							if (P_dominates_Q(a[j].Ml.m, a[j].Ml.b, P[i].m, P[i].b) == 1) {
								a[j].put_data(P[i]);
								a[j].tmax = tmax;
								a[j].tmin = tmin;
								a[j].sub = sub;
								a[j].SE = SE;
								a[j].SW = SW;
								a[j].SN = SN;
								a[j].SS = SS;

								G[i].g[p] = a[j];

								break A;
							}
						}
						System.out.println("bound region error");
						break A;
					}

				}

			}
		}

		 plot(position,P,p);
		for (i = 0; i < a.length; i++) {
			if (a[i].n != 0) {
				a[i].error_process(sub, tmin, tmax);
			}

		}

		return a;
	}

	void plot(double[][] position, data[] P, int p) throws IOException {
		FileWriter fw = new FileWriter("data\\p=" + p + ",1.txt");
		FileWriter fw2 = new FileWriter("data\\p=" + p + ",2.txt");

		int i, j;
		for (i = 0; i < position.length; i++) {
			fw.write(String.valueOf(position[i][0]));
			fw.write(" ");
			fw.write(String.valueOf(position[i][1]));
			fw.write("\r\n");
		}

		for (i = 0; i < P.length; i++) {
			fw2.write(String.valueOf(P[i].m));
			fw2.write(" ");
			fw2.write(String.valueOf(P[i].b));
			fw2.write("\r\n");
		}
		fw.close();
		fw2.close();

		if (p == vector.length - 1) {
			FileWriter fw3 = new FileWriter("data\\vector.txt");
			for (i = 0; i < vector.length; i++) {
				for (j = 0; j < 3; j++) {
					fw3.write(String.valueOf(vector[i][j]));
					if (j != 2)
						fw3.write(" ");
				}

				fw3.write("\r\n");
			}
			fw3.close();

			System.out.println("complete");
			System.exit(1);
		}
	}

	double Euclid(double[] x, double[] y) {
		return (Math.pow(Math.pow(x[0] - y[0], 2) + Math.pow(x[1] - y[1], 2), 0.5));

	}

	int P_dominates_Q(double m_P, double b_P, double m_Q, double b_Q) {

		if (Boolean.valueOf((b_Q >= (-tmin) * m_Q + m_P * tmin + b_P || Math.floor(b_Q * 1000000)
				/ 1000000 == Math.floor((((-tmin) * m_Q + m_P * tmin + b_P) * 1000000)) / 1000000)))
			if (Boolean.valueOf((b_Q >= (-tmax) * m_Q + m_P * tmax + b_P || Math.floor(b_Q * 1000000)
					/ 1000000 == Math.floor((((-tmax) * m_Q + m_P * tmax + b_P) * 1000000)) / 1000000)))
				return 1;

		return 0;
	}

	double Errmin_M(region P, double tj) throws IOException {
		int i;
		double value;
		double temp;

		value = Errmin(P.M[0], tj, P.sub_min);
		for (i = 1; i < P.n; i++) {

			temp = Errmin(P.M[i], tj, P.sub_min);
			if (value > temp) {
				value = temp;
			}
		}
		return value;
	}

	double Errmax_M(region P, double tj) throws IOException {
		int i;
		double value;
		double temp;
		value = Errmax(P.M[0], tj, P.sub_max);
		for (i = 1; i < P.n; i++) {

			temp = Errmax(P.M[i], tj, P.sub_max);

			if (value < temp) {
				value = temp;
			}
		}
		return value;
	}

	double[] process_Q_error(data Q, int p) throws IOException {
		int i;
		double tis;

		double[][] sub_Q;
		sub_Q = new double[vector.length][sub_number];
		for (i = 0; i < sub_number; i++) {
			double temp2;
			tis = (sub[i] + sub[i + 1]) / 2;
			temp2 = Err_integration(Q, tis, tmin);
			temp2 = temp2 - Err_integration(Q, tmax, tis);

			sub_Q[p][i] = temp2;

		}
		return sub_Q[p];

	}

	double EMDLowerBound(region a, data_Q Q_class, int p) throws IOException {
		double tis;
		double value = 0;
		tmax = a.tmax;
		tmin = a.tmin;
		sub = a.sub;
		SW = a.SW;
		SE = a.SE;
		SN = a.SN;
		SS = a.SS;
		double[] temp;
		data Mis;
		double[][] sub_Q = Q_class.sub_Q;
		data Q = Q_class.Q[p];
		if (P_dominates_Q(a.Mu.m, a.Mu.b, Q.m, Q.b) == 1 && P_dominates_Q(a.Ml.m, a.Ml.b, Q.m, Q.b) == 1) {
			tis = get_tis(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation);
			if (tmin < tis && tis < tmax)
				System.out.println("tis error");
			value = value + EMD_Normal(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation, tis);
			value = value - a.no_intersection_error_max;
			value = value + Q_class.no_intersection_error[p];
			return (value);
		} else if (P_dominates_Q(Q.m, Q.b, a.Ml.m, a.Ml.b) == 1 && P_dominates_Q(Q.m, Q.b, a.Mu.m, a.Mu.b) == 1) {
			tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
			if (tmin < tis && tis < tmax)
				System.out.println("tis error");
			value = value + EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation, tis);
			value = value - Q_class.no_intersection_error[p];
			value = value + a.no_intersection_error_min;
			return (value);
		} else if (P_dominates_Q(a.Mu.m, a.Mu.b, Q.m, Q.b) == 0 && P_dominates_Q(a.Ml.m, a.Ml.b, Q.m, Q.b) == 1) {

			if (P_dominates_Q(a.gw_m, a.gw_b, Q.m, Q.b) == 1) {
				temp = get_intersection(-tmin, a.Mu.m, a.Mu.b, -tmax, Q.m, Q.b);
				Mis = new data(temp[0], temp[1]);
				tis = get_tis(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation, tis);

				value = value + error(a, a.Mu, Q, p, sub_Q);

				return (value);
			}

			if (P_dominates_Q(a.ge_m, a.ge_b, Q.m, Q.b) == 1) {
				temp = get_intersection(-tmax, a.Mu.m, a.Mu.b, -tmin, Q.m, Q.b);
				Mis = new data(temp[0], temp[1]);
				tis = get_tis(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation, tis);
				value = value + error(a, a.Mu, Q, p, sub_Q);

				return (value);
			}

		} else if (P_dominates_Q(Q.m, Q.b, a.Mu.m, a.Mu.b) == 1 && P_dominates_Q(Q.m, Q.b, a.Ml.m, a.Ml.b) == 0) {
			if (P_dominates_Q(Q.m, Q.b, a.gw_m, a.gw_b) == 1) {
				temp = get_intersection(-tmax, a.gw_m, a.gw_b, -tmin, Q.m, Q.b);
				Mis = new data(temp[0], temp[1]);
				tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation, tis);
				tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + error(a, a.Ml, Q, p, sub_Q);

				return (value);
			}
			if (P_dominates_Q(Q.m, Q.b, a.ge_m, a.ge_b) == 1) {
				temp = get_intersection(-tmin, a.Ml.m, a.Ml.b, -tmax, Q.m, Q.b);
				Mis = new data(temp[0], temp[1]);
				tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Mis.mean, Mis.standard_deviation, tis);
				value = value + error(a, a.Ml, Q, p, sub_Q);

				return (value);
			}

		}
		if (P_dominates_Q(a.Mu.m, a.Mu.b, Q.m, Q.b) == 0 && P_dominates_Q(a.Ml.m, a.Ml.b, Q.m, Q.b) == 0) {
			double[] temp3;
			temp = get_intersection(-tmin, a.gw_m, a.gw_b, -tmax, SW[0], SW[1]);
			temp3 = get_intersection(-tmax, a.ge_m, a.ge_b, -tmin, SE[0], SE[1]);

			if (P_dominates_Q(temp[0], temp[1], Q.m, Q.b) == 1) {
				double temp1, temp2;

				Mis = new data(a.gw_m, a.gw_b);
				tis = get_tis(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation, tis);
				value = value + error(a, a.Mu, Q, p, sub_Q);
				temp1 = value;
				value = 0;

				tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation, tis);
				value = value + error(a, a.Ml, Q, p, sub_Q);
				temp2 = value;

				if (temp1 > temp2)
					return temp2;
				else
					return temp1;

			} else if (P_dominates_Q(temp3[0], temp3[1], Q.m, Q.b) == 1) {

				double temp1, temp2;
				Mis = new data(a.ge_m, a.ge_b);
				tis = get_tis(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Mu.mean, a.Mu.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation);
				value = value
						- 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, a.Mu.mean, a.Mu.standard_deviation, tis);
				value = value + error(a, a.Mu, Q, p, sub_Q);
				temp1 = value;
				value = 0;

				tis = get_tis(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(a.Ml.mean, a.Ml.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation);
				value = value + 0.5 * EMD_Normal(Mis.mean, Mis.standard_deviation, Q.mean, Q.standard_deviation, tis);
				tis = get_tis(Mis.mean, Mis.standard_deviation, a.Ml.mean, a.Ml.standard_deviation);
				value = value + error(a, a.Ml, Q, p, sub_Q);
				temp2 = value;

				if (temp1 > temp2)
					return temp2;
				else
					return temp1;

			}

		}

		return 0;

	}

	double error(region a, data P, data Q, int p, double[][] sub_Q) throws IOException {
		double value = 0;// P代表REGION

		if (Q.standard_deviation < P.standard_deviation) {

			value = value - si_Q_Errmax(sub_Q[p]);
			value = value + si_Errmin(a);

		} else {

			value = value - si_Errmax(a);
			value = value + si_Q_Errmin(sub_Q[p]);

		}

		/*
		 * double tis=0;
		 * tis=get_tis(P.mean,P.standard_deviation,Q.mean,Q.standard_deviation);
		 * if(Q.standard_deviation<P.standard_deviation) {
		 * 
		 * 
		 * 
		 * value=value-Errmax(Q,tis,sub_Q[p]); value=value+Errmin_M(a,tis);
		 * 
		 * 
		 * }else {
		 * 
		 * value=value-Errmax_M(a,tis); value=value+Errmin(Q,tis,sub_Q[p]);
		 * 
		 * 
		 * }
		 * 
		 * 
		 */

		return value;
	}

	double si_Errmin(region P) {
		int i;
		double value;
		double temp;
		for (value = P.sub_min[0], i = 1; i < P.sub_min.length; i++) {
			temp = P.sub_min[i];
			if (value > temp) {
				value = temp;
			}

		}

		return value;
	}

	double si_Errmax(region P) {

		int i;
		double value;
		double temp;
		for (value = P.sub_max[0], i = 1; i < P.sub_max.length; i++) {
			temp = P.sub_max[i];
			if (value < temp) {
				value = temp;
			}

		}

		return value;
	}

	double si_Q_Errmin(double[] sub_Q) {

		int i;
		double value;
		double temp;
		for (value = sub_Q[0], i = 1; i < sub_Q.length; i++) {
			temp = sub_Q[i];
			if (value > temp) {
				value = temp;
			}

		}
		return value;
	}

	double si_Q_Errmax(double[] sub_Q) {

		int i;
		double value;
		double temp;
		for (value = sub_Q[0], i = 1; i < sub_Q.length; i++) {
			temp = sub_Q[i];
			if (value < temp) {
				value = temp;
			}

		}
		return value;
	}

	double[] get_intersection(double m1, double x1, double y1, double m2, double x2, double y2) {
		double c1, c2;
		double[] value = new double[2];
		c1 = y1 - m1 * x1;
		c2 = y2 - m2 * x2;
		value[0] = (c2 - c1) / (m1 - m2);
		value[1] = m1 * value[0] + c1;
		if ((m1 - m2) == 0) {
			System.out.println("No Intersection between the lines");
			System.exit(1);
			return null;
		} else

			return value;

	}

	int EMDJoin(signature_t s1, signature_t s2) throws CloneNotSupportedException {

		double e = 0;
		double projected_EMDLB;
		projected_base_EMD_lowerbound c = new projected_base_EMD_lowerbound();
		dualEMD_class b = new dualEMD_class();

		emd_class hh = new emd_class();

		projected_EMDLB = c.PASUM(s1, s2);
		if (DEBUG_LEVEL == 1)
			System.out.println("projected_EMDLB=" + projected_EMDLB);
		if (projected_EMDLB < threshold) {

			b.dualEMD(s1, s2);
			if (DEBUG_LEVEL == 1)
				System.out.println("EMDdual=" + b.get_dualEMD_value());
			if (b.get_dualEMD_value() < threshold)

				e = hh.emd(s1, s2, null, null);
			number_of_performance_of_EMD = number_of_performance_of_EMD + 1;
			if (DEBUG_LEVEL == 1)
				System.out.println("EMD=" + e);
			if (e < threshold) {

				return 1;
			}
		}
		return 0;

	}

}

class data {
	feature_t[] f;
	double[] w;
	double shift;
	double mean, standard_deviation;
	double m, b;
	int n;
	int name;
	double[] unit;
	double totoal_weight, t_max, t_min;
	int standard_deviation_is_zero;

	public data() {

	}

	public data(double m, double b) {
		this.m = m;
		this.b = b;
		this.standard_deviation = 1 / m;
		this.mean = -this.standard_deviation * this.b;
	}

	void t_change(double t) {
		this.shift = t;

	}

	void get_feature_and_weight() throws IOException {
		int i, j, k;
		int nn = 0;
		LowerBound a = new LowerBound();
		f = new feature_t[a.MAX_SIG_SIZE];
		w = new double[a.MAX_SIG_SIZE];
		float[][][] pixel;
		pixel = a.histogram(name - 1);
		for (i = 0; i < a.array_length_x; i++)
			for (j = 0; j < a.array_length_y; j++)
				for (k = 0; k < a.array_length_z; k++) {
					if (pixel[i][j][k] != 0) {
						f[nn] = new feature_t(i, j, k);
						w[nn] = pixel[i][j][k];

						nn++;
					}

				}

		for (i = 0; i < nn; i++) {

			f[i].x = f[i].x * unit[0] + f[i].y * unit[1] + f[i].z * unit[2];
			f[i].y = 0;
			f[i].z = 0;
			for (j = 0; j < i; j++) {
				if (f[i].x == f[j].x) {

					w[j] = w[j] + w[i];
					for (k = i; k < nn; k++) {
						f[k] = f[k + 1];
						w[k] = w[k + 1];

					}

					nn--;
					i--;
				}
			}
		}
		n = nn;
		for (i = 0; i < nn; i++) {
			w[i] = w[i] / this.totoal_weight;

		}
		for (i = 0; i < nn; i++) {
			f[i].x = f[i].x - (int) shift;
		}

	}

	public data(feature_t[] f, double[] w, double[] unit, int n) {
		int i, j, k;
		this.unit = unit;
		double temp2 = 0;

		for (i = 0; i < n; i++) {
			temp2 = temp2 + w[i];
		}
		this.totoal_weight = temp2;
		for (i = 0; i < n; i++) {

			f[i].x = f[i].x * unit[0] + f[i].y * unit[1] + f[i].z * unit[2];
			f[i].y = 0;
			f[i].z = 0;
			for (j = 0; j < i; j++) {
				if (f[i].x == f[j].x) {

					w[j] = w[j] + w[i];
					for (k = i; k < n; k++) {
						f[k] = f[k + 1];
						w[k] = w[k + 1];

					}

					n--;
					i--;
				}
			}
		}
		for (i = 0; i < n; i++) {
			w[i] = w[i] / totoal_weight;

		}
		temp2 = 0;

		this.n = n;
		f = null;
		w = null;

	}

	void process_data() throws IOException {
		int i;
		double temp = 0;
		get_feature_and_weight();
		for (i = 0; i < n; i++)
			this.mean = this.mean + this.w[i] * (f[i].x) / 1;

		for (i = 0; i < n; i++)
			temp = temp + this.w[i] * Math.pow((f[i].x) - this.mean, 2);
		this.standard_deviation = Math.pow(temp / 1, 0.5);
		if (standard_deviation == 0) {
			this.standard_deviation_is_zero = 1;

		}
		if (this.standard_deviation_is_zero == 0) {
			this.m = 1 / this.standard_deviation;
			this.b = -mean / this.standard_deviation;

			for (i = 1, this.t_max = f[0].x; i < n; i++) {
				if (this.t_max < f[i].x)
					this.t_max = f[i].x;

			}
			for (i = 1, this.t_min = f[0].x; i < n; i++) {
				if (this.t_min > f[i].x)
					this.t_min = f[i].x;

			}
		}
		f = null;
		w = null;
	}
}

class region {
	int n = 0;
	LowerBound a = new LowerBound();
	data[] M = new data[a.data_number];
	data Mu = new data();
	data Ml = new data();
	double gw_m, ge_m;
	double gw_b, ge_b;
	double[] sub_max;
	double[] sub_min;
	double no_intersection_error_max;
	double no_intersection_error_min;
	double[] sub;
	double[] SS = new double[2];
	double[] SN = new double[2];
	double[] SW = new double[2];
	double[] SE = new double[2];

	double tmax, tmin;

	public region(double[] a, double[] b, double[] c, double[] d) {
		Ml = new data(a[0], a[1]);
		Mu = new data(b[0], b[1]);
		ge_m = c[0];
		ge_b = c[1];
		gw_m = d[0];
		gw_b = d[1];

	}

	void put_data(data P) {
		M[n] = P;
		n++;
	}

	void error_process(double[] sub, double tmin, double tmax) throws IOException {
		int i, j;
		double tis;
		double temp;
		sub_max = new double[sub.length - 1];
		sub_min = new double[sub.length - 1];
		LowerBound a = new LowerBound();
		for (i = 0; i < sub.length - 1; i++) {
			sub_max[i] = -Double.MAX_VALUE;
			sub_min[i] = Double.MAX_VALUE;
		}
		no_intersection_error_max = -Double.MAX_VALUE;
		no_intersection_error_min = Double.MAX_VALUE;

		for (i = 0; i < n; i++) {
			M[i].get_feature_and_weight();
			for (j = 0; j < sub.length - 1; j++) {

				tis = (sub[j] + sub[j + 1]) / 2;

				temp = a.Err_integration(M[i], tis, tmin);
				temp = temp - a.Err_integration(M[i], tmax, tis);

				if (sub_min[j] > temp)
					sub_min[j] = temp;
				if (sub_max[j] < temp)
					sub_max[j] = temp;

			}
			temp = a.Err_integration(M[i], tmax, tmin);
			if (no_intersection_error_min > temp)
				no_intersection_error_min = temp;
			if (no_intersection_error_max < temp)
				no_intersection_error_max = temp;

			M[i].f = null;
			M[i].w = null;
		}

	}

}

class Multiplie_region {
	LowerBound a = new LowerBound();

	region[] g = new region[a.vector.length];

}

class data_Q {
	double[][] sub_Q;
	double[] no_intersection_error;
	data[] Q;

}