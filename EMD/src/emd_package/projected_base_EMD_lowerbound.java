package emd_package;

import java.io.IOException;

public class projected_base_EMD_lowerbound {
	double[][] vector = { { 1, 0, 0 } };

	double PASUM(signature_t s1, signature_t s2) throws CloneNotSupportedException {
		int i, j, k, p;
		signature_t s1_temp;
		signature_t s2_temp;
		double tmax, tmin;
		double s1_tmax, s1_tmin;
		double s2_tmax, s2_tmin;
		double value = 0;
		for (p = 0; p < vector.length; p++) {

			s1_temp = (signature_t) s1.clone();
			s2_temp = (signature_t) s2.clone();
			tmax = -Double.MAX_VALUE;
			tmin = Double.MAX_VALUE;
			s1_tmax = -Double.MAX_VALUE;
			s1_tmin = Double.MAX_VALUE;
			s2_tmax = -Double.MAX_VALUE;
			s2_tmin = Double.MAX_VALUE;

			for (i = 0; i < s1_temp.n; i++) {

				s1_temp.Features[i].x = s1_temp.Features[i].x * vector[p][0] + s1_temp.Features[i].y * vector[p][1]
						+ s1_temp.Features[i].z * vector[p][2];
				s1_temp.Features[i].y = 0;
				s1_temp.Features[i].z = 0;
				if (tmax < s1_temp.Features[i].x)
					tmax = s1_temp.Features[i].x;
				if (s1_tmax < s1_temp.Features[i].x)
					s1_tmax = s1_temp.Features[i].x;
				if (tmin > s1_temp.Features[i].x)
					tmin = s1_temp.Features[i].x;
				if (s1_tmin > s1_temp.Features[i].x)
					s1_tmin = s1_temp.Features[i].x;

				for (j = 0; j < i; j++) {
					if (s1_temp.Features[i].x == s1_temp.Features[j].x) {

						s1_temp.Weights[j] = s1_temp.Weights[j] + s1_temp.Weights[i];
						for (k = i; k < s1_temp.n; k++) {
							s1_temp.Features[k] = s1_temp.Features[k + 1];
							s1_temp.Weights[k] = s1_temp.Weights[k + 1];

						}

						s1_temp.n--;
						i--;
					}
				}
			}

			double temp = 0;
			for (i = 0; i < s1_temp.n; i++) {
				temp = temp + s1_temp.Weights[i];
			}

			for (i = 0; i < s2_temp.n; i++) {

				s2_temp.Features[i].x = s2_temp.Features[i].x * vector[p][0] + s2_temp.Features[i].y * vector[p][1]
						+ s2_temp.Features[i].z * vector[p][2];
				s2_temp.Features[i].y = 0;
				s2_temp.Features[i].z = 0;
				if (tmax < s2_temp.Features[i].x)
					tmax = s2_temp.Features[i].x;
				if (s2_tmax < s2_temp.Features[i].x)
					s2_tmax = s2_temp.Features[i].x;
				if (tmin > s2_temp.Features[i].x)
					tmin = s2_temp.Features[i].x;
				if (s2_tmin > s2_temp.Features[i].x)
					s2_tmin = s2_temp.Features[i].x;

				for (j = 0; j < i; j++) {
					if (s2_temp.Features[i].x == s2_temp.Features[j].x) {

						s2_temp.Weights[j] = s2_temp.Weights[j] + s2_temp.Weights[i];
						for (k = i; k < s2_temp.n; k++) {
							s2_temp.Features[k] = s2_temp.Features[k + 1];
							s2_temp.Weights[k] = s2_temp.Weights[k + 1];

						}

						s2_temp.n--;
						i--;
					}
				}

			}

			value = value
					+ Math.abs(CDF(tmin, s1_tmax, s1_temp) - CDF(tmin, s2_tmax, s2_temp)) * ((Math.ceil(tmin) - tmin));
			for (i = (int) Math.ceil(tmin); i < (int) Math.floor(tmax); i++) {

				value = value + Math.abs(CDF(i, s1_tmax, s1_temp) - CDF(i, s2_tmax, s2_temp));
			}

			value = value + Math.abs(CDF(tmax, s1_tmax, s1_temp) - CDF(tmax, s2_tmax, s2_temp))
					* (((tmax - Math.floor(tmax))));

		}

		return value / Math.pow(vector.length, 0.5);
	}

	double CDF(double n, double tmax, signature_t s) {
		int i;
		double value = 0;

		feature_t[] f = s.Features;
		double[] w = s.Weights;

		if (n < tmax) {
			for (i = 0; i < s.n; i++) {

				if (f[i].x <= n)
					value = value + w[i];

			}

		} else
			return 1;
		return (value);
	}

	double get_second_min(double min, signature_t s1, signature_t s2) {
		double second = 0;
		int n = s1.n + s2.n;
		double[] array = new double[n];

		int i;
		for (i = 0; i < s1.n; i++) {
			array[i] = s1.Features[i].x;

		}
		for (i = 0; i < s2.n; i++) {
			array[i + s1.n] = s2.Features[i].x;

		}

		for (i = 0; i < n; i++) {
			second = array[i];
			if (second > min)
				break;

		}
		if (i == n)
			second = min;
		for (i = 0; i < n; i++) {
			if (Boolean.valueOf(array[i] < second && array[i] > min)) {
				second = array[i];
			}
		}

		return second;
	}

}
