package emd_package;

public class dualEMD_class {
	final int MAX_SIG_SIZE = 3000;
	int s = 0;
	data database[] = new data[8];
	double[] phi;
	double[] pai;
	double dualEMD_value;

	void dualEMD(signature_t s1, signature_t s2) {

		double temp = 0, temp2 = 0;
		int i, j, k;

		k = 0;
		algorithm1(s1, s2);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 1;

		algorithm1(s2, s1);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 2;

		algorithm2(s1, s2);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 3;

		algorithm2(s2, s1);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 4;

		algorithm3(s1, s2);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 5;

		algorithm3(s2, s1);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		k = 6;

		algorithm4(s1, s2);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;
		k = 7;

		algorithm4(s2, s1);

		database[k] = new data();
		database[k].phi = phi;
		database[k].pai = pai;
		database[k].dualEMD_value = dualEMD_value;

		for (temp = database[0].dualEMD_value, i = 1; i < database.length; i++) {
			temp2 = database[i].dualEMD_value;
			if (temp < temp2) {
				temp = temp2;
				s = i;

			}
		}

	}

	double Dist(feature_t F1, feature_t F2) {
		double dX = F1.x - F2.x, dY = F1.y - F2.y, dZ = F1.z - F2.z;
		return Math.pow(dX * dX + dY * dY + dZ * dZ, 0.5);
	}

	class data {
		double[] phi;
		double[] pai;

		double dualEMD_value;

	}

	double get_dualEMD_value() {
		return (database[s].dualEMD_value);
	}

	void algorithm1(signature_t s1, signature_t s2) {

		phi = new double[s1.n];
		pai = new double[s2.n];
		double temp = 0, temp2 = 0, total_weights = 1;
		int i, j, z;

		////////////
		for (i = 0; i < s1.n; i++) {
			for (z = 0, temp = Dist(s1.Features[i], s2.Features[0]) * (-1)
					* (s1.Weights[i] - s2.Weights[0]), j = 1; j < s2.n; j++) {
				temp2 = Dist(s1.Features[i], s2.Features[j]) * (-1) * (s1.Weights[i] - s2.Weights[j]);
				if (temp < temp2) {

					temp = temp2;
					z = j;
				}

			}

			phi[i] = Dist(s1.Features[i], s2.Features[z]);

		}

		/////////

		for (i = 0; i < s2.n; i++) {

			for (temp = Dist(s2.Features[i], s1.Features[0]) - phi[0], j = 1; j < s1.n; j++) {

				temp2 = Dist(s2.Features[i], s1.Features[j]) - phi[j];
				if (temp > temp2) {
					temp = temp2;
				}

			}

			pai[i] = temp;

		}
		////////////
		temp = 0;
		for (i = 0; i < s1.n; i++) {
			temp = temp + phi[i] * s1.Weights[i];

		}
		for (i = 0; i < s2.n; i++) {
			temp = temp + pai[i] * s2.Weights[i];

		}

		dualEMD_value = temp / (total_weights);

	}

	void algorithm2(signature_t s1, signature_t s2) {

		phi = new double[s1.n];
		pai = new double[s2.n];
		double temp = 0, temp2 = 0, total_weights = 1;
		int i, j, z;

		////////////
		for (i = 0; i < s1.n; i++) {
			for (z = 0, temp = Dist(s1.Features[i], s2.Features[0]) * (1)
					* (s1.Weights[i] - s2.Weights[0]), j = 1; j < s2.n; j++) {
				temp2 = Dist(s1.Features[i], s2.Features[j]) * (1) * (s1.Weights[i] - s2.Weights[j]);
				if (temp < temp2) {

					temp = temp2;
					z = j;
				}

			}

			phi[i] = Dist(s1.Features[i], s2.Features[z]);

		}

		/////////

		for (i = 0; i < s2.n; i++) {

			for (temp = Dist(s2.Features[i], s1.Features[0]) - phi[0], j = 1; j < s1.n; j++) {

				temp2 = Dist(s2.Features[i], s1.Features[j]) - phi[j];
				if (temp > temp2) {
					temp = temp2;
				}

			}

			pai[i] = temp;

		}
		////////////
		temp = 0;
		for (i = 0; i < s1.n; i++) {
			temp = temp + phi[i] * s1.Weights[i];

		}
		for (i = 0; i < s2.n; i++) {
			temp = temp + pai[i] * s2.Weights[i];

		}

		dualEMD_value = temp / (total_weights);

	}

	void algorithm3(signature_t s1, signature_t s2) {

		phi = new double[s1.n];
		pai = new double[s2.n];
		double temp = 0, temp2 = 0, total_weights = 1;
		int i, j, z;

		////////////
		for (i = 0; i < s1.n; i++) {
			for (z = 0, temp = Dist(s1.Features[i], s2.Features[0]) * (1)
					* Math.abs((s1.Weights[i] - s2.Weights[0])), j = 1; j < s2.n; j++) {
				temp2 = Dist(s1.Features[i], s2.Features[j]) * (1) * Math.abs((s1.Weights[i] - s2.Weights[j]));
				if (temp < temp2) {

					temp = temp2;
					z = j;
				}

			}

			phi[i] = Dist(s1.Features[i], s2.Features[z]);

		}

		/////////

		for (i = 0; i < s2.n; i++) {

			for (temp = Dist(s2.Features[i], s1.Features[0]) - phi[0], j = 1; j < s1.n; j++) {

				temp2 = Dist(s2.Features[i], s1.Features[j]) - phi[j];
				if (temp > temp2) {
					temp = temp2;
				}

			}

			pai[i] = temp;

		}
		////////////
		temp = 0;
		for (i = 0; i < s1.n; i++) {
			temp = temp + phi[i] * s1.Weights[i];

		}
		for (i = 0; i < s2.n; i++) {
			temp = temp + pai[i] * s2.Weights[i];

		}

		dualEMD_value = temp / (total_weights);

	}

	void algorithm4(signature_t s1, signature_t s2) {

		phi = new double[s1.n];
		pai = new double[s2.n];
		double temp = 0, temp2 = 0, total_weights = 1;
		int i, j, z;

		////////////
		for (i = 0; i < s1.n; i++) {
			for (z = 0, temp = Dist(s1.Features[i], s2.Features[0]) * (1)
					* (s1.Weights[i] - s2.Weights[0]), j = 1; j < s2.n; j++) {
				if (temp < 0)
					temp = 0;
				temp2 = Dist(s1.Features[i], s2.Features[j]) * (1) * (s1.Weights[i] - s2.Weights[j]);
				if (temp2 < 0)
					temp2 = 0;
				if (temp < temp2) {

					temp = temp2;
					z = j;
				}

			}

			phi[i] = Dist(s1.Features[i], s2.Features[z]);

		}

		/////////

		for (i = 0; i < s2.n; i++) {

			for (temp = Dist(s2.Features[i], s1.Features[0]) - phi[0], j = 1; j < s1.n; j++) {

				temp2 = Dist(s2.Features[i], s1.Features[j]) - phi[j];
				if (temp > temp2) {
					temp = temp2;
				}

			}

			pai[i] = temp;

		}
		////////////
		temp = 0;
		for (i = 0; i < s1.n; i++) {
			temp = temp + phi[i] * s1.Weights[i];

		}
		for (i = 0; i < s2.n; i++) {
			temp = temp + pai[i] * s2.Weights[i];

		}

		dualEMD_value = temp / (total_weights);

	}

}
