#include "Eigen/Dense"
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES 
#include <cmath>
using namespace std;
using Vektor = Eigen::Vector3d;

const double e = 2.71828182845904523536;

Vektor Krive_u_Dekartove(Vektor a)
{
	double x = a.x();
	double y = (a.y() - x / 2) * (2 / sqrt(3));
	double z = (3 / sqrt(6))*(a.z() - x / 2 - y * sqrt(3) / 6);
	return { x,y,z };
}

double SkalarniProizvod(Vektor u, Vektor v)
{
	Vektor v1 = Krive_u_Dekartove(u);
	Vektor v2 = Krive_u_Dekartove(v);
	return v1.dot(v2);
}


double Intenzitet(Vektor v)
{
	Vektor a = Krive_u_Dekartove(v);
	return a.norm();
}

double Rastojanje(Vektor u, Vektor v)
{
	return Intenzitet(u - v);
}

random_device r;
default_random_engine rng(r());

int NextInt(int low, int high) {
	auto dist = uniform_int_distribution<int>(low, high - 1);
	return dist(rng);
}

bool FlipCoin(double probTrue) {
	auto dist = bernoulli_distribution(probTrue);
	return dist(rng);
}

struct Pirohlor
{
	struct Celija //unit cell sa vectorom spinova
	{
		vector<Vektor> spin;
		Celija(double s0, double s1, double s2, double s3) //prosledjujemo znak spina za 0. 1. 2. 3., plus ka unutra
		{

			spin.push_back(s0 * Vektor(1.0, 1.0, 1.0));
			spin.push_back(s1 * Vektor(-1., 1, 1));
			spin.push_back(s2 * Vektor(-1, -1, 1));
			spin.push_back(s3 * Vektor(1, 1, -1));
		}

		Celija() = default;

		void Okreni(int broj) //prosledjujemo indeks 0/1/2/3 i okrecemo ga
		{
			spin[broj] = -1 * spin[broj];
		}
	};

	const static int N = 10;
	double A;
	int n = 0;
	double hamiltonijan = 0;
	Celija resetka[N][N][N];
	double J;
	double kb;

	Vektor ukupna_mag = Vektor(); //zbir sin(om*t)*M

	Vektor zbir_spinova = Vektor(0, 0, 0); //treba nam kod dela kad racunamo polje, pa cemo tda apdejtujemo uvek

	double hamiltonijanNN_pocetni(double J) //nirest nejbr hamiltonijan (taj deo)
	{
		double nnHam = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					for (int o = 0; o < 4; o++)
					{
						vector<Vektor> sus = NadjiSusede(i, j, k, o);
						for (Vektor v : sus)
						{
							nnHam += SkalarniProizvod(resetka[i][j][k].spin[o], v);
						}
					}

				}
			}
		}
		return nnHam * J / 6.0; //To je ono jer racunam po dva puta (/2) i sto mi spinovi nisu sqrt(3)/3 nego 1
	}

	double hamiltonijanDD_pocetni(double A, double D) //dipol dipol deo hamiltoijana
	{
		double ham_DD = 0;

		for (int i1 = 0; i1 < n; i1++)
		{
			for (int j1 = 0; j1 < n; j1++)
			{
				for (int k1 = 0; k1 < n; k1++)
				{
					for (int i2 = i1 - 1; i2 < i1 + 2; i2++)
					{
						for (int j2 = j1 - 1; j2 < j1 + 2; j2++)
						{
							for (int k2 = k1 - 1; k2 < k1 + 2; k2++)
							{
								if (i2 > 0 && i2 < n && j2 > 0 && j2 < n && k2 > 0 && k2 < n)
								{
									for (int br1 = 0; br1 < 4; br1++)
									{
										for (int br2 = 0; br2 < 4; br2++)
										{
											Vektor v1 = Vektor(2 * A * i1, 2 * A * j1, 2 * A * k1);
											Vektor v2 = Vektor(2 * A * i2, 2 * A * j2, 2 * A * k2);
											if (br1 == 1) v1 += Vektor(A, 0, 0);
											else if (br1 == 2) v1 += Vektor(0, A, 0);
											else if (br1 == 3) v1 += Vektor(0, 0, A);
											if (br2 == 1) v2 += Vektor(A, 0, 0);
											else if (br2 == 2) v2 += Vektor(0, A, 0);
											else if (br2 == 3) v2 += Vektor(0, 0, A);

											double d = Rastojanje(v1, v2);
											Vektor razlika = v2 - v1;

											if (d > 0)
											{
												Vektor s1 = resetka[i1][j1][k1].spin[br1];
												Vektor s2 = resetka[i2][j2][k2].spin[br2];
												ham_DD = ham_DD + SkalarniProizvod(s1, s2) / pow(d, 3) - (3 * SkalarniProizvod(s1, razlika) * SkalarniProizvod(s2, razlika)) / pow(d, 5);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		ham_DD = (D * A) * ham_DD;
		return 0.5 * ham_DD;
	}

	Pirohlor(int n, double J, double kb, double A, double D)
	{
		this->J = J;
		this->n = n;
		this->kb = kb;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int z = 0; z < n; z++)
				{
					
					int poms0 = FlipCoin(0.5);
					if (poms0 == 0) poms0 = -1;
					int poms1 = FlipCoin(0.5);
					if (poms1 == 0) poms1 = -1;
					int poms2 = FlipCoin(0.5);
					if (poms2 == 0) poms2 = -1;
					int poms3 = FlipCoin(0.5);
					if (poms3 == 0) poms3 = -1;
					this->resetka[i][j][z] = Celija(poms0, poms1, poms2, poms3); //generisemo resetku

					zbir_spinova += resetka[i][j][z].spin[0];
					zbir_spinova += resetka[i][j][z].spin[1];
					zbir_spinova += resetka[i][j][z].spin[2];
					zbir_spinova += resetka[i][j][z].spin[3];
				}
			}
		}

		this->A = A;
		//Console.WriteLine(hamiltonijanNN_pocetni(J));
		//Console.WriteLine(hamiltonijanDD_pocetni(A, D));
		this->hamiltonijan = hamiltonijanNN_pocetni(J) + hamiltonijanDD_pocetni(A, D); //bez polja jer je od polja na pocetku nula
	}

	vector<Vektor> NadjiSusede(int x, int y, int z, int broj) //zadajemo poziciju i broj u celiji i trazimo susede
	{
		vector<Vektor> susedi;

		if (broj == 0)
		{
			susedi.push_back(resetka[x][y][z].spin[1]);
			susedi.push_back(resetka[x][y][z].spin[2]);
			susedi.push_back(resetka[x][y][z].spin[3]);
			if (x > 0) susedi.push_back(resetka[x - 1][y][z].spin[1]);
			if (y > 0) susedi.push_back(resetka[x][y - 1][z].spin[2]);
			if (z > 0) susedi.push_back(resetka[x][y][z - 1].spin[3]);
		}
		else if (broj == 1)
		{
			susedi.push_back(resetka[x][y][z].spin[0]);
			susedi.push_back(resetka[x][y][z].spin[2]);
			susedi.push_back(resetka[x][y][z].spin[3]);
			if (x < n - 1) susedi.push_back(resetka[x + 1][y][z].spin[0]);
			if (x < n - 1 && y > 0) susedi.push_back(resetka[x + 1][y - 1][z].spin[2]);
			if (x < n - 1 && z > 0) susedi.push_back(resetka[x + 1][y][z - 1].spin[3]);
		}
		else if (broj == 2)
		{
			susedi.push_back(resetka[x][y][z].spin[0]);
			susedi.push_back(resetka[x][y][z].spin[1]);
			susedi.push_back(resetka[x][y][z].spin[3]);
			if (y < n - 1) susedi.push_back(resetka[x][y + 1][z].spin[0]);
			if (x > 0 && y < n - 1) susedi.push_back(resetka[x - 1][y + 1][z].spin[1]);
			if (y < n - 1 && z > 0) susedi.push_back(resetka[x][y + 1][z - 1].spin[3]);
		}
		else if (broj == 3)
		{
			susedi.push_back(resetka[x][y][z].spin[0]);
			susedi.push_back(resetka[x][y][z].spin[1]);
			susedi.push_back(resetka[x][y][z].spin[2]);
			if (z < n - 1) susedi.push_back(resetka[x][y][z + 1].spin[0]);
			if (x > 0 && z < n - 1) susedi.push_back(resetka[x - 1][y][z + 1].spin[1]);
			if (y > 0 && z < n - 1) susedi.push_back(resetka[x][y - 1][z + 1].spin[2]);
		}

		return susedi;
	}

	Vektor Magnetizacija()
	{
		return zbir_spinova / (double)(n*n*n * 4);
	}

	double DoprinosHamuNN(int x, int y, int z, int broj) //koliko jedan onaj doprinosi HamuNN
	{
		double dop = 0;
		vector<Vektor> sus = NadjiSusede(x, y, z, broj);
		for (Vektor v : sus)
		{
			dop += SkalarniProizvod(v, resetka[x][y][z].spin[broj]);
		}
		return dop;
	}

	double DoprinosHamuDD(int x, int y, int z, int broj)
	{
		double doprinos_DD = 0;

		for (int i = x - 1; i < x + 2; i++)
		{
			for (int j = y - 1; j < y + 2; j++)
			{
				for (int k = z - 1; k < z + 2; k++)
				{
					for (int br0 = 0; br0 < 4; br0++)
					{
						if (i > 0 && j > 0 && k > 0 && i < n && j < n && k < n)
						{
							Vektor v1 = Vektor(2 * A * x, 2 * A * y, 2 * A * z);
							Vektor v2 = Vektor(2 * A * i, 2 * A * j, 2 * A * k);
							if (broj == 1) v1 += Vektor(A, 0, 0);
							else if (broj == 2) v1 += Vektor(0, A, 0);
							else if (broj == 3) v1 += Vektor(0, 0, A);
							if (br0 == 1) v2 += Vektor(A, 0, 0);
							else if (br0 == 2) v2 += Vektor(0, A, 0);
							else if (br0 == 3) v2 += Vektor(0, 0, A);

							double d = Rastojanje(v1, v2);
							Vektor razlika = v2 - v1;

							if (d > 0)
							{
								Vektor s1 = resetka[x][y][z].spin[broj];
								Vektor s2 = resetka[i][j][k].spin[br0];
								doprinos_DD = doprinos_DD + SkalarniProizvod(s1, s2) / pow(d, 3) - (3 * SkalarniProizvod(s1, razlika) * SkalarniProizvod(s2, razlika)) / pow(d, 5);
							}
						}
					}
				}
			}
		}
		return doprinos_DD;
	}

	void KorakMetropolis(double T, double t /*vreme*/, double dt, Vektor H0, double omega, double mief) //sibamo jedan korak u metropolisu i menjamo razliku hamiltonijana pocetnog i novog
	{
		static double prev_intensity = 0, cur_intensity = 0;
		
		int x = NextInt(0, n);
		int y = NextInt(0, n);
		int z = NextInt(0, n);
		int broj = NextInt(0, 4);

		double hamiltonijan0 = hamiltonijan;

		double doprinosHamuNN = DoprinosHamuNN(x, y, z, broj);
		double doprinosHamuDD = DoprinosHamuDD(x, y, z, broj);
		Vektor zbir_spinova_potencijalno = zbir_spinova + ((-2.0) * resetka[x][y][z].spin[broj]);
		double hamiltonijan_potencijalno = hamiltonijan 
			- 2 * (doprinosHamuNN + doprinosHamuDD) 
			- mief * sin(omega * t) * (SkalarniProizvod(H0, zbir_spinova)) 
			+ mief * sin(omega * (t + dt)) * (SkalarniProizvod(H0, zbir_spinova_potencijalno));

		double dE = hamiltonijan_potencijalno - hamiltonijan;

		if (dE < 0)
		{
			resetka[x][y][z].spin[broj] *= -1;
			hamiltonijan = hamiltonijan_potencijalno;
			zbir_spinova = zbir_spinova_potencijalno;
		}
		else
		{
			double verov = exp((-1 * dE) / (kb * T));
			//cout << verov << endl;
			//Console.WriteLine(verov);
			//Console.ReadLine();
			if (FlipCoin(verov))
			{
				resetka[x][y][z].spin[broj] *= -1;;
				hamiltonijan = hamiltonijan_potencijalno;
				zbir_spinova = zbir_spinova_potencijalno;
			}
		}
		prev_intensity = cur_intensity;
		cur_intensity = Intenzitet(zbir_spinova);
		/*if (prev_intensity > cur_intensity)
		{
			cout << cur_intensity << " " << prev_intensity - cur_intensity << endl;
		}*/
		//cout << cur_intensity << endl;

		ukupna_mag += (sin(omega * (t + dt)) * Magnetizacija());
	}
};


int main()
{
	double J = 3.41;
	double kb = 1;
	double T = 0.7;
	int n = 10;
	double dt = 0.0001;
	int br_koraka = 50000;
	Vektor H0 = Vektor(1, 1, 1);
	double mief = -1;
	double omega = 2;
	double A = 7.191;
	double D = 1.32;

	ofstream sws("sus.txt");
	ofstream swf("frek.txt");

	for (int i1 = 0; i1 < 1; i1++)
	{

		Pirohlor p(n, J, kb, A, D);

		if (i1 % 10 == 0)
		{
			cout << (i1) << endl;
		}

		string sh = ("ham" + to_string(i1)) + ".txt";
		string sm = ("mag" + to_string(i1)) + ".txt";
		ofstream sw1(sh);
		ofstream sw2(sm);


		for (int i = 0; i < br_koraka; i++)
		{
			sw1 << (p.hamiltonijan) << endl;
			sw2 << (Intenzitet(p.Magnetizacija())) << endl;
			p.KorakMetropolis(T, dt * i, dt, H0, omega, mief);
		}

		sws << (Intenzitet(p.ukupna_mag) / (Intenzitet(H0)*br_koraka)) << endl;

		/*        sw1.Close();
		sw2.Close();
		*/
		swf << (omega) << endl;
		omega += 0.10;
	}
	return 0;
}
