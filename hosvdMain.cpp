#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <complex> 
#include <algorithm>
#include <functional>
#define MKL_Complex16 std::complex<double> //Definition to use the proper options
#include "mkl.h" /* To compile with the option icc -mkl 
				  * and with the intel module loaded 
				  */

using namespace std;

void decoupe(const string& chaine, char sep, vector<double>& sortie){
	/* Takes as input a string with spaces inside
	 * and transform it as a vector of numbers.
	 * This is especially useful to read delimited data
	 */
	int i = 0;
	int j = chaine.find(sep);

	while(j >= 0){
		sortie.push_back(stod(chaine.substr(i, j - i)));
		i = ++j;
		j = chaine.find(sep, j);
	}
	if(j < 0){
		sortie.push_back(stod(chaine.substr(i, chaine.length()))); 
	}
}

void loadInputTensor(string path, vector<complex<double>>& tableau){ //Vectors are not given as reference by default
	ifstream fichier;
	fichier.open(path);

	if(!fichier){
		cerr << "Please indicate a name for the files" << endl;
	} else {
		string entree;
		vector<double> entreeVec;
		complex<double> nombre;

		while(getline(fichier, entree)){
			decoupe(entree, ' ', entreeVec);
			nombre = complex<double>(entreeVec[0], entreeVec[1]); //Not elegant, but it works

			tableau.push_back(nombre);
		}
	}
	fichier.close();
}

void conversionVecteurArray(vector<complex<double>>& entree, complex<double> sortie[]){
	for(int i = 0; i < entree.size(); i++){
		sortie[i] = entree[i];
	}
}//Convert a vector to an array of the same format, to use with mkl

void head(double entree[]){ //Affiche le début de la matrice
	for(int i = 0; i < 5; i++){	
		cout << entree[i] << endl;
	}
}

void head(complex<double> entree[]){ //Affiche le début de la matrice
	for(int i = 0; i < 5; i++){	
		cout << entree[i] << endl;
	}
}

void svdPartie(int m, int n, complex<double> entree[], complex<double> matU[]){
	double diag[m];
	complex<double> matV[m * n];
	double superb[m - 1];

	LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'N', m, n, entree, m, diag, matU, m, matV, m, superb);
}

void tuckerProduct(complex<double> tensor[], complex<double> matrix[], complex<double> sortie[],
		int index, int dim[]){
	limit = dim[index];
	for(int i = 0; i < limite; i++){
			
	}
}

int main(){
	setprecision(12); //requires <iomanip>
	
	vector<complex<double>> table;

	//Loading of all the wavefunctions
	for(int i = 0; i < 9000; i++){
		loadInputTensor("wavefunction/vector" + to_string(i) + ".txt", table); //g,e,t
	}
	cout << table.size() << endl;

	complex<double> matA[table.size()]; //t * e_max * g_max + e * g_max + g
	conversionVecteurArray(table, matA);

	//Other layouts
	complex<double> matB[table.size()]; //e * g_max * t_max + g * t_max + t
	complex<double> matC[table.size()]; //g * t_max * e_max + t * e_max + e

	for(int i = 0; i < table.size(); i++){
		int t = i / 3577;
		int e = (i % 3577) / 511;
		int g = (i % 3577) % 511;

		matB[e * 511 * 9000 + g * 9000 + t] = table[i];
		matC[g * 9000 * 7 + t * 7 + e] = table[i];
	}

	complex<double> decomp1[511 * 511];
	svdPartie(511, 9000 * 7, matA, decomp1);
	head(decomp1);
	cout << endl;

	complex<double> decomp2[9000 * 3577];
	svdPartie(9000, 511 * 7, matB, decomp2);
	head(decomp2);
	cout << endl;

	complex<double> decomp3[7 * 7];
	svdPartie(7, 511 * 9000, matC, decomp3);
	head(decomp3);
	cout << endl;

	//sort(diag, diag + 3577, greater<double>());
	//head(diag); //Already sorted

	return 1;
}
