#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <complex> 
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

int main(){
	setprecision(12); //requires <iomanip>
	
	vector<complex<double>> table;

	//Loading of all the wavefunctions
	for(int i = 0; i < 9000; i++){
		loadInputTensor("wavefunction/vector" + to_string(i) + ".txt", table);
	}
	cout << table.size() << endl;

	complex<double> matA[table.size()];
	conversionVecteurArray(table, matA);

	double diag[];
	double offDiag[];
	complex<double> taup[];
	complex<double> tauq[];

	LAPACKE_zgebrd(LAPACK_ROW_MAJOR, 3577, 9000, matA, 9000, diag, offDiag,taup, tauq);

	cout << matA[0] << endl;

	return 1;
}
