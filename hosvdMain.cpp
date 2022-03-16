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

/* Loads the time, electronic states and grid tensor
 * and store it into a simple vector
 */
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
			entreeVec.clear(); //Otherwise I only take into account the first entry
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

void head(complex<double> entree[], int positionInitiale = 0){ //Affiche le début de la matrice
	for(int i = 0; i < 5; i++){	
		cout << entree[positionInitiale + i] << endl;
	}
}

void svdPartie(int m, int n, complex<double> entree[], complex<double> matU[]){
	double diag[m];
	complex<double> matV[m * n];
	double superb[m - 1];

	LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'N', m, n, entree, m, diag, matU, m, matV, m, superb);
}

void produitKronecker(complex<double> entree1[], complex<double> entree2[], complex<double> sortie[],
		int dim1, int dim2){
	for(int i = 0; i < dim1; i++){
		for(int j = 0; j < dim2; j++){
			sortie[i * dim2 + j] = entree1[i] * entree2[j];
		}
	}
}

void conjugueComplexe(complex<double> entree[], complex<double> sortie[], int nrow, int ncol){
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			sortie[j * nrow + i] = conj(entree[i * nrow + j]);
		}
	}
}

void getModulus(complex<double> entree[], double sortie[], int dim){
	for(int i = 0; i < dim; i++){
		sortie[i] = abs(entree[i]);
	}
}

int main(){
	setprecision(12); //requires <iomanip>
	
	vector<complex<double>> table;

	cout << "Loading all the wavefunctions" << endl;
	for(int i = 0; i < 9000; i++){
		loadInputTensor("wavefunction/vector" + to_string(i) + ".txt", table); //g,e,t
	}

	cout << table.size() << endl;
	
	for(int i = 0; i < 5; i++){
		cout << table[i] << endl;
	}

	complex<double> matA[table.size()]; //t * e_max * g_max + e * g_max + g
	conversionVecteurArray(table, matA);

	//Other layouts
	complex<double> matB[table.size()]; //e * g_max * t_max + g * t_max + t
	complex<double> matC[table.size()]; //g * t_max * e_max + t * e_max + e

	cout << "Conversion of the format of the arrays" << endl;

	for(int i = 0; i < table.size(); i++){
		int t = i / 3577;
		int e = (i % 3577) / 511;
		int g = (i % 3577) % 511;

		//Utilisation de permutations non cycliques

		matB[g * 7 * 9000 + e * 9000 + t] = table[i]; //t,e,g
		matC[t * 511 * 7 + g * 7 + e] = table[i]; //e,g,t
	}

	cout << "Test layout matB:" << endl;

	//head(matA);
	//cout << endl;
	head(matB);
	cout << endl;
	head(matB, 9000);
	cout << endl;
	head(matC);
	cout << endl;
	head(matC, 7);

	cout << "Computations of the SVD" << endl;

	complex<double> decomp1[511 * 511];
	svdPartie(511, 9000 * 7, matA, decomp1);
	head(decomp1, 511);
	cout << endl;

	complex<double> decomp2[9000 * 3577];
	svdPartie(9000, 511 * 7, matB, decomp2); //all the signs seem inverted...
	head(decomp2, 9000);
	cout << endl;

	complex<double> decomp3[7 * 7];
	svdPartie(7, 511 * 9000, matC, decomp3);
	//head(decomp3);
	cout << endl;

	for(int i = 0; i < 7; i++){
		for(int j = 0; j < 7; j++){
			cout << decomp3[i * 7 + j];
		}
		cout << endl;
	}

	cout << "Matrix products in preparation" << endl;

	int dimKron = 9000 * 3577 * 7 * 7;
	complex<double> *kronecker = (complex<double>*) malloc(dimKron * sizeof(complex<double>));//Dynamic allocation due to the size of the array
	produitKronecker(decomp2, decomp3, kronecker, 9000 * 3577, 7 * 7); 

	cout << "Kronecker product: " << endl;
	head(kronecker, 63000);

	complex<double> hermit[511 * 511];
	conjugueComplexe(decomp1, hermit, 511, 511);

	cout << "Hermitian matrix:" << endl;
	head(hermit);

	//And now the matrix products
	complex<double> premierResultat[511 * 63000];
	complex<double> alpha = 1.;
	complex<double> beta = 0.;
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 511, 63000, 511, &alpha, hermit, 511, matA, 511, &beta, premierResultat, 511);

	complex<double> *deuxiemeResultat = (complex<double>*) malloc(511 * 25039 * sizeof(complex<double>));
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 511, 25039, 63000, &alpha, 
			premierResultat, 511, kronecker, 63000, &beta, deuxiemeResultat, 511);

	cout << "Non-sorted result:" << endl;
	head(deuxiemeResultat);
	
	cout << "Sorting of the results" << endl;

	double resultats[511 * 25039];
	getModulus(deuxiemeResultat, resultats, 511 * 25039);
	sort(resultats, resultats + 511 * 25039, greater<double>());

	//cout << endl;
	cout << "Here are the results" << endl;
	head(resultats);

	//Now we compute the modulus and sort the list
	
	//sort(diag, diag + 3577, greater<double>());
	//head(diag); //Already sorted

	return 0;
}
