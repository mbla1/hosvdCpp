#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <complex> 
#include <algorithm>
#include <functional>
#define EIGEN_USE_MKL_ALL
#include <Eigen>

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

void head(double entree[]){ //Affiche le début d'une matrice
	for(int i = 0; i < 5; i++){	
		cout << entree[i] << endl;
	}
}

void head(complex<double> entree[], int positionInitiale = 0){ //Affiche le début d'une matrice
	for(int i = 0; i < 5; i++){	
		cout << entree[positionInitiale + i] << endl;
	}
}

void produitKronecker(const Eigen::MatrixXcd &entree1, const Eigen::MatrixXcd &entree2, Eigen::MatrixXcd &sortie){
	// You need to pass it by reference idiot!
	int lignes = entree1.rows();
	int colonnes = entree1.cols();
	for(int i = 0; i < lignes; i++){
		for(int j = 0; j < colonnes; j++){
			sortie.block(i * entree2.rows(), j * entree2.cols(), 
					entree2.rows(), entree2.cols()) = entree1(i, j) * entree2; // Produit bloc par bloc
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

	Eigen::MatrixXcd matA(511, 9000 * 7); // g,e,t	
	Eigen::MatrixXcd matB(9000, 511 * 7); // g,e,t	
	Eigen::MatrixXcd matC(7, 511 * 9000); // g,e,t	

	cout << "Conversion of the format of the arrays" << endl;

	for(int i = 0; i < table.size(); i++){
		int t = i / 3577;
		int e = (i % 3577) / 511;
		int g = (i % 3577) % 511;

		//Utilisation de permutations non cycliques

		matA(g, t * 7 + e) = table[i]; // g,e,t
		matB(t, e * 511 + g) = table[i]; // t,e,g
		matC(e, g * 9000 + t) = table[i]; // e,g,t
	}

	cout << "Test layout matB:" << endl;

	cout << matA.topLeftCorner(5, 5) << endl << endl;
	cout << matB.topLeftCorner(5, 5) << endl << endl;
	cout << matC.topLeftCorner(5, 5) << endl << endl;


	cout << "Computations of the SVD" << endl;

	Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(matA, Eigen::ComputeThinU);
	cout << "First SVD done!\n";
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd2(matB, Eigen::ComputeThinU);
	cout << "Second SVD done!\n";
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd3(matC, Eigen::ComputeThinU);
	cout << "Third SVD done!\n";


	cout << "Matrix products in preparation" << endl;

	Eigen::MatrixXcd kronecker(9000 * 7, 3577 * 7);
	produitKronecker(svd2.matrixU(), svd3.matrixU(), kronecker);

	cout << kronecker.topLeftCorner(5, 5) << "\n\n";

	Eigen::MatrixXcd resultat;
	resultat = svd1.matrixU().adjoint() * matA * kronecker;

	cout << resultat.topLeftCorner(5, 5) << "\n\n";

	vector<double> entrees;

	for(int i = 0; i < resultat.size(); i++){
		entrees.push_back(abs(resultat(i)));
	}

	sort(entrees.begin(), entrees.end(), greater<double>());

	for(int i = 0; i < 5; i++){
		cout << entrees[i] << "\n";
	}

	return 0;
}
