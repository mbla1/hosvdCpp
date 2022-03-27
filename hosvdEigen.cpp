#include <iostream>
#include <iomanip> // setprecision of cout
#include <fstream>
#include <string>
#include <vector>
#include <complex> 
#include <algorithm> // sort function
//#include <functional> // max functional
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
 * and stores it into a simple vector
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

// Kronecker product between two matrices

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

int main(){
	setprecision(12); //requires <iomanip>

	constexpr int gridSize = 511;
	constexpr int nElec = 7;	
	constexpr int initialTime = 700; // in au
	constexpr int finalTime = 8999;
	constexpr int timeSize = finalTime - initialTime + 1;
	constexpr int outputSize = 10000; // Number of values outputted
	string path{"resultats/resultatsHosvdAfterPulse.txt"};

	vector<complex<double>> table;

	cout << "Loading all the wavefunctions" << endl;
	for(int i = initialTime; i < (finalTime + 1); i++){
		loadInputTensor("wavefunction/vector" + to_string(i) + ".txt", table); //g,e,t
	}

	Eigen::MatrixXcd matA(gridSize, timeSize * nElec); 	
	Eigen::MatrixXcd matB(timeSize, gridSize * nElec); 	
	Eigen::MatrixXcd matC(nElec, gridSize * timeSize); 	

	cout << "Conversion of the format of the arrays" << endl;

	for(int i = 0; i < table.size(); i++){
		int t = i / (gridSize * nElec);
		int e = (i % (gridSize * nElec)) / gridSize;
		int g = (i % (gridSize * nElec)) % gridSize;

		//Utilisation de permutations non cycliques

		matA(g, t * nElec + e) = table[i]; // g,e,t
		matB(t, e * gridSize + g) = table[i]; // t,e,g
		matC(e, g * timeSize + t) = table[i]; // e,g,t
	}

	cout << "Computations of the SVD" << endl;

	Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(matA, Eigen::ComputeThinU); // ThinSVD for only the u vectors
	cout << "First SVD done!\n";
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd2(matB, Eigen::ComputeThinU);
	cout << "Second SVD done!\n";
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd3(matC, Eigen::ComputeThinU);
	cout << "Third SVD done!\n";


	cout << "Matrix products in preparation" << endl;

	Eigen::MatrixXcd kronecker(timeSize * nElec, gridSize * nElec * nElec);
	produitKronecker(svd2.matrixU(), svd3.matrixU(), kronecker);

	Eigen::MatrixXcd resultat;
	resultat = svd1.matrixU().adjoint() * matA * kronecker;

	cout << "Outputting the results\n";

	vector<double> entrees;
	vector<int> index;

	// stores the modulus of the value
	for(int i = 0; i < resultat.size(); i++){
		entrees.push_back(abs(resultat(i)));
		index.push_back(i);
	}

	//sort(entrees.begin(), entrees.end(), greater<double>());
	sort(index.begin(), index.end(), // Personnalized sort on indices
			[&](int i, int j){return entrees[i] > entrees[j];}); // lambda function

	ofstream sortieFichier;
	sortieFichier.open(path);

	if(!sortieFichier){
		cerr << "Wrong filename!\n";
		return -1;
	}

	sortieFichier << "valeur\tindex\tsomme_carre\n";

	double somme = 0.;

	for(double d : entrees){
		somme += d * d;
	}

	for(int i = 0; i < outputSize; i++){
		sortieFichier << entrees[index[i]] << "\t" << index[i] << "\t" << somme << "\n";
	}

	sortieFichier.close();
	
	return 0;
}
