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

	constexpr int gridSize = 511;
	constexpr int nElec = 7;	
	constexpr int initialTime = 0;
	constexpr int finalTime = 8999;
	constexpr int timeSize = finalTime - initialTime + 1;
	constexpr int outputSize = 10000;
	string path{"code/resultatsHosvd.txt"};

	vector<complex<double>> table;

	cout << "Loading all the wavefunctions" << endl;
	for(int i = initialTime; i < (finalTime + 1); i++){
		loadInputTensor("wavefunction/vector" + to_string(i) + ".txt", table); //g,e,t
	}

	Eigen::MatrixXcd matA(gridSize, timeSize * nElec); // g,e,t	
	Eigen::MatrixXcd matB(timeSize, gridSize * nElec); // g,e,t	
	Eigen::MatrixXcd matC(nElec, gridSize * timeSize); // g,e,t	

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

	Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(matA, Eigen::ComputeThinU);
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

	sortieFichier << "valeur\tindex\n";

	for(int i = 0; i < outputSize; i++){
		sortieFichier << entrees[index[i]] << "\t" << index[i] << "\n";
	}

	sortieFichier.close();
	
	double somme = 0.;

	for(double d : entrees){
		somme += d * d;
	}

	cout << "The sum of all values is: " << somme << "\n";

	return 0;
}
