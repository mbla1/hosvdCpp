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
#include <omp.h>

#define NUM_THREADS 4 // Number of threads to use

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
 * and appends it into a simple vector for storage
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

// Storage for the index of an element of the tensor

struct Index{
	int e;
	int q1;
	int q2;
};

Index conversion_index(int idx, int coord1Size, int coord2Size){
	
	Index index_temp;
	
	index_temp.e = idx / (coord1Size * coord2Size);
	index_temp.q1 = (idx % (coord1Size * coord2Size)) / coord2Size;
	index_temp.q2 = (idx % (coord1Size * coord2Size)) % coord2Size;

	return index_temp;
}

int main(){
	setprecision(12); //requires <iomanip>
	omp_set_num_threads(NUM_THREADS); // requires <omp.h>

	constexpr int coord1Size = 146;
	constexpr int coord2Size = 184;
	constexpr int nElec = 3;	
	constexpr int outputSize = 100; // Number of values outputted
	constexpr int time_limit = 1800; // AU
	constexpr int time_step = 50; // AU
	constexpr double tolerance = 0.99;
	string path{"resultats/"};

	ofstream sortieFichier;
	sortieFichier.open(path + "resultatsHosvd2.txt");

	if(!sortieFichier){
		cerr << "Wrong filename!\n";
		return -1;
	}

	sortieFichier << "temps\tvaleur\te\tq1\tq2\tsomme_carre\n";

	ofstream sortieU;
	sortieU.open(path + "matU2.txt");

	ofstream sortieV;
	sortieV.open(path + "matV2.txt");

	ofstream sortieW;
	sortieW.open(path + "matW2.txt");

	ofstream sortieNombre;
	sortieNombre.open(path + "nombre2.txt");

	#pragma omp parallel for ordered schedule(static, 1) // Parallel computing
	for(int time = 0; time < time_limit; time += time_step){ // Main loop
		vector<complex<double>> table;

		cout << "Loading a wavefunction" << endl;
		cout << "Time: " << time << "\n";
		loadInputTensor("waveFunction6O1/vector" + to_string(time) + ".txt", table); // q2,q1,e

		Eigen::MatrixXcd matA(coord2Size, nElec * coord1Size); 	
		Eigen::MatrixXcd matB(nElec, coord2Size * coord1Size); 	
		Eigen::MatrixXcd matC(coord1Size, coord2Size * nElec); 	

		cout << "Conversion of the format of the arrays" << endl;

		for(int i = 0; i < table.size(); i++){
			int e = i / (coord1Size * coord2Size);
			int q1 = (i % (coord1Size * coord2Size)) / coord2Size;
			int q2 = (i % (coord1Size * coord2Size)) % coord2Size;

			// Utilisation de permutations non cycliques

			matA(q2, e * coord1Size + q1) = table[i]; // q2,q1,e
			matB(e, q2 * coord1Size + q1) = table[i]; // e,q1,q2
			matC(q1, e * coord2Size + q2) = table[i]; // q1,q2,e
		}

		cout << "Computations of the SVD" << endl;

		/* Comparing with the Ramesh article:
		 * matA = U
		 * matB = V
		 * matC = W
		 */

		Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(matA, Eigen::ComputeThinU); // ThinSVD for only the u vectors
		cout << "First SVD done!\n";
		Eigen::JacobiSVD<Eigen::MatrixXcd> svd2(matB, Eigen::ComputeThinU);
		cout << "Second SVD done!\n";
		Eigen::JacobiSVD<Eigen::MatrixXcd> svd3(matC, Eigen::ComputeThinU);
		cout << "Third SVD done!\n";


		cout << "Matrix products in preparation" << endl;

		Eigen::MatrixXcd kronecker(nElec * coord1Size, coord2Size * coord1Size * coord1Size);
		produitKronecker(svd2.matrixU(), svd3.matrixU(), kronecker);

		Eigen::MatrixXcd resultat;
		resultat = svd1.matrixU().adjoint() * matA * kronecker;

		cout << "Outputting the results\n";

		vector<double> entrees;
		vector<int> index;
		vector<double> cumul;

		// stores the modulus of the value
		for(int i = 0; i < resultat.size(); i++){
			entrees.push_back(abs(resultat(i)));
			index.push_back(i);
		}

		//sort(entrees.begin(), entrees.end(), greater<double>());
		sort(index.begin(), index.end(), // Personnalized sort on indices
				[&](int i, int j){return entrees[i] > entrees[j];}); // lambda function


		cout << "Performing calculations\n";

		double somme = 0.;

		cout << "Cumul : \n";

		for(int i = 0; i < entrees.size(); ++i){ // Lecture dans l'ordre décroissant
			somme += entrees[index[i]] * entrees[index[i]];
			cumul.push_back(somme);
		}

		int compteur = 0;

		for(double d : cumul)
			if(d < tolerance){
				//cout << d << "\n";
				++compteur;
			}

		++compteur; // One more to get to 0.99

		cumul.clear();

		Index indices;

		#pragma omp ordered
		{
			for( int i = 0; i < outputSize; ++i){
				indices = conversion_index(index[i], coord1Size, coord2Size);
				//sortieFichier << time << "\t" << entrees[index[i]] << "\t" 
				//	<< indices.e << "\t" << indices.q1 << "\t"
				//	<< indices.q2 << "\t" << somme << "\n";
			}

			for(int i = 0; i < svd1.matrixU().size(); ++i)
				//sortieU << svd1.matrixU()(i) << "\t"; // column major by default in Eigen
			sortieU << "\n";

			for(int i = 0; i < svd2.matrixU().size(); ++i)
				//sortieV << svd2.matrixU()(i) << "\t"; // column major by default in Eigen
			sortieV << "\n";

			for(int i = 0; i < svd3.matrixU().size(); ++i)
				//sortieW << svd3.matrixU()(i) << "\t"; // column major by default in Eigen
			sortieW << "\n";

			sortieNombre << time << "\t" << compteur << "\n";
		}
	}
	sortieFichier.close();
	sortieU.close();
	sortieV.close();
	sortieW.close();
	sortieNombre.close();
	
	return 0;
}
