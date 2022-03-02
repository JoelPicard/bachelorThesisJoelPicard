#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include <string>
#include <gsl/gsl_rng.h>
using namespace std;

const double T_c = 1 / (2 * asinh(pow(1 - pow(2 * 0.5 - 1, 8), -0.25)));

//tests if a number has a useful value or NaN
bool IsNumber(double x)
{
	// This looks like it should always be true, 
	// but it's false if x is a NaN.
	return (x == x);
}

int sgn(int value) {
	return (0 < value) - (value < 0);
}

void fillArray(double* a, int length, double value = 0) {
	for (int i = 0; i < length; i++) {
		a[i] = value;
	}
}

double meanValue(double* liste, int start, int stop) {
	double x = 0;
	for (int i = start; i < stop; i++) {
		x += liste[i];
	}
	return x / (stop - start);
}

//prints relevant stuff
void printLattice(int* lattice, int xdim, int ydim) {
	int cStanding = 0;
	int cLying = 0;
	for (int y = 0; y < ydim; y++) {
		for (int x = 0; x < xdim; x++) {
			cout << lattice[y * xdim + x] << " ";
			cStanding += (lattice[y * xdim + x] == 2);
			cLying += (lattice[y * xdim + x] == 1);
		}
		cout << endl;
	}
	cout << cStanding << " are standing" << endl;
	cout << cLying << " are lying" << endl;
}
void printNMolecules(int* nMolecules) {
	cout << "nMolecules: " << nMolecules[0] << " " << nMolecules[1] << " " << nMolecules[2] << endl;
}
void printEnergy(double energy) {
	cout << energy << endl;
}
void printX(double T, double coverage, double meanX) {
	cout << T << " " << coverage << " " << meanX << endl;
}
void printFraction(int xdim, double* fraction) {
	for (int x = 0; x < xdim; x++) {
		cout << fraction[x] << ";";
	}
	cout << endl;
}

//returns the position of the {right, left, above, below} neighbor under periodic boundary conditions
int positionRight(int xdim, int ydim, int position) {
	if (position % xdim == xdim - 1) return position - xdim + 1;
	else return position + 1;
}
int positionLeft(int xdim, int ydim, int position) {
	if (position % xdim == 0) return position + xdim - 1;
	else return position - 1;
}
int positionAbove(int xdim, int ydim, int position) {
	if (position < xdim) return position + xdim * (ydim - 1);
	else return position - xdim;
}
int positionBelow(int xdim, int ydim, int position) {
	if (position >= xdim * (ydim - 1)) return position - xdim * (ydim - 1);
	else return position + xdim;
}

// fills the lattice with 2s and 0s and sets initial values for n1 and n2
void fillLatticeRandom2(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int* nMolecules) {
	for (int i = 0; i < xdim * ydim; i++) {
		lattice[i] = 0;
	}
	int k;
	int c = 0;
	do {
		k = gsl_rng_uniform_int(ranuge, xdim * ydim);
		if (lattice[k] != 0) continue;
		lattice[k] = 2;
		nMolecules[2]++;
		c++;
	} while (c < nMolecules[0]);
}
void fillLatticeSlab2(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int* nMolecules) {
	for (int i = 0; i < xdim * ydim; i++) {
		lattice[i] = 0;
	}
	int c = 0;
	int xMiddle = xdim / 2;
	for (int x = 0; x < xMiddle; x++) {
		for (int y = 0; y < ydim; y++) {
			lattice[xMiddle + x + xdim * y] = 2;
			nMolecules[2]++;
			c++;
			if (c == nMolecules[0]) { break; }
			lattice[xMiddle + -x + xdim * y] = 2;
			if (x != 0) {
				nMolecules[2]++;
				c++;
			}
			if (c == nMolecules[0]) { break; }
		}
		if (c == nMolecules[0]) { break; }
	}
	printNMolecules(nMolecules);
}
//shifts the whole lattice to the right if direction > 0 and to the left if direction < 0
void shiftLattice(int* lattice, int xdim, int ydim, double translation) {
	int* tempLattice = new int[xdim * ydim];
	if (translation < 0) {
		for (int position = 0; position < xdim * ydim; position++) {
			tempLattice[positionRight(xdim, ydim, position)] = lattice[position];
		}
	}
	else {
		for (int position = 0; position < xdim * ydim; position++) {
			tempLattice[positionLeft(xdim, ydim, position)] = lattice[position];
		}
	}

	for (int position = 0; position < xdim * ydim; position++) {
		lattice[position] = tempLattice[position];
	}
	delete[] tempLattice;
}
//returns differns between the average position of Molecules and the middle of the lattice
double slabOfset(double* fraction, int xdim) {
	double integral1 = 0;
	double integral2 = 0;
	for (int i = 0; i < xdim;i++) {
		integral1 += fraction[i] * i;
		integral2 += fraction[i];
	}
	integral1 /= integral2;
	//cout << integral1 << " - " << xdim/2 << " = " << integral1 - xdim / 2 << endl;
	return integral1 - xdim / 2;
}
//Energy difference h(T) = f2 - f1 
double h(double T, double epsilon);

//counts the connections to other standing molecules considering PBC
int countNeighbors(int* lattice, int xdim, int ydim, int position) {
	int neighbors = 0;
	//above
	neighbors += (lattice[positionAbove(xdim, ydim, position)] == 2);
	//below
	neighbors += (lattice[positionBelow(xdim, ydim, position)] == 2);
	//left
	neighbors += (lattice[positionLeft(xdim, ydim, position)] == 2);
	//right
	neighbors += (lattice[positionRight(xdim, ydim, position)] == 2);
	return neighbors;
}

//energy
void initEnergy(double& energy, int* lattice, int* nMolecules, int xdim, int ydim, double T, double hT) {
	int neighbors = 0;
	for (int i = 0; i < xdim * ydim; i++) {
		neighbors += (lattice[i] == 2) * countNeighbors(lattice, xdim, ydim, i);
	}
	energy = -neighbors / 2.0 + hT * nMolecules[2];//2.0 because bonds are counted twice
}

//Chance if a molecule with the given number of neighbors lays down
double chanceLayDown(int nNeighbors, double T, double hT) {
	double dE = nNeighbors - hT;
	return exp(-dE / T);
}
//Chance if a standing molecule moves
double chanceMoveStanding(int dNeighbors, double T) {
	double dE = -dNeighbors;
	return exp(-dE / T);
}

//tries to lay down a molecule at position
void layDown(gsl_rng* ranuge, int* lattice, int* nMolecules, int xdim, int ydim, int position, double* rLayDown, double& energy, double T, double hT) {
	int nNeighbors = countNeighbors(lattice, xdim, ydim, position);
	if (gsl_rng_uniform(ranuge) < rLayDown[nNeighbors]) {
		lattice[position] = 1;
		nMolecules[1]++;
		nMolecules[2]--;
		energy += (nNeighbors - hT);
	}
}

//tries to stand up a molecule at position
void standUp(gsl_rng* ranuge, int* lattice, int* nMolecules, int xdim, int ydim, int position, double* rStandUp, double& energy, double T, double hT) {
	int nNeighbors = countNeighbors(lattice, xdim, ydim, position);
	if (gsl_rng_uniform(ranuge) < rStandUp[nNeighbors]) {
		lattice[position] = 2;
		nMolecules[1]--;
		nMolecules[2]++;
		energy -= (nNeighbors - hT);
	}
}

//tests if the neighboring lattice site is free; returns True if it is free
bool testRight(int* lattice, int xdim, int ydim, int position) {
	return !(lattice[positionRight(xdim, ydim, position)]);
}
bool testLeft(int* lattice, int xdim, int ydim, int position) {
	return !(lattice[positionLeft(xdim, ydim, position)]);
}
bool testBelow(int* lattice, int xdim, int ydim, int position) {
	return !(lattice[positionBelow(xdim, ydim, position)]);
}
bool testAbove(int* lattice, int xdim, int ydim, int position) {
	return !(lattice[positionAbove(xdim, ydim, position)]);
}

//moves a lying molecule to the neighboring lattice if it is free
//ausbaufähig durch multiplikation statt if else
void moveLyingRight(int* lattice, int xdim, int ydim, int position) {
	int x = positionRight(xdim, ydim, position);
	if (lattice[x] == 0) {
		lattice[x] = 1;
		lattice[position] = 0;
	}
}
void moveLyingLeft(int* lattice, int xdim, int ydim, int position) {
	int x = positionLeft(xdim, ydim, position);
	if (lattice[x] == 0) {
		lattice[x] = 1;
		lattice[position] = 0;
	}
}
void moveLyingAbove(int* lattice, int xdim, int ydim, int position) {
	int x = positionAbove(xdim, ydim, position);
	if (lattice[x] == 0) {
		lattice[x] = 1;
		lattice[position] = 0;
	}
}
void moveLyingBelow(int* lattice, int xdim, int ydim, int position) {
	int x = positionBelow(xdim, ydim, position);
	if (lattice[x] == 0) {
		lattice[x] = 1;
		lattice[position] = 0;
	}
}

//tries to move a lying molecule
void moveLying(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position) {
	int direction = gsl_rng_uniform_int(ranuge, 4);
	if (direction == 0) {
		moveLyingAbove(lattice, xdim, ydim, position);
	}
	if (direction == 1) {
		moveLyingRight(lattice, xdim, ydim, position);
	}
	if (direction == 2) {
		moveLyingBelow(lattice, xdim, ydim, position);
	}
	if (direction == 3) {
		moveLyingLeft(lattice, xdim, ydim, position);
	}
}

//moves a standing molecule to the neighboring lattice
//ausbaufähig durch multiplikation statt if else
void moveStandingAbove(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position, double* rMoveStanding, double& energy) {
	int x = positionAbove(xdim, ydim, position);
	int dNeighbors = countNeighbors(lattice, xdim, ydim, x) - 1 - countNeighbors(lattice, xdim, ydim, position);// -1 because otherwise I count the Molecule which is moved
	if (gsl_rng_uniform(ranuge) < rMoveStanding[dNeighbors + 3]) {
		lattice[x] = 2;
		lattice[position] = 0;
		energy -= dNeighbors;
	}
}
void moveStandingBelow(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position, double* rMoveStanding, double& energy) {
	int x = positionBelow(xdim, ydim, position);
	int dNeighbors = countNeighbors(lattice, xdim, ydim, x) - 1 - countNeighbors(lattice, xdim, ydim, position);// -1 because otherwise I count the Molecule which is moved
	if (gsl_rng_uniform(ranuge) < rMoveStanding[dNeighbors + 3]) {
		lattice[x] = 2;
		lattice[position] = 0;
		energy -= dNeighbors;
	}
}
void moveStandingLeft(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position, double* rMoveStanding, double& energy) {
	int x = positionLeft(xdim, ydim, position);
	int dNeighbors = countNeighbors(lattice, xdim, ydim, x) - 1 - countNeighbors(lattice, xdim, ydim, position);// -1 because otherwise I count the Molecule which is moved
	if (gsl_rng_uniform(ranuge) < rMoveStanding[dNeighbors + 3]) {
		lattice[x] = 2;
		lattice[position] = 0;
		energy -= dNeighbors;
	}
}
void moveStandingRight(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position, double* rMoveStanding, double& energy) {
	int x = positionRight(xdim, ydim, position);
	int dNeighbors = countNeighbors(lattice, xdim, ydim, x) - 1 - countNeighbors(lattice, xdim, ydim, position);// -1 because otherwise I count the Molecule which is moved
	if (gsl_rng_uniform(ranuge) < rMoveStanding[dNeighbors + 3]) {
		lattice[x] = 2;
		lattice[position] = 0;
		energy -= dNeighbors;
	}
}

//tries to move a standing molecule
void moveStanding(gsl_rng* ranuge, int* lattice, int xdim, int ydim, int position, double* rMoveStanding, double& energy) {
	int direction = gsl_rng_uniform_int(ranuge, 4);
	if ((direction == 0) && testAbove(lattice, xdim, ydim, position)) {
		moveStandingAbove(ranuge, lattice, xdim, ydim, position, rMoveStanding, energy);
	}
	if ((direction == 1) && testRight(lattice, xdim, ydim, position)) {
		moveStandingRight(ranuge, lattice, xdim, ydim, position, rMoveStanding, energy);
	}
	if ((direction == 2) && testBelow(lattice, xdim, ydim, position)) {
		moveStandingBelow(ranuge, lattice, xdim, ydim, position, rMoveStanding, energy);
	}
	if ((direction == 3) && testLeft(lattice, xdim, ydim, position)) {
		moveStandingLeft(ranuge, lattice, xdim, ydim, position, rMoveStanding, energy);
	}
}

//reduces 2D lattice into 1D mean value of colums
void columFraction(int* lattice, int xdim, int ydim, int mType, double* fraction) {
	for (int x = 0; x < xdim; x++) {
		int c = 0;
		for (int y = 0; y < ydim; y++) {
			if (lattice[x + y * xdim] == mType) { c++; }
		}
		fraction[x] = 1.0 * c / ydim;
	}
}
//for summig up the average fractions x1 and x2
void sumFraction(double* meanX, int xdim, double* fraction) {
	for (int x = 0; x < xdim; x++) {
		meanX[x] += fraction[x];
	}
}
void divideFraction(double* meanX, int xdim, int tMax) {
	for (int x = 0; x < xdim; x++) {
		meanX[x] /= tMax;
	}
}

//writghs the current state of a lattice in a file
void writeLattice(ofstream& stream, int* lattice, int nLattice) {
	for (int i = 0; i < nLattice - 1; i++) {
		stream << lattice[i] << ";";
	}
	stream << lattice[nLattice - 1] << endl;
}
//writghs the current state of nMolecules in a file
void writeNMolecules(ofstream& stream, int* nMolecules) {
	stream << nMolecules[0] << ";" << nMolecules[1] << ";" << nMolecules[2] << endl;
}
//writghs the current state of energy in a file
void writeEnergy(ofstream& stream, double& energy) {
	stream << energy << endl;
}
//writghs the fraction x of a simulation in a file
void writeMeanX(ofstream& stream, int xdim, double T, double coverage, double* meanX1, double* meanX2) {
	stream << xdim << ";" << T << ";" << coverage;
	for (int x = 0; x < xdim; x++) {
		stream << ";" << meanX1[x];
	}
	for (int x = 0; x < xdim; x++) {
		stream << ";" << meanX2[x];
	}
	stream << endl;
}
void writeFraction(ofstream& stream, int xdim, double* fraction) {
	for (int x = 0; x < xdim - 1; x++) {
		stream << fraction[x] << ";";
	}
	stream << fraction[xdim - 1] << endl;
}

//attempts nMolecules times to alter the system
void sweep(gsl_rng* ranuge, int* lattice, int xdim, int ydim, double* rLayDown, double* rStandUp, double* rMoveStanding, double rTryFlip, int* nMolecules, double& energy, double T, double hT) {
	int c = 0;
	int position;
	double r;
	do {
		//which latticesite undergoes change
		position = gsl_rng_uniform_int(ranuge, xdim * ydim);
		//change
		if (lattice[position] == 0) continue;
		r = gsl_rng_uniform(ranuge);
		if (lattice[position] == 1) {
			if (r < rTryFlip) {
				standUp(ranuge, lattice, nMolecules, xdim, ydim, position, rStandUp, energy, T, hT);
			}
			else {
				moveLying(ranuge, lattice, xdim, ydim, position);
			}
			c++;
			continue; //continue is needed, because otherwise after standUp the lattice[position] could be 2
		}
		if (lattice[position] == 2) {
			if (r < rTryFlip) {
				layDown(ranuge, lattice, nMolecules, xdim, ydim, position, rLayDown, energy, T, hT);
			}
			else {
				moveStanding(ranuge, lattice, xdim, ydim, position, rMoveStanding, energy);
			}
			c++;
			continue;
		}
	} while (c < nMolecules[0]);
}

void oneSimulation(ofstream& streamX, int xdim, int ydim, double T, double coverage, double epsilon, int tMax, int nSnapshots, double rTryFlip, int& sim, bool saveAll) {
	//system constants and memory reservation
	int nLattice = xdim * ydim;
	int* nMolecules = new int[3];//{N_ges, N_1, N_2}
	nMolecules[0] = round(nLattice * coverage);
	nMolecules[1] = 0;
	nMolecules[2] = 0;
	int* lattice = new int[nLattice];
	double energy = 0;
	double* fraction1 = new double[xdim];
	double* fraction2 = new double[xdim];
	double* meanX1 = new double[xdim];
	double* meanX2 = new double[xdim];
	fillArray(meanX1, xdim);
	fillArray(meanX2, xdim);
	int tRelax = round(tMax * 0.3);
	int tWrite = tMax / nSnapshots;
	double translation = 0;

	//set h(T) so it is not calculated every time
	double hT = h(T,epsilon);
	//Calculate transition chances
	int nR = 5;
	//Lay down
	double* rLayDown = new double[nR];
	for (int i = 0; i < nR; i++) {
		rLayDown[i] = chanceLayDown(i, T, hT);
		//cout << chanceLayDown(i, T) << endl;
	}
	//Stand up
	double* rStandUp = new double[5];
	for (int i = 0; i < nR; i++) {
		rStandUp[i] = 1.0 / rLayDown[i];
		//cout << 1.0 / rLayDown[i] << endl;
	}
	//Move standing
	double* rMoveStanding = new double[7];
	for (int i = 0; i < 7; i++) {
		rMoveStanding[i] = chanceMoveStanding(i - 3, T);
		//cout << rMoveStanding[i] << endl;
	}
	//random number generator
	const gsl_rng_type* rngTyp;
	rngTyp = gsl_rng_mt19937;
	gsl_rng* ranuge;
	gsl_rng_env_setup();
	ranuge = gsl_rng_alloc(rngTyp);
	gsl_rng_set(ranuge, time(NULL));//seed
	//printf("generator type: %s\n", gsl_rng_name(ranuge));

	//initialize system
	fillLatticeSlab2(ranuge, lattice, xdim, ydim, nMolecules);
	//fillLatticeRandom2(ranuge, lattice, xdim, ydim, nMolecules);
	initEnergy(energy, lattice, nMolecules, xdim, ydim, T, hT);
	//open filestreams
	string stringLattice = "lattice.csv";
	string stringNMolecules = "nMolecules.csv";
	string stringEnergy = "energy.csv";
	string stringFraction1 = "fraction1.csv";
	string stringFraction2 = "fraction2.csv";
	if (saveAll) {
		stringLattice = to_string(sim) + stringLattice;
		stringNMolecules = to_string(sim) + stringNMolecules;
		stringEnergy = to_string(sim) + stringEnergy;
		stringFraction1 = to_string(sim) + stringFraction1;
		stringFraction2 = to_string(sim) + stringFraction2;

	}
	ofstream streamLattice(stringLattice);
	ofstream streamNMolecules(stringNMolecules);
	ofstream streamEnergy(stringEnergy);
	ofstream streamFraction1(stringFraction1);
	ofstream streamFraction2(stringFraction2);



	//wright comment lines in files
	streamLattice << "timesteps " << tMax << ", xdim = " << xdim << ", ydim = " << ydim << ", T = " << T << ", coverage = " << coverage << endl;
	streamEnergy << "timesteps " << tMax << ", xdim = " << xdim << ", ydim = " << ydim << ", T = " << T << ", coverage = " << coverage << endl;
	streamNMolecules << "timesteps " << tMax << ", xdim = " << xdim << ", ydim = " << ydim << ", T = " << T << ", coverage = " << coverage << endl;
	streamFraction1 << "timesteps " << tMax << ", xdim = " << xdim << ", ydim = " << ydim << ", T = " << T << ", coverage = " << coverage << endl;
	streamFraction2 << "timesteps " << tMax << ", xdim = " << xdim << ", ydim = " << ydim << ", T = " << T << ", coverage = " << coverage << endl;
	//relax system
	for (int t = 0; t < tRelax; t++) {
		//evolve Lattice
		sweep(ranuge, lattice, xdim, ydim, rLayDown, rStandUp, rMoveStanding, rTryFlip, nMolecules, energy, T, hT);
		//recenter
		if ((t % 10) == 0) {
			translation = slabOfset(fraction2, xdim);
			while (abs(translation) > 0.501) {//0.501 to prevent getting stuck between -0.5 and 0.5
				//cout << translation << " => shift!" << endl;
				shiftLattice(lattice, xdim, ydim, translation);
				columFraction(lattice, xdim, ydim, 2, fraction2);
				translation = slabOfset(fraction2, xdim);
				//cout << translation << " => shift?" << endl;
			}
		}
		//needed to correct ofset
		columFraction(lattice, xdim, ydim, 1, fraction1);
		columFraction(lattice, xdim, ydim, 2, fraction2);
		//extract data
		if (t % tWrite == 0) {
			writeLattice(streamLattice, lattice, nLattice);
			writeFraction(streamFraction1, xdim, fraction1);
			writeFraction(streamFraction2, xdim, fraction2);
			cout << "relax " << 100.0 * t / tRelax << "% von " << tRelax << endl;
		}

		if (t % 10 == 0) {
			writeNMolecules(streamNMolecules, nMolecules);
			writeEnergy(streamEnergy, energy);
		}
	}
	//time evolution
	for (int t = 0; t < tMax; t++) {
		//evolve Lattice
		sweep(ranuge, lattice, xdim, ydim, rLayDown, rStandUp, rMoveStanding, rTryFlip, nMolecules, energy, T, hT);
		//recenter
		if ((t % 10) == 0) {
			translation = slabOfset(fraction2, xdim);
			while (abs(translation) > 0.501) {//0.501 to prevent getting stuck between -0.5 and 0.5
				//cout << translation << " => shift!" << endl;
				shiftLattice(lattice, xdim, ydim, translation);
				columFraction(lattice, xdim, ydim, 2, fraction2);
				translation = slabOfset(fraction2, xdim);
				//cout << translation << " => shift?" << endl;
			}
		}

		//needed to correct ofset
		columFraction(lattice, xdim, ydim, 1, fraction1);
		columFraction(lattice, xdim, ydim, 2, fraction2);
		//extract data
		sumFraction(meanX1, xdim, fraction1);
		sumFraction(meanX2, xdim, fraction2);

		if (t % tWrite == 0) {
			//for animations
			writeLattice(streamLattice, lattice, nLattice);
			writeFraction(streamFraction1, xdim, fraction1);
			writeFraction(streamFraction2, xdim, fraction2);
			cout << 100.0 * t / tMax << "% von " << tMax << endl;
		}

		if (t % 10 == 0) {
			writeNMolecules(streamNMolecules, nMolecules);
			writeEnergy(streamEnergy, energy);
		}
	}
	//save final state
	writeLattice(streamLattice, lattice, nLattice);
	writeFraction(streamFraction1, xdim, fraction1);
	writeFraction(streamFraction2, xdim, fraction2);
	writeNMolecules(streamNMolecules, nMolecules);
	writeEnergy(streamEnergy, energy);
	//save single simulation result
	divideFraction(meanX1, xdim, tMax);
	divideFraction(meanX2, xdim, tMax);
	writeMeanX(streamX, xdim, T, coverage, meanX1, meanX2);
	//print final state
	//printLattice(lattice, xdim, ydim);
	printFraction(xdim, fraction1);
	printFraction(xdim, fraction2);
	//close filestreams
	streamLattice.close();
	streamNMolecules.close();
	streamEnergy.close();
	streamFraction1.close();
	streamFraction2.close();

	//free memory
	gsl_rng_free(ranuge); // löscht den random number generator
	delete[] nMolecules;
	delete[] lattice;
	delete[] rLayDown;
	delete[] rStandUp;
	delete[] rMoveStanding;
	delete[] fraction1;
	delete[] fraction2;
	delete[] meanX1;
	delete[] meanX2;

	sim++;
}

int main() {
	chrono::high_resolution_clock::time_point tStart = chrono::high_resolution_clock::now();
	cout << "Final 7" << endl;
	//open filestream
	ofstream streamX("xFractions.csv");
	streamX << "xdim;T;coverage;x1[::];x2[::]" << endl;

	//set physical variables
	int ydim = 40;
	int xdim = 3 * ydim;
	double T = 0.2;
	double coverage = 0.5;
	double epsilon = 0.6;
	//set simulation variables
	int tMax = 3e6;
	int nSnapshots = 5;
	double rTryFlip = 0.1;
	bool saveAll = false;

	int sim = 0;
	//simulate
	for (T = 0.1; T < 0.61; T += 0.02) {
		cout << "Simulating for T = " << T << endl;
		oneSimulation(streamX, xdim, ydim, T, coverage, epsilon, tMax, nSnapshots, rTryFlip, sim, saveAll);
		cout << "Simulations completed " << sim << endl;
	}

	//close filestreams
	streamX.close();

	chrono::high_resolution_clock::time_point tEnd = chrono::high_resolution_clock::now();
	chrono::high_resolution_clock::duration tDiff = tEnd - tStart;
	int runtime = chrono::duration_cast<chrono::milliseconds>(tDiff).count();
	cout << runtime << " Laufzeit in ms = " << round(runtime / (1000 * 60)) << "min" << endl;
	ofstream streamRunTime("runtime.txt");
	streamRunTime << runtime << " Laufzeit in ms = " << round(runtime / (1000 * 60)) << "min" << endl;
	streamRunTime.close();
	return 42;
}

//Energy difference h(T) = f2 - f1 
double h(double T, double epsilon) {
	return epsilon;
	//return epsilon - T * 0.4;
	//return epsilon - T * 6 * (T - T_c);
	//return epsilon - T * 12 * (T - T_c);
	//return epsilon - T * 3 / (2 * T_c) * (T - T_c);
}