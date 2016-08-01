#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



using namespace std;



char randomNucleotide(){
	switch (rand() % 4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}



string randomSequence(const uint length){
	string result(length, 'A');
	for(uint i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}



string mutateSequence(const string& referenceSequence, double mutRate=4.0, vector <double> ratioMutation={0.06,0.73,0.21}){
	double substitutionRate(mutRate * ratioMutation[0]);
	double insertionRate(mutRate * ratioMutation[1]);
	double deletionRate(mutRate * ratioMutation[2]);
	string result;
	result.reserve(2 * referenceSequence.size());
	for(uint i(0); i < referenceSequence.size(); ++i){
		uint dice(rand() % 100);
		if(dice<substitutionRate){
			//SUBSTITUTION
			char newNucleotide(randomNucleotide());
			while(newNucleotide == referenceSequence[i]){
				newNucleotide = randomNucleotide();
			}
			result.push_back(newNucleotide);
			continue;
		}
		if(dice < deletionRate+substitutionRate){
			//DELETION
			
			continue;
		}
		if(dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randomNucleotide());
			result.push_back(newNucleotide);
			--i;
			continue;
		}
		//NO ERROR
		result.push_back(referenceSequence[i]);
	}
	return result;
}



vector<string> generateAlternativeTranscriptReferences(uint transcriptNumber=2, uint totalExonNumber=5, uint exonNumber=3, uint sizeExons=100){
	vector<string> result;
	vector<string> exonList;
	for(uint i(0); i < totalExonNumber; ++i){
		exonList.push_back(randomSequence(sizeExons));
	}
	string transcript;
	transcript.reserve(exonNumber*sizeExons);
	unordered_set<uint> selectedExons;
	for(uint i(0); i < transcriptNumber; ++i){
		transcript = "";
		selectedExons = {};
		while(selectedExons.size() != exonNumber){
			selectedExons.insert(rand() % totalExonNumber);
		}
		for(uint ii(0); ii < totalExonNumber; ++ii){
			if(selectedExons.count(ii) == 1){
				transcript += exonList[ii];
				
			}
		}
		result.push_back(transcript);
	}
	return result;
}




void generateReads(uint numberReads, uint referencesNumber=2, const string& outFileName="simulatedReads.fa", const string& refFileName="RefFile"){
	ofstream out(outFileName);
	ofstream outRef(refFileName);
	vector<vector<string>> referenceList;
	for(uint i(0);i < referencesNumber; ++i){
		referenceList.push_back(generateAlternativeTranscriptReferences());
		for(uint ii(0); ii<referencesNumber; ++ii){
			outRef << ">referenceNumber:" << i << " alternativeNumber" << ii << endl;
			outRef << referenceList[i][ii] << endl;
		}
	}
	
	string refRead,realRead;
	for(uint i(0); i < numberReads; ++i){
		uint dice1(rand() % referencesNumber);
		uint dice2(rand() % referenceList[dice1].size());
		refRead = referenceList[dice1][dice2];
		realRead = mutateSequence(refRead);
		out << ">referenceNumber:" << dice1 << " alternativeNumber" << dice2 << endl;
		out << realRead << endl;
	}
}



int main(int argc, char ** argv){
	srand (time(NULL));
	auto startChrono = chrono::system_clock::now();
	generateReads(1000);
	auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
	cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	return 0;
}
