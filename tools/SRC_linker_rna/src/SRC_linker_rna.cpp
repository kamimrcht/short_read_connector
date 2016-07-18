#include <SRC_linker_rna.hpp>


using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_URI_BANK_INPUT = "-bank";
//~ static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_GAMMA = "-gamma";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_WINDOW = "-window_size";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";


SRC_linker_rna::SRC_linker_rna ()  : Tool ("SRC_linker_rna"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	//~ getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal percentage of shared kmer in a region for considering 2 reads in a same group.",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_WINDOW, "Size of a region (putative exon).",    false, "80"));
	getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value",    false, "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
}


void SRC_linker_rna::create_quasi_dictionary (int fingerprint_size, int nbCores){
	const int display = getInput()->getInt (STR_VERBOSE);
	// We get a handle on the HDF5 storage object.
	// Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	// We get the group for dsk
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
	if(nbSolidKmers==0){
		cout<<"No solid kmers in bank -- exit"<<endl;
		exit(0);
	}
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	quasiDico = quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper, u_int32_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size, nbCores, gamma_value);
}


struct FunctorIndexer{
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int kmer_size;

	FunctorIndexer(quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size)  :  quasiDico(quasiDico), kmer_size(kmer_size) {
	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq,kmer_size)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
		u_int32_t read_id = static_cast<u_int32_t>(seq.getIndex()+1);
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// Adding the read id to the list of ids associated to this kmer.note that the kmer may not exist in the dictionary if it was under the solidity threshold.in this case, nothing is done
			quasiDico.set_value((itKmer)->value().getVal(), read_id);
		}
	}
};



void SRC_linker_rna::fill_quasi_dictionary (const int nbCores, const string& bankName){
	bool exists;
	IBank* bank = Bank::open (bankName);
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size));
}



class FunctorQuerySpanKmers // FunctorQuery used after claires discussion: number of positions covered by a shared kmer.
{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	int kmer_size;
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t>* quasiDico;
	std::unordered_map<uint64_t, vector<uint32_t>> reads_sharing_kmer_2_positions;  // store the position where a k-mer is seen in a read that can be potentially recruited
	std::unordered_map<uint64_t, vector<readGrouped>> read_group;  // for a read, get all reads sharing at least a window
	uint threshold;
	uint size_window;
	vector<uint32_t> associated_read_ids;
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;
    
	FunctorQuerySpanKmers(const FunctorQuerySpanKmers& lol) // used by the dispatcher
	{
		size_window=lol.size_window;
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
		quasiDico=lol.quasiDico;
		threshold=lol.threshold;
		associated_read_ids=lol.associated_read_ids;
		reads_sharing_kmer_2_positions = lol.reads_sharing_kmer_2_positions;
		read_group =  lol.read_group;
		model=lol.model;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >* quasiDico, const int threshold, const uint size_window, std::unordered_map<uint64_t, vector<uint>>& reads_sharing_kmer_2_positions, std::unordered_map<uint64_t, vector<readGrouped>> read_group)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), quasiDico(quasiDico), threshold(threshold),  size_window(size_window), reads_sharing_kmer_2_positions(reads_sharing_kmer_2_positions), read_group(read_group){
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
	}
    
	FunctorQuerySpanKmers (){
	}
    
	void operator() (Sequence& seq){
	    if (not valid_sequence(seq, kmer_size)){return;}
		bool exists;
		associated_read_ids = {}; // list of the ids of reads from the bank where a kmer occurs
		reads_sharing_kmer_2_positions = {};  // store the position where a k-mer is seen in a read that can be potentially recruited
		//~ read_group = {};
		itKmer->setData (seq.getData());
		uint i(0); // position on the read
		uint64_t seqIndex(seq.getIndex() + 1);
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
		    quasiDico->get_value((*itKmer)->value().getVal(), exists, associated_read_ids);
		    if(!exists) {++i; continue;}
		    for (uint r(0); r < associated_read_ids.size(); ++r){
			if (reads_sharing_kmer_2_positions.count(associated_read_ids[r])){
			    reads_sharing_kmer_2_positions[associated_read_ids[r]].push_back(i);
			} else {
			    if (associated_read_ids[r] != seq.getIndex() + 1){  // we dont want to store the info about a read similar to itself
				//~ cout << associated_read_ids[r] << endl;
				reads_sharing_kmer_2_positions[associated_read_ids[r]].push_back(i);
			    }
			}
		    }
		    ++i;
		}
		for (auto r(reads_sharing_kmer_2_positions.begin()); r != reads_sharing_kmer_2_positions.end(); ++r){
		    size_t lenseq = seq.getDataSize();
		    vector<uint> presence(uint(lenseq) - kmer_size + 1, 0);
		    uint count(0);
		    bool found(false);
		    for (uint j(0); j < r->second.size(); ++j){
			presence[r->second[j]] = 1;
		    }
		    for (uint w(0); w < presence.size(); ++w){
			if (w < size_window){
			    if (presence[w] == 1){
				++count;
			    }
			} else {
			    uint start(w - size_window + 1);
			    if (presence[start - 1] == 1){
				--count;
			    }
			    if (presence[w] == 1){
				++count;
			    }
			}
			if (count / (size_window - kmer_size + 1) * 100 >= threshold){
			    found = true;
			    break;
			}
			
		    }
		    if (found){
			bool confirm(false);
			if (read_group.count(r->first)){
			    for (uint rr(0); rr < read_group[r->first].size(); ++rr){
				if (read_group[r->first][rr].index == seqIndex){ // the two reads have been found matching
				    read_group[r->first][rr].confirmed = true;
				    confirm = true;
				    break;
				}
			    }
			}
			cout << confirm << endl;
			if (read_group.count(seqIndex)){
			    read_group[seqIndex].push_back({r->first, confirm});
			} else {
			    readGrouped rg({r->first, confirm});
			    vector <readGrouped> v({rg});
			    read_group[seqIndex] = {v};
			}
		    }
		}
	    string toPrint;
	    bool read_id_printed = false; // Print (and sync file) only if the read is similar to something.
	    if (not read_group[seqIndex].empty()){
		if (not read_id_printed){
		    read_id_printed = true;
		    toPrint = to_string(seqIndex) + ":";
		}
		for (uint i(0); i < read_group[seqIndex].size(); ++i){
		    //~ if (read_group[seqIndex][i].confirmed){
			toPrint += to_string(read_group[seqIndex][i].index) + " ";
		    //~ }
		}
		if (read_id_printed){
		    synchro->lock();
		    toPrint += "\n";
		    fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
		    synchro->unlock();
		}
	    }
	}
};






void SRC_linker_rna::parse_query_sequences(int threshold, uint size_window, const int nbCores, const string& bankName){
    IBank* bank = Bank::open(bankName);
    cout<<"Query "<<kmer_size<<"-mers from bank "<< bankName <<endl;
    FILE * outFile;
    outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n");
    fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    LOCAL (bank);
    ProgressIterator<Sequence> itSeq (*bank);
    ISynchronizer* synchro = System::thread().newSynchronizer();
    Dispatcher dispatcher (1, 10000);
    //~ Dispatcher dispatcher (nbCores, 10000);
    std::unordered_map<uint64_t, vector<uint>> reads_sharing_kmer_2_positions;
    std::unordered_map<uint64_t, vector<readGrouped>> read_group;
    dispatcher.iterate(itSeq, FunctorQuerySpanKmers(synchro, outFile, kmer_size, &quasiDico, threshold, size_window, reads_sharing_kmer_2_positions, read_group));
    fclose (outFile);
    delete synchro;
}


void SRC_linker_rna::execute (){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
	gamma_value = getInput()->getInt(STR_GAMMA);
	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (ven with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently us a vector<u_int32_t> for storing read ids associated to a kmer.

	// We need a non null finger print because of non solid non indexed kmers
	//	if (getInput()->getStr(STR_URI_BANK_INPUT).compare(getInput()->getStr(STR_URI_QUERY_INPUT))==0)
	//		fingerprint_size=0;
	cout<<"fingerprint = "<<fingerprint_size<<endl;
	create_quasi_dictionary(fingerprint_size, nbCores);
	string bankName(getInput()->getStr(STR_URI_BANK_INPUT));
	fill_quasi_dictionary(nbCores, bankName);

	int threshold = getInput()->getInt(STR_THRESHOLD);
	uint size_window =  getInput()->getInt(STR_WINDOW);
	parse_query_sequences(threshold, size_window, nbCores, bankName);
	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Sequences bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	//~ getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
	getInfo()->add (2, "gamma",  "%d",  gamma_value);
	getInfo()->add (2, "Minimal kmer span percentage",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}
