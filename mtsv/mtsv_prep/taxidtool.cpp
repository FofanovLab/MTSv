#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <thread>
#include <mutex>
#include <regex>
#include <iterator>
#include <unordered_map>
#include <forward_list>

using namespace std;
mutex source_lock, GI_lock, key_words_lock, file_lock;
forward_list<unsigned long> gi_list;
forward_list<string> keywords_list, source_words_list;
set<string> key_words, source_words, empty_set;

/*
Helper function to acquire the size of files in input lists.
*/
ifstream::streampos get_capacity(string file_name){
    ifstream file(file_name.c_str(), ios::in|ios::ate );
    ifstream::streampos file_size = file.tellg();
    file.close();
    return file_size;
}

/*
Modified first fit algorithm that is used to divide work load between threads
*/
std::vector<std::set<std::string>> first_fit(std::vector<std::pair<std::string, ifstream::streampos>> size_list, ifstream::streampos threshold, unsigned long max_bins){
	std::vector<set<std::string>> chunks;
	std::pair<std::string, ifstream::streampos> key_value;
	std::string tax_id;
	ifstream::streampos tax_id_size;
	bool chunk_found;
	ifstream::streampos sum = 0;
	std::string sum_str = "total";
	std::vector<std::unordered_map<std::string, ifstream::streampos>> totals;
	for (unsigned int i = 0; i < size_list.size(); i++) {
		key_value = size_list[i];
		chunk_found = false;
		tax_id_size = key_value.second;
		tax_id = key_value.first;
		for (unsigned  int j = 0; j < chunks.size() ; j++) {
			std::set<std::string> temp;

			temp = chunks[j];

			std::unordered_map<std::string, ifstream::streampos> tmp;
			tmp = totals[j];
			sum = tmp[sum_str];
			if (tax_id_size <= (threshold - sum)) {
				chunks[j].insert(tax_id);
				sum += tax_id_size;
				totals[j][sum_str] = sum;
				chunk_found = true;
				break;
			}
		}
		if (!chunk_found) {
            if(chunks.size() < max_bins ){
                std::set<std::string> chunk;
                std::unordered_map<std::string, ifstream::streampos> tmp;
                tmp.insert({ sum_str, tax_id_size });
                totals.push_back(tmp);
                chunk.insert(tax_id);
                chunks.push_back(chunk);
            }
            else{
                unsigned int min_ind = 0;
                std::unordered_map<std::string, ifstream::streampos> tmp = totals[0];
                ifstream::streampos tot =  tmp[sum_str];

                for(unsigned int ind = 0; ind < totals.size(); ++ind){
                    tmp = totals[ind];
                    if(tmp[sum_str] < tot){
                        tot = tmp[sum_str];
                        min_ind = ind;
                    }
                }
                chunks[min_ind].insert(tax_id);
                totals[min_ind][sum_str] = tot + tax_id_size;
            }
		}

}
	return std::move(chunks);
}

/*
Builds a string vector of all files in input file list
*/
vector<pair<string,ifstream::streampos>> get_filesize_list(char* file_name){
    vector<pair<string, ifstream::streampos>> file_list;
    ifstream input;
    input.open(file_name );
    string buffer;
    getline(input, buffer);
    while(true){
        if (buffer.length()){
            file_list.push_back(make_pair(buffer, get_capacity(buffer)));
        }
        else break;
        getline(input, buffer);
    }
    sort(file_list.begin(), file_list.end(), [](const std::pair<string,ifstream::streampos> &left, const std::pair<string,ifstream::streampos> &right) {
    return left.second > right.second;});

    return file_list;
}

/*
This function takes a string of whitespace delimited words and inserts each word
into an set of strings, the delimited words are split again and each word added
i.e alpha-beta becomes alpha beta. The set is returned
*/
set<string> split(string words_str) {

	transform(words_str.begin(), words_str.end(), words_str.begin(),::tolower);
	stringstream sentence(words_str);
	set<string> words;
	string cur_word, delimited_words;

	regex alnum("[^a-zA-Z0-9]");


	while (sentence >> cur_word) {
		words.insert(cur_word);
		delimited_words = regex_replace(cur_word,alnum, " ");
		if (delimited_words.compare(cur_word) != 0) {
			stringstream sub_stream(delimited_words);
			while (sub_stream >> delimited_words)
				words.insert(delimited_words);
		}
	}
	return words;
}

/*
This function takes a string of whitespace delimited words and inserts each word
into an set of strings. Only alphanumerics are allowed.
*/
set<string> source_split(string words_str) {
    regex al_num("[^a-zA-Z0-9]");
    words_str = regex_replace(words_str, al_num, " ");
	transform(words_str.begin(), words_str.end(), words_str.begin(),::tolower);

	stringstream sentence(words_str);
	set<string> words;
	string cur_word;

	while (sentence >> cur_word) {
		words.insert(cur_word);
	}
	return words;
}

/*
This function produces a two column tab delimited database file with keywords, GIs.
The GIs in each row are ordered by the set container property and white space delimited
The rows are sorted by keyword using the set container property.
*/
void output_db(char* outfile, vector<string> &word_list ){
    vector<string> GI_sets(word_list.size());
    string line;
    ofstream output;
    forward_list<unsigned long>::iterator cur_gi = gi_list.begin();

    while(!keywords_list.empty()){
        stringstream words(keywords_list.front());
        string word;
        while (words >> word) {
            auto idx = lower_bound(word_list.begin(), word_list.end(),word) - word_list.begin();
            GI_sets[idx] += to_string(*cur_gi) + " ";
        }
        keywords_list.pop_front();
        ++cur_gi;
    }
    output.open(outfile);
    unsigned long idx = 0;
    for(vector<string>::iterator word = word_list.begin(); word != word_list.end(); ++word ){
        line += *word + "\t";
        line += GI_sets[idx];
//        for(set<unsigned long>::iterator gi = GI_sets[idx].begin(); gi != GI_sets[idx].end(); ++gi){
//            line += to_string(*gi) + " ";
//        }
        output << line.substr(0,line.length()-1) << endl;
        line.clear();
        idx++;
    }
    output.close();
}

/*
This function produces a two column tab delimited database file with Source words, GIs.
The GIs in each row are ordered by the set container property and white space delimited
The rows are sorted by Source word using the set container property.
*/
void output_source(char* outfile, vector<string> &word_list){
    vector<string> GI_sets(word_list.size());
    string line;
    ofstream output;
    forward_list<unsigned long>::iterator cur_gi = gi_list.begin();

    while(!source_words_list.empty()){
        stringstream words(source_words_list.front());
        string word;
        while (words >> word) {


//        for(set<string>::iterator word = source_words_list.front().begin(); word!= source_words_list.front().end(); ++word){
            auto idx = lower_bound(word_list.begin(), word_list.end(),word) - word_list.begin();
            GI_sets[idx] += to_string(*cur_gi) + " ";
        }
        source_words_list.pop_front();
        ++cur_gi;
    }
    output.open(outfile);
    unsigned long idx = 0;

    for(vector<string>::iterator word = word_list.begin(); word != word_list.end(); ++word ){
        line += *word + "\t";
        line += GI_sets[idx] + " ";
//        for(set<unsigned long>::iterator gi = GI_sets[idx].begin(); gi != GI_sets[idx].end(); ++gi){
//            line += to_string(*gi) + " ";
//        }
        output << line.substr(0,line.length()-1) << endl;
        line.clear();
        idx++;
    }
    output.close();
}

/*
This function produces a Three column tab delimited database file with GIs, Keywords, Source words.
The Keywords and source words in each row are ordered by the set container property and white space delimited.
The rows are unsorted but are later sorted with the unix sort command.
*/
void output_gi_key_src(char* outfile){
    string line ="";
    ofstream output;
    forward_list<unsigned long>::iterator cur_gi = gi_list.begin();
    forward_list<string>::iterator cur_key = keywords_list.begin();
    forward_list<string>::iterator cur_src = source_words_list.begin();

    output.open(outfile);
    while(cur_gi != gi_list.end()){
        line += to_string(*cur_gi) + "\t";
        line += *cur_key + "\t";
        line += *cur_src;
        output << line << endl;
        line.clear();
        ++cur_key;
        ++cur_src;
        ++cur_gi;
    }
    output.close();
}

/*
This function parses out various information from a set of flat files and writes the sequence data found is written
into the fasta format used by the taxclipper pipeline. The information parsed are GIs, Keywords, Source
words maintaining association by building forward_list structures at the same time. The words themselves
are initially inserted into sets to acquire the ordered and unique properties that container has. When
the words are inserted into the global keyword set a space delimited string is built to use for final processing.
*/
void producer_seqs(set<std::string> file_list, char* seqs_out, char* position_out){

    string command, buffer, line, temp, temp_key ="", temp_src = "", sequences ="", positions="";
    const int BUFFER_SIZE = 2500000;
    char in_buffer[BUFFER_SIZE];
    sequences.reserve(BUFFER_SIZE);
    FILE* in;
    ofstream output;
    regex allowed("[^a-zA-Z0-9:'_. *-]");
    regex alpha("[^a-zA-Z]");
    unsigned long  gi ;

    for(set<string>::iterator file = file_list.begin(); file != file_list.end(); ++file){
            buffer.clear();
            command = "cat ";
            command += *file;
            command += " | tr -d '\\r'";

            in = popen(command.c_str() , "r");

            while(fgets(in_buffer, sizeof(in_buffer), in) != NULL){
                buffer += in_buffer;
            }
            pclose(in);

//        }
        if (buffer.length() != 0){
            stringstream input(buffer);
            set<string> definitions, keywords;
            string version, source;
            getline(input,line);
            while(!input.eof()){
                if(line.length())

                    line = line.substr(line.find_first_not_of(" \t"), line.length());
                    temp = line.substr(0, line.find_first_of(" ") );
                    if(temp == "DEFINITION"){
                        command.clear();
                        temp = line.substr(line.find_first_of(" "), line.length());

                        getline(input,line);
                        command = line.substr(0, line.find_first_of(" ") );
                        while(command != "ACCESSION"){
                            temp += " " + line.substr(line.find_first_of(" "), line.length());
                            getline(input,line);
                            command = line.substr(0, line.find_first_of(" ") );
                        }

                        temp = regex_replace( temp , allowed, " " );
                        definitions = split(temp);
                    }
                    else if(temp == "VERSION"){
                        line = line.substr(line.find_first_of(" "), line.length() );
                        line = line.substr(line.find_first_not_of(" "), line.length() );
                        version =  regex_replace( line, regex("  +"), " " );
                        positions += version;
                    }
                    else if(temp == "SOURCE"){
                        line = line.substr(line.find_first_of(" "), line.length() );
                        source = line.substr(line.find_first_not_of(" "), line.length() );
                        keywords = source_split(source);
                    }
                    else if(temp == "ORIGIN"){
                        sequences += ">" + version + " " + source + "\n";
                        getline(input, temp);
                        temp = regex_replace(temp, alpha, "");
                        transform(temp.begin(), temp.end(), temp.begin(),::toupper);
                        while(temp.length()){
                            sequences += temp +"\n";
                            getline(input, temp);
                            temp = regex_replace(temp, alpha, "");
                            transform(temp.begin(), temp.end(), temp.begin(),::toupper);
                        }
                        version = version.substr(version.find_first_of(":")+1, version.length()  );
                        gi = strtoul(version.c_str() , NULL, 0);
                        for(set<string>::iterator word = keywords.begin(); word != keywords.end(); ++word ){
                            source_lock.lock();
                            source_words.insert(*word);
                            source_lock.unlock();
                            temp_src += *word + " ";
                        }
                        temp_src.replace(temp_src.length()-1, 1, "");
                        for(set<string>::iterator word = definitions.begin(); word != definitions.end(); ++word ){
                            key_words_lock.lock();
                            key_words.insert(*word);
                            key_words_lock.unlock();
                            temp_key += *word +" ";
                        }
                        temp_key.replace(temp_key.length()-1, 1, "");
                        GI_lock.lock();
                        gi_list.push_front(gi);
                        source_words_list.push_front(temp_src);
                        keywords_list.push_front(temp_key);
                        GI_lock.unlock();
                        temp_src.clear();
                        temp_key.clear();
                    }
                    else if(temp == "gene"){
//                        temp = line.substr(line.find_first_of(" \t"), line.length());
//                        temp = temp.substr(temp.find_first_not_of(" \t"), temp.length());
//                        getline(input,line);
//                        temp = line.substr(line.find_first_of('"'), line.find_last_of('"')) + "|" + temp;
//                        getline(input,line);
//                        temp = line.substr(line.find_first_of('"'), line.find_last_of('"')) + "|" + temp;
//                        positions += "\t" + temp;
                    }
                    getline(input,line);
                }


                if (sequences.length() > 0){
                    file_lock.lock();
                    output.open(seqs_out, ofstream::app);
                    output <<sequences;
                    output.close();
//                    output.open(position_out, ofstream::app);
//                    output << positions;
//                    output << "\n";
//                    output.close();
                    file_lock.unlock();
                    string().swap(sequences);
                    string().swap(positions);
                }
//                else if( file_lock.try_lock()){
//                    output.open(seqs_out, ofstream::app);
//                    output <<sequences;
//                    output.close();
//                    file_lock.unlock();
//                    sequences.clear();
//                    }

                cout << *file << endl;

            }
        }

        if (sequences.length()){
            file_lock.lock();
            output.open(seqs_out, ofstream::app);
            output <<sequences;
            output.close();
//            output.open(position_out, ofstream::app);
//            output <<positions;
//            output << "\n";
//            output.close();
            file_lock.unlock();
            sequences.clear();
            positions.clear();
        }

    }

/*
Thread function that will have a thread compress a file using gzip
*/
void compress(char* file){
            string buffer, command;
            FILE* in;
            const int BUFFER_SIZE=128;
            char in_buffer[BUFFER_SIZE];
            command = "gzip ";
            command += file;

            in = popen(command.c_str() , "r");


            while(fgets(in_buffer, sizeof(in_buffer), in) != NULL){
                buffer += in_buffer;
            }
            pclose(in);
            if(buffer.length()){
                cout << "Compression Failed\n";
            }

}

/*
Thread function that will have a thread sort the gi to keyword/source word file
*/
void sort_file(char* file){
            string buffer, command;
            FILE* in;
            const int BUFFER_SIZE=128;
            char in_buffer[BUFFER_SIZE];
            command = "sort -gk1 ";
            command += file;
            command += " -o ";
            command += file;
            in = popen(command.c_str() , "r");


            while(fgets(in_buffer, sizeof(in_buffer), in) != NULL){
                buffer += in_buffer;
            }
            pclose(in);
            if(buffer.length()){
                cout << "Sort Failed\n";
            }

}

/*  Main Function outline
1) Get total threads to spawn from parameters and build a vector pairing filenames and file sizes from the file list parameter
2) Use modified first fit algorithm to distribute work load between threads
3) Call producer function for passing each thread a set of files to work on.
3) Upon threads return if the 6th parameter is present use it to build a gi,keyword,source word source word file tab delimited file
4) A thread is then dispatched to build the Keyword to GIs database file consuming the keywords data structure in the process
5) A thread is then dispatched to build the source words to GIs database file consuming the source words data structure
6) Threads return and main complete
*/
int main(int argc, char** argv){
        vector<thread> thread_vector;
        unsigned long threads = atol(argv[4]);

        vector<pair<string,ifstream::streampos>> file_pairs =  get_filesize_list(argv[1]);

        ifstream::streampos files_size = 0;
        for(unsigned int i =0 ; i< file_pairs.size(); i++ )
            files_size += file_pairs[i].second;

        vector<set<string>> file_sets = first_fit(file_pairs, files_size/threads, threads );
        ofstream out;
        out.open(argv[2]);
        out.close();
        out.open(argv[3]);
        out.close();

        while(file_sets.size()){
            thread_vector.push_back(thread(producer_seqs, file_sets.back(), argv[2], argv[3] ));
            file_sets.pop_back();
        }
        for (auto& thread : thread_vector) {
            thread.join();
        }
        thread_vector.clear();

    return 0;
}
