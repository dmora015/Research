#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

/* data structures used to store ISOLATE, LOCI, and CONTAINS data */
struct iso_entry
{
    string id;
    string site;
    string latitude;
    string longitude;
    string year;
    string name;
    string plant;
    string nodule;
};

struct loci_entry
{
    char type;
    string num;
    string sequence;
};

struct contains_entry
{
    string id;
    char type;
    string num;
};

/* file operatore wrappers */
bool openfile(char* filename,  ifstream* fs)
{
    fs->open(filename);
    return fs->is_open();
}

bool openfile(char* filename,  ofstream* fs)
{
    fs->open(filename);
    return fs->is_open();
}

/* simple grouping operators */
inline bool is_num(char c)
{
    return c >= '0' && c <= '9';
}

inline bool is_alpha(char c)
{
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

/* Returns the length of the number starting at str
    Ex: get_num_count("345hello") returns 3 */
int get_num_count(const char* str)
{
    int i;
    for(i = 0; str[i]; ++i)
        if(!is_num(str[i]))
            return i;
    return i;
}

/* Returns the length of the string of alpha chars (a-z, A-Z) starting at str
    Ex: get_char_count("hello345") returns 5 */
int get_char_count(const char* str)
{
    int i;
    for(i = 0; str[i]; ++i)
        if(!is_alpha(str[i]))
            return i;
    return i;
}

/* Returns true if parsed correctly. Fullname should be of the form "05LocS2_4" */
bool parse_iso(string fullname, string& year, string& name, string& plant, string& nodule)
{
    const char* cstr = fullname.c_str();
    int i = 0;
    int j;
    
    /* get year */
    j = get_num_count(cstr+i);
    year = fullname.substr(i, j);
    i = i + j;
    
    /* get name */
    j = get_char_count(cstr+i);
    name = fullname.substr(i, j);
    i = i +j;
    
    /* get plant */
    j = get_num_count(cstr+i);
    plant = fullname.substr(i, j);
    i = i + j + 1; /* to skip over "_" */
    
    /* get nodule */
    j = get_num_count(cstr+i);
    nodule = fullname.substr(i, j);
}

bool is_valid_loci_type(char c)
{
    return c == 'K' || c == 'G' || c == 'I' || c == 'R'||
           c == 'D' || c == 'N' || c == 'Z' || c == 'L';
}

int get_next_word(string& line, int start, string& word)
{
    const char* cstr = line.c_str();
    int j = 0;
    
    /* skip to next char if first one is a ',' */
    if(line[start] == ',')
        ++start;
    
    for(int i = start; cstr[i]; ++i, ++j)
        if(cstr[i] == ',')
            break;
            
    word = line.substr(start,j);
    return start+j;
}

void remove_degree_sign(string& latitude)
{
    latitude = latitude.substr(0, latitude.size() - 3);
}

void get_iso(string& line, iso_entry& iso)
{
    string fullname;
    int i = 0;
    
    i = get_next_word(line, i, iso.id);
    i = get_next_word(line, i, iso.site);
    i = get_next_word(line, i, iso.latitude);
    i = get_next_word(line, i, iso.longitude);
    i = get_next_word(line, i, fullname);
    
    remove_degree_sign(iso.latitude);
    remove_degree_sign(iso.longitude);
    parse_iso(fullname, iso.year, iso.name, iso.plant, iso.nodule);
}

/* Parses string LINE and stores the type, num, and sequence in LOCI.
    LINE must have the form "XXXXXXXXX_K01". Returns true if the loci is valid,
    false otherwise. */
bool get_loci(string& loci_name, loci_entry& loci)
{
    /* if LOCI_NAME is null or has the form XXXXXX_UNIQUE */
    if(loci_name.size() == 0 || loci_name[loci_name.size()-4] != '_')
        return false;
    
    loci.type = loci_name[loci_name.size()-3];
    loci.num = loci_name.substr(loci_name.size()-2,2);
}

void print_iso(ostream& out, iso_entry& iso)
{
    out << iso.id << ',';
    out << iso.site << ',';
    out << iso.latitude << ',';
    out << iso.longitude << ',';
    out << iso.year << ',';
    out << iso.name << ',';
    out << iso.plant << ',';
    out << iso.nodule;
}

void print_loci(ostream& out, loci_entry& loci)
{
    out << loci.type << ',';
    out << loci.num << ',';
    out << loci.sequence;
}


ostream& operator<<(ostream& out, iso_entry& iso)
{
    print_iso(out, iso);
    return out;
}

ostream& operator<<(ostream& out, loci_entry& loci)
{
    print_loci(out, loci);
    return out;
}

void process_next_iso_line(ifstream& in, iso_entry& iso)
{
    string line;
    getline(in, line);
    get_iso(line, iso);
}

bool is_in_vec(loci_entry& loci, vector<loci_entry>& loci_vec)
{
    for(int i = 0; i < loci_vec.size(); ++i)
        if(loci.type == loci_vec.at(i).type && loci.num == loci_vec.at(i).num)
        {
            return true;
        }
    return false;
}

/* Reads in a line from IN and processes it. If the loci is not in LOCI_VEC, then
    a the new loci is saved in LOCI and added to LOCI_VEC. If it is found inside LOCI_VEC,
    then the entry is discarded. */
bool process_next_loci_lines(ifstream& in, loci_entry& loci, vector<loci_entry>& loci_vec)
{
    /* read in first line */
    string loci_name;
    getline(in, loci_name);
    
    /* check if the loci is in LOCI_VEC. If it is not found, add it to LOCI_VEC */
    get_loci(loci_name, loci);
    loci.sequence = "";
    if(!is_in_vec(loci, loci_vec))
    {
        /* the LOCI was not found in LOCI_VEC */
        loci_vec.push_back(loci);
        
        /* read in the nucleotide sequence and store it in LOCI.SEQUENCE */
        getline(in, loci.sequence);
        
        return true;
    }
    else
    {
        /* skip next sequence */
        string junk;
        getline(in, junk);
        return false;
    }
}

void process_iso_file(ifstream& in, ofstream& out)
{
    iso_entry iso;
    while(!in.eof())
    {
        process_next_iso_line(in, iso);
        out << iso << endl;
    }
}



void process_loci_file(ifstream& in, ostream& out)
{
    loci_entry loci;
    vector<loci_entry> loci_vec;
    while(!in.eof())
    {
        if(process_next_loci_lines(in, loci, loci_vec))
            out << loci << endl;
    }
}

bool open_all_files(char* argv[],  ifstream& in_iso_fs, ifstream& in_loci_K_fs, ifstream& in_loci_G_fs, ifstream& in_loci_I_fs,
    ofstream& out_iso_fs, ofstream& out_loci_fs)
{
    if(!openfile(argv[1], &in_iso_fs))
    {
        cout << "Input iso file could not be opened." << endl;
        return false;
    }
    if(!openfile(argv[2], &in_loci_K_fs))
    {
        cout << "Input loci K file could not be opened." << endl;
        return false;
    }
		if(!openfile(argv[3], &in_loci_G_fs))
		{
			cout << "Input loci G file could not be opened." << endl;
			return false;
		}
		if(!openfile(argv[4], &in_loci_I_fs))
		{
			cout << "Input loci I file could not be opened." << endl;
			return false;
		}
    if(!openfile(argv[5], &out_iso_fs))
    {
        cout << "Output iso file could not be opened." << endl;
        return false;
    }
    if(!openfile(argv[6], &out_loci_fs))
    {
        cout << "Output loci file could not be opened." << endl;
        return false;
    }
    return true;
}

void close_all_files(ifstream& in_iso_fs, ifstream& in_loci_K_fs,
    ofstream& out_iso_fs, ofstream& out_loci_fs)

{
    in_iso_fs.close();
    in_loci_K_fs.close();
    out_iso_fs.close();
    out_loci_fs.close();
}

int main(int argc, char* argv[])
{
    /* make sure input and output filenames are passed in */
    if(argc != 7)
    {
        cout << "USAGE: csv_create in_iso_file in_loci_K_file in_loci_G_file in_loci_I_file" /*in_loci_R_file*/
			 << /*in_loci_D_file in_loci_N_file in_loci_Z_file in_loci_L_file*/" out_iso_file out_loci_file"<< endl;
        return 0;
    }
    
    /* open files */
    ifstream in_iso_fs, in_loci_K_fs, in_loci_G_fs, in_loci_I_fs, in_loci_R,
						in_loci_D_fs, in_loci_N_fs, in_loci_Z_fs, in_loci_L_fs;
    ofstream out_iso_fs, out_loci_fs;
    if(! open_all_files(argv, in_iso_fs, in_loci_K_fs, in_loci_G_fs, in_loci_I_fs, out_iso_fs, out_loci_fs) )
        return 0;
    
    process_iso_file(in_iso_fs, out_iso_fs);
    process_loci_file(in_loci_K_fs, out_loci_fs);
		process_loci_file(in_loci_G_fs, out_loci_fs);
	//	process_loci_file(in_loci_I_fs, out_loci_fs);
   
    /* close files */
    close_all_files(in_iso_fs, in_loci_K_fs, out_iso_fs, out_loci_fs);
}
