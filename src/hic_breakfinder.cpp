#include "api/BamReader.h"
#include "api/BamAux.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <getopt.h>
#include <unordered_map>
#include <math.h>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace BamTools;
using namespace std;

void usage () {

  std::cerr << "./hic_breakfinder\n\n";
  std::cerr << "Required options:\n\n";
  std::cerr << "--bam-file [input bam file]\n";
  std::cerr << "--exp-file-inter [Inter-chromosomal 1Mb expectation file]\n";
  std::cerr << "--exp-file-intra [Intra-chromosomal 100kb expectation file]\n";
  std::cerr << "--name [output file name prefix, will append with *.super_matrix.txt and *.SR.txt]\n\n";  
  
}

void split_string (std::vector<std::string> &output_vector,const std::string &input_string, const char &input_delimiter) {
  
  std::stringstream iss0(input_string);
  std::string token0;
  int count0 = 0;
  while (getline(iss0,token0,input_delimiter)) {
    output_vector.push_back(token0);
  }

}

void get_super_matrix (char* bam_file,\
		       const string name,\
		       const int bin_size,\
		       const int clean_flag,\
		       const int filter_flag,\
		       const float upper,\
		       const float lower) {
  
  unordered_map<string,int> chr_name;
  for (long long int i = 1; i <= 22; i++) {
    string out_chr_1 = "chr" + to_string(i);
    string out_chr = to_string(i);
    chr_name[out_chr_1] = 1;
    chr_name[out_chr] = 1;
  }

  string chrX_1 = "chrX";
  string chrY_1 = "chrY";
  string chrX = "X";
  string chrY = "Y";
  chr_name[chrX_1] = 1;
  chr_name[chrY_1] = 1;
  chr_name[chrX] = 1;
  chr_name[chrY] = 1;

  BamReader reader;

  if ( !reader.Open(bam_file) ) {
    std::cerr << "Could not open input BAM file.\n";
    return;
  }

  reader.Open(bam_file);

  RefVector refs = reader.GetReferenceData();

  unordered_map<string,unordered_map<string,int> > data_map;
  unordered_map<string,int> data_sum;
  unordered_map<string,int> SR_sum;
  unordered_map<string,int> chr_hash;

  BamAlignment bam;

  string last_chr;

  while ( reader.GetNextAlignmentCore(bam) ) {
    
    if (bam.RefID >= 0) {
      
      if (refs.at(bam.RefID).RefName != last_chr) {
	cerr << "Going through " << refs.at(bam.RefID).RefName << " now\n";
	last_chr = refs.at(bam.RefID).RefName;
      }

      long long int Bin = bin_size*floor(bam.Position/bin_size);
      string tag1 = refs.at(bam.RefID).RefName + ":" + to_string(Bin) + "-" + to_string(Bin + bin_size);
      chr_hash[tag1] = bam.RefID;

      if (bam.MateRefID >= 0) {

	long long int MateBin = bin_size*floor(bam.MatePosition/bin_size);
	string tag2 = refs.at(bam.MateRefID).RefName + ":" + to_string(MateBin) + "-" + to_string(MateBin + bin_size);

	if (bam.RefID == bam.MateRefID) {

	  if (Bin < MateBin) {

	    if (clean_flag == 1) {
	         
	      if ((chr_name.count(refs.at(bam.RefID).RefName) > 0) &&
		  (chr_name.count(refs.at(bam.MateRefID).RefName)) > 0) {
		  
		data_map[tag1][tag2]++;
		data_sum[tag1]++;
		data_sum[tag2]++;

	      }

	    } else {

	      data_map[tag1][tag2]++;
	      data_sum[tag1]++;
	      data_sum[tag2]++;

	    }

	  }
	  
	} else {

	  if (bam.RefID < bam.MateRefID) {

	    if (clean_flag == 1) {
	          
	      if ((chr_name.count(refs.at(bam.RefID).RefName) > 0) &&
                  (chr_name.count(refs.at(bam.MateRefID).RefName)) > 0) {

		data_map[tag1][tag2]++;
		data_sum[tag1]++;
                data_sum[tag2]++;

              } 

	    } else {

	      data_map[tag1][tag2]++;
	      data_sum[tag1]++;
	      data_sum[tag2]++;

	    }

	  }
	    
	}

      } else {
      
	if (clean_flag == 1) {
	    
	  if (chr_name.count(refs.at(bam.RefID).RefName) > 0) {

	    SR_sum[tag1]++;
	    data_sum[tag1]++;

	  }

	} else {

	  SR_sum[tag1]++;
	  data_sum[tag1]++;
      
	}

      }

    }
  }

  reader.Close();

  unordered_map<int,int> sum_hist;
  unordered_map<string,int>::iterator iter1;
  int total_N = 0;
  for (iter1 = data_sum.begin(); iter1 != data_sum.end(); iter1++) {
    sum_hist[iter1->second]++;
    total_N++;
  }
  
  vector<int> count_vector;
  unordered_map<int,int>::iterator iter2;
  for (iter2 = sum_hist.begin(); iter2 != sum_hist.end(); iter2++) {
    count_vector.push_back(iter2->first);
  }
  
  sort(count_vector.begin(),count_vector.end());
  
  int lower_N;
  int upper_N;
  int lower_check = 0;
  
  int run_sum = 0;
  
  for (int i=0; i < count_vector.size(); i++) {
    
    if ((run_sum >= lower*total_N) &&
	(lower_check == 0)) {
      lower_N = count_vector[i];
      lower_check = 1;
    }
    
    if (upper*total_N >= run_sum) {
      upper_N = count_vector[i];
    }
    
    run_sum += sum_hist[count_vector[i]];
    
  }
  
  ofstream super_matrix_output;
  ofstream SR_file_output;
  
  string SR_file_name = name + ".SR.txt";
  string super_matrix_file_name = name + ".super_matrix.txt";
  
  SR_file_output.open(SR_file_name);
  super_matrix_output.open(super_matrix_file_name);
    
  unordered_map<string,unordered_map<string,int> >::iterator out_iter1;
  for (out_iter1 = data_map.begin(); out_iter1 != data_map.end(); out_iter1++) {
    unordered_map<string,int>::iterator out_iter2;
    
    if ((data_sum[out_iter1->first] >= lower_N) &&
	(data_sum[out_iter1->first] <= upper_N)) {
      
      SR_file_output << out_iter1->first << "\t" << SR_sum[out_iter1->first] << "\n";

    }
    
    for (out_iter2 = data_map[out_iter1->first].begin(); out_iter2 != data_map[out_iter1->first].end(); out_iter2++) {
      
      if ((data_sum[out_iter1->first] >= lower_N) &&
	  (data_sum[out_iter1->first] <= upper_N) &&
	  (data_sum[out_iter2->first] >= lower_N) &&
	  (data_sum[out_iter2->first] <= upper_N)) {
	  
	super_matrix_output << out_iter1->first << "\t" << out_iter2->first << "\t" << data_map[out_iter1->first][out_iter2->first] << "\n";
	
      }
      
    }
    
  }
  
  SR_file_output.close();
  super_matrix_output.close();
  
}

int string_to_int (const string &input_string) {

  int value;
  stringstream convert(input_string);
  if ( !(convert >> value) );
  return value;

}

void read_super_matrix_iterative_correction (const string super_matrix_file, \
					     vector<double>& data_vector, \
					     vector<string>& tag1_vector, \
					     vector<string>& tag2_vector, \
					     unordered_map<string,double>& bias_map, \
					     int inter_flag) {

  std::time_t result = std::time(0);

  cerr << "Reading in super matrix file now... ";

  string line;
  ifstream myfile (super_matrix_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    vector<string> loc1;
    split_string(loc1,array[0],':');
    
    vector<string> loc2;
    split_string(loc2,array[1],':');

    double n = stod(array[2]);

    if (inter_flag == 1) {
    
      if (loc1[0] != loc2[0]) {
	
	data_vector.push_back(n);
	tag1_vector.push_back(array[0]);
	tag2_vector.push_back(array[1]);

	bias_map[array[0]] = 1;
	bias_map[array[1]] = 1;
	
      }

    } else {

      data_vector.push_back(n);
      tag1_vector.push_back(array[0]);
      tag2_vector.push_back(array[1]);
      
      bias_map[array[0]] = 1;
      bias_map[array[1]] = 1;
      
    }

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void read_SR_file (const string SR_file_name,			\
		   vector<double>& SR_vector,			\
		   vector<string>& SR_tag_vector,		\
		   unordered_map<string,double>& bias_map) {
  
  cerr << "Reading in SR file now... ";
  
  string line;
  ifstream myfile (SR_file_name);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    double n = stod(array[1]);
    
    SR_vector.push_back(n);
    SR_tag_vector.push_back(array[0]);
    bias_map[array[0]] = 1;

  }

  cerr << "Done\n";

}

void iterative_correction(const string super_matrix,	\
			  const string SR_file,		\
			  const string name,		\
			  const int iterations) {
  
  vector<double> data_vector;
  vector<string> tag1_vector;
  vector<string> tag2_vector;
  unordered_map<string,double> bias_map;
  read_super_matrix_iterative_correction(super_matrix,data_vector,tag1_vector,tag2_vector,bias_map,0);
  
  vector<double> SR_vector;
  vector<string> SR_tag_vector;
  
  read_SR_file(SR_file,SR_vector,SR_tag_vector,bias_map);
  
  double mean_b = 1;
  
  int count = 1;
  
  while (count <= iterations) {
    
    cerr << "Going through iteration " << count << " now... ";
    
    //Step 1: Calculate the sum of W over all rows
    
    unordered_map<string,double> sum_hash;
    
    for (int i=0; i < data_vector.size(); i++) {
      sum_hash[tag1_vector[i]] += data_vector[i];
      sum_hash[tag2_vector[i]] += data_vector[i];
    }
    
    //Step 2: Calculate the vector of corrected SS reads
    
    unordered_map<string,double> SR_map;
    
    for (int i=0; i < SR_vector.size(); i++) {
      double val = SR_vector[i]/(mean_b*bias_map[SR_tag_vector[i]]);
      SR_map[SR_tag_vector[i]] = val;
    }
    
    //Step 3: Calculate the additional vector of biases
    
    unordered_map<string,double> bias_vector;
    double bias_sum = 0;
    int bias_N = 0;
    
    unordered_map<string,double>::iterator iter4;
    for (iter4 = bias_map.begin(); iter4 != bias_map.end(); iter4++) {
      bias_vector[iter4->first] = sum_hash[iter4->first] + SR_map[iter4->first];
      if ((sum_hash[iter4->first] + SR_map[iter4->first]) > 0) {
	bias_sum += (sum_hash[iter4->first] + SR_map[iter4->first]);
	bias_N++;
      }
    }
    
    double bias_mean = bias_sum/bias_N;
    
    //Step 4: Renormalize bias vector by its mean value to avoid instabilities
    //Step 5: Set zero values of bias vector to 1 to avoid 0/0;
    
    double var_sum = 0;
    double var_N = 0;
    
    unordered_map<string,double>::iterator iter5;
    for (iter5 = bias_vector.begin(); iter5 != bias_vector.end(); iter5++) {
      double temp = bias_vector[iter5->first];
      double new_val;
      if (temp == 0) {
	new_val = 1;
      } else {
	new_val = temp/bias_mean;
	var_sum += (temp - bias_mean)*(temp - bias_mean);
	var_N++;
      }
      bias_vector[iter5->first] = new_val;      
    }
    
    double variance = var_sum/var_N;
    
    cerr << "variance is " << variance << "\n";
    
    //Step 6: Divide all Wij by product of bias vectors\n";
    
    for (int i=0; i < data_vector.size(); i++) {
      double b1 = bias_vector[tag1_vector[i]];
      double b2 = bias_vector[tag2_vector[i]];
      double val = data_vector[i]/(b1*b2);
      data_vector[i] = val;
    }
    
    //Step 7: Multiple the total vector of biases by additional biases
    
    double sum_b = 0;
    double N_b = 0;
    
    unordered_map<string,double>::iterator iter8;
    for (iter8 = bias_map.begin(); iter8 != bias_map.end(); iter8++) {
      double temp = bias_map[iter8->first]*bias_vector[iter8->first];
      bias_map[iter8->first] = temp;
      sum_b += temp;
      N_b++;
    }
    
    mean_b = sum_b/N_b;
    
    count++;
  }
  
  ofstream super_matrix_output;
  
  string super_matrix_file_name = name + ".super_matrix.norm.txt";
  
  super_matrix_output.open(super_matrix_file_name);
  
  for (int i=0; i < data_vector.size(); i++) {
    super_matrix_output << tag1_vector[i] << "\t" << tag2_vector[i] << "\t" <<  data_vector[i] << "\n";
  }
  
  super_matrix_output.close();
  
  ofstream bias_output;

  string bias_file_name = name + ".bias_vector.txt";

  bias_output.open(bias_file_name);

  unordered_map<string,double>::iterator iter11;
  for (iter11 = bias_map.begin(); iter11 != bias_map.end(); iter11++) {
    bias_output << iter11->first << "\t" << bias_map[iter11->first] << "\n";
  }

  bias_output.close();

}

void iterative_correction_no_SR (const string super_matrix,\
				 const string name,\
				 const int iterations,\
				 int inter_flag) {
  
  vector<double> data_vector;
  vector<string> tag1_vector;
  vector<string> tag2_vector;
  unordered_map<string,double> bias_map;
  read_super_matrix_iterative_correction(super_matrix,data_vector,tag1_vector,tag2_vector,bias_map,inter_flag);
  
  vector<double> SR_vector;
  vector<string> SR_tag_vector;
  
  unordered_map<string,double>::iterator it;
  for (it = bias_map.begin(); it != bias_map.end(); it++) {
    SR_vector.push_back(0);
    SR_tag_vector.push_back(it->first);
  }
  
  double mean_b = 1;

  int count = 1;

  while (count <= iterations) {
  
    cerr << "Going through iteration " << count << " now... ";

    //Step 1: Calculate the sum of W over all rows

    unordered_map<string,double> sum_hash;

    for (int i=0; i < data_vector.size(); i++) {
      sum_hash[tag1_vector[i]] += data_vector[i];
      sum_hash[tag2_vector[i]] += data_vector[i];
    }
   
    //Step 2: Calculate the vector of corrected SS reads

    unordered_map<string,double> SR_map;
    
    for (int i=0; i < SR_vector.size(); i++) {
      double val = SR_vector[i]/(mean_b*bias_map[SR_tag_vector[i]]);
      SR_map[SR_tag_vector[i]] = val;
    }

    //Step 3: Calculate the additional vector of biases

    unordered_map<string,double> bias_vector;
    double bias_sum = 0;
    int bias_N = 0;

    unordered_map<string,double>::iterator iter4;
    for (iter4 = bias_map.begin(); iter4 != bias_map.end(); iter4++) {
      bias_vector[iter4->first] = sum_hash[iter4->first] + SR_map[iter4->first];
      if ((sum_hash[iter4->first] + SR_map[iter4->first]) > 0) {
	bias_sum += (sum_hash[iter4->first] + SR_map[iter4->first]);
	bias_N++;
      }
    }

    double bias_mean = bias_sum/bias_N;

    //Step 4: Renormalize bias vector by its mean value to avoid instabilities
    //Step 5: Set zero values of bias vector to 1 to avoid 0/0;

    double var_sum = 0;
    double var_N = 0;

    unordered_map<string,double>::iterator iter5;
    for (iter5 = bias_vector.begin(); iter5 != bias_vector.end(); iter5++) {
      double temp = bias_vector[iter5->first];
      double new_val;
      if (temp == 0) {
	new_val = 1;
      } else {
	new_val = temp/bias_mean;
	var_sum += (temp - bias_mean)*(temp - bias_mean);
	var_N++;
      }
      bias_vector[iter5->first] = new_val;      
    }

    double variance = var_sum/var_N;

    cerr << "variance is " << variance << "\n";

    //Step 6: Divide all Wij by product of bias vectors\n";

    for (int i=0; i < data_vector.size(); i++) {
      double b1 = bias_vector[tag1_vector[i]];
      double b2 = bias_vector[tag2_vector[i]];
      double val = data_vector[i]/(b1*b2);
      data_vector[i] = val;
    }

    //Step 7: Multiple the total vector of biases by additional biases

    double sum_b = 0;
    double N_b = 0;

    unordered_map<string,double>::iterator iter8;
    for (iter8 = bias_map.begin(); iter8 != bias_map.end(); iter8++) {
      double temp = bias_map[iter8->first]*bias_vector[iter8->first];
      bias_map[iter8->first] = temp;
      sum_b += temp;
      N_b++;
    }

    mean_b = sum_b/N_b;

    count++;
  }
  
  ofstream super_matrix_output;
  
  string super_matrix_file_name = name + ".super_matrix.norm.txt";

  super_matrix_output.open(super_matrix_file_name);

  for (int i=0; i < data_vector.size(); i++) {
    super_matrix_output << tag1_vector[i] << "\t" << tag2_vector[i] << "\t" <<  data_vector[i] << "\n";
  }

  super_matrix_output.close();

  ofstream bias_output;

  string bias_file_name = name + ".bias_vector.txt";

  bias_output.open(bias_file_name);

  unordered_map<string,double>::iterator iter11;
  for (iter11 = bias_map.begin(); iter11 != bias_map.end(); iter11++) {
    bias_output << iter11->first << "\t" << bias_map[iter11->first] << "\n";
  }

  bias_output.close();

}

void read_super_matrix_eigen (const string super_matrix_file,		\
			      unordered_map<string,unordered_map<string,double> >& ref_map, \
			      unordered_map<string,int>& def_map) {
  
  cerr << "Reading in super matrix file now... ";

  string line;
  ifstream myfile (super_matrix_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    double n = stod(array[2]);

    ref_map[array[0]][array[1]] = n;

    def_map[array[0]] = 1;
    def_map[array[1]] = 1;

  }

  cerr << "Done\n";

}

void get_eigen (const string super_matrix_file,\
		const string name) {
  
  unordered_map<string,unordered_map<string,double> > data_map;
  unordered_map<string,int> def_map;
  read_super_matrix_eigen(super_matrix_file,data_map,def_map);

  cerr << "Construction matrix now... ";

  vector<string> def_vector;
  unordered_map<string,int>::iterator iter;
  for (iter = def_map.begin(); iter != def_map.end(); iter++) {
    def_vector.push_back(iter->first);
  }

  Eigen::MatrixXd m(def_vector.size(),def_vector.size());
  
  for (int i=0; i < def_vector.size(); i++) {
    for (int j=0; j < def_vector.size(); j++) {
      double val;
      if (data_map.count(def_vector[i]) > 0) {
	if (data_map[def_vector[i]].count(def_vector[j]) > 0) {
	  val = data_map[def_vector[i]][def_vector[j]];
	} else {
	  if (data_map.count(def_vector[j]) > 0) {
	    if (data_map[def_vector[j]].count(def_vector[i]) > 0) {
	      val = data_map[def_vector[j]][def_vector[i]];
	    } else {
	      val = 0;
	    }
	  } else {
	    val = 0;
	  }
	}
      } else {
	if (data_map.count(def_vector[j]) > 0) {
	  if (data_map[def_vector[j]].count(def_vector[i]) > 0) {
	    val = data_map[def_vector[j]][def_vector[i]];
	  } else {
	    val = 0;
	  }
	} else {
	  val = 0;
	}
      }
  
      m(i,j) = val;
      
    }
  }

  cerr << " done\n";

  cerr << "Creating covariance matrix now... ";
  Eigen::MatrixXd cen = m.rowwise() - m.colwise().mean();
  Eigen::MatrixXd cov = (cen.adjoint() * cen) / double(m.rows() - 1);
  cerr << " done\n";

  cerr << "Finding eigenvectors now... ";    
  Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
  cerr << " done\n";

  Eigen::MatrixXd d = m*es.eigenvectors().col(0).real()*es.eigenvectors().col(0).real().transpose(); 

  ofstream eig_index_output;
 
  string eig_file_name = name + ".super_matrix.norm.eig.txt";

  eig_index_output.open(eig_file_name);

  for (int i=0; i < (def_vector.size() - 1); i++) {
    for (int j=(i + 1); j < def_vector.size(); j++) {
      eig_index_output << def_vector[i] << "\t" << def_vector[j] << "\t" << d(i,j) << "\n";
    }
  }

  eig_index_output.close();

}

void read_bias_vector (string bias_vector_file,				\
		       unordered_map<string,vector<int> >& ref_map,	\
                       unordered_map<string,unordered_map<int,double> >& bias_map) {

  std::time_t result = std::time(0);

  cerr << "Reading in bias vector file now... ";

  unordered_map<string,vector<int> > all_data;
  unordered_map<string,unordered_map<int,double> > data;

  string line;
  ifstream myfile (bias_vector_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');
    
    vector<string> loc1;
    split_string(loc1,array[0],':');
    
    vector<string> range1;
    split_string(range1,loc1[1],'-');
    
    int start = string_to_int(range1[0]);
    all_data[loc1[0]].push_back(start);

    double n = stod(array[1]);
    data[loc1[0]][start] = n;
  }

  unordered_map<string,vector<int> >::iterator it;

  for (it = all_data.begin(); it != all_data.end(); it++) {

    vector<int> n_vector = it->second;
    sort(n_vector.begin(),n_vector.end());

    ref_map[it->first] = n_vector;

    for (int i = 0; i < n_vector.size(); i++) {
      double n_val = data[it->first][n_vector[i]];
      bias_map[it->first][n_vector[i]] = n_val;
    }

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void read_super_matrix (string super_matrix_file,   \
			unordered_map<string,unordered_map<string,int> >& ref_map) {

  std::time_t result = std::time(0);

  cerr << "Reading in super matrix file now... ";

  string line;
  ifstream myfile (super_matrix_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    vector<string> loc1;
    split_string(loc1,array[0],':');
    
    vector<string> loc2;
    split_string(loc2,array[1],':');

    vector<string> range1;
    split_string(range1,loc1[1],'-');

    vector<string> range2;
    split_string(range2,loc2[1],'-');

    string tag1 = loc1[0] + "_" + range1[0];
    string tag2 = loc2[0] + "_" + range2[0];

    int n = string_to_int(array[2]);

    ref_map[tag1][tag2] = n;
    ref_map[tag2][tag1] = n;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void read_index_file (string index_file,				\
                      unordered_map<string,unordered_map<string,double> >& index_map) {
  
  std::time_t result = std::time(0);
  
  cerr << "Reading in exp file now... ";

  string line;
  ifstream myfile (index_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    vector<string> loc1;
    split_string(loc1,array[0],':');

    vector<string> loc2;
    split_string(loc2,array[1],':');

    vector<string> range1;
    split_string(range1,loc1[1],'-');

    vector<string> range2;
    split_string(range2,loc2[1],'-');

    string tag1 = loc1[0] + "_" + range1[0];
    string tag2 = loc2[0] + "_" + range2[0];

    double n = stod(array[2]);

    index_map[tag1][tag2] = n;
    index_map[tag2][tag1] = n;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void find_parameters_1Mb (unordered_map<int,double>& m_hash,		\
			  unordered_map<int,double>& r_hash,		\
			  unordered_map<int,double>& var_hash,		\
			  unordered_map<int,int>& N_hash,		\
			  unordered_map<string,unordered_map<string,int> >& data_map, \
			  unordered_map<string,vector<int> >& ref_map,	\
			  unordered_map<string,unordered_map<int,double> >& bias_map, \
			  unordered_map<string,unordered_map<string,double> >& exp_map, \
			  unordered_map<string,unordered_map<string,double> >& pc1_map, \
			  int bin_size,					\
			  int max_size,					\
			  double &inter_m,				\
			  double &inter_r,				\
			  double &submatrix,				\
			  int &total_weight) {

  std::time_t result = std::time(0);

  cerr << "Finding parameters now... ";

  unordered_map<string,vector<int> >::iterator it;

  vector<string> chr_array;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {
    string chr = it->first;
    chr_array.push_back(chr);
  }

  unordered_map<int,double> sum_hash;
  unordered_map<int,double> var_sum_hash;

  double inter_sum = 0;
  double inter_var_sum = 0;
  double inter_N = 0;

  int max_dist = 0;

  for (int u = 0; u < chr_array.size(); u++) {

    string chr1 = chr_array[u];
    vector<int> loc_array1 = ref_map[chr1];

    for (int v = u; v < chr_array.size(); v++) {

      if (v == u) {

        for (int i=0; i < (loc_array1.size() - 1); i++) {
          for (int j=(i + 1); j < loc_array1.size(); j++) {
            string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[i]) )->str();
            string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[j]) )->str();

            double b1 = bias_map[chr1][loc_array1[i]];
            double b2 = bias_map[chr1][loc_array1[j]];

            string tag1 = chr1 + "_" + loc1;
            string tag2 = chr1 + "_" + loc2;

            double n;

            if (data_map.count(tag1) == 1) {
              if (data_map[tag1].count(tag2) == 1) {
                n = data_map[tag1][tag2];
              } else {
                n = 0;
              }
            } else {
              n = 0;
            }

            int dist = loc_array1[j] - loc_array1[i];

            if (dist > max_dist) {
              max_dist = dist;
            }

            double norm = n/(b1*b2);

            sum_hash[dist] += norm;
            N_hash[dist]++;
            if (dist <= max_size) {
              total_weight++;
            }

          }
        }

      } else {

        string chr2 = chr_array[v];
        vector<int> loc_array2 = ref_map[chr2];

        double v1 = loc_array1.size();
        double v2 = loc_array2.size();

        double temp_N = v1*(v1 + 1)*v2*(v2 + 1)/4;

        submatrix += temp_N;

        for (int i=0; i < loc_array1.size(); i++) {
          for (int j=0; j < loc_array2.size(); j++) {

            string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[i]) )->str();
            string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array2[j]) )->str();

            double b1 = bias_map[chr1][loc_array1[i]];
            double b2 = bias_map[chr2][loc_array2[j]];

            string tag1 = chr1 + "_" + loc1;
            string tag2 = chr2 + "_" + loc2;

            double n;

            if (data_map.count(tag1) == 1) {
              if (data_map[tag1].count(tag2) == 1) {
                n = data_map[tag1][tag2];
              } else {
                n = 0;
              }
            } else {
              n = 0;
            }

	    double pc1;

            if (pc1_map.count(tag1) == 1) {
              if (pc1_map[tag1].count(tag2) == 1) {
                pc1 = pc1_map[tag1][tag2];
              } else {
                pc1 = 0;
              }
            } else {
              pc1 = 0;
            }


            double expect;

            if (exp_map.count(tag1) == 1) {
              if (exp_map[tag1].count(tag2) == 1) {
                if (exp_map[tag1][tag2] > 0) {
                  expect = exp_map[tag1][tag2];
                } else {
                  expect = 1;
                }
              } else {
                expect = 1;
              }
            } else {
              expect = 1;
            }

            double norm = (n/(b1*b2) - pc1)/expect;

            inter_sum += norm;
            inter_N++;

          }
        }
      }
    }
  }

  for (int u = 0; u < chr_array.size(); u++) {

    string chr1 = chr_array[u];
    vector<int> loc_array1 = ref_map[chr1];

    for (int v = u; v < chr_array.size(); v++) {

      if (v == u) {

	for (int i=0; i < (loc_array1.size() - 1); i++) {
          for (int j=(i + 1); j < loc_array1.size(); j++) {
            string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[i]) )->str();
            string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[j]) )->str();

            double b1 = bias_map[chr1][loc_array1[i]];
            double b2 = bias_map[chr1][loc_array1[j]];

            string tag1 = chr1 + "_" + loc1;
            string tag2 = chr1 + "_" + loc2;

            double n;

            if (data_map.count(tag1) == 1) {
              if (data_map[tag1].count(tag2) == 1) {
                n = data_map[tag1][tag2];
              } else {
                n = 0;
              }
            } else {
              n = 0;
            }

            int dist = loc_array1[j] - loc_array1[i];

            double norm = n/(b1*b2);

            double val = (norm - sum_hash[dist]/N_hash[dist])*(norm - sum_hash[dist]/N_hash[dist]);

            var_sum_hash[dist] += val;

          }
        }

      }

      else {

        string chr2 = chr_array[v];
        vector<int> loc_array2 = ref_map[chr2];

        for (int i=0; i < loc_array1.size(); i++) {
          for (int j=0; j < loc_array2.size(); j++) {

            string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[i]) )->str();
            string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array2[j]) )->str();

            double b1 = bias_map[chr1][loc_array1[i]];
            double b2 = bias_map[chr2][loc_array2[j]];

            string tag1 = chr1 + "_" + loc1;
            string tag2 = chr2 + "_" + loc2;

            double n;

            if (data_map.count(tag1) == 1) {
              if (data_map[tag1].count(tag2) == 1) {
                n = data_map[tag1][tag2];
              } else {
                n = 0;
              }
            } else {
              n = 0;
            }

	    double pc1;

            if (pc1_map.count(tag1) == 1) {
              if (pc1_map[tag1].count(tag2) == 1) {
                pc1 = pc1_map[tag1][tag2];
              } else {
                pc1 = 0;
              }
            } else {
              pc1 = 0;
            }

            double expect;

            if (exp_map.count(tag1) == 1) {
              if (exp_map[tag1].count(tag2) == 1) {
                if (exp_map[tag1][tag2] > 0) {
                  expect = exp_map[tag1][tag2];
                } else {
                  expect = 1;
                }
              } else {
                expect = 1;
              }
            } else {
              expect = 1;
            }

            double norm = (n/(b1*b2) - pc1)/expect;

            double val = (norm - inter_sum/inter_N)*(norm - inter_sum/inter_N);

            inter_var_sum += val;

          }
        }
      }
    }
  }

  inter_m = inter_sum/inter_N;
  double inter_var = inter_var_sum/inter_N;
  inter_r = inter_m*inter_m/(inter_var - inter_m);

  if (inter_r > 5) {
    inter_r = 5;
  }

  unordered_map<int,double>::iterator it3;

  for (it3 = sum_hash.begin(); it3 != sum_hash.end(); it3++) {

    int dist = it3->first;

    double sum_val;
    int N_val;
    if (sum_hash[it3->first] > 0) {
      sum_val = sum_hash[it3->first];
      N_val =  N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      sum_val = sum_hash[temp_dist];
      N_val = N_hash[temp_dist];
    }
    double var_sum_val;
    int var_sum_N;
    if (var_sum_hash[it3->first] > 0) {
      var_sum_val = var_sum_hash[it3->first];
      var_sum_N = N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (var_sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      var_sum_val = var_sum_hash[temp_dist];
      var_sum_N = N_hash[temp_dist];
    }

    double m_val = sum_val/(double)N_val;
    double var = var_sum_val/(double)var_sum_N;
    double r_val;
    if (var > m_val) {
      r_val = m_val*m_val/(var - m_val);
    } else {
      r_val = 20;
    }

    m_hash[dist] = m_val;
    r_hash[dist] = r_val;
    var_hash[dist] = var;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

double get_corr_map (unordered_map<int,double>& data_map,       \
                     unordered_map<int,int>& N_map) {

  vector<double> val_vector;
  vector<int> N_vector;

  unordered_map<int,double>::iterator iter1;
  for (iter1 = data_map.begin(); iter1 != data_map.end(); iter1++) {
    double mean_val = data_map[iter1->first]/N_map[iter1->first];
    val_vector.push_back(mean_val);
    N_vector.push_back(iter1->first);
  }

  vector<double> sort_val_vector = val_vector;
  vector<int> sort_N_vector = N_vector;

  sort(sort_val_vector.begin(),sort_val_vector.end());
  sort(sort_N_vector.begin(),sort_N_vector.end());

  unordered_map<double,int> val_sort_map;
  unordered_map<double,int> N_sort_map;

  int last = 0;
  for (int i=0; i < sort_val_vector.size(); i++) {
    if (i > 0) {
      if (sort_val_vector[i] == sort_val_vector[i - 1]) {
        val_sort_map[sort_val_vector[i]] = last;
      } else {
        val_sort_map[sort_val_vector[i]] = i;
        last = i;
      }
    } else {
      val_sort_map[sort_val_vector[i]] = i;
    }
  }

  int new_last = 0;
  for (int i=0;i < sort_N_vector.size(); i++) {
    if (i > 0) {
      if (sort_N_vector[i] == sort_N_vector[i - 1]) {
	N_sort_map[sort_N_vector[i]] = new_last;
      } else {
	N_sort_map[sort_N_vector[i]] = i;
        new_last = i;
      }
    } else {
      N_sort_map[sort_N_vector[i]] = i;
    }
  }

  vector<int> x_vector;

  for (int i=0; i < N_vector.size(); i++) {
    x_vector.push_back(N_sort_map[N_vector[i]]);
  }

  vector<int> y_vector;
  for (int i=0; i < val_vector.size(); i++) {
    y_vector.push_back(val_sort_map[val_vector[i]]);
  }

  double x_sum = 0;
  double y_sum = 0;
  for (int i=0; i < y_vector.size(); i++) {
    y_sum += y_vector[i];
    x_sum += x_vector[i];
  }

  double x_mean = x_sum/x_vector.size();
  double y_mean = y_sum/y_vector.size();

  double xy_sum = 0;
  double xx_sum = 0;
  double yy_sum = 0;

  for (int i=0; i < x_vector.size(); i++) {
    xy_sum += (x_vector[i] - x_mean)*(y_vector[i] - y_mean);
    yy_sum += (y_vector[i] - y_mean)*(y_vector[i] - y_mean);
    xx_sum += (x_vector[i] - x_mean)*(x_vector[i] - x_mean);
  }

  double cor;
  if ((xx_sum > 0) && (yy_sum > 0)) {
    cor = xy_sum/(pow(yy_sum,0.5)*pow(xx_sum,0.5));
  } else {
    cor = 0;
  }

  return(cor);

}

double get_corr_vec (vector<double>& val_vector,        \
                     vector<int>& N_vector) {

  vector<double> sort_val_vector = val_vector;
  vector<int> sort_N_vector = N_vector;

  sort(sort_val_vector.begin(),sort_val_vector.end());
  sort(sort_N_vector.begin(),sort_N_vector.end());

  unordered_map<double,int> val_sort_map;
  unordered_map<double,int> N_sort_map;

  int last = 0;
  for (int i=0; i < sort_val_vector.size(); i++) {
    if (i > 0) {
      if (sort_val_vector[i] == sort_val_vector[i - 1]) {
        val_sort_map[sort_val_vector[i]] = last;
      } else {
        val_sort_map[sort_val_vector[i]] = i;
        last = i;
      }
    } else {
      val_sort_map[sort_val_vector[i]] = i;
    }
  }

  int new_last = 0;
  for (int i=0;i < sort_N_vector.size(); i++) {
    if (i > 0) {
      if (sort_N_vector[i] == sort_N_vector[i - 1]) {
        N_sort_map[sort_N_vector[i]] = new_last;
      } else {
        N_sort_map[sort_N_vector[i]] = i;
        new_last = i;
      }
    } else {
      N_sort_map[sort_N_vector[i]] = i;
    }
  }

  vector<int> x_vector;

  for (int i=0; i < N_vector.size(); i++) {
    x_vector.push_back(N_sort_map[N_vector[i]]);
  }

  vector<int> y_vector;
  for (int i=0; i < val_vector.size(); i++) {
    y_vector.push_back(val_sort_map[val_vector[i]]);
  }

  double x_sum = 0;
  double y_sum = 0;
  for (int i=0; i < y_vector.size(); i++) {
    y_sum += y_vector[i];
    x_sum += x_vector[i];
  }

  double x_mean = x_sum/x_vector.size();
  double y_mean = y_sum/y_vector.size();

  double xy_sum = 0;
  double xx_sum = 0;
  double yy_sum = 0;

  for (int i=0; i < x_vector.size(); i++) {
    xy_sum += (x_vector[i] - x_mean)*(y_vector[i] - y_mean);
    yy_sum += (y_vector[i] - y_mean)*(y_vector[i] - y_mean);
    xx_sum += (x_vector[i] - x_mean)*(x_vector[i] - x_mean);
  }

  double cor;
  if ((xx_sum > 0) && (yy_sum > 0)) {
    cor = xy_sum/(pow(yy_sum,0.5)*pow(xx_sum,0.5));
  } else {
    cor = 0;
  }

  return(cor);

}

void find_breaks_1Mb (unordered_map<int,double>& m_hash,		\
		      unordered_map<int,double>& r_hash,		\
		      unordered_map<int,double>& var_hash,		\
		      unordered_map<int,int>& N_hash,			\
		      double &inter_m,					\
		      double &inter_r,					\
		      unordered_map<string,unordered_map<string,int> >& data_map, \
		      unordered_map<string,vector<int> >& ref_map,	\
		      unordered_map<string,unordered_map<int,double> >& bias_map, \
		      unordered_map<string,unordered_map<string,double> >& exp_map, \
		      unordered_map<string,unordered_map<string,double> >& pc1_map, \
		      int bin_size,					\
		      int max_size,					\
		      int total_weight,					\
		      double thresh,					\
		      string name) {
  
  double inter_var = inter_m*inter_m/inter_r + inter_m;

  unordered_map<int,double> fact_hash;
  fact_hash[0] = 0;
  fact_hash[1] = 0;

  std::time_t result = std::time(0);

  cerr << "Finding breaks now...\n";

  ofstream break_output;
  string break_output_file_name = name + ".breaks.txt";
  
  break_output.open(break_output_file_name);

  unordered_map<string,vector<int> >::iterator it;

  vector<string> chr_array;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {
    string chr = it->first;
    chr_array.push_back(chr);
  }

  for (int u = 0; u < (chr_array.size() - 1); u++) {

    string chr1 = chr_array[u];
    vector<int> loc_array1 = ref_map[chr1];

    cerr << chr1 << "\n";

    for (int v = (u + 1); v < chr_array.size(); v++) {

      string chr2 = chr_array[v];
      vector<int> loc_array2 = ref_map[chr2];

      vector<vector<double> > odds_matrix(loc_array1.size(),vector<double>(loc_array2.size()));
      vector<vector<double> > ori_odds_matrix(loc_array1.size(),vector<double>(loc_array2.size()));
      vector<vector<double> > corr_matrix(loc_array1.size(),vector<double>(loc_array2.size(),0));
      vector<vector<double> > norm_matrix(loc_array1.size(),vector<double>(loc_array2.size(),0));

      for (int i=0; i < loc_array1.size(); i++) {
        for (int j=0; j < loc_array2.size(); j++) {

          string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array1[i]) )->str();
          string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array2[j]) )->str();

          double b1 = bias_map[chr1][loc_array1[i]];
          double b2 = bias_map[chr2][loc_array2[j]];

          string tag1 = chr1 + "_" + loc1;
          string tag2 = chr2 + "_" + loc2;

          double n;

          if (data_map.count(tag1) == 1) {
            if (data_map[tag1].count(tag2) == 1) {
              n = data_map[tag1][tag2];
            } else {
              n = 0;
            }
          } else {
            n = 0;
          }

          if (fact_hash.count(n) == 0) {
            int temp_n = n;
            while (fact_hash.count(temp_n) == 0) {
              temp_n--;
            }
            for (int k = (temp_n + 1); k <= n; k++) {
              double t_val = fact_hash[k - 1] + log(k);
              fact_hash[k] = t_val;
            }
          }

	  double pc1;

          if (pc1_map.count(tag1) == 1) {
            if (pc1_map[tag1].count(tag2) == 1) {
              pc1 = pc1_map[tag1][tag2];
            } else {
              pc1 = 0;
            }
          } else {
            pc1 = 0;
          }

          double expect;

          if (exp_map.count(tag1) == 1) {
            if (exp_map[tag1].count(tag2) == 1) {
              if (exp_map[tag1][tag2] > 0) {
                //              expect = 1;                                                                                                                                                                                      
                expect = exp_map[tag1][tag2];
              } else {
                expect = 1;
                //              expect = 1;                                                                                                                                                                                      
              }
            } else {
              expect = 1;
            }
          } else {
            expect = 1;
          }

          double test_m;

          if (b1*b2*((inter_m + corr_matrix[i][j])*expect + pc1) < 0.1*inter_m) {
            test_m = 0.1*inter_m;
          } else {
            test_m = b1*b2*((inter_m + corr_matrix[i][j])*expect + pc1);
          }

          double test_r;

          double term1 = inter_r*log(inter_r/(inter_r + test_m));

          double term2 = 0;

	  for (int k=1; k <= n; k++) {
            double temp_val = log(n + inter_r - k);
            term2 += temp_val;
          }

          term2 += -1*fact_hash[n];

          double term3 = n*log(test_m/(test_m + inter_r));

          double inter_p = term1 + term2 + term3;

          double intra_p = 0;

          unordered_map<int,double>::iterator it0;

          int current_weight = 0;

	  for (it0 = m_hash.begin(); it0 != m_hash.end(); it0++) {

            int dist = it0->first;

            double t_m = m_hash[dist]*b1*b2;

            double t_r;
            if (t_m > var_hash[dist]*b1*b1*b2*b2) {
              t_r = 20;
            } else {
              t_r = t_m*t_m/(var_hash[dist]*b1*b1*b2*b2 - t_m);
            }

            if ((dist <= max_size) &&
                (test_m < t_m)) {

              double t1 = r_hash[dist]*log(r_hash[dist]/(r_hash[dist] + t_m));

              double t2 = 0;

              for (int l=1; l <= n; l++) {
                double temp_val = log(n + r_hash[dist] - l);

                t2 += temp_val;
              }

              t2 += -1*fact_hash[n];

              double t3 = n*log(t_m/(t_m + r_hash[dist]));

              double new_weight = (max_size - dist)/bin_size + 1;

              double temp_p = new_weight*exp(t1 + t2 + t3);
              current_weight += new_weight;

              intra_p += temp_p;

            }

          }

          double final_intra_p = intra_p/current_weight;
          double odds = log(final_intra_p) - inter_p;
          double norm_val = n/(b1*b2);

          if ((exp_map.count(tag1) == 1) &&
              (exp_map[tag1].count(tag2) == 1) &&
              (exp_map[tag1][tag2] > 0)) {

            odds_matrix[i][j] = odds;
            ori_odds_matrix[i][j] = odds;
            norm_matrix[i][j] = norm_val;

          } else {

            odds_matrix[i][j] = 0;
            ori_odds_matrix[i][j] = 0;
            norm_matrix[i][j] = 0;

          }

        }
      }

      int check = 1;
      vector<string> top_hits;
      vector<double> top_sums;
      double t_thresh = thresh;

      while (check == 1) {

        double max_sum = 0;

        int k1;
        int k2;
        int j1;
        int j2;

        for (int i=0; i < loc_array2.size(); i++) {

          double init_val = 0;

          vector<double> temp(loc_array1.size(),init_val);

          for (int j=i; j < loc_array2.size(); j++) {

            for (int k=0; k < loc_array1.size(); k++) {

              temp[k] += odds_matrix[k][j];

            }

            double temp_sum = 0;
            double temp_max = 0;

            int t_k1 = 0;
            int t_k2 = 0;

            int f_k1;
            int f_k2;

	    for (int k=0; k< loc_array1.size(); k++) {

              double new_sum = temp_sum + temp[k];

              if (new_sum < 0) {
                temp_sum = 0;
                t_k1 = k + 1;
              } else {
                temp_sum = new_sum;
                t_k2 = k;
                if (temp_sum > temp_max) {
                  temp_max = temp_sum;
                  f_k1 = t_k1;
                  f_k2 = t_k2;
                }
              }
            }

            if (temp_max > max_sum) {
              max_sum = temp_max;
              j1 = i;
              j2 = j;
              k1 = f_k1;
              k2 = f_k2;
            }

          }

        }

	if (max_sum > t_thresh) {

          t_thresh += thresh;

          for (int k = k1; k <= k2; k++) {

            for (int j = j1; j <= j2; j++) {

              double change = -5;
              odds_matrix[k][j] = change;

            }
          }

          string j_strand = "s";
          string k_strand = "s";

          string s1 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
          string e1 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();
          string s2 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
          string e2 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();

          string new_string = chr1 + "_" + s1 + "_" + e1 + "_" + chr2 + "_" + s2 + "_" + e2;
          top_hits.push_back(new_string);
          top_sums.push_back(max_sum);

        } else {
          break;
        }

      }

      int overlap_check = 0;
      if (top_hits.size() > 1) {

        for (int i=0; i < (top_hits.size() - 1); i++) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[4]);
          int e2 = stoi(new_array1[5]);
          for (int j=(i + 1); j < top_hits.size(); j++) {
            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');
            int s3 = stoi(new_array2[1]);
            int e3 = stoi(new_array2[2]);
            int s4 = stoi(new_array2[4]);
            int e4 = stoi(new_array2[5]);

            int ot1 = 0;
            if (s1 < s3) {
              if (e1 >= s3) {
                ot1 = 1;
              }
            } else {
              if (s1 <= e3) {
                ot1 = 1;
              }
            }

            int ot2 = 0;
            if (s2 < s4) {
              if (e2 >= s4) {
                ot2 = 1;
              }
            } else {
              if (s2 <= e4) {
                ot2 = 1;
              }
            }

	    if ((ot1 == 1) &&
                (ot2 == 1)) {
              overlap_check = 1;
            }

          }
        }

      }

      if (overlap_check == 1) {

        int set = 0;

        while ((set + 1) < top_hits.size()) {

          vector<vector<double> > new_odds_matrix = ori_odds_matrix;
          vector<vector<double> > test_odds_matrix = ori_odds_matrix;

          vector<string> new_top_hits;
          vector<double> new_top_sums;

          if (set > 0) {
            for (int i = 0; i < set; i++) {
              new_top_hits.push_back(top_hits[i]);
              new_top_sums.push_back(top_sums[i]);

              vector<std::string> new_array1;
	      split_string(new_array1,top_hits[i],'_');

              int s1 = stoi(new_array1[1]);
              int e1 = stoi(new_array1[2]);
              int s2 = stoi(new_array1[4]);
              int e2 = stoi(new_array1[5]);

              for (int k = s1; k <= e1; k++) {
                for (int j = s2; j <= e2; j++) {

                  double change = -5;
                  new_odds_matrix[k][j] = change;
                  test_odds_matrix[k][j] = change;

                }
              }

            }
          }

	  for (int i=0; i < top_hits.size(); i++) {

            if (i == (set + 1)) {

              vector<std::string> new_array1;
	      split_string(new_array1,top_hits[i],'_');
              int s1 = stoi(new_array1[1]);
              int e1 = stoi(new_array1[2]);
              int s2 = stoi(new_array1[4]);
              int e2 = stoi(new_array1[5]);

              for (int k = s1; k <= e1; k++) {

                for (int j = s2; j <= e2; j++) {

                  double change = -5;

                  new_odds_matrix[k][j] = change;

                }
              }
            }

          }

          int new_check = 1;
          int pass = 0;

          double t_thresh = thresh;

	  while (new_check == 1) {

            double max_sum = 0;

            int k1;
            int k2;
            int j1;
            int j2;

            for (int i=0; i < loc_array2.size(); i++) {

              double init_val = 0;

              vector<double> temp(loc_array1.size(),init_val);

              for (int j=i; j < loc_array2.size(); j++) {

                for (int k=0; k < loc_array1.size(); k++) {

                  if (pass == 0) {

                    temp[k] += new_odds_matrix[k][j];

                  } else {

                    temp[k] += test_odds_matrix[k][j];

                  }

                }

                double temp_sum = 0;
                double temp_max = 0;

                int t_k1 = 0;
                int t_k2 = 0;

                int f_k1;
                int f_k2;

                for (int k=0; k< loc_array1.size(); k++) {

                  double new_sum = temp_sum + temp[k];

                  if (new_sum < 0) {
                    temp_sum = 0;
                    t_k1 = k + 1;
                  } else {
                    temp_sum = new_sum;
                    t_k2 = k;
                    if (temp_sum > temp_max) {
                      temp_max = temp_sum;
                      f_k1 = t_k1;
                      f_k2 = t_k2;
                    }
                  }
                }

		if (temp_max > max_sum) {
                  max_sum = temp_max;
                  j1 = i;
                  j2 = j;
                  k1 = f_k1;
                  k2 = f_k2;
                }

              }
            }



            if (max_sum > t_thresh) {

              t_thresh += thresh;

              for (int k = k1; k <= k2; k++) {

                for (int j = j1; j <= j2; j++) {

                  double change = -5;
                  test_odds_matrix[k][j] = change;

                }
              }

              string j_strand = "s";
              string k_strand = "s";

              string s1 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
              string e1 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();
              string s2 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
              string e2 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();

              string new_string = chr1 + "_" + s1 + "_" + e1 + "_" + chr2 + "_" + s2 + "_" + e2;
              new_top_hits.push_back(new_string);
              new_top_sums.push_back(max_sum);

	      pass = 1;

            } else {
              break;
            }

          }

	  double all_sum = 0;
          double all_new_sum = 0;

          for (int i = 0; i < top_sums.size(); i++) {
            all_sum += top_sums[i];
          }

          for (int i=0; i < new_top_sums.size(); i++) {
            all_new_sum += new_top_sums[i];
          }

          if (all_new_sum > all_sum) {
            top_hits = new_top_hits;
            top_sums = new_top_sums;
          } else {
            set++;
          }

	  int final_overlap_check = 0;
          for (int i=set; i < (top_hits.size() - 1); i++) {
            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);
            for (int j=(i + 1); j < top_hits.size(); j++) {
              vector<std::string> new_array2;
	      split_string(new_array2,top_hits[j],'_');
              int s3 = stoi(new_array2[1]);
              int e3 = stoi(new_array2[2]);
              int s4 = stoi(new_array2[4]);
              int e4 = stoi(new_array2[5]);

              int ot1 = 0;
              if (s1 < s3) {
                if (e1 >= s3) {
                  ot1 = 1;
                }
              } else {
                if (s1 <= e3) {
                  ot1 = 1;
                }
              }

              int ot2 = 0;
              if (s2 < s4) {
                if (e2 >= s4) {
                  ot2 = 1;
                }
              } else {
                if (s2 <= e4) {
                  ot2 = 1;
                }
              }

              if ((ot1 == 1) &&
                  (ot2 == 1)) {
                final_overlap_check = 1;
              }
            }

          }

          if (final_overlap_check == 0) {
            break;
          }

        }


      }

      int stop_test = 0;

      while (stop_test == 0) {

        if (top_hits.size() > 1) {
          int over_check = 0;
          double over_max = 0;
          int over_i;
          int over_j;

          double all_sum = 0;
          vector<vector<double> > all_odds_matrix = ori_odds_matrix;

          for (int i=0; i < top_hits.size(); i++) {
            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);
            for (int k = s1; k <= e1; k++) {
              for (int l = s2; l <= e2; l++) {
                all_sum += all_odds_matrix[k][l];
                all_odds_matrix[k][l] = -5;
              }
            }
          }

	  for (int i=0; i < (top_hits.size() - 1); i++) {
            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);
            for (int j=(i + 1); j < top_hits.size(); j++) {
              vector<std::string> new_array2;
	      split_string(new_array2,top_hits[j],'_');
              int s3 = stoi(new_array2[1]);
              int e3 = stoi(new_array2[2]);
              int s4 = stoi(new_array2[4]);
              int e4 = stoi(new_array2[5]);
              int new_start1;
              if (s1 < s3) {
                new_start1 = s1;
              } else {
                new_start1 = s3;
              }
              int new_end1;
              if (e1 > e3) {
                new_end1 = e1;
              } else {
                new_end1 = e3;
              }
              int new_start2;
              if (s2 < s4) {
                new_start2 = s2;
              } else {
                new_start2 = s4;
              }
              int new_end2;
              if (e2 > e4) {
                new_end2 = e2;
              } else {
                new_end2 = e4;
              }

              vector<vector<double> > all_new_odds_matrix = ori_odds_matrix;
              double all_new_sum = 0;
              double odds_sum = 0;

              for (int k = new_start1; k <= new_end1; k++) {
                for (int l = new_start2; l <= new_end2; l++) {
                  all_new_sum += all_new_odds_matrix[k][l];
                  all_new_odds_matrix[k][l] = 0;
                  odds_sum += ori_odds_matrix[k][l];
                }
              }

              for (int k=0; k < top_hits.size(); k++) {
                if ((k == i) ||
                    (k == j)) {
                  continue;
                }

                vector<std::string> new_array3;
		split_string(new_array3,top_hits[k],'_');
                int s5 = stoi(new_array3[1]);
                int e5 = stoi(new_array3[2]);
                int s6 = stoi(new_array3[4]);
                int e6 = stoi(new_array3[5]);
                for (int l = s5; l <= e5; l++) {
                  for (int m = s6; m <= e6; m++) {
                    all_new_sum += all_new_odds_matrix[l][m];
                    all_new_odds_matrix[l][m] = 0;
                  }
                }
              }

	      if ((all_new_sum > over_max) &&
                  (odds_sum >= thresh)) {
                over_max = all_new_sum;
                over_i = i;
                over_j = j;
                over_check++;
              }
            }
          }

	  if (over_check > 0) {
            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[over_i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);
            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[over_j],'_');
            int s3 = stoi(new_array2[1]);
            int e3 = stoi(new_array2[2]);
            int s4 = stoi(new_array2[4]);
            int e4 = stoi(new_array2[5]);
            string new_start1;
            if (s1 <= s3) {
              new_start1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
            } else {
              new_start1 = static_cast<ostringstream*>( &(ostringstream() << s3) )->str();
            }
            string new_end1;
            if (e1 >= e3) {
              new_end1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
            } else {
              new_end1 = static_cast<ostringstream*>( &(ostringstream() << e3) )->str();
            }
            string new_start2;
            if (s2 <= s4) {
              new_start2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
            } else {
              new_start2 = static_cast<ostringstream*>( &(ostringstream() << s4) )->str();
            }
            string new_end2;
            if (e2 >= e4) {
              new_end2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();
            } else {
              new_end2 = static_cast<ostringstream*>( &(ostringstream() << e4) )->str();
            }
            string new_string = new_array1[0] + "_" + new_start1 + "_" + new_end1 + "_" + new_array1[3] + "_" + new_start2 + "_" + new_end2;
            vector<string> new_top_hits;
            new_top_hits.push_back(new_string);

	    for (int i=0; i < top_hits.size(); i++) {
              if ((i != over_i) &&
                  (i != over_j)) {
                new_top_hits.push_back(top_hits[i]);
              }
            }
            top_hits = new_top_hits;
          } else {
            stop_test = 1;
          }
        } else {
          stop_test = 1;
        }
      }

      for (int i=0; i < top_hits.size(); i++) {
        vector<std::string> new_array1;
	split_string(new_array1,top_hits[i],'_');
        int s1 = stoi(new_array1[1]);
        int e1 = stoi(new_array1[2]);
        int s2 = stoi(new_array1[4]);
        int e2 = stoi(new_array1[5]);

        unordered_map<int,double> i_sum;
        unordered_map<int,double> j_sum;
        unordered_map<int,int> i_N;
        unordered_map<int,int> j_N;
        vector<double> all_diag1;
        vector<double> all_diag2;
        vector<int> all_diag1_N;
        vector<int> all_diag2_N;
        unordered_map<int,double> diag1_sum;
        unordered_map<int,double> diag2_sum;
        unordered_map<int,int> diag1_N;
        unordered_map<int,int> diag2_N;
        unordered_map<int,double> i_log_sum;
        unordered_map<int,double> j_log_sum;
        unordered_map<int,double> diag1_log_sum;
        unordered_map<int,double> diag2_log_sum;
        unordered_map<int,double> diag1_bin_sum;
        unordered_map<int,double> diag2_bin_sum;

        double odds_sum = 0;
        for (int k=s1; k <= e1; k++) {

          int i_bin = k - s1;
          i_N[i_bin]++;

          for (int j=s2; j <= e2; j++) {
            odds_sum += ori_odds_matrix[k][j];

            int j_bin = j - s2;
            int rev_j_bin = e2 - j;
            j_N[j_bin]++;

            int diag1_bin = i_bin + j_bin;
            int diag2_bin = i_bin + rev_j_bin;

            diag1_N[diag1_bin]++;
            diag2_N[diag2_bin]++;

            all_diag1_N.push_back(diag1_bin);
            all_diag2_N.push_back(diag2_bin);

            i_sum[i_bin] += norm_matrix[k][j];
            j_sum[j_bin] += norm_matrix[k][j];
            i_log_sum[i_bin] += log(norm_matrix[k][j] + 0.1);
            j_log_sum[j_bin] += log(norm_matrix[k][j] + 0.1);
            diag1_sum[diag1_bin] += norm_matrix[k][j];
            diag2_sum[diag2_bin] += norm_matrix[k][j];
            diag1_log_sum[diag1_bin] += log(norm_matrix[k][j] + 0.1);
            diag2_log_sum[diag2_bin] += log(norm_matrix[k][j] + 0.1);
            diag1_bin_sum[diag1_bin] += 1;
            diag2_bin_sum[diag2_bin] += 1;
            all_diag1.push_back(norm_matrix[k][j]);
            all_diag2.push_back(norm_matrix[k][j]);

          }
        }

	double i_cor = get_corr_map(i_sum,i_N);
        double j_cor = get_corr_map(j_sum,j_N);
        double i_log_cor = get_corr_map(i_log_sum,i_N);
        double j_log_cor = get_corr_map(j_log_sum,j_N);
        double diag1_cor = get_corr_map(diag1_sum,diag1_N);
        double diag2_cor = get_corr_map(diag2_sum,diag2_N);
        double diag1_log_cor = get_corr_map(diag1_log_sum,diag1_N);
        double diag2_log_cor = get_corr_map(diag2_log_sum,diag2_N);
        double diag1_bin_cor = get_corr_map(diag1_bin_sum,diag1_N);
        double diag2_bin_cor = get_corr_map(diag2_bin_sum,diag2_N);
        double all_diag1_cor = get_corr_vec(all_diag1,all_diag1_N);
        double all_diag2_cor = get_corr_vec(all_diag2,all_diag2_N);

        double final_i = all_diag1_cor + all_diag2_cor + diag1_bin_cor + diag2_bin_cor + diag1_cor + diag2_cor + diag1_log_cor + diag2_log_cor + i_cor + i_log_cor;
        double final_j = all_diag1_cor + -1*all_diag2_cor + diag1_bin_cor + -1*diag2_bin_cor + diag1_cor + -1*diag2_cor + diag1_log_cor + -1*diag2_log_cor + j_cor + j_log_cor;

        string strand1;
        string strand2;

        if (final_i > 0) {
          strand1 = "+";
        } else {
          strand1 = "-";
        }

        if (final_j > 0) {
          strand2 = "+";
        } else {
          strand2 = "-";
        }

        break_output << odds_sum << "\t" << chr1 << "\t" << loc_array1[s1] << "\t" << (loc_array1[e1] + bin_size) << "\t" << strand1 << "\t" << chr2 << "\t" << loc_array2[s2] << "\t" << (loc_array2[e2] + bin_size) << "\t" << strand2 << "\n";
	
      }

    }
  }

  break_output.close();

}

void get_1Mb_odds_ratio_inter (string super_matrix_file,\
			       string bias_vector_file,\
			       char* exp_file,	       \
			       string pc1_file,	       \
			       int bin_size,\
			       int max_size,\
			       string name) {
  
  unordered_map<string,vector<int> > ref_map;
  unordered_map<string,unordered_map<int,double> > bias_map;
  read_bias_vector(bias_vector_file,ref_map,bias_map);

  unordered_map<string,unordered_map<string,int> > data_map;
  read_super_matrix(super_matrix_file,data_map);
  
  unordered_map<string,unordered_map<string,double> > pc1_map;
  read_index_file(pc1_file,pc1_map);

  unordered_map<string,unordered_map<string,double> > exp_map;
  read_index_file(exp_file,exp_map);

  unordered_map<int,double> m_hash;
  unordered_map<int,double> r_hash;
  unordered_map<int,double> var_hash;
  unordered_map<int,int> N_hash;
  double inter_m;
  double inter_r;
  double submatrix;
  int total_weight;
  find_parameters_1Mb(m_hash,r_hash,var_hash,N_hash,data_map,ref_map,bias_map,exp_map,pc1_map,bin_size,max_size,inter_m,inter_r,submatrix,total_weight);

  double prior = 0.000001;
  double corr = log(prior/(1 - prior));
  double thresh = -1*(log(0.05/submatrix) + corr);

  cerr << "Thresh is " << thresh << "\n";

  find_breaks_1Mb(m_hash,r_hash,var_hash,N_hash,inter_m,inter_r,data_map,ref_map,bias_map,exp_map,pc1_map,bin_size,max_size,total_weight,thresh,name);

}

vector<string> read_blocks (string block_file,   \
			    int max_size) {
  
  vector<string> temp_array;
  
  string line;
  ifstream myfile (block_file);
  while ( getline(myfile,line) ) {
    temp_array.push_back(line);
  }
  
  int i=0;
  
  while (i < temp_array.size()) {

    vector<string> array1;
    split_string(array1,temp_array[i],'\t');

    string chr1 = array1[1];
    string chr2 = array1[5];
    
    int start1;
    if ((stoi(array1[2]) - max_size) > 0) {
      start1 = stoi(array1[2]) - max_size;
    } else {
      start1 = 0;
    }
    int end1 = stoi(array1[3]) + max_size;
    int start2;
    if ((stoi(array1[6]) - max_size) > 0) {
      start2 = stoi(array1[6]) - max_size;
    } else {
      start2 = 0;
    }
    int end2 = stoi(array1[7]) + max_size;

    int test = 0;

    for (int j=(i + 1); j < temp_array.size(); j++) {

      vector<string> array2;
      split_string(array2,temp_array[j],'\t');

      string chr3 = array2[1];
      string chr4 = array2[5];

      int start3;
      if ((stoi(array2[2]) - max_size) > 0) {
	start3 = stoi(array2[2]) - max_size;
      } else {
	start3 = 0;
      }
      int end3 = stoi(array2[3]) + max_size;
      int start4;
      if ((stoi(array2[6]) - max_size) > 0) {
	start4 = stoi(array2[6]) - max_size;
      } else {
	start4 = 0;
      }
      int end4 = stoi(array2[7]) + max_size;

      int over_check = 0;

      if ((chr1 == chr3) &&
          (chr2 == chr4)) {

        if (start1 < start3) {
          if (end1 > start3) {
            if (start2 < start4) {
              if (end2 > start4) {
                over_check = 1;
              }
            } else {
              if (start2 < end4) {
		over_check = 1;
              }
            }
          }
        } else {
          if (start1 < end3) {
            if (start2 < start4) {
              if (end2 > start4) {
                over_check = 1;
              }
            } else {
              if (start2 < end4) {
                over_check = 1;
              }
            }
          }
        }
	
	if (over_check == 1) {
          test = 1;
	  
          vector<string> new_temp_array;
	  
          for (int k=0; k < temp_array.size(); k++) {
            if ((k != i) &&
                (k != j)) {
              new_temp_array.push_back(temp_array[k]);
            }
          }
          string new_start1;
          if (start1 < start3) {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << start1) )->str();
          } else {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << start3) )->str();
          }
	  
          string new_end1;
          if (end1 > end3) {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << end1) )->str();
          } else {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << end3) )->str();
          }
	  
	  string new_start2;
          if (start2 < start4) {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << start2) )->str();
          } else {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << start4) )->str();
          }
	  
          string new_end2;
          if (end2 > end4) {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << end2) )->str();
          } else {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << end4) )->str();
          }
	  
          string new_string = "-1\t" + chr1 + "\t" + new_start1 + "\t" + new_end1 + "\t-1\t" + chr2 + "\t" + new_start2 + "\t" + new_end2 + "\t-1";
	  
          new_temp_array.push_back(new_string);
          temp_array = new_temp_array;
          break;
        }
      }
      
    }
    
    if (test == 0) {
      i++;
    } else {
      i=0;
    }
    
  }

  return(temp_array);

}

void find_breaks_100kb_inter (vector<string>& block_data,		\
			      unordered_map<int,double>& m_hash,	\
			      unordered_map<int,double>& r_hash,	\
			      unordered_map<int,int>& N_hash,		\
			      double &inter_m,				\
			      double &inter_r,				\
			      unordered_map<string,unordered_map<string,int> >& data_map, \
			      unordered_map<string,vector<int> >& ref_map, \
			      unordered_map<string,unordered_map<int,double> >& bias_map, \
			      int bin_size,				\
			      int max_size,				\
			      int total_weight,				\
			      double thresh,				\
			      string name) {
  
  unordered_map<int,double> fact_hash;
  fact_hash[0] = 0;
  fact_hash[1] = 0;

  std::time_t result = std::time(0);

  cerr << "Finding breaks now...\n";

  ofstream break_output;
  string break_output_file_name = name + ".breaks.txt";

  break_output.open(break_output_file_name);

  for (int pod = 0; pod < block_data.size(); pod++) {

    string line = block_data[pod];

    cerr << "Going through\t" << line << "\n";
    vector<string> array;
    split_string(array,line,'\t');

    string chr1 = array[1];
    string chr2 = array[5];

    vector<int> loc_array1 = ref_map[chr1];
    vector<int> loc_array2 = ref_map[chr2];

    int start1 = stoi(array[2]);
    int end1 = stoi(array[3]);
    int start2 = stoi(array[6]);
    int end2 = stoi(array[7]);

    int i_size = (end1 - start1)/bin_size;
    int j_size = (end2 - start2)/bin_size;

    double submatrix = i_size*(i_size + 1)*j_size*(j_size + 1)/4;
    double prior = 0.000001;
    double corr = log(prior/(1 - prior));

    vector<vector<double> > odds_matrix(i_size,vector<double>(j_size));
    vector<vector<double> > ori_odds_matrix(i_size,vector<double>(j_size));
    vector<vector<double> > norm_matrix(i_size,vector<double>(j_size));

    for (int i=start1; i < end1; i += bin_size) {
      for (int j=start2; j < end2; j += bin_size) {

        string loc1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        string loc2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();

        double b1;
        if (bias_map[chr1].count(i) > 0) {
          b1 = bias_map[chr1][i];
        } else {
          b1 = 1;
        }
        double b2;
        if (bias_map[chr2].count(j) > 0) {
          b2 = bias_map[chr2][j];
        } else {
          b2 = 1;
        }

        string tag1 = chr1 + "_" + loc1;
        string tag2 = chr2 + "_" + loc2;

        double n;

	if (data_map.count(tag1) == 1) {
          if (data_map[tag1].count(tag2) == 1) {
            n = data_map[tag1][tag2];
          } else {
            n = 0;
          }
        } else {
          n = 0;
        }

        double norm = n/(b1*b2);

        if (fact_hash.count(n) == 0) {
          int temp_n = n;
          while (fact_hash.count(temp_n) == 0) {
            temp_n--;
          }
          for (int k = (temp_n + 1); k <= n; k++) {
            double t_val = fact_hash[k - 1] + log(k);
            fact_hash[k] = t_val;
          }
        }

        double test_m;

        if (b1*b2 < 0.1) {
          test_m = 0.1*inter_m;
        } else {
          test_m = b1*b2*inter_m;
        }

	double term1 = inter_r*log(inter_r/(inter_r + test_m));

        double term2 = 0;

        for (int k=1; k <= n; k++) {
          double temp_val = log(n + inter_r - k);
          term2 += temp_val;
        }

        term2 += -1*fact_hash[n];

        double term3 = n*log(test_m/(test_m + inter_r));

        double inter_p = term1 + term2 + term3;

        double intra_p = 0;

        unordered_map<int,double>::iterator it0;

        int current_weight = 0;
 
	for (it0 = m_hash.begin(); it0 != m_hash.end(); it0++) {

          int dist = it0->first;

          double t_m = m_hash[dist]*b1*b2;

          if ((dist <= max_size) &&
              (test_m < t_m)) {

            double t1 = r_hash[dist]*log(r_hash[dist]/(r_hash[dist] + t_m));

            double t2 = 0;

            for (int l=1; l <= n; l++) {
              double temp_val = log(n + r_hash[dist] - l);
              t2 += temp_val;
            }

            t2 += -1*fact_hash[n];

            double t3 = n*log(t_m/(t_m + r_hash[dist]));

            double temp_p = N_hash[dist]*exp(t1 + t2 + t3);

            current_weight += N_hash[dist];

            intra_p += temp_p;

          }

	}

        double final_intra_p = intra_p/current_weight;

        double odds = log(final_intra_p) - inter_p;

        int i_bin = (i - start1)/bin_size;
        int j_bin = (j - start2)/bin_size;

        if ((bias_map[chr1].count(i) > 0) &&
            (bias_map[chr2].count(j) > 0)) {

          odds_matrix[i_bin][j_bin] = odds;
          ori_odds_matrix[i_bin][j_bin] = odds;
          norm_matrix[i_bin][j_bin] = norm;

        } else {

          odds_matrix[i_bin][j_bin] = 0;
          ori_odds_matrix[i_bin][j_bin] = 0;
          norm_matrix[i_bin][j_bin] = 0;

        }

      }
    }

    vector<string> top_hits;
    vector<double> top_sums;
    int check = 1;
    double t_thresh = thresh;
    while (check == 1) {

      double max_sum = 0;

      int k1;
      int k2;
      int j1;
      int j2;

      for (int i=0; i < odds_matrix.size(); i++) {

        double init_val = 0;

        vector<double> temp(odds_matrix[0].size(),init_val);

	for (int j=i; j < odds_matrix.size(); j++) {

          for (int k=0; k < odds_matrix[0].size(); k++) {

            temp[k] += odds_matrix[j][k];

          }

          double temp_sum = 0;
          double temp_max = 0;

          int t_k1 = 0;
          int t_k2 = 0;

          int f_k1;
          int f_k2;

	  for (int k=0; k < odds_matrix[0].size(); k++) {

            double new_sum = temp_sum + temp[k];

            if (new_sum < 0) {
              temp_sum = 0;
              t_k1 = k + 1;
            } else {
              temp_sum = new_sum;
              t_k2 = k;
              if (temp_sum > temp_max) {
                temp_max = temp_sum;
                f_k1 = t_k1;
                f_k2 = t_k2;
              }
            }
          }

          if (temp_max > max_sum) {
            max_sum = temp_max;
            j1 = i;
            j2 = j;
            k1 = f_k1;
            k2 = f_k2;

          }

        }

      }

      if (max_sum > t_thresh) {

        t_thresh += thresh;

        for (int k = k1; k <= k2; k++) {

          for (int j = j1; j <= j2; j++) {

            double change = -5;
            odds_matrix[j][k] = change;

          }
        }

        string s1 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
        string e1 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();
        string s2 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
        string e2 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();

        string new_string = chr1 + "_" + s1 + "_" + e1 + "_" + chr2 + "_" + s2 + "_" + e2;

        top_hits.push_back(new_string);
        top_sums.push_back(max_sum);

        string j_strand = "s";
        string k_strand = "s";

      } else {
        break;
      }

    }

    int overlap_check = 0;
    if (top_hits.size() > 1) {

      for (int i=0; i < (top_hits.size() - 1); i++) {
        vector<std::string> new_array1;
	split_string(new_array1,top_hits[i],'_');
        int s1 = stoi(new_array1[1]);
        int e1 = stoi(new_array1[2]);
        int s2 = stoi(new_array1[4]);
        int e2 = stoi(new_array1[5]);
        for (int j=(i + 1); j < top_hits.size(); j++) {
          vector<std::string> new_array2;
	  split_string(new_array2,top_hits[j],'_');
          int s3 = stoi(new_array2[1]);
          int e3 = stoi(new_array2[2]);
          int s4 = stoi(new_array2[4]);
          int e4 = stoi(new_array2[5]);

	  int ot1 = 0;
          if (s1 < s3) {
            if (e1 >= s3) {
              ot1 = 1;
            }
          } else {
            if (s1 <= e3) {
              ot1 = 1;
            }
          }

          int ot2 = 0;
          if (s2 < s4) {
            if (e2 >= s4) {
              ot2 = 1;
            }
          } else {
            if (s2 <= e4) {
              ot2 = 1;
            }
          }

          if ((ot1 == 1) &&
              (ot2 == 1)) {
            overlap_check = 1;
          }

        }

      }

    }
    
    if (overlap_check == 1) {

      int set = 0;

      while ((set + 1) < top_hits.size()) {

        vector<vector<double> > new_odds_matrix = ori_odds_matrix;
        vector<vector<double> > test_odds_matrix = ori_odds_matrix;

        vector<string> new_top_hits;
        vector<double> new_top_sums;

        if (set > 0) {
          for (int i = 0; i < set; i++) {
            new_top_hits.push_back(top_hits[i]);
            new_top_sums.push_back(top_sums[i]);

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');

            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);

	    for (int k = s2; k <= e2; k++) {
              for (int j = s1; j <= e1; j++) {

                double change = -5;
                new_odds_matrix[j][k] = change;
                test_odds_matrix[j][k] = change;

              }
            }

          }
        }

        for (int i=0; i < top_hits.size(); i++) {

          if (i == (set + 1)) {

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
	    
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);

            for (int k = s2; k <= e2; k++) {
              for (int j = s1; j <= e1; j++) {

                double change = -5;
                new_odds_matrix[j][k] = change;

              }
            }
          }
        }

	int new_check = 1;
	int pass = 0;

        while (new_check == 1) {

          double max_sum = 0;

          int k1;
          int k2;
          int j1;
          int j2;

          for (int i=0; i < odds_matrix.size(); i++) {

            double init_val = 0;

            vector<double> temp(odds_matrix[0].size(),init_val);

	    for (int j=i; j < odds_matrix.size(); j++) {

	      for (int k=0; k < odds_matrix[0].size(); k++) {

                if (pass == 0) {

                  temp[k] += new_odds_matrix[j][k];

                } else {

                  temp[k] += test_odds_matrix[j][k];

                }

              }

              double temp_sum = 0;
              double temp_max = 0;

              int t_k1 = 0;
              int t_k2 = 0;

              int f_k1;
              int f_k2;

	      for (int k=0; k < odds_matrix[0].size(); k++) {

                double new_sum = temp_sum + temp[k];

                if (new_sum < 0) {
                  temp_sum = 0;
                  t_k1 = k + 1;
                } else {
                  temp_sum = new_sum;
                  t_k2 = k;
                  if (temp_sum > temp_max) {
                    temp_max = temp_sum;
                    f_k1 = t_k1;
                    f_k2 = t_k2;
                  }
                }
              }

              if (temp_max > max_sum) {
                max_sum = temp_max;
                j1 = i;
                j2 = j;
                k1 = f_k1;
                k2 = f_k2;

              }
            }
          }

	  if (max_sum > t_thresh) {

            t_thresh += thresh;

            for (int k = k1; k <= k2; k++) {

              for (int j = j1; j <= j2; j++) {

                double change = -5;
                test_odds_matrix[j][k] = change;

              }
            }

            string s1 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
            string e1 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();
            string s2 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
            string e2 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();

            string new_string = chr1 + "_" + s1 + "_" + e1 + "_" + chr2 + "_" + s2 + "_" + e2;

            string j_strand = "s";
            string k_strand = "s";

	    new_top_hits.push_back(new_string);
            new_top_sums.push_back(max_sum);

            pass = 1;

          } else {
            break;
          }
	}

	double all_sum = 0;
	double all_new_sum = 0;
	
	for (int i = 0; i < top_sums.size(); i++) {
	  all_sum += top_sums[i];
	}
	
	for (int i=0; i < new_top_sums.size(); i++) {
	  all_new_sum += new_top_sums[i];
	}
	
	if (all_new_sum > all_sum) {
	  top_hits = new_top_hits;
	  top_sums = new_top_sums;
	  
	} else {
	  set++;
	}
	
	int final_overlap_check = 0;
	for (int i=set; i < (top_hits.size() - 1); i++) {
	  vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');

	  int s1 = stoi(new_array1[1]);
	  int e1 = stoi(new_array1[2]);
	  int s2 = stoi(new_array1[4]);
	  int e2 = stoi(new_array1[5]);
	  for (int j=(i + 1); j < top_hits.size(); j++) {
	    vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');

	    int s3 = stoi(new_array2[1]);
	    int e3 = stoi(new_array2[2]);
	    int s4 = stoi(new_array2[4]);
	    int e4 = stoi(new_array2[5]);
	    
	    int ot1 = 0;
	    if (s1 < s3) {
	      if (e1 >= s3) {
		ot1 = 1;
	      }
	    } else {
	      if (s1 <= e3) {
		ot1 = 1;
	      }
	    }
	    
	    int ot2 = 0;
	    if (s2 < s4) {
	      if (e2 >= s4) {
		ot2 = 1;
	      }
	    } else {
	      if (s2 <= e4) {
		ot2 = 1;
	      }
	    }
	    
	    if ((ot1 == 1) &&
		(ot2 == 1)) {
	      final_overlap_check = 1;
	    }
	  }
	}
	
	if (final_overlap_check == 0) {
	  break;
	}
	
      }
      
    }
    
    int stop_test = 0;
    
    while (stop_test == 0) {
      if (top_hits.size() > 1) {
	int over_check = 0;
	double over_max = 0;
	double over_min = -10000000;
	int over_i;
	int over_j;

	double all_sum = 0;
	vector<vector<double> > all_odds_matrix = ori_odds_matrix;
	
	for (int i=0; i < top_hits.size(); i++) {
	  vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');

	  int s1 = stoi(new_array1[1]);
	  int e1 = stoi(new_array1[2]);
	  int s2 = stoi(new_array1[4]);
	  int e2 = stoi(new_array1[5]);
	  for (int k = s1; k <= e1; k++) {
	    for (int l = s2; l <= e2; l++) {
	      all_sum += all_odds_matrix[k][l];
	      all_odds_matrix[k][l] = -5;
	    }
	  }
	}
	
	for (int i=0; i < (top_hits.size() - 1); i++) {
	  vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');

	  int s1 = stoi(new_array1[1]);
	  int e1 = stoi(new_array1[2]);
	  int s2 = stoi(new_array1[4]);
	  int e2 = stoi(new_array1[5]);

	  for (int j=(i + 1); j < top_hits.size(); j++) {
	    vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');

	    int s3 = stoi(new_array2[1]);
	    int e3 = stoi(new_array2[2]);
	    int s4 = stoi(new_array2[4]);
	    int e4 = stoi(new_array2[5]);
	    int new_start1;
	    if (s1 < s3) {
	      new_start1 = s1;
	    } else {
	      new_start1 = s3;
	    }
	    int new_end1;
	    if (e1 > e3) {
	      new_end1 = e1;
	    } else {
	      new_end1 = e3;
	    }
	    int new_start2;
	    if (s2 < s4) {
	      new_start2 = s2;
	    } else {
	      new_start2 = s4;
	    }
	    
	    int new_end2;
	    if (e2 > e4) {
	      new_end2 = e2;
	    } else {
	      new_end2 = e4;
	    }
	    
	    vector<vector<double> > all_new_odds_matrix = ori_odds_matrix;
	    
	    double all_new_sum = 0;
	    double odds_sum = 0;
	    for (int k = new_start1; k <= new_end1; k++) {
	      for (int l = new_start2; l <= new_end2; l++) {
		all_new_sum += all_new_odds_matrix[k][l];
		all_new_odds_matrix[k][l] = 0;
		odds_sum += ori_odds_matrix[k][l];
	      }
	    }
	    
	    for (int k=0; k < top_hits.size(); k++) {
	      if ((k == i) ||
		  (k == j)) {
		continue;
	      }
	      vector<std::string> new_array3;
	      split_string(new_array3,top_hits[k],'_');

	      int s5 = stoi(new_array3[1]);
	      int e5 = stoi(new_array3[2]);
	      int s6 = stoi(new_array3[4]);
	      int e6 = stoi(new_array3[5]);
	      for (int l = s5; l <= e5; l++) {
		for (int m = s6; m <= e6; m++) {
		  all_new_sum += all_new_odds_matrix[l][m];
		  all_new_odds_matrix[l][m] = 0;
		}
	      }
	    }
	    
	    double loss = all_new_sum - all_sum;
	    
	    if ((all_new_sum >= over_max) &&
		(odds_sum >= thresh)) {
	      over_min = loss;
	      over_max = all_new_sum;
	      over_i = i;
	      over_j = j;
	      over_check++;
	    }
	  }
	}
	
	if (over_check > 0) {
	  vector<std::string> new_array1;
	  split_string(new_array1,top_hits[over_i],'_');
	  int s1 = stoi(new_array1[1]);
	  int e1 = stoi(new_array1[2]);
	  int s2 = stoi(new_array1[4]);
	  int e2 = stoi(new_array1[5]);
	  vector<std::string> new_array2;
	  split_string(new_array2,top_hits[over_j],'_');
	  int s3 = stoi(new_array2[1]);
	  int e3 = stoi(new_array2[2]);
	  int s4 = stoi(new_array2[4]);
	  int e4 = stoi(new_array2[5]);
	  int ns1;
	  string new_start1;
	  if (s1 <= s3) {
	    new_start1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
	    ns1 = s1;
	  } else {
	    new_start1 = static_cast<ostringstream*>( &(ostringstream() << s3) )->str();
	    ns1 = s3;
	  }
	  int ne1;
	  string new_end1;
	  if (e1 >= e3) {
	    new_end1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
	    ne1 = e1;
	  } else {
	    new_end1 = static_cast<ostringstream*>( &(ostringstream() << e3) )->str();
	    ne1 = e3;
	  }
	  int ns2;
	  string new_start2;
	  if (s2 <= s4) {
	    new_start2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
	    ns2 = s2;
	  } else {
	    new_start2 = static_cast<ostringstream*>( &(ostringstream() << s4) )->str();
	    ns2 = s4;
	  }
	  int ne2;
	  string new_end2;
	  if (e2 >= e4) {
	    new_end2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();
	    ne2 = e2;
	  } else {
	    new_end2 = static_cast<ostringstream*>( &(ostringstream() << e4) )->str();
	    ne2 = e4;
	  }
	  double new_sum = 0;
	  for (int k = ns1; k <= ne1; k++) {
	    for (int l = ns2; l <= ne2; l++) {
	      odds_matrix[k][l] = -5;
	      new_sum += ori_odds_matrix[k][l];
	    }
	  }
	  string new_string = new_array1[0] + "_" + new_start1 + "_" + new_end1 + "_" + new_array1[3] + "_" + new_start2 + "_" + new_end2;
	  
	  string j_strand = "s";
	  string k_strand = "s";
	  
	  vector<string> new_top_hits;
	  vector<double> new_top_sums;
	  new_top_hits.push_back(new_string);
	  new_top_sums.push_back(new_sum);
	  for (int i=0; i < top_hits.size(); i++) {
	    if ((i != over_i) &&
		(i != over_j)) {
	      new_top_hits.push_back(top_hits[i]);
	      new_top_sums.push_back(top_sums[i]);
	    }
	  }
	  top_hits = new_top_hits;
	  top_sums = new_top_sums;
	} else {
	  stop_test = 1;
	}
      } else {
	stop_test = 1;
      }
    }
    for (int i=0; i < top_hits.size(); i++) {
      vector<std::string> new_array1;
      split_string(new_array1,top_hits[i],'_');
      int s1 = stoi(new_array1[1]);
      int e1 = stoi(new_array1[2]);
      int s2 = stoi(new_array1[4]);
      int e2 = stoi(new_array1[5]);
      
      unordered_map<int,double> i_sum;
      unordered_map<int,double> j_sum;
      unordered_map<int,int> i_N;
      unordered_map<int,int> j_N;
      vector<double> all_diag1;
      vector<double> all_diag2;
      vector<int> all_diag1_N;
      vector<int> all_diag2_N;
      unordered_map<int,double> diag1_sum;
      unordered_map<int,double> diag2_sum;
      unordered_map<int,int> diag1_N;
      unordered_map<int,int> diag2_N;
      unordered_map<int,double> i_log_sum;
      unordered_map<int,double> j_log_sum;
      unordered_map<int,double> diag1_log_sum;
      unordered_map<int,double> diag2_log_sum;
      unordered_map<int,double> diag1_bin_sum;
      unordered_map<int,double> diag2_bin_sum;
      
      double odds_sum = 0;
      for (int k=s1; k <= e1; k++) {
	
	int i_bin = k - s1;
	i_N[i_bin]++;
	
	for (int j=s2; j <= e2; j++) {
	  odds_sum += ori_odds_matrix[k][j];
	  
	  int j_bin = j - s2;
	  int rev_j_bin = e2 - j;
	  j_N[j_bin]++;
	  
	  int diag1_bin = i_bin + j_bin;
	  int diag2_bin = i_bin + rev_j_bin;
	  
	  diag1_N[diag1_bin]++;
	  diag2_N[diag2_bin]++;
	  
	  all_diag1_N.push_back(diag1_bin);
	  all_diag2_N.push_back(diag2_bin);
	  
	  i_sum[i_bin] += norm_matrix[k][j];
	  j_sum[j_bin] += norm_matrix[k][j];
	  i_log_sum[i_bin] += log(norm_matrix[k][j] + 0.1);
	  j_log_sum[j_bin] += log(norm_matrix[k][j] + 0.1);
	  diag1_sum[diag1_bin] += norm_matrix[k][j];
	  diag2_sum[diag2_bin] += norm_matrix[k][j];
	  diag1_log_sum[diag1_bin] += log(norm_matrix[k][j] + 0.1);
	  diag2_log_sum[diag2_bin] += log(norm_matrix[k][j] + 0.1);
	  diag1_bin_sum[diag1_bin] += 1;
	  diag2_bin_sum[diag2_bin] += 1;
	  all_diag1.push_back(norm_matrix[k][j]);
	  all_diag2.push_back(norm_matrix[k][j]);
	  
	  
	}
      }
      
      double i_cor = get_corr_map(i_sum,i_N);
      double j_cor = get_corr_map(j_sum,j_N);
      double i_log_cor = get_corr_map(i_log_sum,i_N);
      double j_log_cor = get_corr_map(j_log_sum,j_N);
      double diag1_cor = get_corr_map(diag1_sum,diag1_N);
      double diag2_cor = get_corr_map(diag2_sum,diag2_N);
      double diag1_log_cor = get_corr_map(diag1_log_sum,diag1_N);
      double diag2_log_cor = get_corr_map(diag2_log_sum,diag2_N);
      double diag1_bin_cor = get_corr_map(diag1_bin_sum,diag1_N);
      double diag2_bin_cor = get_corr_map(diag2_bin_sum,diag2_N);
      double all_diag1_cor = get_corr_vec(all_diag1,all_diag1_N);
      double all_diag2_cor = get_corr_vec(all_diag2,all_diag2_N);
      
      double final_i = all_diag1_cor + all_diag2_cor + diag1_bin_cor + diag2_bin_cor + diag1_cor + diag2_cor + diag1_log_cor + diag2_log_cor + i_cor + i_log_cor;
      double final_j = all_diag1_cor + -1*all_diag2_cor + diag1_bin_cor + -1*diag2_bin_cor + diag1_cor + -1*diag2_cor + diag1_log_cor + -1*diag2_log_cor + j_cor + j_log_cor;
      
      string strand1;
      string strand2;
      
      if (final_i > 0) {
        strand1 = "+";
      } else {
        strand1 = "-";
      }
      
      if (final_j > 0) {
        strand2 = "+";
      } else {
        strand2 = "-";
      }
      
      break_output << odds_sum << "\t" << new_array1[0] << "\t" << (start1 + s1*bin_size) << "\t" << (start1 + e1*bin_size + bin_size) << "\t" << strand1 << "\t" << new_array1[3] << "\t" << (start2 + s2*bin_size) << "\t" << (start2 + e2*bin_size + bin_size) << "\t" << strand2 << "\n";

      
      
    }
    
  }

  break_output.close();
  
}

void find_parameters_100kb_inter (unordered_map<int,double>& m_hash,                \
				  unordered_map<int,double>& r_hash,	\
				  unordered_map<int,int>& N_hash,	\
				  unordered_map<string,unordered_map<string,int> >& data_map, \
				  unordered_map<string,vector<int> >& ref_map, \
				  unordered_map<string,unordered_map<int,double> >& bias_map, \
				  int bin_size,				\
				  int max_size,				\
				  double &inter_m,			\
				  double &inter_r,			\
				  double &submatrix,	\
				  int &total_weight) {

  std::time_t result = std::time(0);

  cerr << "Finding parameters now... ";

  double inter_N = 0;

  unordered_map<string,vector<int> >::iterator it;

  vector<string> chr_array;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {
    string chr = it->first;
    chr_array.push_back(chr);
  }

  for (int u = 0; u < chr_array.size(); u++) {

    string chr1 = chr_array[u];
    vector<int> loc_array1 = ref_map[chr1];

    for (int v = u; v < chr_array.size(); v++) {

      if (v == u) {

        int max_d = loc_array1[loc_array1.size() - 1] - loc_array1[0];

        for (int i=bin_size; i < max_d; i += bin_size) {
          int tally = (max_d - i)/bin_size;
          N_hash[i] += tally;
          if (i <= max_size) {
            total_weight++;
          }
        }

      } else {

	string chr2 = chr_array[v];
        vector<int> loc_array2 = ref_map[chr2];

        double v1 = loc_array1.size();
        double v2 = loc_array2.size();

        double temp_N = v1*(v1 + 1)*v2*(v2 + 1)/4;

        submatrix += temp_N;

        inter_N += v1*v2;

      }
    }
  }

  unordered_map<int,double> sum_hash;
  unordered_map<int,double> var_sum_hash;

  double inter_sum = 0;
  double inter_var_sum = 0;

  unordered_map<string,unordered_map<string,int> >::iterator it1;

  for (it1 = data_map.begin(); it1 != data_map.end(); it1++) {
    string tag1 = it1->first;
    unordered_map<string,int>::iterator it2;
    vector<string> array1;
    split_string(array1,tag1,'_');

    string chr1 = array1[0];
    vector<int> loc_array1 = ref_map[chr1];
    int loc1 = stoi(array1[1]);

    for (it2 = data_map[tag1].begin(); it2 != data_map[tag1].end(); it2++) {
      string tag2 = it2->first;

      vector<string> array2;
      split_string(array2,tag2,'_');

      string chr2 = array2[0];
      vector<int> loc_array2 = ref_map[chr2];
      int loc2 = stoi(array2[1]);

      double b1 = bias_map[chr1][loc1];
      double b2 = bias_map[chr2][loc2];

      double n;

      if (data_map.count(tag1) == 1) {
        if (data_map[tag1].count(tag2) == 1) {
          n = data_map[tag1][tag2];
        } else {
          n = 0;
        }
      } else {
        n = 0;
      }

      double norm = n/(b1*b2);

      if (chr1 == chr2) {
	if (loc2 > loc1) {
	  int dist = loc2 - loc1;
	  sum_hash[dist] += norm;
	}
      } else {

        inter_sum += norm;

      }

    }
  }

  for (it1 = data_map.begin(); it1 != data_map.end(); it1++) {
    string tag1 = it1->first;
    unordered_map<string,int>::iterator it2;

    vector<string> array1;
    split_string(array1,tag1,'_');

    string chr1 = array1[0];
    vector<int> loc_array1 = ref_map[chr1];
    int loc1 = stoi(array1[1]);

    for (it2 = data_map[tag1].begin(); it2 != data_map[tag1].end(); it2++) {
      string tag2 = it2->first;

      vector<string> array2;
      split_string(array2,tag2,'_');

      string chr2 = array2[0];
      vector<int> loc_array2 = ref_map[chr2];
      int loc2 = stoi(array2[1]);

      double b1 = bias_map[chr1][loc1];
      double b2 = bias_map[chr2][loc2];

      double n;

      if (data_map.count(tag1) == 1) {
        if (data_map[tag1].count(tag2) == 1) {
          n = data_map[tag1][tag2];
        } else {
          n = 0;
        }
      } else {
        n = 0;
      }

      double norm = n/(b1*b2);

      if (chr1 == chr2) {
	if (loc2 > loc1) {

	  int dist = loc2 - loc1;

	  double val = (norm - sum_hash[dist]/N_hash[dist])*(norm - sum_hash[dist]/N_hash[dist]);

	  var_sum_hash[dist] += val;

	}

      } else {

        double val = (norm - inter_sum/inter_N)*(norm - inter_sum/inter_N);

        inter_var_sum += val;

      }

    }

  }

  inter_m = inter_sum/inter_N;
  double inter_var = inter_var_sum/inter_N;
  inter_r = inter_m*inter_m/(inter_var - inter_m);

  if ((inter_r > 5) ||
      (inter_r < 0)) {
    inter_r = 5;
  }

  unordered_map<int,int>::iterator it3;

  for (it3 = N_hash.begin(); it3 != N_hash.end(); it3++) {

    int dist = it3->first;

    double sum_val;
    int N_val;

    if (sum_hash[it3->first] > 0) {
      sum_val = sum_hash[it3->first];
      N_val =  N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      sum_val = sum_hash[temp_dist];
      N_val = N_hash[temp_dist];
    }
    double var_sum_val;
    int var_sum_N;
    if (var_sum_hash[it3->first] > 0) {
      var_sum_val = var_sum_hash[it3->first];
      var_sum_N = N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (var_sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      var_sum_val = var_sum_hash[temp_dist];
      var_sum_N = N_hash[temp_dist];
    }

    double m_val = sum_val/(double)N_val;
    double var = var_sum_val/(double)var_sum_N;
    double r_val;
    if (var > m_val) {
      r_val = m_val*m_val/(var - m_val);
    } else {
      r_val = 20;
    }

    m_hash[dist] = m_val;
    r_hash[dist] = r_val;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

  cerr << inter_m << "\t" << inter_var << "\t" << inter_r << "\t" << inter_sum << "\t" << inter_var_sum << "\t" << inter_N << "\n";

}

void get_100kb_odds_ratio_inter (string super_matrix_file,	\
				 string bias_vector_file,	\
				 string block_file,		\
				 int bin_size,			\
				 int max_size,			\
				 string name) {
  
  vector<string> block_data = read_blocks(block_file,max_size);

  unordered_map<string,vector<int> > ref_map;
  unordered_map<string,unordered_map<int,double> > bias_map;
  read_bias_vector(bias_vector_file,ref_map,bias_map);

  unordered_map<string,unordered_map<string,int> > data_map;
  read_super_matrix(super_matrix_file,data_map);

  unordered_map<int,double> m_hash;
  unordered_map<int,double> r_hash;
  unordered_map<int,int> N_hash;
  double inter_m;
  double inter_r;
  double submatrix = 0;
  int total_weight;

  find_parameters_100kb_inter(m_hash,r_hash,N_hash,data_map,ref_map,bias_map,bin_size,max_size,inter_m,inter_r,submatrix,total_weight);

  double prior = 0.000001;
  double corr = log(prior/(1 - prior));
  double thresh = -1*(log(0.05/submatrix) + corr);

  cerr << "threshold is\t" << thresh << "\n";

  find_breaks_100kb_inter(block_data,m_hash,r_hash,N_hash,inter_m,inter_r,data_map,ref_map,bias_map,bin_size,max_size,total_weight,thresh,name);

}

void find_parameters_100kb_intra (unordered_map<int,double>& m_hash,                \
				  unordered_map<int,double>& r_hash,	\
				  unordered_map<int,int>& N_hash,	\
				  unordered_map<string,int>& data_map,	\
				  unordered_map<string,vector<int> >& ref_map, \
				  unordered_map<string,unordered_map<int,double> >& bias_map, \
				  unordered_map<string,double>& exp_map, \
				  int bin_size,				\
				  int max_size,				\
				  int& max_weight,			\
				  double& submatrix) {
  
  std::time_t result = std::time(0);

  cerr << "Finding parameters now... ";

  unordered_map<int,double> sum_hash;
  unordered_map<int,double> var_sum_hash;

  unordered_map<string,vector<int> >::iterator it;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {

    vector<int> loc_array = ref_map[it->first];
    string chr = it->first;

    double big_N = loc_array.size() - (double)1;

    for (int i=1; i < big_N; i++) {
      submatrix += (big_N + 1 - i)*(big_N + 2 - i)*i/2;
    }

    for (int i=0; i < (loc_array.size() - 1); i++) {
      for (int j=(i + 1); j < loc_array.size(); j++) {
        string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array[i]) )->str();
        string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array[j]) )->str();
        double b1 = bias_map[chr][loc_array[i]];
        double b2 = bias_map[chr][loc_array[j]];
        string tag = chr + "_" + loc1 + "_" + loc2;
        int dist = loc_array[j] - loc_array[i];
        double n;

        if (data_map.count(tag) == 1) {
          n = data_map[tag];
        } else {
          n = 0;
        }

        double expect;
        if (exp_map.count(tag) == 1) {
          expect = exp_map[tag];
        } else {
          expect = 1;
        }

        double norm = n/(b1*b2);

        sum_hash[dist] += norm;

        N_hash[dist]++;

        if (dist <= max_size) {
          max_weight++;
        }
      }
    }

  }

  unordered_map<string,vector<int> >::iterator it2;

  for (it2 = ref_map.begin(); it2 != ref_map.end(); it2++) {

    vector<int>loc_array = ref_map[it2->first];
    string chr = it2->first;

    for (int i=0; i < (loc_array.size() - 1); i++) {
      for (int j=(i + 1); j < loc_array.size();j++) {
        string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array[i]) )->str();
        string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array[j]) )->str();
        double b1 = bias_map[chr][loc_array[i]];
        double b2 = bias_map[chr][loc_array[j]];
        string tag = chr + "_" + loc1 + "_" + loc2;
        string tag1 = chr + "_" + loc1;
        string tag2 = chr + "_" + loc2;
        int dist = loc_array[j] - loc_array[i];
        double n;
        if (data_map.count(tag) == 1) {
          n = data_map[tag];
        } else {
          n = 0;
        }

        double norm = n/(b1*b2);

        double m_val = sum_hash[dist]/(double)N_hash[dist];

        double var_val = (norm - m_val)*(norm - m_val);

        var_sum_hash[dist] += var_val;

      }
    }

  }

  unordered_map<int,double> r_mod_hash;

  unordered_map<int,double>::iterator it3;

  for (it3 = sum_hash.begin(); it3 != sum_hash.end(); it3++) {
    int dist = it3->first;
    double sum_val;
    int N_val;
    if (sum_hash[it3->first] > 0) {
      sum_val = sum_hash[it3->first];
      N_val =  N_hash[it3->first];

    } else {
      int temp_dist = dist;
      while (sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      sum_val = sum_hash[temp_dist];
      N_val = N_hash[temp_dist];
    }
    double var_sum_val;
    int var_sum_N;
    if (var_sum_hash[it3->first] > 0) {
      var_sum_val = var_sum_hash[it3->first];
      var_sum_N = N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (var_sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      var_sum_val = var_sum_hash[temp_dist];
      var_sum_N = N_hash[temp_dist];
    }
    double m_val = sum_val/(double)N_val;
    double var = var_sum_val/(double)var_sum_N;
    double var_mod = var*0.95;

    double r_val;
    if (var > m_val) {
      r_val = m_val*m_val/(var - m_val);
    } else {
      r_val = 20;
    }

    m_hash[dist] = m_val;
    r_hash[dist] = r_val;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";
  
}

void find_breaks_100kb_intra (unordered_map<int,double>& m_hash,	\
			      unordered_map<int,double>& r_hash,	\
			      unordered_map<int,int>& N_hash,		\
			      unordered_map<string,int>& data_map,	\
			      unordered_map<string,vector<int> >& ref_map, \
			      unordered_map<string,unordered_map<int,double> >& bias_map, \
			      unordered_map<string,double>& exp_map,	\
			      int bin_size,				\
			      int max_size,				\
			      int max_weight,				\
			      double thresh,				\
			      string name) {
  
  unordered_map<int,double> fact_hash;
  fact_hash[0] = 0;
  fact_hash[1] = 0;

  std::time_t result = std::time(0);

  cerr << "Finding breaks now...\n";

  ofstream break_output;
  string break_output_file_name = name + ".breaks.txt";

  break_output.open(break_output_file_name);

  unordered_map<string,vector<int> >::iterator it;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {

    vector<int> loc_array = ref_map[it->first];
    string chr = it->first;

    cerr << chr << "\n";

    if (chr != "chr19") {
      //     continue;
    }

    vector<vector<double> > odds_matrix(loc_array.size(),vector<double>(loc_array.size()));
    vector<vector<double> > ori_odds_matrix(loc_array.size(),vector<double>(loc_array.size()));
    vector<vector<double> > norm_matrix(loc_array.size(),vector<double>(loc_array.size()));

    for (int i=0; i < (loc_array.size() - 1); i++) {
      for (int j=(i + 1); j < loc_array.size(); j++) {
	string loc1 = static_cast<ostringstream*>( &(ostringstream() << loc_array[i]) )->str();
	string loc2 = static_cast<ostringstream*>( &(ostringstream() << loc_array[j]) )->str();
	double b1 = bias_map[chr][loc_array[i]];
	double b2 = bias_map[chr][loc_array[j]];
        string tag = chr + "_" + loc1 + "_" + loc2;
        string tag1 = chr + "_" + loc1;
        string tag2 = chr + "_" + loc2;
        int n;
        if (data_map.count(tag) == 1) {
          n = data_map[tag];
        } else {
          n = 0;
        }

        if (fact_hash.count(n) == 0) {
          int temp_n = n;
          while (fact_hash.count(temp_n) == 0) {
            temp_n--;
          }
          for (int k = (temp_n + 1); k <= n; k++) {
            double t_val = fact_hash[k - 1] + log(k);
            fact_hash[k] = t_val;
          }
        }

        double norm = n/(b1*b2);

        double expect;
        if (exp_map.count(tag) > 0) {
          if (exp_map[tag] > 1) {
            expect = exp_map[tag];
          } else {
            expect = 1;
          }
        } else {
          expect = 1;
        }

        int dist = loc_array[j] - loc_array[i];

        double test_m = m_hash[dist]*b1*b2*expect;

        if (test_m < 0.5*m_hash[dist]) {
          test_m = 0.5*m_hash[dist];
        }

	double term1 = r_hash[dist]*log(r_hash[dist]/(r_hash[dist] + test_m));

        double term2 = 0;

        for (int k=1; k <= n; k++) {
          double temp_val = log(n + r_hash[dist] - k);
          term2 += temp_val;
        }

        term2 += -1*fact_hash[n];

        double term3 = n*log(test_m/(test_m + r_hash[dist]));

        double act_p = term1 + term2 + term3;

        double low_p = 0;

        double weight_sum = 0;
        double p_sum = 0;

        for (int k = bin_size; k <= max_size; k += bin_size) {

          double t_m = m_hash[k]*b1*b2;

          if (t_m > test_m) {

            double t1 = r_hash[k]*log(r_hash[k]/(r_hash[k] + t_m));

            double t2 = 0;

            for (int l=1; l <= n; l++) {
              double temp_val = log(n + r_hash[k] - l);
              t2 += temp_val;
            }

            t2 += -1*fact_hash[n];

            double t3 = n*log(t_m/(t_m + r_hash[k]));

            double temp_p = N_hash[k]*exp(t1 + t2 + t3)/max_weight;

            low_p += temp_p;

            double new_weight = (max_size - k)/bin_size + 1;
            p_sum += new_weight*exp(t1 + t2 + t3);
            weight_sum += new_weight;

          }

        }

	double final_p = p_sum/weight_sum;

        double odds = log(final_p) - act_p;

        if (dist <= max_size) {

          odds_matrix[i][j] = -1;
          ori_odds_matrix[i][j] = -1;
          norm_matrix[i][j] = norm;
	  
        } else {
	  
	  if ((exp_map.count(tag) > 0) &&
	      (exp_map[tag] > 0)) {
		
	    odds_matrix[i][j] = odds;
	    ori_odds_matrix[i][j] = odds;
	    norm_matrix[i][j] = norm;
	      
	  } else {
		
	    odds_matrix[i][j] = 0;
	    ori_odds_matrix[i][j] = 0;
	    norm_matrix[i][j] = norm;
	    
	  }
	  
	}

      }
    }

    int check = 1;
    vector<string> top_hits;
    vector<double> top_sums;

    while (check == 1) {

      double max_sum = 0;

      int k1;
      int k2;
      int j1;
      int j2;

      for (int i=0; i < loc_array.size(); i++) {

        vector<double> temp(i,0);

        for (int j=i; j < loc_array.size(); j++) {

          for (int k=0; k < i; k++) {

            temp[k] += odds_matrix[k][j];

          }

          double temp_sum = 0;
          double temp_max = 0;

          int t_k1 = 0;
          int t_k2 = 0;

          int f_k1;
          int f_k2;

          for (int k=0; k< i; k++) {

            double new_sum = temp_sum + temp[k];

            if (new_sum < 0) {
              temp_sum = 0;
              t_k1 = k + 1;
            } else {
              temp_sum = new_sum;
              t_k2 = k;
              if (temp_sum > temp_max) {
                temp_max = temp_sum;
                f_k1 = t_k1;
                f_k2 = t_k2;
              }
            }
          }

          if (temp_max > max_sum) {
            max_sum = temp_max;
            j1 = i;
            j2 = j;
            k1 = f_k1;
            k2 = f_k2;
          }

        }
      }

      if (max_sum > thresh) {

        for (int k = k1; k <= k2; k++) {
          for (int j = j1; j <= j2; j++) {

            double change = -5;
            odds_matrix[k][j] = change;

          }
        }

        string s1 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
        string e1 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();
        string s2 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
        string e2 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();

        string new_string = chr + "_" + s1 + "_" + e1 + "_" + chr + "_" + s2 + "_" + e2;

        top_hits.push_back(new_string);
        top_sums.push_back(max_sum);

      } else {
        check = 0;
        break;
      }
    }

    int overlap_check = 0;
    if (top_hits.size() > 1) {

      for (int i=0; i < (top_hits.size() - 1); i++) {
        vector<std::string> new_array1;
	split_string(new_array1,top_hits[i],'_');
        int s1 = stoi(new_array1[1]);
        int e1 = stoi(new_array1[2]);
        int s2 = stoi(new_array1[4]);
        int e2 = stoi(new_array1[5]);
        for (int j=(i + 1); j < top_hits.size(); j++) {
          vector<std::string> new_array2;
	  split_string(new_array2,top_hits[j],'_');
          int s3 = stoi(new_array2[1]);
          int e3 = stoi(new_array2[2]);
          int s4 = stoi(new_array2[4]);
          int e4 = stoi(new_array2[5]);

          int ot1 = 0;
          if (s1 < s3) {
            if (e1 >= s3) {
              ot1 = 1;
            }
          } else {
            if (s1 <= e3) {
              ot1 = 1;
            }
          }

          int ot2 = 0;
          if (s2 < s4) {
            if (e2 >= s4) {
              ot2 = 1;
            }
          } else {
            if (s2 <= e4) {
              ot2 = 1;
            }
          }

          if ((ot1 == 1) &&
              (ot2 == 1)) {
            overlap_check = 1;
          }

        }

      }

    }

    if (overlap_check == 1) {

      int set = 0;

      while ((set + 1) < top_hits.size()) {
	
        vector<vector<double> > new_odds_matrix = ori_odds_matrix;
        vector<vector<double> > test_odds_matrix = ori_odds_matrix;

        vector<string> new_top_hits;
        vector<double> new_top_sums;
	
        if (set > 0) {
          for (int i = 0; i < set; i++) {
            new_top_hits.push_back(top_hits[i]);
            new_top_sums.push_back(top_sums[i]);

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);

            for (int k = s2; k <= e2; k++) {
              for (int j = s1; j <= e1; j++) {

                double change = -5;
                new_odds_matrix[j][k] = change;
                test_odds_matrix[j][k] = change;

              }
            }

          }
        }

        for (int i=0; i < top_hits.size(); i++) {

          if (i == (set + 1)) {

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = stoi(new_array1[1]);
            int e1 = stoi(new_array1[2]);
            int s2 = stoi(new_array1[4]);
            int e2 = stoi(new_array1[5]);

            for (int k = s2; k <= e2; k++) {
              for (int j = s1; j <= e1; j++) {

                double change = -5;
                new_odds_matrix[j][k] = change;

              }
            }
          }
        }

	int new_check = 1;
        int pass = 0;
	
        while (new_check == 1) {

          double max_sum = 0;
	  
          int k1;
          int k2;
          int j1;
          int j2;

	    
          for (int i=0; i < loc_array.size(); i++) {

            vector<double> temp(i,0);
            vector<double> temp_raw(i,0);

	    

            for (int j=i; j < loc_array.size(); j++) {

	      

              for (int k=0; k < i; k++) {

                if (pass == 0) {

                  temp[k] += new_odds_matrix[k][j];

                } else {

                  temp[k] += test_odds_matrix[k][j];

                }

              }

	      

              double temp_sum = 0;
              double temp_max = 0;

              int t_k1 = 0;
              int t_k2 = 0;

              int f_k1;
              int f_k2;

	      for (int k=0; k< i; k++) {

                double new_sum = temp_sum + temp[k];

		if (new_sum < 0) {
			    
                  temp_sum = 0;
                  t_k1 = k + 1;
		
                } else {

                  temp_sum = new_sum;
                  t_k2 = k;

                  if (temp_sum > temp_max) {

		  
		    temp_max = temp_sum;
                    f_k1 = t_k1;
                    f_k2 = t_k2;
                  }

                }

              }

              if (temp_max > max_sum) {
		max_sum = temp_max;
                j1 = i;
                j2 = j;
                k1 = f_k1;
                k2 = f_k2;

              }
	      

            }

	    

	  }

          if (max_sum > thresh) {

            for (int k = k1; k <= k2; k++) {

              for (int j = j1; j <= j2; j++) {

		double change = -5;
		test_odds_matrix[k][j] = change;
              }
            }
	      
	    string s1 = static_cast<ostringstream*>( &(ostringstream() << k1) )->str();
            string e1 = static_cast<ostringstream*>( &(ostringstream() << k2) )->str();
            string s2 = static_cast<ostringstream*>( &(ostringstream() << j1) )->str();
            string e2 = static_cast<ostringstream*>( &(ostringstream() << j2) )->str();
	    
            string new_string = chr + "_" + s1 + "_" + e1 + "_" + chr + "_" + s2 + "_" + e2;

            string j_strand = "s";
            string k_strand = "s";
	    
            new_top_hits.push_back(new_string);
            new_top_sums.push_back(max_sum);
	    	    	    
	    pass = 1;
	  	    
          } else {
            break;
          }
	  
	  
	}

        double all_sum = 0;
        double all_new_sum = 0;
	
        for (int i = 0; i < top_sums.size(); i++) {
          all_sum += top_sums[i];
        }

        for (int i=0; i < new_top_sums.size(); i++) {
          all_new_sum += new_top_sums[i];
        }

        if (all_new_sum > all_sum) {
          top_hits = new_top_hits;
          top_sums = new_top_sums;

        } else {
          set++;
        }
	
	int final_overlap_check = 0;
		
        for (int i=set; i < (top_hits.size() - 1); i++) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[4]);
          int e2 = stoi(new_array1[5]);
          for (int j=(i + 1); j < top_hits.size(); j++) {
            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');
            int s3 = stoi(new_array2[1]);
            int e3 = stoi(new_array2[2]);
            int s4 = stoi(new_array2[4]);
            int e4 = stoi(new_array2[5]);

            int ot1 = 0;
            if (s1 < s3) {
              if (e1 >= s3) {
                ot1 = 1;
              }
            } else {
              if (s1 <= e3) {
                ot1 = 1;
              }
            }

            int ot2 = 0;
            if (s2 < s4) {
              if (e2 >= s4) {
                ot2 = 1;
              }
            } else {
              if (s2 <= e4) {
                ot2 = 1;
              }
            }

            if ((ot1 == 1) &&
                (ot2 == 1)) {
              final_overlap_check = 1;
            }
          }
        }

        if (final_overlap_check == 0) {
	  break;
        }

      }

    }

    int stop_test = 0;
    while (stop_test == 0) {
      if (top_hits.size() > 1) {
        int over_check = 0;
        double over_max = 0;
        double over_min = -10000000;
        int over_i;
        int over_j;

        double all_sum = 0;
        vector<vector<double> > all_odds_matrix = ori_odds_matrix;

        for (int i=0; i < top_hits.size(); i++) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[4]);
          int e2 = stoi(new_array1[5]);
          for (int k = s1; k <= e1; k++) {
            for (int l = s2; l <= e2; l++) {
              all_sum += all_odds_matrix[k][l];
              all_odds_matrix[k][l] = -5;
            }
          }
        }

	for (int i=0; i < (top_hits.size() - 1); i++) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[4]);
          int e2 = stoi(new_array1[5]);
          for (int j=(i + 1); j < top_hits.size(); j++) {
            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');
            int s3 = stoi(new_array2[1]);
            int e3 = stoi(new_array2[2]);
            int s4 = stoi(new_array2[4]);
            int e4 = stoi(new_array2[5]);

            int new_start1;
            if (s1 < s3) {
              new_start1 = s1;
            } else {
              new_start1 = s3;
            }
            int new_end1;
            if (e1 > e3) {
              new_end1 = e1;
            } else {
              new_end1 = e3;
            }
            int new_start2;
            if (s2 < s4) {
              new_start2 = s2;
            } else {
              new_start2 = s4;
            }
            int new_end2;
            if (e2 > e4) {
              new_end2 = e2;
            } else {
              new_end2 = e4;
            }

            vector<vector<double> > all_new_odds_matrix = ori_odds_matrix;

            double all_new_sum = 0;
            double odds_sum = 0;
            for (int k = new_start1; k <= new_end1; k++) {
              for (int l = new_start2; l <= new_end2; l++) {
                all_new_sum += all_new_odds_matrix[k][l];
                all_new_odds_matrix[k][l] = 0;
                odds_sum += ori_odds_matrix[k][l];
              }
            }

            for (int k=0; k < top_hits.size(); k++) {
              if ((k == i) ||
                  (k == j)) {
                continue;
              }
              vector<std::string> new_array3;
	      split_string(new_array3,top_hits[k],'_');
              int s5 = stoi(new_array3[1]);
              int e5 = stoi(new_array3[2]);
              int s6 = stoi(new_array3[4]);
              int e6 = stoi(new_array3[5]);
              for (int l = s5; l <= e5; l++) {
                for (int m = s6; m <= e6; m++) {
                  all_new_sum += all_new_odds_matrix[l][m];
                  all_new_odds_matrix[l][m] = 0;
                }
              }
            }

	    double loss = all_new_sum - all_sum;

            if ((all_new_sum >= over_max) &&
                (odds_sum >= thresh)) {
              over_min = loss;
              over_max = all_new_sum;
              over_i = i;
              over_j = j;
              over_check++;
            }
          }
        }

        if (over_check > 0) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[over_i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[4]);
          int e2 = stoi(new_array1[5]);
          vector<std::string> new_array2;
	  split_string(new_array2,top_hits[over_j],'_');
          int s3 = stoi(new_array2[1]);
          int e3 = stoi(new_array2[2]);
          int s4 = stoi(new_array2[4]);
          int e4 = stoi(new_array2[5]);
          int ns1;
          string new_start1;
          if (s1 <= s3) {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
            ns1 = s1;
          } else {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << s3) )->str();
            ns1 = s3;
          }
          int ne1;
          string new_end1;
          if (e1 >= e3) {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
            ne1 = e1;
          } else {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << e3) )->str();
            ne1 = e3;
          }
	  int ns2;
          string new_start2;
          if (s2 <= s4) {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
            ns2 = s2;
          } else {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << s4) )->str();
            ns2 = s4;
          }
          int ne2;
          string new_end2;
          if (e2 >= e4) {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();
            ne2 = e2;
          } else {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << e4) )->str();
            ne2 = e4;
          }
          double new_sum = 0;
          for (int k = ns1; k <= ne1; k++) {
            for (int l = ns2; l <= ne2; l++) {
              odds_matrix[k][l] = -5;
              new_sum += ori_odds_matrix[k][l];
            }
          }
          string new_string = new_array1[0] + "_" + new_start1 + "_" + new_end1 + "_" + new_array1[3] + "_" + new_start2 + "_" + new_end2;

          string j_strand = "s";
          string k_strand = "s";

          vector<string> new_top_hits;
          vector<double> new_top_sums;
          new_top_hits.push_back(new_string);
          new_top_sums.push_back(new_sum);
          for (int i=0; i < top_hits.size(); i++) {
            if ((i != over_i) &&
                (i != over_j)) {
              new_top_hits.push_back(top_hits[i]);
              new_top_sums.push_back(top_sums[i]);
            }
          }
	  top_hits = new_top_hits;
          top_sums = new_top_sums;
        } else {
          stop_test = 1;
        }
      } else {
        stop_test = 1;
      }
    }

    for (int i=0; i < top_hits.size(); i++) {

      vector<std::string> new_array1;
      split_string(new_array1,top_hits[i],'_');
      int s1 = stoi(new_array1[1]);
      int e1 = stoi(new_array1[2]);
      int s2 = stoi(new_array1[4]);
      int e2 = stoi(new_array1[5]);

      unordered_map<int,double> i_sum;
      unordered_map<int,double> j_sum;
      unordered_map<int,int> i_N;
      unordered_map<int,int> j_N;
      vector<double> all_diag1;
      vector<double> all_diag2;
      vector<int> all_diag1_N;
      vector<int> all_diag2_N;
      unordered_map<int,double> diag1_sum;
      unordered_map<int,double> diag2_sum;
      unordered_map<int,int> diag1_N;
      unordered_map<int,int> diag2_N;
      unordered_map<int,double> i_log_sum;
      unordered_map<int,double> j_log_sum;
      unordered_map<int,double> diag1_log_sum;
      unordered_map<int,double> diag2_log_sum;
      unordered_map<int,double> diag1_bin_sum;
      unordered_map<int,double> diag2_bin_sum;
      double odds_sum = 0;
      for (int k=s1; k <= e1; k++) {

        int i_bin = k - s1;
        i_N[i_bin]++;

        for (int j=s2; j <= e2; j++) {
          odds_sum += ori_odds_matrix[k][j];

          int j_bin = j - s2;
          int rev_j_bin = e2 - j;
          j_N[j_bin]++;

          int diag1_bin = i_bin + j_bin;
          int diag2_bin = i_bin + rev_j_bin;

          diag1_N[diag1_bin]++;
          diag2_N[diag2_bin]++;

          all_diag1_N.push_back(diag1_bin);
          all_diag2_N.push_back(diag2_bin);

          i_sum[i_bin] += norm_matrix[k][j];
          j_sum[j_bin] += norm_matrix[k][j];
          i_log_sum[i_bin] += log(norm_matrix[k][j] + 0.1);
          j_log_sum[j_bin] += log(norm_matrix[k][j] + 0.1);
          diag1_sum[diag1_bin] += norm_matrix[k][j];
          diag2_sum[diag2_bin] += norm_matrix[k][j];
          diag1_log_sum[diag1_bin] += log(norm_matrix[k][j] + 0.1);
          diag2_log_sum[diag2_bin] += log(norm_matrix[k][j] + 0.1);
          diag1_bin_sum[diag1_bin] += 1;
          diag2_bin_sum[diag2_bin] += 1;
          all_diag1.push_back(norm_matrix[k][j]);
          all_diag2.push_back(norm_matrix[k][j]);


        }
      }

      double i_cor = get_corr_map(i_sum,i_N);
      double j_cor = get_corr_map(j_sum,j_N);
      double i_log_cor = get_corr_map(i_log_sum,i_N);
      double j_log_cor = get_corr_map(j_log_sum,j_N);
      double diag1_cor = get_corr_map(diag1_sum,diag1_N);
      double diag2_cor = get_corr_map(diag2_sum,diag2_N);
      double diag1_log_cor = get_corr_map(diag1_log_sum,diag1_N);
      double diag2_log_cor = get_corr_map(diag2_log_sum,diag2_N);
      double diag1_bin_cor = get_corr_map(diag1_bin_sum,diag1_N);
      double diag2_bin_cor = get_corr_map(diag2_bin_sum,diag2_N);
      double all_diag1_cor = get_corr_vec(all_diag1,all_diag1_N);
      double all_diag2_cor = get_corr_vec(all_diag2,all_diag2_N);

      double final_i = all_diag1_cor + all_diag2_cor + diag1_bin_cor + diag2_bin_cor + diag1_cor + diag2_cor + diag1_log_cor + diag2_log_cor + i_cor + i_log_cor;
      double final_j = all_diag1_cor + -1*all_diag2_cor + diag1_bin_cor + -1*diag2_bin_cor + diag1_cor + -1*diag2_cor + diag1_log_cor + -1*diag2_log_cor + j_cor + j_log_cor;

      string strand1;
      string strand2;

      if (final_i > 0) {
        strand1 = "+";
      } else {
        strand1 = "-";
      }

      if (final_j > 0) {
        strand2 = "+";
      } else {
        strand2 = "-";
      }

      break_output << odds_sum << "\t" << chr << "\t" << loc_array[s1] << "\t" << (loc_array[e1] + bin_size) << "\t" << strand1 << "\t" << chr << "\t" << loc_array[s2] << "\t" << (loc_array[e2] + bin_size) << "\t" << strand2 << "\n";


    }
  }

  break_output.close();

}
    
void read_super_matrix_intra (string super_matrix_file,	\
				    unordered_map<string,int>& ref_map) {
  
  std::time_t result = std::time(0);
  
  cerr << "Reading in super matrix file now... ";
  
  string line;
  ifstream myfile (super_matrix_file);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    vector<string> loc1;
    split_string(loc1,array[0],':');
    
    vector<string> loc2;
    split_string(loc2,array[1],':');

    if (loc1[0] == loc2[0]) {
      vector<string> range1;
      split_string(range1,loc1[1],'-');

      vector<string> range2;
      split_string(range2,loc2[1],'-');

      string tag = loc1[0] + "_" + range1[0] + "_" + range2[0];

      int n = string_to_int(array[2]);

      ref_map[tag] = n;

    }

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void read_index_file_100kb_intra (string index_file_name,		\
				  unordered_map<string,double>& index_map) {

  std::time_t result = std::time(0);

  cerr << "Reading in exp file now... ";

  string line;
  ifstream myfile (index_file_name);
  while ( getline(myfile,line) ) {
    vector<string> array;
    split_string(array,line,'\t');

    vector<string> loc1;
    split_string(loc1,array[0],':');

    vector<string> loc2;
    split_string(loc2,array[1],':');

    if (loc1[0] == loc2[0]) {
      vector<string> range1;
      split_string(range1,loc1[1],'-');

      vector<string> range2;
      split_string(range2,loc2[1],'-');

      string tag = loc1[0] + "_" + range1[0] + "_" + range2[0];

      double n = stod(array[2]);

      index_map[tag] = n;

    }

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void get_100kb_odds_ratio_intra (string super_matrix_file,	\
				 string bias_vector_file,	\
				 char* exp_file_intra,		\
				 int bin_size,			\
				 int max_size,			\
				 string name) {

  unordered_map<string,vector<int> > ref_map;
  unordered_map<string,unordered_map<int,double> > bias_map;
  read_bias_vector(bias_vector_file,ref_map,bias_map);

  unordered_map<string,int> data_map;
  read_super_matrix_intra(super_matrix_file,data_map);

  unordered_map<string,double> exp_map;
  read_index_file_100kb_intra(exp_file_intra,exp_map);

  unordered_map<int,double> m_hash;
  unordered_map<int,double> r_hash;
  unordered_map<int,int> N_hash;
  int max_weight = 0;
  double submatrix = 0;
  find_parameters_100kb_intra(m_hash,r_hash,N_hash,data_map,ref_map,bias_map,exp_map,bin_size,max_size,max_weight,submatrix);

  double prior = 0.000001;
  double corr = log(prior/(1 - prior));
  double thresh = -1*(log(0.05/submatrix) + corr);

  find_breaks_100kb_intra(m_hash,r_hash,N_hash,data_map,ref_map,bias_map,exp_map,bin_size,max_size,max_weight,thresh,name);
  
}

void find_parameters_10kb_intra (unordered_map<int,double>& m_hash,                \
				 unordered_map<int,double>& r_hash,	\
				 unordered_map<int,int>& N_hash,	\
				 unordered_map<string,int>& data_map,	\
				 unordered_map<string,vector<int> >& ref_map, \
				 unordered_map<string,unordered_map<int,double> >& bias_map, \
				 int bin_size,				\
				 int max_size,				\
				 int& max_weight,			\
				 double& submatrix) {
  
  std::time_t result = std::time(0);
  cerr << "Finding parameters now... ";
  
  unordered_map<string,vector<int> >::iterator it;
  
  vector<string> chr_array;

  for (it = ref_map.begin(); it != ref_map.end(); it++) {
    string chr = it->first;
    chr_array.push_back(chr);
  }

  for (int u = 0; u < chr_array.size(); u++) {

    string chr = chr_array[u];
    vector<int> loc_array = ref_map[chr];

    double big_N = loc_array.size() - (double)1;

    for (int i=1; i < big_N; i++) {
      submatrix += (big_N + 1 - i)*(big_N + 2 - i)*i/2;
    }

    int max_d = loc_array[loc_array.size() - 1] - loc_array[0];

    for (int i=bin_size; i < max_d; i += bin_size) {
      int tally = (max_d - i)/bin_size;
      N_hash[i] += tally;
      if (i <= max_size) {
        max_weight += tally;
      }
    }
  }

  unordered_map<int,double> sum_hash;
  unordered_map<int,double> var_sum_hash;

  unordered_map<string,int>::iterator it1;

  for (it1 = data_map.begin(); it1 != data_map.end(); it1++) {

    string tag = it1->first;

    vector<string> array;
    split_string(array,tag,'_');

    string chr = array[0];
    int loc1 = stoi(array[1]);
    int loc2 = stoi(array[2]);
    int dist = loc2 - loc1;

    double b1 = bias_map[chr][loc1];
    double b2 = bias_map[chr][loc2];

    double n;
    if (data_map.count(tag) == 1) {
      n = data_map[tag];
    } else {
      n = 0;
    }

    double norm = n/(b1*b2);

    sum_hash[dist] += norm;

  }

  for (it1 = data_map.begin(); it1 != data_map.end(); it1++) {

    string tag = it1->first;

    vector<string> array;
    split_string(array,tag,'_');
    
    string chr = array[0];
    int loc1 = stoi(array[1]);
    int loc2 = stoi(array[2]);
    int dist = loc2 - loc1;

    double b1 = bias_map[chr][loc1];
    double b2 = bias_map[chr][loc2];

    double n;
    if (data_map.count(tag) == 1) {
      n = data_map[tag];
    } else {
      n = 0;
    }

    double norm = n/(b1*b2);

    double val = (norm - sum_hash[dist]/N_hash[dist])*(norm - sum_hash[dist]/N_hash[dist]);

    var_sum_hash[dist] += val;

  }

  unordered_map<int,double>::iterator it3;

  for (it3 = sum_hash.begin(); it3 != sum_hash.end(); it3++) {
    int dist = it3->first;
    double sum_val;
    int N_val;
    if (sum_hash[it3->first] > 0) {
      sum_val = sum_hash[it3->first];
      N_val =  N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      sum_val = sum_hash[temp_dist];
      N_val = N_hash[temp_dist];
    }
    double var_sum_val;
    int var_sum_N;
    if (var_sum_hash[it3->first] > 0) {
      var_sum_val = var_sum_hash[it3->first];
      var_sum_N = N_hash[it3->first];
    } else {
      int temp_dist = dist;
      while (var_sum_hash[temp_dist] == 0) {
        temp_dist += -1*bin_size;
      }
      var_sum_val = var_sum_hash[temp_dist];
      var_sum_N = N_hash[temp_dist];
    }
    double m_val = sum_val/(double)N_val;
    double var = var_sum_val/(double)var_sum_N;

    double r_val;
    if (var > m_val) {
      r_val = m_val*m_val/(var - m_val);
    } else {
      r_val = 20;
    }

    m_hash[dist] = m_val;
    r_hash[dist] = r_val;

  }

  std::time_t result2 = std::time(0);
  int diff = result2 - result;
  cerr << "Done... Required " << diff << " seconds\n";

}

void find_breaks_10kb_intra (vector<string>& block_data,                           \
			     unordered_map<int,double>& m_hash,		\
			     unordered_map<int,double>& r_hash,		\
			     unordered_map<int,int>& N_hash,		\
			     unordered_map<string,int>& data_map,	\
			     unordered_map<string,vector<int> >& ref_map, \
			     unordered_map<string,unordered_map<int,double> >& bias_map, \
			     int bin_size,				\
			     int max_size,				\
			     int max_weight,				\
			     double thresh,\
			     string name) {
  
  unordered_map<int,double> fact_hash;
  fact_hash[0] = 0;
  fact_hash[1] = 0;

  std::time_t result = std::time(0);

  cerr << "Finding breaks now...\n";

  ofstream break_output;
  string break_output_file_name = name + ".breaks.txt";

  break_output.open(break_output_file_name);

  for (int pod = 0; pod < block_data.size(); pod++) {

    string line = block_data[pod];

    cerr << "Going through\t" << line << "\n";

    vector<string> array;
    split_string(array,line,'\t');
    string chr = array[1];
    vector<int> loc_array = ref_map[chr];

    int start1 = stoi(array[2]);
    int end1 = stoi(array[3]);
    int start2 = stoi(array[6]);
    int end2 = stoi(array[7]);

    int i_size = (end1 - start1)/bin_size;
    int j_size = (end2 - start2)/bin_size;

    vector<vector<double> > odds_matrix(i_size,vector<double>(j_size));
    vector<vector<double> > ori_odds_matrix(i_size,vector<double>(j_size));
    vector<vector<double> > norm_matrix(i_size,vector<double>(j_size));

    for (int i=0; i < odds_matrix.size(); i++) {
      for (int j=0; j < odds_matrix[0].size(); j++) {

        int l1 = i*bin_size + start1;
        int l2 = j*bin_size + start2;

        string loc1 = static_cast<ostringstream*>( &(ostringstream() << l1) )->str();
        string loc2 = static_cast<ostringstream*>( &(ostringstream() << l2) )->str();
        double b1;
        if (bias_map[chr].count(l1) > 0) {
          b1 = bias_map[chr][l1];
        } else {
          b1 = 1;
        }
        double b2;
        if (bias_map[chr].count(l2) > 0) {
          b2 = bias_map[chr][l2];
        } else {
          b2 = 1;
        }
        string tag = chr + "_" + loc1 + "_" + loc2;
        int n;
        if (data_map.count(tag) == 1) {
          n = data_map[tag];
        } else {
          n = 0;
        }

        double norm = n/(b1*b2);

        if (fact_hash.count(n) == 0) {
          int temp_n = n;
          while (fact_hash.count(temp_n) == 0) {
            temp_n--;
          }
          for (int k = (temp_n + 1); k <= n; k++) {
            double t_val = fact_hash[k - 1] + log(k);
            fact_hash[k] = t_val;
          }
        }

	int dist = l2 - l1;

        double test_m;

        if (b1*b2 > 0.5) {

          test_m = m_hash[dist]*b1*b2;

        } else {

          test_m = m_hash[dist]*0.5;

        }

        double term1 = r_hash[dist]*log(r_hash[dist]/(r_hash[dist] + test_m));

        double term2 = 0;

        for (int k=1; k <= n; k++) {
          double temp_val = log(n + r_hash[dist] - k);
          term2 += temp_val;
        }

        term2 += -1*fact_hash[n];

        double term3 = n*log(test_m/(test_m + r_hash[dist]));

        double act_p = term1 + term2 + term3;

        double low_p = 0;

        for (int k = bin_size; k <= max_size; k += bin_size) {

          double t_m = m_hash[k]*b1*b2;

          double t1 = r_hash[k]*log(r_hash[k]/(r_hash[k] + t_m));

          double t2 = 0;

          for (int l=1; l <= n; l++) {
            double temp_val = log(n + r_hash[k] - l);
            t2 += temp_val;
          }

          t2 += -1*fact_hash[n];

          double t3 = n*log(t_m/(t_m + r_hash[k]));

          double temp_p = N_hash[k]*exp(t1 + t2 + t3)/max_weight;

	  low_p += temp_p;

        }

        double odds = log(low_p) - act_p;

        if (dist <= max_size) {

          odds_matrix[i][j] = -1;
          ori_odds_matrix[i][j] = -1;
          norm_matrix[i][j] = norm;

        } else {

          odds_matrix[i][j] = odds;
          ori_odds_matrix[i][j] = odds;
          norm_matrix[i][j] = norm;

        }

      }
    }

    vector<string> top_hits;
    vector<double> top_sums;
    int check = 1;
    double t_thresh = thresh;

    while (check == 1) {

      double max_sum = 0;

      int k1;
      int k2;
      int j1;
      int j2;

      for (int i=0; i < odds_matrix[0].size(); i++) {

        vector<double> temp(odds_matrix.size(),0);

        for (int j=i; j < odds_matrix[0].size(); j++) {

          for (int k=0; k < odds_matrix.size(); k++) {

            temp[k] += odds_matrix[k][j];

          }

          double temp_sum = 0;
          double temp_max = 0;

          int t_k1 = 0;
          int t_k2 = 0;

          int f_k1;
          int f_k2;

          for (int k=0; k < odds_matrix.size(); k++) {

            double new_sum = temp_sum + temp[k];

            if (new_sum < 0) {
              temp_sum = 0;
              t_k1 = k + 1;
            } else {
              temp_sum = new_sum;
              t_k2 = k;
              if (temp_sum > temp_max) {
                temp_max = temp_sum;
                f_k1 = t_k1;
                f_k2 = t_k2;
              }
            }
          }

          if (temp_max > max_sum) {
            max_sum = temp_max;
            j1 = i;
            j2 = j;
            k1 = f_k1;
            k2 = f_k2;
          }

        }
      }

      if (max_sum > t_thresh) {

        t_thresh += thresh;

        for (int k = k1; k <= k2; k++) {
          for (int j = j1; j <= j2; j++) {
            double change = -5;
            odds_matrix[k][j] = change;
          }
        }

        string strand3;
        string strand4;

        int s1 = start1 + k1*bin_size;
        string st1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
        int e1 = start1 + k2*bin_size;
        string en1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
        int s2 = start2 + j1*bin_size;
        string st2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
        int e2 = start2 + j2*bin_size;
        string en2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();

        string out_string = chr + "_" + st1 + "_" + en1 + "_" + st2 + "_" + en2;

        top_hits.push_back(out_string);
        top_sums.push_back(max_sum);

      } else {
        check = 0;
        break;
      }
    }

    int stop_test = 0;

    int overlap_check = 0;
    if (top_hits.size() > 1) {
      for (int i=0; i < (top_hits.size() - 1); i++) {
        vector<std::string> new_array1;
	split_string(new_array1,top_hits[i],'_');
        int s1 = (stoi(new_array1[1]) - start1)/bin_size;
        int e1 = (stoi(new_array1[2]) - start1)/bin_size;
        int s2 = (stoi(new_array1[3]) - start2)/bin_size;
        int e2 = (stoi(new_array1[4]) - start2)/bin_size;

        for (int j=(i + 1); j < top_hits.size(); j++) {
          vector<std::string> new_array2;
	  split_string(new_array2,top_hits[j],'_');
          int s3 = (stoi(new_array2[1]) - start1)/bin_size;
          int e3 = (stoi(new_array2[2]) - start1)/bin_size;
          int s4 = (stoi(new_array2[3]) - start2)/bin_size;
          int e4 = (stoi(new_array2[4]) - start2)/bin_size;

          int ot1 = 0;
          if (s1 < s3) {
            if (e1 >= s3) {
              ot1 = 1;
            }
          } else {
            if (s1 <= e3) {
              ot1 = 1;
            }
          }

          int ot2 = 0;
          if (s2 < s4) {
            if (e2 >= s4) {
              ot2 = 1;
	    }
          } else {
            if (s2 <= e4) {
              ot2 = 1;
            }
          }

          if ((ot1 == 1) &&
              (ot2 == 1)) {
            overlap_check = 1;
          }

        }
      }

    }

    if (overlap_check == 1) {

      int set = 0;

      while ((set + 1) < top_hits.size()) {

        vector<vector<double> > new_odds_matrix = ori_odds_matrix;
        vector<vector<double> > test_odds_matrix = ori_odds_matrix;

        vector<string> new_top_hits;
        vector<double> new_top_sums;

        if (set > 0) {
          for (int i = 0; i < set; i++) {
            new_top_hits.push_back(top_hits[i]);
            new_top_sums.push_back(top_sums[i]);

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = (stoi(new_array1[1]) - start1)/bin_size;
            int e1 = (stoi(new_array1[2]) - start1)/bin_size;
            int s2 = (stoi(new_array1[3]) - start2)/bin_size;
            int e2 = (stoi(new_array1[4]) - start2)/bin_size;

            for (int k = s1; k <= e1; k++) {
              for (int j = s2; j <= e2; j++) {

                double change = -5;
                new_odds_matrix[k][j] = change;
                test_odds_matrix[k][j] = change;

              }
            }
          }
        }

	for (int i=0; i < top_hits.size(); i++) {

          if (i == (set + 1)) {

            vector<std::string> new_array1;
	    split_string(new_array1,top_hits[i],'_');
            int s1 = (stoi(new_array1[1]) - start1)/bin_size;
            int e1 = (stoi(new_array1[2]) - start1)/bin_size;
            int s2 = (stoi(new_array1[3]) - start2)/bin_size;
            int e2 = (stoi(new_array1[4]) - start2)/bin_size;

            for (int k = s1; k <= e1; k++) {
              for (int j = s2; j <= e2; j++) {

                double change = -5;
                new_odds_matrix[k][j] = change;

              }
            }
          }
        }

        int new_check = 1;
        int pass = 0;
        double t_thresh = thresh;
        while (new_check == 1) {

          double max_sum = 0;

          int k1;
          int k2;
          int j1;
          int j2;

          for (int i=0; i < odds_matrix[0].size(); i++) {

            vector<double> temp(odds_matrix.size(),0);

            for (int j=i; j < odds_matrix[0].size(); j++) {

              for (int k=0; k < odds_matrix.size(); k++) {

                if (pass == 0) {

                  temp[k] += new_odds_matrix[k][j];

                } else {

                  temp[k] += test_odds_matrix[k][j];

                }

              }

              double temp_sum = 0;
              double temp_max = 0;

              int t_k1 = 0;
              int t_k2 = 0;
	      
	      int f_k1;
              int f_k2;

              for (int k=0; k < odds_matrix.size(); k++) {

                double new_sum = temp_sum + temp[k];

                if (new_sum < 0) {
                  temp_sum = 0;
                  t_k1 = k + 1;
                } else {
                  temp_sum = new_sum;
                  t_k2 = k;
                  if (temp_sum > temp_max) {
                    temp_max = temp_sum;
                    f_k1 = t_k1;
                    f_k2 = t_k2;
                  }
                }
              }

              if (temp_max > max_sum) {
                max_sum = temp_max;
                j1 = i;
                j2 = j;
                k1 = f_k1;
                k2 = f_k2;
              }

            }

          }
	   
	  if (max_sum > t_thresh) {

            t_thresh += thresh;

            for (int k = k1; k <= k2; k++) {
              for (int j = j1; j <= j2; j++) {
                double change = -5;
                test_odds_matrix[k][j] = change;
              }
            }

            int s1 = start1 + k1*bin_size;
            string st1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
            int e1 = start1 + k2*bin_size;
            string en1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
            int s2 = start2 + j1*bin_size;
            string st2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
            int e2 = start2 + j2*bin_size;
            string en2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();

            string out_string = chr + "_" + st1 + "_" + en1 + "_" + st2 + "_" + en2;

            new_top_hits.push_back(out_string);
            new_top_sums.push_back(max_sum);

            pass = 1;

          } else {
            break;
          }

        }

        double all_sum = 0;
        double all_new_sum = 0;

        for (int i = 0; i < top_sums.size(); i++) {
          all_sum += top_sums[i];
        }

        for (int i=0; i < new_top_sums.size(); i++) {
          all_new_sum += new_top_sums[i];
        }

        if (all_new_sum > all_sum) {
          top_hits = new_top_hits;
          top_sums = new_top_sums;

        } else {
          set++;
        }

	int final_overlap_check = 0;
        for (int i=set; i < (top_hits.size() - 1); i++) {
          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = (stoi(new_array1[1]) - start1)/bin_size;
          int e1 = (stoi(new_array1[2]) - start1)/bin_size;
          int s2 = (stoi(new_array1[3]) - start2)/bin_size;
          int e2 = (stoi(new_array1[4]) - start2)/bin_size;

          for (int j=(i + 1); j < top_hits.size(); j++) {
            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');
            int s3 = (stoi(new_array2[1]) - start1)/bin_size;
            int e3 = (stoi(new_array2[2]) - start1)/bin_size;
            int s4 = (stoi(new_array2[3]) - start2)/bin_size;
            int e4 = (stoi(new_array2[4]) - start2)/bin_size;

            int ot1 = 0;
            if (s1 < s3) {
              if (e1 >= s3) {
                ot1 = 1;
              }
            } else {
              if (s1 <= e3) {
                ot1 = 1;
              }
            }

            int ot2 = 0;
            if (s2 < s4) {
              if (e2 >= s4) {
                ot2 = 1;
              }
            } else {
              if (s2 <= e4) {
                ot2 = 1;
              }
            }

            if ((ot1 == 1) &&
                (ot2 == 1)) {
              final_overlap_check = 1;
            }
          }
        }

	if (final_overlap_check == 0) {
          break;
        }

	

      }
    }

    while (stop_test == 0) {

      if (top_hits.size() > 1) {

        int over_check = 0;
        int over_max = 0;
        int over_i = 0;
        int over_j = 0;

        double all_sum = 0;
        vector<vector<double> > all_odds_matrix = ori_odds_matrix;

        for (int i=0; i < (top_hits.size() - 1); i++) {

          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = (stoi(new_array1[1]) - start1)/bin_size;
          int e1 = (stoi(new_array1[2]) - start1)/bin_size;
          int s2 = (stoi(new_array1[3]) - start2)/bin_size;
          int e2 = (stoi(new_array1[4]) - start2)/bin_size;

          for (int k = s1; k <= e1; k++) {
            for (int j = s2; j <= e2; j++) {
              all_sum += all_odds_matrix[k][j];
              all_odds_matrix[k][j] = -5;
            }
          }
        }

	for (int i=0; i < (top_hits.size() - 1); i++) {

          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[i],'_');
          int s1 = (stoi(new_array1[1]) - start1)/bin_size;
          int e1 = (stoi(new_array1[2]) - start1)/bin_size;
          int s2 = (stoi(new_array1[3]) - start2)/bin_size;
          int e2 = (stoi(new_array1[4]) - start2)/bin_size;

          for (int j=(i + 1); j < top_hits.size(); j++) {

            vector<std::string> new_array2;
	    split_string(new_array2,top_hits[j],'_');
            int s3 = (stoi(new_array2[1]) - start1)/bin_size;
            int e3 = (stoi(new_array2[2]) - start1)/bin_size;
            int s4 = (stoi(new_array2[3]) - start2)/bin_size;
            int e4 = (stoi(new_array2[4]) - start2)/bin_size;

            int new_start1;
            if (s1 < s3) {
              new_start1 = s1;
            } else {
              new_start1 = s3;
            }
            int new_end1;
            if (e1 > e3) {
              new_end1 = e1;
            } else {
              new_end1 = e3;
            }
            int new_start2;
            if (s2 < s4) {
              new_start2 = s2;
            } else {
              new_start2 = s4;
            }
            int new_end2;
            if (e2 > e4) {
              new_end2 = e2;
            } else {
              new_end2 = e4;
            }

	    vector<vector<double> > all_new_odds_matrix = ori_odds_matrix;

            double all_new_sum = 0;
            double odds_sum = 0;
            for (int k = new_start1; k <= new_end1; k++) {
              for (int l = new_start2; l <= new_end2; l++) {
                all_new_sum += all_new_odds_matrix[k][l];
                all_new_odds_matrix[k][l] = 0;
                odds_sum += ori_odds_matrix[k][l];
              }
            }

            for (int k=0; k < top_hits.size(); k++) {
              if ((k == i) ||
                  (k == j)) {
                continue;
              }
              vector<std::string> new_array3;
	      split_string(new_array3,top_hits[k],'_');
              int s5 = (stoi(new_array3[1]) - start1)/bin_size; //seems to be an error here?
              int e5 = (stoi(new_array3[2]) - start1)/bin_size;
              int s6 = (stoi(new_array3[3]) - start2)/bin_size;
              int e6 = (stoi(new_array3[4]) - start2)/bin_size;

              for (int l = s5; l <= e5; l++) {
                for (int m = s6; m <= e6; m++) {
                  all_new_sum += all_new_odds_matrix[l][m];
                  all_new_odds_matrix[l][m] = 0;
                }
              }
            }

            if ((all_new_sum >= over_max) &&
                (odds_sum >= thresh)) {
              over_max = all_new_sum;
              over_i = i;
              over_j = j;
              over_check++;
            }

          }
        }

	if (over_check > 0) {

          vector<std::string> new_array1;
	  split_string(new_array1,top_hits[over_i],'_');
          int s1 = stoi(new_array1[1]);
          int e1 = stoi(new_array1[2]);
          int s2 = stoi(new_array1[3]);
          int e2 = stoi(new_array1[4]);
          vector<std::string> new_array2;
	  split_string(new_array2,top_hits[over_j],'_');
          int s3 = stoi(new_array2[1]);
          int e3 = stoi(new_array2[2]);
          int s4 = stoi(new_array2[3]);
          int e4 = stoi(new_array2[4]);
          int ns1;
          string new_start1;
          if (s1 <= s3) {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << s1) )->str();
            ns1 = (s1 - start1)/bin_size;
          } else {
            new_start1 = static_cast<ostringstream*>( &(ostringstream() << s3) )->str();
            ns1 = (s3 - start1)/bin_size;
          }
          int ne1;
          string new_end1;
          if (e1 >= e3) {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << e1) )->str();
            ne1 = (e1 - start1)/bin_size;
          } else {
            new_end1 = static_cast<ostringstream*>( &(ostringstream() << e3) )->str();
            ne1 = (e3 - start1)/bin_size;
          }
          int ns2;
          string new_start2;
          if (s2 <= s4) {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << s2) )->str();
            ns2 = (s2 - start2)/bin_size;
          } else {
            new_start2 = static_cast<ostringstream*>( &(ostringstream() << s4) )->str();
            ns2 = (s4 - start2)/bin_size;
          }
          int ne2;
          string new_end2;
          if (e2 >= e4) {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << e2) )->str();
            ne2 = (e2 - start2)/bin_size;
          } else {
            new_end2 = static_cast<ostringstream*>( &(ostringstream() << e4) )->str();
            ne2 = (e4 - start2)/bin_size;
          }

	  double new_sum = 0;
          for (int k = ns1; k < ne1; k++) {
            for (int l = ns2; l < ne2; l++) {
              new_sum += ori_odds_matrix[k][l];
              odds_matrix[k][l] = -5;

            }
          }

          string new_string = new_array1[0] + "_" + new_start1 + "_" + new_end1 + "_" + new_start2 + "_" + new_end2;

          vector<string> new_top_hits;
          vector<double> new_top_sums;

          new_top_hits.push_back(new_string);
          new_top_sums.push_back(new_sum);
          for (int i=0; i < top_hits.size(); i++) {
            if ((i != over_i) &&
                (i != over_j)) {
              new_top_hits.push_back(top_hits[i]);
              new_top_sums.push_back(top_sums[i]);
            }
          }
          top_hits = new_top_hits;
          top_sums = new_top_sums;

        } else {
          stop_test = 1;
        }
      } else {
        stop_test = 1;
      }

    }

    for (int i=0; i < top_hits.size(); i++) {

      vector<std::string> new_array;
      split_string(new_array,top_hits[i],'_');
      int s1 = (stoi(new_array[1]) - start1)/bin_size;
      int e1 = (stoi(new_array[2]) - start1)/bin_size;
      int s2 = (stoi(new_array[3]) - start2)/bin_size;
      int e2 = (stoi(new_array[4]) - start2)/bin_size;

      unordered_map<int,double> i_sum;
      unordered_map<int,double> j_sum;
      unordered_map<int,int> i_N;
      unordered_map<int,int> j_N;
      vector<double> all_diag1;
      vector<double> all_diag2;
      vector<int> all_diag1_N;
      vector<int> all_diag2_N;
      unordered_map<int,double> diag1_sum;
      unordered_map<int,double> diag2_sum;
      unordered_map<int,int> diag1_N;
      unordered_map<int,int> diag2_N;
      unordered_map<int,double> i_log_sum;
      unordered_map<int,double> j_log_sum;
      unordered_map<int,double> diag1_log_sum;
      unordered_map<int,double> diag2_log_sum;
      unordered_map<int,double> diag1_bin_sum;
      unordered_map<int,double> diag2_bin_sum;

      double new_odds_sum = 0;
      for (int k = s1; k <= e1; k++) {

        int i_bin = k - s1;
        i_N[i_bin]++;

        for (int j = s2; j <= e2; j++) {
          new_odds_sum += ori_odds_matrix[k][j];

          int j_bin = j - s2;
          int rev_j_bin = e2 - j;
          j_N[j_bin]++;

          int diag1_bin = i_bin + j_bin;
          int diag2_bin = i_bin + rev_j_bin;

          diag1_N[diag1_bin]++;
          diag2_N[diag2_bin]++;

          all_diag1_N.push_back(diag1_bin);
          all_diag2_N.push_back(diag2_bin);

          i_sum[i_bin] += norm_matrix[k][j];
          j_sum[j_bin] += norm_matrix[k][j];
          i_log_sum[i_bin] += log(norm_matrix[k][j] + 0.1);
          j_log_sum[j_bin] += log(norm_matrix[k][j] + 0.1);
          diag1_sum[diag1_bin] += norm_matrix[k][j];
          diag2_sum[diag2_bin] += norm_matrix[k][j];
          diag1_log_sum[diag1_bin] += log(norm_matrix[k][j] + 0.1);
          diag2_log_sum[diag2_bin] += log(norm_matrix[k][j] + 0.1);
          diag1_bin_sum[diag1_bin] += 1;
          diag2_bin_sum[diag2_bin] += 1;
          all_diag1.push_back(norm_matrix[k][j]);
          all_diag2.push_back(norm_matrix[k][j]);

        }
      }

      double i_cor = get_corr_map(i_sum,i_N);
      double j_cor = get_corr_map(j_sum,j_N);
      double i_log_cor = get_corr_map(i_log_sum,i_N);
      double j_log_cor = get_corr_map(j_log_sum,j_N);
      double diag1_cor = get_corr_map(diag1_sum,diag1_N);
      double diag2_cor = get_corr_map(diag2_sum,diag2_N);
      double diag1_log_cor = get_corr_map(diag1_log_sum,diag1_N);
      double diag2_log_cor = get_corr_map(diag2_log_sum,diag2_N);
      double diag1_bin_cor = get_corr_map(diag1_bin_sum,diag1_N);
      double diag2_bin_cor = get_corr_map(diag2_bin_sum,diag2_N);
      double all_diag1_cor = get_corr_vec(all_diag1,all_diag1_N);
      double all_diag2_cor = get_corr_vec(all_diag2,all_diag2_N);

      double final_i = all_diag1_cor + all_diag2_cor + diag1_bin_cor + diag2_bin_cor + diag1_cor + diag2_cor + diag1_log_cor + diag2_log_cor + i_cor + i_log_cor;
      double final_j = all_diag1_cor + -1*all_diag2_cor + diag1_bin_cor + -1*diag2_bin_cor + diag1_cor + -1*diag2_cor + diag1_log_cor + -1*diag2_log_cor + j_cor + j_log_cor;

      string strand1;
      string strand2;

      if (final_i > 0) {
        strand1 = "+";
      } else {
        strand1 = "-";
      }

      if (final_j > 0) {
        strand2 = "+";
      } else {
        strand2 = "-";
      }

      break_output << new_odds_sum << "\t" << new_array[0] << "\t" << new_array[1] << "\t" << (stoi(new_array[2]) + bin_size) << "\t" << strand1 << "\t" << new_array[0] << "\t" << new_array[3] << "\t" << (stoi(new_array[4]) + bin_size) << "\t" << strand2 << "\n";

    }

  }

  break_output.close();

}
	
void get_10kb_odds_ratio_intra (string super_matrix_file,       \
                                string bias_vector_file,        \
                                string block_file,              \
                                int bin_size,                   \
                                int max_size,                   \
                                string name) {

  vector<string> block_data = read_blocks(block_file,bin_size);

  unordered_map<string,vector<int> > ref_map;
  unordered_map<string,unordered_map<int,double> > bias_map;
  read_bias_vector(bias_vector_file,ref_map,bias_map);

  unordered_map<string,int> data_map;
  read_super_matrix_intra(super_matrix_file,data_map);

  unordered_map<int,double> m_hash;
  unordered_map<int,double> r_hash;
  unordered_map<int,int> N_hash;
  int max_weight = 0;
  double submatrix = 0;
  find_parameters_10kb_intra(m_hash,r_hash,N_hash,data_map,ref_map,bias_map,bin_size,max_size,max_weight,submatrix);

  double prior = 0.000001;
  double corr = log(prior/(1 - prior));
  double thresh = -1*(log(0.05/submatrix) + corr);

  cerr << thresh << "\n";

  find_breaks_10kb_intra(block_data,m_hash,r_hash,N_hash,data_map,ref_map,bias_map,bin_size,max_size,max_weight,thresh,name);

}

void read_break_file (string break_file,\
                      vector<string>& break_vector,\
                      vector<string>& size_vector,\
		      string size) {

  cerr << "Reading in break file now... ";

  string line;
  ifstream myfile (break_file);
  while ( getline(myfile,line) ) {
    break_vector.push_back(line);
    size_vector.push_back(size);
  }

  cerr << "Done\n";;

}

void merge_break_file (string break_file,		\
                       int bin_size,                    \
                       vector<string>& break_vector,    \
                       vector<string>& size_vector,
                       string size) {
  
  cerr << "Reading in break file now... ";

  vector<string> final_break_vector;
  vector<string> final_size_vector;
  vector<string> new_break_vector;

  string line;
  ifstream myfile (break_file);
  while ( getline(myfile,line) ) {
    final_break_vector.push_back(line);
    new_break_vector.push_back(line);
    final_size_vector.push_back(size);
  }

  int i=0;

  vector<string> ori_break_vector = break_vector;

  while (i < break_vector.size()) {

    vector<string> array1;
    split_string(array1,break_vector[i],'\t'); 

    string chr1 = array1[1];
    string chr2 = array1[5];
    
    int start1 = stoi(array1[2]) - bin_size;
    int end1 = stoi(array1[3]) + bin_size;
    int start2 = stoi(array1[6]) - bin_size;
    int end2 = stoi(array1[7]) + bin_size;

    int test = 0;
    
    for (int j=(i+1); j < break_vector.size(); j++) {

      vector<string> array2;
      split_string(array2,break_vector[j],'\t');

      string chr3 = array2[1];
      string chr4 = array2[5];

      int start3 = stoi(array2[2]);
      int end3 = stoi(array2[3]);
      int start4 = stoi(array2[6]);
      int end4 = stoi(array2[7]);

      int over_check = 0;

      if ((chr1 == chr3) &&
          (chr2 == chr4)) {
	
	if (start1 < start3) {
	  if (end1 > start3) {
	    if (start2 < start4) {
	      if (end2 > start4) {
		over_check = 1;
	      }
	    } else {
	      if (start2 < end4) {
		over_check = 1;
	      }
	    }
	  }
	} else {
	  if (start1 < end3) {
	    if (start2 < start4) {
	      if (end2 > start4) {
		over_check = 1;
	      }
	    } else {
	      if (start2 < end4) {
		over_check = 1;
	      }
	    }
	  }
	}
	
	if (over_check == 1) {
	  test = 1;
	  
	  vector<string> new_temp_array;
	  
	  for (int k=0; k < break_vector.size(); k++) {
	    if ((k != i) &&
		(k != j)) {
	      new_temp_array.push_back(break_vector[k]);
	    }
	  }
	  
	  string new_start1;
	  if (start1 < start3) {
	    new_start1 = static_cast<ostringstream*>( &(ostringstream() << start1) )->str();
	  } else {
	    new_start1 = static_cast<ostringstream*>( &(ostringstream() << start3) )->str();
	  }
	  
	  string new_end1;
	  if (end1 > end3) {
	    new_end1 = static_cast<ostringstream*>( &(ostringstream() << end1) )->str();
	  } else {
	    new_end1 = static_cast<ostringstream*>( &(ostringstream() << end3) )->str();
	  }
	  
	  string new_start2;
	  if (start2 < start4) {
	    new_start2 = static_cast<ostringstream*>( &(ostringstream() << start2) )->str();
	  } else {
	    new_start2 = static_cast<ostringstream*>( &(ostringstream() << start4) )->str();
	  }
	  
	  string new_end2;
	  if (end2 > end4) {
	    new_end2 = static_cast<ostringstream*>( &(ostringstream() << end2) )->str();
	  } else {
	    new_end2 = static_cast<ostringstream*>( &(ostringstream() << end4) )->str();
	  }
	  
	  string new_string = "-1\t" + chr1 + "\t" + new_start1 + "\t" + new_end1 + "\t-1\t" + chr2 + "\t" + new_start2 + "\t" + new_end2 + "\t-1";
	  
	  new_temp_array.push_back(new_string);
	  break_vector = new_temp_array;
	  break;
	}
      }
      
    }
    
    if (test == 0) {
      i++;
    } else {
      i=0;
    }

  }

  for (int i=0; i < break_vector.size(); i++) {
    int check = 0;
    vector<string> array1;
    split_string(array1,break_vector[i],'\t');

    string chr1 = array1[1];
    string chr2 = array1[5];

    int start1 = stoi(array1[2]);
    int end1 = stoi(array1[3]);
    int start2 = stoi(array1[6]);
    int end2 = stoi(array1[7]);

    for (int j=0; j < new_break_vector.size(); j++) {
      vector<string> array2;
      split_string(array2,new_break_vector[j],'\t');

      string chr3 = array2[1];
      string chr4 = array2[5];

      int start3 = stoi(array2[2]);
      int end3 = stoi(array2[3]);
      int start4 = stoi(array2[6]);
      int end4 = stoi(array2[7]);
      
      if ((chr1 == chr3) &&
	  (chr2 == chr4)) {
	if ((array1[1] == array2[1]) &&
	    (array1[5] == array2[5])) {
	  if ((start3 >= start1) &&
	      (end3 <= end1) &&
	      (start4 >= start2) &&
	      (end4 <= end2)) {
	    check = 1;
	  }
	}
      }
    }
    if (check == 0) {
      for (int k = 0; k < ori_break_vector.size(); k++) {
	
	vector<string> array2;
	split_string(array2,ori_break_vector[k],'\t');

	string chr3 = array2[1];
	string chr4 = array2[5];

	int start3 = stoi(array2[2]);
	int end3 = stoi(array2[3]);
	int start4 = stoi(array2[6]);
	int end4 = stoi(array2[7]);

	if ((chr1== chr3) &&
	    (chr2== chr4)) {

	  if ((array1[1] == array2[1]) &&
	      (array1[5] == array2[5])) {
	    if ((start3 >= start1) &&
		(end3 <= end1) &&
		(start4 >= start2) &&
		(end4 <= end2)) {
	  
	      final_break_vector.push_back(ori_break_vector[k]);
	      final_size_vector.push_back(size_vector[k]);

	    }
	  }
	}
      }
    }
  }

  break_vector = final_break_vector;
  size_vector = final_size_vector;

  cerr << "Done\n";
  
}

void merge_output (string name,\
		   int min_1kb_flag) {

  vector<string> inter_break_vector;
  vector<string> inter_size_vector;

  string inter_1Mb_file = name + "_1Mb.breaks.txt";
  string size_1Mb = "1Mb";

  read_break_file(inter_1Mb_file,inter_break_vector,inter_size_vector,size_1Mb);
  
  string inter_100kb_file = name + "_100kb.breaks.txt";
  string size_100kb = "100kb";

  merge_break_file(inter_100kb_file,1000000,inter_break_vector,inter_size_vector,size_100kb);

  string inter_10kb_file = name + "_10kb.breaks.txt";
  string size_10kb = "10kb";

  merge_break_file(inter_10kb_file,100000,inter_break_vector,inter_size_vector,size_10kb);

  if (min_1kb_flag == 1) {
    string inter_1kb_file = name + "_1kb.breaks.txt";
    string size_1kb = "1kb";
    merge_break_file(inter_1kb_file,10000,inter_break_vector,inter_size_vector,size_1kb);
  }

  vector<string> intra_break_vector;
  vector<string> intra_size_vector;

  string intra_100kb_file = name + "_100kb_intra.breaks.txt";

  read_break_file(intra_100kb_file,intra_break_vector,intra_size_vector,size_100kb);

  string intra_10kb_file = name + "_10kb_intra.breaks.txt";

  merge_break_file(intra_10kb_file,100000,intra_break_vector,intra_size_vector,size_10kb);

  if (min_1kb_flag == 1) {
    string intra_1kb_file = name + "_1kb_intra.breaks.txt";
    string size_1kb = "1kb";
    merge_break_file(intra_1kb_file,10000,intra_break_vector,intra_size_vector,size_1kb);
  }

  ofstream break_output;
  string break_output_file_name = name + ".breaks.txt";

  break_output.open(break_output_file_name);

  for (int i=0; i < inter_break_vector.size(); i++) {
    break_output << inter_break_vector[i] << "\t" << inter_size_vector[i] << "\n";
  }

  for (int i=0; i < intra_break_vector.size(); i++) {
    break_output << intra_break_vector[i] << "\t" << intra_size_vector[i] << "\n";
  }

  break_output.close();
  
}

int main(int argc, char **argv) {

  int c;
  char* bam_file = NULL;
  char* out_name = NULL;
  char* exp_file_inter = NULL;
  char* exp_file_intra = NULL;
  int min_flag = 0;

  if (argc == 1) {
    usage();
    return(0);
  }

  while (1) {

    static struct option long_options[] = {
      /* These options set a flag. */
      {"bam-file",required_argument, 0, 'b'},
      {"name",required_argument, 0, 'q'},
      {"exp-file-inter",required_argument, 0, 'e'},
      {"exp-file-intra",required_argument, 0, 'p'},
      {"min-1kb",no_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
    
    int option_index = 0;

    c = getopt_long (argc, argv, "abc:d:f:",
		     long_options, &option_index);
    
    if (c == -1) {
      break;
    }

    switch (c)
      {
	
      case 'b':
	bam_file = optarg;
	break;
	
      case 'q':
        out_name = optarg;
        break;

      case 'e':
        exp_file_inter = optarg;
        break;

      case 'p':
        exp_file_intra = optarg;
        break;

      case 'm':
	min_flag = 1;
        break;

      case '?':
	break;

      default:
	abort ();
      }
  }

  if ((bam_file == NULL) ||
      (out_name == NULL) ||
      (exp_file_inter == NULL) ||
      (exp_file_intra == NULL)) {
    usage();
    return(0);
  }

  string name(out_name);
  string name_1Mb = name + "_1Mb";
  string name_100kb = name + "_100kb";
  string name_100kb_intra = name + "_100kb_intra";
  string name_10kb = name + "_10kb";
  string name_10kb_intra = name + "_10kb_intra";
  string name_1kb = name + "_1kb";
  string name_1kb_intra = name + "_1kb_intra";

  get_super_matrix(bam_file,name_1Mb,1000000,1,1,1,0.05);
  get_super_matrix(bam_file,name_100kb,100000,1,1,1,0.05);
  get_super_matrix(bam_file,name_10kb,10000,1,1,1,0.05);

  if (min_flag == 1) {
    get_super_matrix(bam_file,name_1kb,1000,1,1,1,0.05);
  }

  string super_matrix_1Mb = name_1Mb + ".super_matrix.txt";
  string super_matrix_100kb = name_100kb + ".super_matrix.txt";
  string super_matrix_10kb = name_10kb + ".super_matrix.txt";
  string super_matrix_1kb = name_1kb + ".super_matrix.txt";
  string SR_100kb = name_100kb + ".SR.txt";
  string SR_10kb = name_10kb + ".SR.txt";
  string SR_1kb = name_1kb + ".SR.txt";

  iterative_correction_no_SR(super_matrix_1Mb,name_1Mb,50,1);
  iterative_correction(super_matrix_100kb,SR_100kb,name_100kb,50);
  iterative_correction(super_matrix_10kb,SR_10kb,name_10kb,50);
  
  if (min_flag == 1) {
    iterative_correction(super_matrix_1kb,SR_1kb,name_1kb,50);
  }

  string super_matrix_1Mb_norm = name_1Mb + ".super_matrix.norm.txt";
  
  get_eigen(super_matrix_1Mb_norm,name_1Mb);

  string bias_vector_1Mb = name_1Mb + ".bias_vector.txt";
  string pc1_file = name_1Mb + ".super_matrix.norm.eig.txt";

  get_1Mb_odds_ratio_inter(super_matrix_1Mb,bias_vector_1Mb,exp_file_inter,pc1_file,1000000,10000000,name_1Mb);

  string bias_vector_100kb = name_100kb + ".bias_vector.txt";
  string block_name_1Mb = name_1Mb + ".breaks.txt";

  get_100kb_odds_ratio_inter(super_matrix_100kb,bias_vector_100kb,block_name_1Mb,100000,1000000,name_100kb);

  string bias_vector_10kb = name_10kb + ".bias_vector.txt";
  string block_name_100kb = name_100kb + ".breaks.txt";

  get_100kb_odds_ratio_inter(super_matrix_10kb,bias_vector_10kb,block_name_100kb,10000,100000,name_10kb);

  if (min_flag == 1) {
    
    string bias_vector_1kb = name_1kb + ".bias_vector.txt";
    string block_name_10kb = name_10kb + ".breaks.txt";

    get_100kb_odds_ratio_inter(super_matrix_1kb,bias_vector_1kb,block_name_10kb,1000,10000,name_1kb);

  }
  
  get_100kb_odds_ratio_intra(super_matrix_100kb,bias_vector_100kb,exp_file_intra,100000,1000000,name_100kb_intra);

  string block_name_100kb_intra = name_100kb_intra + ".breaks.txt";

  get_10kb_odds_ratio_intra(super_matrix_10kb,bias_vector_10kb,block_name_100kb_intra,10000,100000,name_10kb_intra);

  if (min_flag == 1) {

    string bias_vector_1kb = name_1kb + ".bias_vector.txt";
    string block_name_10kb_intra = name_10kb_intra + ".breaks.txt";

    get_10kb_odds_ratio_intra(super_matrix_1kb,bias_vector_1kb,block_name_10kb_intra,1000,10000,name_1kb_intra);
    
  }

  merge_output(name,min_flag);

}
