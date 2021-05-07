/*
Commands to compile and run this code:

Compile: g++ -o tv2.o -fopenmp -std=gnu++11 openmp_version2.cpp
Run: ./tv2.o

*/

#include<bits/stdc++.h>
#include<fstream>
#include<omp.h>
using namespace std;

typedef basic_string<unsigned char> ustring;

const int bits_unit = 8;

struct node {
    int frequency;
    unsigned char in;
    node * left;
    node * right;
    node(int freq, unsigned char i)
    {
        frequency = freq;
        in = i;
        left = NULL;
        right = NULL;
    }
   
};

struct comparator{

    bool operator()(node* temp1, node* temp2)
    {
        if(temp1->frequency > temp2->frequency)
            return true;
        else 
            return false;
    }
};

void obtain_huffmancodes(node* root, int byte, pair<int, int>* encode, int cur_bit_id = 0)
{
    if(root->left == NULL && root->right == NULL)
    {
        encode[root->in].first = byte;
        encode[root->in].second = cur_bit_id;
        assert(cur_bit_id <= 32);
        return;
    }
    int nxt_id = cur_bit_id + 1;
    int nxt_byte = byte | (1 << cur_bit_id);
    obtain_huffmancodes(root->left, byte, encode, nxt_id);
    obtain_huffmancodes(root->right, nxt_byte, encode, nxt_id);
}

int main()
{
    cout << "Enter the file name with extension (e.g. bonsai.raw)" << endl;
    string filename;
    cin >> filename;
    string decompressed_file_name = "output_" + filename;

    // Stores the no. of times we will repeat the whole procedure to measure the time accurately by taking average in the end.
    int runs = 10;
    long double freq_table_time, huffman_tree_time, hashtable_time, compression_time, decompression_time, whole_time;
    freq_table_time = huffman_tree_time = hashtable_time = compression_time = decompression_time = whole_time = 0;
    bool printed = false;
    
    // Experiments for different chunk sizes
    for(int chunk_size = 1 << 8; chunk_size < 1 << 20; chunk_size <<= 2) {
        cout << "Chunk Size: " << chunk_size << endl;
        
        // Experiments for different number of processors
        for(int processors = 1; processors <= 16; processors *= 2) {
            cout << "Processors: " << processors << endl;

            // Reset parameters
            freq_table_time = huffman_tree_time = hashtable_time = compression_time = decompression_time = whole_time = 0;
            printed = false;
            for(int run = 0; run < runs; run++){

                struct timespec whole_start, whole_end;
                clock_gettime(CLOCK_MONOTONIC, &whole_start);
                
                ifstream file(filename, ios::in | ios::binary);
                file.seekg (0, file.end);
                int file_size = file.tellg();
                assert(file_size <= (1 << 29));
                file.seekg (0, file.beg);
                unsigned char* buffer = new unsigned char[file_size];
                file.read((char*)buffer, file_size);
                file.close();
                
                // file_size must be divisible by the no. of processors and chunk_size
                assert(file_size % processors == 0 && file_size % chunk_size == 0);
                int part_sz = file_size / processors;
                omp_set_num_threads(processors);
                
                struct timespec start, end;
                // Create frequency table --------------------------------------------------START
                clock_gettime(CLOCK_MONOTONIC, &start);
            
                int freq[256] = {0};
                #pragma omp parallel
                {
                    unsigned char a;
                    int freq_temp[256];
                    memset(freq_temp, 0, sizeof freq_temp);
                    int thread_num = omp_get_thread_num();
                    int thread_chunk_offset = part_sz * thread_num;
                    for(int i = 0; i < part_sz; i++){
                        a = buffer[thread_chunk_offset++];
                        freq_temp[a]++;
                    }
                    #pragma omp critical
                    {
                        for(int i = 0; i < 256; i++)
                            freq[i] += freq_temp[i];
                    }
                }
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time_taken;
                time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                freq_table_time += time_taken;
                // Create frequency table -----------------------------------------------------------------END
    
                // Create huffman tree -------------------------------------------------------------------START
                clock_gettime(CLOCK_MONOTONIC, &start);
            
                priority_queue<node*, vector<node*>, comparator> container;
            
                for(int i = 0; i < 256; i++)
                {
                    if(freq[i] > 0) {
                        node* temp = new node(freq[i], i);
                        container.push(temp);
                    }
                }
                
                while(container.size() > 1)
                {
                    node* left = container.top();
                    container.pop();
            
                    node * right = container.top();
                    container.pop();
            
                    node* curr = new node(left->frequency + right->frequency, '0');
                    curr->left = left;
                    curr->right = right;
                    container.push(curr);
                }
                
                assert(container.size() == 1);
                node* root = container.top();
                container.pop();
                node* root_store = root;
                
                clock_gettime(CLOCK_MONOTONIC, &end);
                time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                huffman_tree_time += time_taken;
                // Create huffman tree----------------------------------------------------------------------------END
        
                // Create Hashtable of encodes------------------------------------------------------------------START
                clock_gettime(CLOCK_MONOTONIC, &start);
            
                int ans = 0;
                pair<int, int> encode[256];
                obtain_huffmancodes(root, ans, encode);

                clock_gettime(CLOCK_MONOTONIC, &end);
                time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                hashtable_time += time_taken;
                // Create Hashtable of encodes-----------------------------------------------------------------------END
    
                // Compression -----------------------------------------------------------------------------------------START
                clock_gettime(CLOCK_MONOTONIC, &start);
                int no_of_chunks = file_size / chunk_size;
                vector<ustring> comp_chunks(no_of_chunks);
                vector<unsigned int> comp_chunk_size(no_of_chunks);
                
                #pragma omp parallel
                {
                    #pragma omp for schedule(dynamic, 1)
                    for(int k = 0; k < no_of_chunks; k++) {
                        int idx = 7;
                        unsigned char a, c;
                        unsigned int total_bits = 0;
                        int thread_offset = chunk_size * k;
                        int cur_byte = 0;
                        c = 0;
                        ustring comp_data;  
                        for(int i = 0; i < chunk_size; i++) {
                            a = buffer[thread_offset++];
                            int byte = encode[a].first;
                            int len = encode[a].second;
                            total_bits += len;
                            for(int j = 0; j < len; j++) {
                                bool bit = byte & (1 << j);
                                c |= (bit << idx);
                                idx--;
                                if(idx < 0) {
                                    comp_data += c;
                                    idx = 7;
                                    // refresh c
                                    c &= 0;
                                }
                            }
                        }
                        if(idx != 7) {
                            comp_data += c;
                        }
                        comp_chunks[k] = comp_data;
                        comp_chunk_size[k] = total_bits;
                    }
                }

                clock_gettime(CLOCK_MONOTONIC, &end);
                time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                compression_time += time_taken;
                
                // Size of comp_chunk_size array, which has to be stored
                // File can be stored and retrieved in this format
                // No of total_cores used + comp_chunk_size array + comp_chunks array + original file size
                if(!printed) {
                    int compressed_size = 4 + 4;
                    for(int i = 0; i < no_of_chunks; i++)
                        compressed_size += comp_chunks[i].length() + 4;
                    cout << "Compressed file size in bytes: " << compressed_size << endl;
                    printed = true;
                }
                // free the buffer
                delete[] buffer;
                // Compression ------------------------------------------------------------------------------------END
          
                // Decompression ---------------------------------------------------------------------------------START
                clock_gettime(CLOCK_MONOTONIC, &start);
                buffer = new unsigned char[file_size];
                #pragma omp parallel private(root)
                {
                    #pragma omp for schedule(dynamic, 1)
                    for(int k = 0; k < no_of_chunks; k++) {
                        unsigned char a;
                        root = root_store;    

                        ustring cur_chunk = comp_chunks[k];
                        unsigned int total_bits = comp_chunk_size[k];
    
                        int thread_offset = chunk_size * k;
                        int idx = 7;
                        int bytes = cur_chunk.length();
                        
                        for(int i = 0; i < bytes && total_bits > 0; i++) {
                            a = cur_chunk[i];
                            idx = 7;
                            while(idx >= 0 && total_bits > 0) {
                                if(a & (1 << idx))
                                    root = root->right;
                                else
                                    root = root->left;
                                if(!root->right && !root->left) {
                                    unsigned char cur = (root->in);
                                    buffer[thread_offset++] = cur;
                                    root = root_store;
                                }
                                idx--;
                                total_bits--;
                            }
                        }
                    }
                }
    
                clock_gettime(CLOCK_MONOTONIC, &end);
                time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                decompression_time += time_taken;

                std::ofstream of_c(decompressed_file_name, std::ios_base::binary);
                of_c.write((char *)buffer, file_size);
                of_c.close();
                delete[] buffer;
                // Decompression---------------------------------------------------------------------------------END

                clock_gettime(CLOCK_MONOTONIC, &whole_end);
                time_taken = (whole_end.tv_sec - whole_start.tv_sec) * 1e9;
                time_taken = (time_taken + (whole_end.tv_nsec - whole_start.tv_nsec)) * 1e-9;
                whole_time += time_taken;
            }
        
            freq_table_time /= double(runs);
            huffman_tree_time /= double(runs);
            hashtable_time /= double(runs);
            compression_time /= double(runs);
            decompression_time /= double(runs);
            whole_time /= double(runs);
            
            cout << fixed << setprecision(20) << "Creating frequency table: " << freq_table_time << endl;
            cout << fixed << setprecision(20) << "Creating Huffman Tree: " << huffman_tree_time << endl;
            cout << fixed << setprecision(20) << "Creating hashtable of encodes: " << hashtable_time << endl;
            cout << fixed << setprecision(20) << "Compression: " << compression_time << endl;
            cout << fixed << setprecision(20) << "Decompression: " << decompression_time << endl;
            cout << fixed << setprecision(20) << "Total time: " << whole_time << endl;
        
        }
    }
    return 0;
}