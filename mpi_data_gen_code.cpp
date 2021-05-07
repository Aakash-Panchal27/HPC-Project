/*
Commands to compile and run this code:

Compile: mpicxx -std=gnu++11 parallel_mpi.cpp
Run: mpirun -machinefile machines -np {no_of_cores} ./a.out

Note: Replace {no_of_cores} with the number of cores you want to run the code with.
e.g. mpirun -machinefile machines -np 16 ./a.out

Please replace "filename" variable at the very start of the main function to the one you want to run the code on.
e.g filename = "aneurism.raw"
*/

#include<bits/stdc++.h>
#include<fstream>
#include<time.h>
#include<mpi.h>
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

int main(int argc, char** argv)
{
    string filename = "bonsai.raw";

    long double freq_table_time, huffman_tree_time, hashtable_time, compression_time, decompression_time, whole_time;
    freq_table_time = huffman_tree_time = hashtable_time = compression_time = decompression_time = whole_time = 0;
    
    struct timespec whole_start, whole_end;
    clock_gettime(CLOCK_MONOTONIC, &whole_start);

    string decompressed_file_name = "output_" + filename;

    int proc_rank, total_cores;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_cores);
    
    MPI_Barrier(MPI_COMM_WORLD);

    ifstream file(filename, ios::in | ios::binary);
    file.seekg (0, file.end);
    int file_size = file.tellg();
    assert(file_size <= (1 << 29));
    file.seekg (0, file.beg);
    unsigned char* buffer = new unsigned char[file_size];
    file.read((char*)buffer, file_size);
    file.close();
    
    int chunk_size = file_size / total_cores;
    
    auto end = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    // Create frequency table in Parallel--------------------------------------------------START
    if(proc_rank == 0)
        start = chrono::high_resolution_clock::now();
    
    unsigned char a;
    int freq[256];
    memset(freq, 0, sizeof freq);
    int thread_chunk_offset = chunk_size * proc_rank;
    for(int i = 0; i < chunk_size; i++){
        a = buffer[thread_chunk_offset++];
        freq[a]++;
    }
    
    if(proc_rank != 0)
        MPI_Send(freq, 256, MPI_INT, 0, 0, MPI_COMM_WORLD);
    
    MPI_Status status;
    
    if(proc_rank == 0) {
        for(int i = 1; i < total_cores; i++) {
            int temp_buff[256];
            MPI_Recv(temp_buff, 256, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            for(int i = 0; i < 256; i++)
                freq[i] += temp_buff[i];
        }
        for(int i = 1; i < total_cores; i++)
            MPI_Send(freq, 256, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
    
    if(proc_rank != 0)
        MPI_Recv(freq, 256, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    
    if(proc_rank == 0) {
        end = chrono::high_resolution_clock::now();
        double total_time =  
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
        freq_table_time += total_time;
    }
    // Create frequency table--------------------------------------------------END
    
    // Create huffman tree-----------------------------------------------------START
    if(proc_rank == 0)
        start = chrono::high_resolution_clock::now();
    
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
    
    if(proc_rank == 0) {
        end = chrono::high_resolution_clock::now();
        double total_time =  
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
        huffman_tree_time += total_time;
    }
    // Create huffman tree--------------------------------------------------END
    
    // Creating Hashtable of encodes----------------------------------------START
    if(proc_rank == 0)
        start = chrono::high_resolution_clock::now();
    
    int ans = 0;
    pair<int, int> encode[256];
    obtain_huffmancodes(root, ans, encode);
    //for(auto it:encode)cout << it.second << endl;
    if(proc_rank == 0) {
        end = chrono::high_resolution_clock::now();
        double total_time =  
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
        hashtable_time += total_time;
    }
    // Creating Hashtable of encodes----------------------------------------END
    
    // Compress ----------------------------------------------------------------START
    if(proc_rank == 0)
        start = chrono::high_resolution_clock::now();
    int no_of_chunks = file_size / chunk_size;
    //vector<ustring> comp_chunks;
    unsigned char **comp_chunks;
    vector<unsigned int> comp_chunk_size;
    
    if(proc_rank == 0) {
        comp_chunk_size.assign(total_cores, 0);
        comp_chunks = new unsigned char*[total_cores];
        for(int i = 0; i < total_cores; i++)
            comp_chunks[i] = new unsigned char[chunk_size];
    }
        
    int idx = 7;
    unsigned char c;
    unsigned int total_bits = 0;
    int thread_offset = chunk_size * proc_rank;
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
    int cur_byte = (total_bits + 7) / 8;
    if(proc_rank != 0) {
        MPI_Send(&total_bits, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(comp_data.c_str(), cur_byte, MPI_UNSIGNED_CHAR, 0, 3, MPI_COMM_WORLD);
    }
    
    if(proc_rank == 0) {
        comp_chunk_size[proc_rank] = total_bits;
        // Receive compressed data from all the cores
        for(int i = 1; i < total_cores; i++) {
            unsigned int temp_total_bits;
            MPI_Recv(&temp_total_bits, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
            int size = (temp_total_bits + 7) / 8;
            MPI_Recv(comp_chunks[i], size, MPI_UNSIGNED_CHAR, i, 3, MPI_COMM_WORLD, &status);
            comp_chunk_size[i] = temp_total_bits;
        }
        // As it is a just data generator code.
        /*int compressed_size = 4 + 4;
        for(int i = 0; i < total_cores; i++) {
            compressed_size += (comp_chunk_size[i] + 7)/8 + 4;
        }
        cout << "Compressed file size in bytes: " << compressed_size << endl;
        */
    }
    
    if(proc_rank == 0) {
        end = chrono::high_resolution_clock::now();
        double total_time =  
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
        compression_time += total_time;
    }
    
    // free the buffer
    delete[] buffer;
    
    // Compressed ----------------------------------------------------------------END
    
    
    // Decompress---------------------------------------------------------------START
    if(proc_rank == 0) {
        start = chrono::high_resolution_clock::now();
        buffer = new unsigned char[file_size];
    }

    root = root_store;
    idx = 7;
    if(proc_rank != 0) {
        unsigned char* decomp_buffer = new unsigned char[chunk_size];
        int byte_no = 0;
        int bytes = (total_bits + 7) / 8;
        for(int i = 0; i < bytes && total_bits > 0; i++) {
            a = comp_data[i];
            idx = 7;
            while(idx >= 0 && total_bits > 0) {
                if(a & (1 << idx))
                    root = root->right;
                else
                    root = root->left;
                if(!root->right && !root->left) {
                    unsigned char cur = (root->in);
                    decomp_buffer[byte_no++] = cur;
                    root = root_store;
                }
                idx--;
                total_bits--;
            }
        }
        // send decompressed
        MPI_Send(decomp_buffer, byte_no, MPI_UNSIGNED_CHAR, 0, 6, MPI_COMM_WORLD);
    }
    else {
        int bytes = (total_bits + 7) / 8;
        int byte_no = 0;
        for(int i = 0; i < bytes && total_bits > 0; i++) {
            a = comp_data[i];
            idx = 7;
            while(idx >= 0 && total_bits > 0) {
                if(a & (1 << idx))
                    root = root->right;
                else
                    root = root->left;
                if(!root->right && !root->left) {
                    unsigned char cur = (root->in);
                    buffer[byte_no++] = cur;
                    root = root_store;
                }
                idx--;
                total_bits--;
            }
        }
        for(int i = 1; i < total_cores; i++) {
            int offset = i * chunk_size;
            MPI_Recv(buffer + offset, chunk_size, MPI_UNSIGNED_CHAR, i, 6, MPI_COMM_WORLD, &status);
        }
    }
    
    if(proc_rank == 0) {
        end = chrono::high_resolution_clock::now();
        double total_time =  
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
        decompression_time += total_time;

        std::ofstream of_c(decompressed_file_name, std::ios_base::binary);
        of_c.write((char *)buffer, file_size);
        of_c.close();
        delete[] buffer;
    }
    
    MPI_Finalize();
    // Decompress-----------------------------------------------------------END

    if(proc_rank == 0) {
        clock_gettime(CLOCK_MONOTONIC, &whole_end);
        double time_taken;
        time_taken = (whole_end.tv_sec - whole_start.tv_sec) * 1e9;
        time_taken = (time_taken + (whole_end.tv_nsec - whole_start.tv_nsec)) * 1e-9;
        whole_time += time_taken;

        cout << freq_table_time << endl << huffman_tree_time << endl << hashtable_time << endl;
        cout << compression_time << endl << decompression_time << endl << whole_time << endl;
        /*
        cout << "Time taken by different steps in seconds:" << endl;
        cout << fixed << setprecision(20) << "Create frequency table: " << freq_table_time << endl;
        cout << fixed << setprecision(20) << "Create Huffman Tree: " << huffman_tree_time << endl;
        cout << fixed << setprecision(20) << "Create hashtable of encodes: " << hashtable_time << endl;
        cout << fixed << setprecision(20) << "Compression: " << compression_time << endl;
        cout << fixed << setprecision(20) << "Decompression: " << decompression_time << endl;
        cout << fixed << setprecision(20) << "Total Time: " << whole_time << endl;
        */
    }
    return 0;
}